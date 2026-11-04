"""Headless integration regressions for the UI's worker and data pipelines.

Tk widgets and file dialogs are replaced only at the UI boundary. The tests use
the real analyzers, table I/O, BLAST parser and license validator.
"""

from contextlib import contextmanager
import copy
import importlib.util
import os
from pathlib import Path
import queue
import shutil
import tempfile
import threading
import unittest
from unittest import mock

import pandas as pd

from blast_support import BlastError
from probe_io import read_probe_table, write_probe_table


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('fish_ui_pipeline', ROOT / '0C9.py')
APP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(APP)


def worker_ui():
    """Construct worker state without creating a Tk window or Tcl interpreter."""
    ui = APP.DNAProbeAnalyzerUI.__new__(APP.DNAProbeAnalyzerUI)
    ui.ui_update_queue = queue.Queue()
    ui.pause_flag = threading.Event()
    ui.cancel_flag = threading.Event()
    ui.design_cancel_flag = threading.Event()
    ui.design_results = []
    ui.root = object()
    return ui


def drain_messages(ui):
    messages = []
    while not ui.ui_update_queue.empty():
        messages.append(ui.ui_update_queue.get_nowait())
    return messages


def write_report(path, query_subjects):
    """Write a complete standard outfmt-7 fixture with 20-nt exact alignments."""
    lines = []
    for query_id, subjects in query_subjects.items():
        lines.extend(['# BLASTN 2.17.0+', f'# Query: {query_id}',
                      '# Database: regression_fixture', f'# {len(subjects)} hits found'])
        for subject in subjects:
            lines.append('\t'.join([
                query_id, subject, '100', '20', '0', '0', '1', '20',
                '101', '120', '1e-10', '40', '20', '100',
            ]))
    lines.append(f'# BLAST processed {len(query_subjects)} queries')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class MergePipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_merge_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.ui = worker_ui()
        sequence = 'ACGTACGTACGTACGTACGT'
        self.frame = pd.DataFrame({
            'id': [1, 'sample A', 'duplicate', 'duplicate', '001', 'short', 'N', 'zero'],
            'sequence': [sequence] * 5 + ['ACGT', 'ACGTNACGTACGTACGTACG', sequence],
            'valid_sequence': [True] * 6 + [False, True],
        })
        # Deliberately query rows out of input order. IDs and row labels cannot
        # substitute for this explicit mapping, particularly for duplicate IDs.
        self.query_map = {
            'q000001': 0, 'q000002': 2, 'q000003': 3,
            'q000004': 4, 'q000005': 7, 'q000006': 1,
        }
        self.report = self.directory / 'blast.txt'
        write_report(self.report, {
            query_id: ([] if position == 7 else [f'subject_{position}'])
            for query_id, position in self.query_map.items()
        })

    def assert_attribution(self, result):
        self.assertEqual(result['id'].tolist(),
                         ['1', 'sample A', 'duplicate', 'duplicate', '001', 'short', 'N', 'zero'])
        self.assertEqual(result.loc[:4, 'blast_hit_1'].tolist(),
                         [f'subject_{position}' for position in range(5)])
        self.assertEqual(result.loc[:4, 'blast_hits_count'].map(int).tolist(), [1] * 5)
        # A completed query with zero hits is different from a skipped row.
        self.assertEqual(int(result.loc[7, 'blast_hits_count']), 0)
        self.assertEqual(result.loc[7, 'specificity'], '无匹配（当前参数）')
        self.assertEqual(result.loc[7, 'blast_query_id'], 'q000005')
        self.assertEqual(result.loc[5, 'blast_hits_count'], '')
        self.assertEqual(result.loc[6, 'blast_hits_count'], '')
        self.assertIn('未检索', result.loc[5, 'specificity'])
        self.assertIn('长度不足', result.loc[5, 'specificity'])
        self.assertIn('无效序列', result.loc[6, 'specificity'])
        self.assertEqual(result.loc[[5, 6], 'blast_query_id'].tolist(), ['', ''])

    def test_csv_tsv_and_excel_merge_by_source_row_not_user_id(self):
        for extension in ('.csv', '.tsv', '.xlsx'):
            with self.subTest(extension=extension):
                source = self.directory / ('analysis' + extension)
                destination = self.directory / ('merged' + extension)
                write_probe_table(self.frame, source)
                before = source.read_bytes()
                self.assertTrue(self.ui.merge_blast_results(
                    source, self.report, destination, query_map=self.query_map,
                ))
                self.assert_attribution(read_probe_table(destination))
                self.assertEqual(source.read_bytes(), before)

    def test_merge_uses_dataframe_snapshot_without_reopening_analysis_file(self):
        snapshot = self.frame.copy(deep=True)
        snapshot.index = [101, 205, 309, 413, 517, 621, 725, 829]
        before = snapshot.copy(deep=True)
        destination = self.directory / 'from_memory.tsv'
        with mock.patch.object(APP, 'read_probe_table', side_effect=AssertionError('unexpected reread')):
            self.ui.merge_blast_results(
                self.directory / 'does_not_exist.xlsx', self.report, destination,
                query_map=self.query_map, results_df=snapshot,
            )
        pd.testing.assert_frame_equal(snapshot, before)
        self.assert_attribution(read_probe_table(destination))

    def test_parse_failure_does_not_overwrite_previous_merged_output(self):
        destination = self.directory / 'merged.csv'
        destination.write_text('previous complete output', encoding='utf-8')
        self.report.write_text('# BLASTN 2.17.0+\n# Query: q000001\n', encoding='utf-8')
        before = self.frame.copy(deep=True)
        with self.assertRaises(BlastError):
            self.ui.merge_blast_results(
                'unused.csv', self.report, destination,
                query_map=self.query_map, results_df=self.frame,
            )
        self.assertEqual(destination.read_text(encoding='utf-8'), 'previous complete output')
        pd.testing.assert_frame_equal(self.frame, before)


class AnalysisWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_analysis_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.input_file = self.directory / 'input.csv'
        self.input_file.write_text(
            'id,sequence\n001,CGTTCCAAAGATGTGGGCATGAGCTTAC\nmissing,\n', encoding='utf-8',
        )

    def test_real_analysis_statistics_and_blast_use_in_memory_results(self):
        original_reader = APP.read_probe_table
        for extension in ('.csv', '.tsv', '.xlsx'):
            with self.subTest(extension=extension):
                ui = worker_ui()
                ui.analyzer = APP.DNAProbeAnalyzer()
                output = self.directory / ('analysis' + extension)
                merged = self.directory / ('analysis_with_blast' + extension)
                ui.run_blast_analysis = mock.Mock(return_value=str(merged))
                config = {
                    'input_file': str(self.input_file), 'output_file': str(output),
                    'tm_method': 'santalucia', 'run_blast': True,
                    'blast_output': str(self.directory / 'blast.txt'),
                }

                def input_only(path, *args, **kwargs):
                    if Path(path) != self.input_file:
                        raise AssertionError('The UI reopened an analysis output instead of using results_df')
                    return original_reader(path, *args, **kwargs)

                with mock.patch.object(APP, 'read_probe_table', side_effect=input_only) as reader:
                    ui.run_analysis(config)
                self.assertEqual(reader.call_count, 1)
                self.assertTrue(output.exists())
                ui.run_blast_analysis.assert_called_once()
                passed_config, passed_frame = ui.run_blast_analysis.call_args.args
                self.assertEqual(passed_config, config)
                self.assertIsNot(passed_frame, ui.analyzer.results_df)
                pd.testing.assert_frame_equal(passed_frame, ui.analyzer.results_df)
                messages = drain_messages(ui)
                statistics = [args for task, args in messages if task == 'update_results']
                self.assertEqual(statistics, [['2', '1', '60.32°C', '50.00%']])
                self.assertIn(('analysis_blast_phase', []), messages)
                self.assertIn(('update_status', ['分析完成']), messages)
                self.assertIn(('enable_buttons', []), messages)
                self.assertFalse(any(task == 'show_message' and args[0] == 'showerror'
                                     for task, args in messages))

    def test_empty_results_show_no_mean_instead_of_zero(self):
        self.input_file.write_text('id,sequence\n', encoding='utf-8')
        ui = worker_ui()
        ui.analyzer = APP.DNAProbeAnalyzer()
        ui.run_blast_analysis = mock.Mock()
        ui.run_analysis({
            'input_file': str(self.input_file),
            'output_file': str(self.directory / 'empty.tsv'),
            'tm_method': 'santalucia', 'run_blast': False,
        })
        statistics = [args for task, args in drain_messages(ui) if task == 'update_results']
        self.assertEqual(statistics, [['0', '0', '—', '—']])
        ui.run_blast_analysis.assert_not_called()


class DesignResultTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_design_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.ui = worker_ui()
        self.probes = [{
            'id': 1, 'sequence': 'ACGTACGTACGTACGTACGT',
            'specificity': 'old status', 'blast_hits_count': 4,
            'blast_obsolete': 'old metadata',
            'blast_hits': [{'subject': 'old_target', 'identity': 100}],
        }]
        self.ui.design_results = self.probes
        self.report = self.directory / 'blast.txt'

    def test_zero_hit_rerun_clears_old_evidence_in_a_new_result(self):
        before = copy.deepcopy(self.probes)
        write_report(self.report, {'q000001': []})
        current = self.ui.parse_blast_results(self.report, self.probes, {'q000001': 0})
        self.assertIsNot(current, self.probes)
        self.assertIsNot(current[0], self.probes[0])
        self.assertEqual(self.probes, before)
        self.assertEqual(current[0]['blast_hits'], [])
        self.assertEqual(current[0]['blast_hits_count'], 0)
        self.assertEqual(current[0]['blast_hsp_count'], 0)
        self.assertEqual(current[0]['specificity'], '无匹配（当前参数）')
        self.assertNotIn('blast_obsolete', current[0])
        self.assertEqual(current[0]['sequence'], before[0]['sequence'])

    def test_malformed_rerun_does_not_mutate_design_results(self):
        before = copy.deepcopy(self.probes)
        self.report.write_text('# Query: q000001\n# 0 hits found\n', encoding='utf-8')
        with self.assertRaises(BlastError):
            self.ui.parse_blast_results(self.report, self.probes, {'q000001': 0})
        self.assertIs(self.ui.design_results, self.probes)
        self.assertEqual(self.ui.design_results, before)

    def test_detail_uses_source_snapshot_and_counts_uracil_pairing(self):
        self.ui.design_results = [{
            'id': 1, 'sequence': 'GCAT', 'rna_fragment': 'AUGC',
            'start': 1, 'end': 4, 'gc_content': 50.0, 'tm': 10.0,
            'complexity': 1.0, 'tm_method': 'rna_dna', 'specificity': '未检查',
        }]
        self.ui.design_source_sequence = 'AUGC'
        self.ui.seq_input = mock.Mock()
        self.ui.seq_input.get.return_value = 'AAAA'  # User edited the input after design.
        self.ui.design_tree = mock.Mock()
        self.ui.design_tree.selection.return_value = ('selected',)
        self.ui.design_tree.item.return_value = ('1',)
        self.ui.probe_detail_text = mock.Mock()
        self.ui.draw_mrna_location_probes = mock.Mock()
        self.ui.on_probe_select(None)
        detail = self.ui.probe_detail_text.insert.call_args.args[1]
        self.assertIn("目标RNA: 5'-AUGC-3'", detail)
        self.assertIn('100.0% (4/4 碱基匹配)', detail)
        self.ui.seq_input.get.assert_not_called()
        self.ui.draw_mrna_location_probes.assert_called_once_with(self.ui.design_results[0], 'AUGC')


class FastaSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_fasta_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.fasta = self.directory / 'targets.fasta'
        self.fasta.write_text('>target_1\nATCG\n>target_2\nGCTA\n', encoding='utf-8')
        self.ui = worker_ui()
        self.ui.seq_input = mock.Mock()

    def test_multiple_records_load_only_the_selected_target(self):
        with mock.patch.object(APP.filedialog, 'askopenfilename', return_value=str(self.fasta)), \
                mock.patch.object(APP.simpledialog, 'askinteger', return_value=2) as selection, \
                mock.patch.object(APP.messagebox, 'showerror') as error:
            self.ui.load_fasta()
        self.ui.seq_input.insert.assert_called_once_with(APP.tk.END, 'GCTA')
        self.assertEqual(selection.call_args.kwargs['maxvalue'], 2)
        error.assert_not_called()

    def test_cancelled_record_selection_keeps_existing_input(self):
        with mock.patch.object(APP.filedialog, 'askopenfilename', return_value=str(self.fasta)), \
                mock.patch.object(APP.simpledialog, 'askinteger', return_value=None), \
                mock.patch.object(APP.messagebox, 'showerror') as error:
            self.ui.load_fasta()
        self.ui.seq_input.delete.assert_not_called()
        self.ui.seq_input.insert.assert_not_called()
        error.assert_not_called()


class ActiveJobPathTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_paths_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        (self.directory / 'nested').mkdir()
        self.ui = worker_ui()
        self.ui.active_job_paths = {}
        self.analysis_output = self.directory / 'analysis.csv'
        self.design_output = self.directory / 'design_blast.txt'

    def test_reserved_file_and_normalized_alias_block_another_job(self):
        self.ui._reserve_job_paths('analysis', [self.analysis_output])
        reserved = copy.deepcopy(self.ui.active_job_paths)
        aliases = [self.analysis_output, self.directory / 'nested' / '..' / 'analysis.csv']
        for path in aliases:
            with self.subTest(path=path):
                with self.assertRaises(ValueError):
                    self.ui._assert_paths_available([path])
                with self.assertRaises(ValueError):
                    self.ui._reserve_job_paths('design_blast', [path])
                self.assertEqual(self.ui.active_job_paths, reserved)

    def test_distinct_files_can_be_reserved_by_concurrent_jobs(self):
        self.ui._reserve_job_paths('analysis', [self.analysis_output])
        self.ui._reserve_job_paths('design_blast', [self.design_output])
        self.assertEqual(set(self.ui.active_job_paths), {'analysis', 'design_blast'})
        self.assertFalse(self.ui.active_job_paths['analysis'] & self.ui.active_job_paths['design_blast'])
        free = self.directory / 'separate_export.xlsx'
        self.assertEqual(self.ui._assert_paths_available([free]),
                         {os.path.normcase(os.path.realpath(free))})

    def test_completion_releases_only_its_own_job_paths(self):
        for button in ('start_button', 'pause_button', 'cancel_button', 'design_btn',
                       'export_btn', 'blast_btn', 'design_cancel_btn'):
            setattr(self.ui, button, mock.Mock())
        self.ui.analysis_busy = self.ui.design_busy = True
        self.ui._reserve_job_paths('analysis', [self.analysis_output])
        self.ui._reserve_job_paths('design_blast', [self.design_output])

        self.ui._dispatch_ui_task('enable_buttons', [])
        self.assertEqual(set(self.ui.active_job_paths), {'design_blast'})
        self.assertFalse(self.ui.analysis_busy)
        self.assertTrue(self.ui.design_busy)
        with self.assertRaises(ValueError):
            self.ui._assert_paths_available([self.design_output])

        # A subsequent analysis may reuse its released file while design BLAST
        # is still running. Design completion must preserve that reservation.
        self.ui._reserve_job_paths('analysis', [self.analysis_output])
        self.ui.analysis_busy = True
        self.ui._dispatch_ui_task('enable_design_buttons', [])
        self.assertEqual(set(self.ui.active_job_paths), {'analysis'})
        self.assertFalse(self.ui.design_busy)
        self.assertTrue(self.ui.analysis_busy)
        self.ui._assert_paths_available([self.design_output])
        with self.assertRaises(ValueError):
            self.ui._assert_paths_available([self.analysis_output])


class LicenseLocationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fish_ui_license_')
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.script_directory = self.directory / 'program'
        self.script_directory.mkdir()
        self.other_directory = self.directory / 'other_working_directory'
        self.other_directory.mkdir()
        self.script = self.script_directory / '0C9.py'
        shutil.copyfile(ROOT / '0C9.py', self.script)
        self.readme = self.script_directory / 'README.md'
        self.original = (ROOT / 'README.md').read_bytes().replace(b'\r\n', b'\n')

    def test_script_relative_license_accepts_lf_and_crlf_from_another_cwd(self):
        for content in (self.original, self.original.replace(b'\n', b'\r\n')):
            with self.subTest(newlines='CRLF' if b'\r\n' in content else 'LF'):
                self.readme.write_bytes(content)
                with mock.patch.object(APP, '__file__', str(self.script)), \
                        working_directory(self.other_directory):
                    self.assertTrue(APP._obfuscated_license_check())

    def test_changed_license_content_is_rejected(self):
        self.readme.write_bytes(self.original + b'\nmodified content\n')
        with mock.patch.object(APP, '__file__', str(self.script)), \
                working_directory(self.other_directory):
            self.assertFalse(APP._obfuscated_license_check())

    def test_readme_in_working_directory_cannot_replace_missing_adjacent_license(self):
        (self.other_directory / 'README.md').write_bytes(self.original)
        with mock.patch.object(APP, '__file__', str(self.script)), \
                working_directory(self.other_directory):
            self.assertFalse(APP._obfuscated_license_check())


if __name__ == '__main__':
    unittest.main()

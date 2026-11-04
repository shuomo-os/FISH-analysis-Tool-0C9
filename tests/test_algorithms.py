"""Regression tests for probe calculations and the analysis file pipeline.

Run from the repository root: python -m unittest discover -s tests -v
These tests use real Biopython calculations; no BLAST executable or GUI is needed.
"""

import importlib.util
from pathlib import Path
import tempfile
import threading
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('fish_algorithms', ROOT / '0C9.py')
APP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(APP)


def parameters(**changes):
    values = {
        'probe_length': 20,
        'min_gc': 0,
        'max_gc': 100,
        'min_tm': -200,
        'max_tm': 200,
        'spacing': 0,
        'min_complexity': 0,
        'filter_repeats': False,
        'tm_method': 'santalucia',
        'max_homopolymer_length': 100,
    }
    values.update(changes)
    return values


class ThermodynamicsTests(unittest.TestCase):
    def setUp(self):
        self.analyzer = APP.DNAProbeAnalyzer()
        self.designer = APP.RNAProbeDesigner()

    def test_santalucia_documented_reference(self):
        # Biopython's reference example is 60.32 C with DNA_NN3.
        # The former DNA_NN1 table produced 72.19 C for this same sequence.
        sequence = 'CGTTCCAAAGATGTGGGCATGAGCTTAC'
        for calculator in (self.analyzer, self.designer):
            with self.subTest(calculator=type(calculator).__name__):
                self.assertEqual(calculator.calculate_tm(sequence, 'santalucia'), 60.32)
                self.assertEqual(calculator.calculate_tm(sequence, 'nn'), 60.32)

    def test_unknown_bases_are_not_deleted_before_tm(self):
        for calculator in (self.analyzer, self.designer):
            for method in calculator.TM_METHODS:
                with self.subTest(calculator=type(calculator).__name__, method=method):
                    self.assertIsNone(calculator.calculate_tm('ACGTNACGT', method))
                    self.assertIsNone(calculator.calculate_tm('NAN', method))

    def test_missing_values_and_blank_sequences_are_invalid(self):
        for value in (None, float('nan'), pd.NA, '', '  ', 123):
            for calculator in (self.analyzer, self.designer):
                with self.subTest(value=repr(value), calculator=type(calculator).__name__):
                    self.assertFalse(calculator.check_sequence_validity(value))
                    self.assertIsNone(calculator.calculate_tm(value))

    def test_gc_uses_complete_known_sequence(self):
        self.assertEqual(self.analyzer.calculate_gc_content('  acGT  '), 50.0)
        self.assertEqual(self.designer.calculate_gc('aCGU'), 50.0)
        self.assertIsNone(self.analyzer.calculate_gc_content('ACGTN'))
        self.assertIsNone(self.designer.calculate_gc('ACGUN'))

    def test_rna_dna_uses_target_orientation(self):
        target = 'AUGCUACGGAUUCAGUGCAG'
        # RNA/DNA Tm for this target is 50.77 C. Feeding the reverse-complement
        # DNA probe as the RNA strand would produce 49.66 C and fail this gate.
        probes = self.designer.design_probes(
            target, parameters(tm_method='rna_dna', min_tm=50.5, max_tm=51.0)
        )
        self.assertEqual(len(probes), 1)
        self.assertEqual(probes[0]['sequence'], 'CTGCACTGAATCCGTAGCAT')
        self.assertEqual(probes[0]['rna_fragment'], target)
        self.assertEqual(probes[0]['tm'], 50.77)
        self.assertEqual(probes[0]['tm_method'], 'rna_dna')

    def test_dna_analyzer_does_not_silently_use_rna_model(self):
        self.assertIsNone(self.analyzer.calculate_tm('ACGTACGTACGT', 'rna_dna'))
        self.assertIsNone(self.analyzer.calculate_tm('ACGTACGTACGT', 'unknown'))


class DesignWindowTests(unittest.TestCase):
    def setUp(self):
        self.designer = APP.RNAProbeDesigner()

    def test_target_exactly_one_probe_long_is_scanned(self):
        progress = []
        probes = self.designer.design_probes(
            'AUGCUACGGAUUCAGUGCAG', parameters(), progress_callback=progress.append
        )
        self.assertEqual(len(probes), 1)
        self.assertEqual((probes[0]['start'], probes[0]['end']), (1, 20))
        self.assertEqual(probes[0]['sequence'], 'CTGCACTGAATCCGTAGCAT')
        self.assertEqual(progress[-1], 100)

    def test_last_window_is_scanned_after_earlier_windows_fail(self):
        probes = self.designer.design_probes(
            'AAACGCG', parameters(probe_length=4, min_gc=100, tm_method='wallace')
        )
        self.assertEqual([(p['start'], p['end']) for p in probes], [(4, 7)])

    def test_spacing_counts_uncovered_bases(self):
        probes = self.designer.design_probes(
            'ACGT' * 7, parameters(probe_length=4, spacing=3, tm_method='wallace')
        )
        self.assertEqual([(p['start'], p['end']) for p in probes],
                         [(1, 4), (8, 11), (15, 18), (22, 25)])
        self.assertTrue(all(right['start'] - left['end'] - 1 == 3
                            for left, right in zip(probes, probes[1:])))

    def test_zero_spacing_allows_adjacent_probes(self):
        probes = self.designer.design_probes(
            'ACGTACGT', parameters(probe_length=4, spacing=0, tm_method='wallace')
        )
        self.assertEqual([(p['start'], p['end']) for p in probes], [(1, 4), (5, 8)])

    def test_overlapping_internal_repeats_are_detected(self):
        sequence = 'ATATATATGCCGTAGCTACG'
        self.assertTrue(self.designer.has_repeats(sequence, sequence))
        self.assertFalse(self.designer.has_repeats('ACGTTGCA', 'ACGTTGCA'))

    def test_overlapping_target_occurrences_are_detected(self):
        self.assertTrue(self.designer.has_repeats('ATATAT', 'ATATATAT', min_repeat_length=8))

    def test_complexity_is_normalized_to_zero_through_one(self):
        self.assertEqual(self.designer.calculate_complexity('ACGU'), 1.0)
        self.assertEqual(self.designer.calculate_complexity('AAAA'), 0.0)

    def test_invalid_parameter_ranges_raise_explanatory_errors(self):
        for changes in (
            {'probe_length': 1}, {'probe_length': 19.5}, {'spacing': -1},
            {'min_complexity': 1.01}, {'min_complexity': -0.1},
            {'min_gc': 70, 'max_gc': 60}, {'max_gc': 101},
            {'min_tm': 70, 'max_tm': 60}, {'min_tm': float('nan')},
            {'max_tm': float('inf')}, {'max_homopolymer_length': 0},
            {'tm_method': 'unknown'}, {'spacing': True},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.designer.validate_parameters(parameters(**changes))

    def test_invalid_or_short_target_is_rejected(self):
        for sequence in ('', 'ACGTNACGT', 'ACGT'):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                self.designer.design_probes(sequence, parameters())

    def test_cancelled_design_does_not_return_partial_or_old_results(self):
        flag = threading.Event()
        self.designer.results = [{'id': 'old'}]
        probes = self.designer.design_probes(
            'AUGCUACGGAUUCAGUGCAG', parameters(),
            progress_callback=lambda value: flag.set(), cancel_flag=flag,
        )
        self.assertEqual(probes, [])
        self.assertEqual(self.designer.results, [])


class AnalysisPipelineTests(unittest.TestCase):
    def setUp(self):
        self.analyzer = APP.DNAProbeAnalyzer()
        self.tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def config(self, input_file, output_file, **changes):
        config = {'input_file': str(input_file), 'output_file': str(output_file),
                  'tm_method': 'santalucia'}
        config.update(changes)
        return config

    def test_empty_and_ambiguous_sequences_stay_invalid_and_ids_survive(self):
        input_file = self.directory / 'input.csv'
        input_file.write_text(
            'id,sequence\n001,ACGTACGT\nsample A,\n1,ACGTNACGT\n1,nan\nNA,aCgT\n',
            encoding='utf-8',
        )
        output_file = self.directory / 'results.csv'
        self.assertTrue(self.analyzer.analyze_probes(self.config(input_file, output_file)))
        result = self.analyzer.results_df
        self.assertEqual(result['id'].tolist(), ['001', 'sample A', '1', '1', 'NA'])
        self.assertEqual(result['valid_sequence'].tolist(), [True, False, False, False, True])
        self.assertEqual(result['alphabet_valid'].tolist(), [True, False, True, True, True])
        self.assertEqual(result['tm_calculable'].tolist(), [True, False, False, False, True])
        self.assertEqual(result.loc[1, 'sequence'], '')
        self.assertEqual(result.loc[4, 'input_sequence'], 'aCgT')
        self.assertTrue(result.loc[[1, 2, 3], ['tm', 'gc_content']].isna().all().all())
        self.assertTrue(all(result.loc[[1, 2, 3], 'sequence_error']))

    def test_csv_tsv_and_excel_analysis_pipelines(self):
        from probe_io import read_probe_table, write_probe_table

        for input_ext, output_ext in (('.csv', '.csv'), ('.tsv', '.xlsx'), ('.xlsx', '.tsv')):
            with self.subTest(input_ext=input_ext, output_ext=output_ext):
                input_file = self.directory / ('input' + input_ext)
                output_file = self.directory / ('result' + output_ext)
                write_probe_table(pd.DataFrame({
                    'id': ['001', 'sample A'],
                    'sequence': ['ACGTACGTACGT', 'CGTTCCAAAGATGTGGGCATGAGCTTAC'],
                }), input_file)
                self.assertTrue(self.analyzer.analyze_probes(self.config(input_file, output_file)))
                result = read_probe_table(output_file)
                self.assertEqual(result['id'].tolist(), ['001', 'sample A'])
                self.assertEqual(result['valid_sequence'].tolist(), [True, True])
                self.assertAlmostEqual(float(result.loc[1, 'tm']), 60.32, places=2)

    def test_empty_table_has_headers(self):
        input_file = self.directory / 'empty.csv'
        input_file.write_text('id,sequence\n', encoding='utf-8')
        output_file = self.directory / 'output.csv'
        self.assertTrue(self.analyzer.analyze_probes(self.config(input_file, output_file)))
        self.assertEqual(len(self.analyzer.results_df), 0)
        self.assertIn('valid_sequence', pd.read_csv(output_file).columns)

    def test_cancel_before_save_leaves_existing_output_and_clears_results(self):
        input_file = self.directory / 'input.csv'
        input_file.write_text('id,sequence\n1,ACGTACGT\n', encoding='utf-8')
        output_file = self.directory / 'output.csv'
        output_file.write_text('existing output', encoding='utf-8')
        self.analyzer.results_df = pd.DataFrame({'id': ['old']})
        self.analyzer.results = [{'id': 'old'}]
        flag = threading.Event()
        success = self.analyzer.analyze_probes(
            self.config(input_file, output_file), cancel_flag=flag,
            progress_callback=lambda value: flag.set(),
        )
        self.assertFalse(success)
        self.assertIsNone(self.analyzer.results_df)
        self.assertEqual(self.analyzer.results, [])
        self.assertEqual(output_file.read_text(encoding='utf-8'), 'existing output')

    def test_missing_sequence_column_fails_without_old_results(self):
        input_file = self.directory / 'invalid.csv'
        input_file.write_text('id,other\n1,ACGT\n', encoding='utf-8')
        output_file = self.directory / 'output.csv'
        self.analyzer.results_df = pd.DataFrame({'id': ['old']})
        self.assertFalse(self.analyzer.analyze_probes(self.config(input_file, output_file)))
        self.assertIsNone(self.analyzer.results_df)
        self.assertFalse(output_file.exists())


if __name__ == '__main__':
    unittest.main()

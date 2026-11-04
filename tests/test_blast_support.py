"""Regressions for query attribution, incomplete reports, and process cancellation."""

from pathlib import Path
import os
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock

from blast_support import (
    BLAST_OUTFMT,
    BlastCancelled,
    BlastError,
    build_blast_command,
    parse_blast_results,
    publish_blast_report,
    run_blast,
    summarize_hits,
    write_query_fasta,
)


def alignment(query="q000001", subject="target", **changes):
    fields = [query, subject, "100.000", "20", "0", "0", "1", "20",
              "101", "120", "1e-5", "40", "20", "100"]
    positions = {
        "identity": 2, "length": 3, "mismatch": 4, "gapopen": 5,
        "qstart": 6, "qend": 7, "sstart": 8, "send": 9,
        "evalue": 10, "bitscore": 11, "qlen": 12, "qcovhsp": 13,
    }
    for name, value in changes.items():
        fields[positions[name]] = str(value)
    return "\t".join(fields)


def report(query_hits):
    text = []
    for query, hits in query_hits:
        text.extend(["# BLASTN 2.17.0+", f"# Query: {query}", "# Database: fixture",
                     f"# {len(hits)} hits found", *hits])
    text.append(f"# BLAST processed {len(query_hits)} queries")
    return "\n".join(text) + "\n"


class BlastSupportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="fish blast test ")
        self.directory = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def parse(self, text, expected=("q000001",)):
        path = self.directory / "blast results.txt"
        path.write_text(text, encoding="utf-8")
        return parse_blast_results(path, expected)

    def test_numeric_spaced_duplicate_and_padded_ids_map_to_original_rows(self):
        # IDs stay in the source table; only row positions identify BLAST queries.
        rows = [
            {"id": 1, "sequence": "ACGTACGTACGTACGTACGT"},
            {"id": "001", "sequence": ""},  # Invalid row skipped before BLAST.
            {"id": "sample A", "sequence": "GCTAGCTAGCTAGCTAGCTA"},
            {"id": "sample A", "sequence": "CGATCGATCGATCGATCGAT"},
            {"id": "001", "sequence": "TCGATCGATCGATCGATCGA"},
        ]
        path = self.directory / "queries.fasta"
        query_rows = write_query_fasta(
            [(position, row["sequence"]) for position, row in enumerate(rows) if row["sequence"]],
            path,
        )
        self.assertEqual(query_rows, {"q000001": 0, "q000002": 2, "q000003": 3, "q000004": 4})
        self.assertNotIn("sample A", path.read_text(encoding="utf-8"))
        parsed = self.parse(report([
            (query, [alignment(query, f"subject_{position}")])
            for query, position in query_rows.items()
        ]), query_rows)
        attributed = {position: parsed[query][0]["subject"] for query, position in query_rows.items()}
        self.assertEqual(attributed, {0: "subject_0", 2: "subject_2", 3: "subject_3", 4: "subject_4"})
        self.assertEqual(rows[4]["id"], "001")

    def test_invalid_queries_do_not_truncate_existing_destination(self):
        path = self.directory / "queries.fasta"
        path.write_text("previous", encoding="utf-8")
        for queries in ([], [(0, "ACGT")], [(0, "ACGTNACGT")], [(0, None)],
                        [(0, "ACGTACGT"), (0, "ACGTACGT")]):
            with self.subTest(queries=queries), self.assertRaises(ValueError):
                write_query_fasta(queries, path)
            self.assertEqual(path.read_text(encoding="utf-8"), "previous")

    def test_command_explicitly_searches_short_sequences_without_masking(self):
        command = build_blast_command("folder with spaces/blastn", "database prefix", "query.fa", "out.txt")
        options = dict(zip(command[1:-1:2], command[2:-1:2]))
        self.assertEqual(command[0], "folder with spaces/blastn")
        self.assertEqual(options["-task"], "blastn-short")
        self.assertEqual(options["-word_size"], "7")
        self.assertEqual(options["-dust"], "no")
        self.assertEqual(options["-soft_masking"], "false")
        self.assertEqual(options["-evalue"], "1000")
        self.assertEqual(options["-outfmt"], BLAST_OUTFMT)
        self.assertEqual(options["-max_target_seqs"], "1000")
        self.assertEqual(options["-db"], '"database prefix"')
        self.assertEqual(command[-1], "-parse_deflines")

    def test_database_prefix_is_quoted_once_for_ncbi_database_list_parser(self):
        for prefix in ('C:\\data folder\\my_db', '"C:\\data folder\\my_db"'):
            with self.subTest(prefix=prefix):
                command = build_blast_command("blastn", prefix, "query.fa", "out.txt")
                self.assertEqual(command[command.index("-db") + 1], '"C:\\data folder\\my_db"')

    def test_zero_hits_is_distinct_from_a_missing_query(self):
        parsed = self.parse(report([
            ("q000001", [alignment()]), ("q000002", []),
        ]), ("q000001", "q000002"))
        self.assertEqual(parsed["q000002"], [])
        self.assertEqual(summarize_hits(parsed["q000002"])["specificity"], "无匹配（当前参数）")
        with self.assertRaises(BlastError):
            self.parse(report([("q000001", [alignment()])]), ("q000001", "q000002"))

    def test_explicit_local_id_spelling_and_query_description_are_supported(self):
        text = report([("lcl|q000001 description", [alignment("lcl|q000001")])])
        self.assertEqual(self.parse(text)["q000001"][0]["subject"], "target")

    def test_empty_truncated_and_corrupt_reports_raise_instead_of_zero_hits(self):
        full_report = report([("q000001", [alignment()])])
        malformed = [
            "", "# BLASTN 2.17.0+\n", full_report.rsplit("# BLAST processed", 1)[0],
            full_report.replace("# 1 hits found", "# 2 hits found"),
            full_report.replace("# 1 hits found\n", ""),
            full_report.replace("# BLAST processed 1 queries", "# BLAST processed 2 queries"),
            full_report.replace("q000001\ttarget", "unexpected\ttarget"),
            full_report.replace("# Query: q000001", "# Query: unexpected"),
            full_report.replace(alignment(), alignment().rsplit("\t", 1)[0]),
            full_report + "# Query: q000001\n# 0 hits found\n",
            report([("q000001", []), ("q000001", [])]),
            alignment() + "\n# BLAST processed 1 queries\n",
        ]
        for text in malformed:
            with self.subTest(text=text), self.assertRaises(BlastError):
                self.parse(text)

    def test_illegal_numbers_cannot_enter_results(self):
        for changes in (
            {"identity": "NaN"}, {"evalue": "inf"}, {"evalue": "bad"},
            {"identity": 101}, {"evalue": -1}, {"bitscore": -1},
            {"qlen": 0}, {"qstart": 21}, {"qcovhsp": 101},
            {"length": 0}, {"mismatch": -1}, {"gapopen": -1},
        ):
            with self.subTest(changes=changes), self.assertRaises(BlastError):
                self.parse(report([("q000001", [alignment(**changes)])]))

    def test_full_matches_require_entire_query_coverage_without_gaps(self):
        lines = [
            alignment(subject="full"),
            alignment(subject="reverse", qstart=20, qend=1),
            alignment(subject="partial", length=10, qend=10, send=110, qcovhsp=50),
            alignment(subject="rounded", length=199, qend=199, send=299, qlen=200, qcovhsp=100),
            alignment(subject="gapped", identity=95.238, length=21, gapopen=1, send=121),
        ]
        # Query lengths must stay constant within one query; test rounded coverage separately.
        rounded = lines.pop(3)
        hits = self.parse(report([("q000001", lines)]))["q000001"]
        self.assertEqual([hit["full_length_exact"] for hit in hits], [True, True, False, False])
        self.assertFalse(self.parse(report([("q000001", [rounded])]))["q000001"][0]["full_length_exact"])

    def test_query_length_must_be_consistent_across_hsps(self):
        with self.assertRaises(BlastError):
            self.parse(report([("q000001", [alignment(), alignment(qlen=21, qcovhsp=95)])]))

    def test_multiple_perfect_subjects_never_receive_high_specificity_label(self):
        hits = self.parse(report([("q000001", [
            alignment(subject="isoform_1"), alignment(subject="isoform_1"),
            alignment(subject="isoform_2"), alignment(subject="other_gene"),
            alignment(subject="another_gene"),
        ])]))["q000001"]
        summary = summarize_hits(hits, max_target_seqs=4)
        self.assertEqual(summary["blast_hits_count"], 4)
        self.assertEqual(summary["blast_hsp_count"], 5)
        self.assertEqual(summary["blast_full_match_count"], 4)
        self.assertTrue(summary["blast_hit_limit_reached"])
        self.assertIn("待核实", summary["specificity"])
        self.assertNotIn("高", summary["specificity"])

    def test_subsequent_zero_hit_parse_returns_no_previous_hits(self):
        initial = self.parse(report([("q000001", [alignment()])]))
        current = self.parse(report([("q000001", [])]))
        self.assertEqual(len(initial["q000001"]), 1)
        self.assertEqual(current, {"q000001": []})

    def staged_reports(self, existing_mapping=True):
        output = self.directory / "published.txt"
        mapping = Path(str(output) + ".queries.csv")
        raw_source = self.directory / "staged.txt"
        mapping_source = self.directory / "staged.csv"
        output.write_text("old report", encoding="utf-8")
        if existing_mapping:
            mapping.write_text("old mapping", encoding="utf-8")
        raw_source.write_text("new report", encoding="utf-8")
        mapping_source.write_text("new mapping", encoding="utf-8")
        return output, mapping, raw_source, mapping_source

    def test_publishing_replaces_both_files_and_removes_backup(self):
        output, mapping, raw_source, mapping_source = self.staged_reports()
        publish_blast_report(raw_source, mapping_source, output)
        self.assertEqual(output.read_text(encoding="utf-8"), "new report")
        self.assertEqual(mapping.read_text(encoding="utf-8"), "new mapping")
        self.assertEqual(list(self.directory.glob("*.backup-*")), [])

    def test_raw_publication_failure_restores_previous_mapping(self):
        output, mapping, raw_source, mapping_source = self.staged_reports()
        real_replace = os.replace

        def replace(source, destination):
            if Path(destination) == output:
                raise PermissionError("raw report is locked")
            return real_replace(source, destination)

        with mock.patch("blast_support.os.replace", side_effect=replace):
            with self.assertRaisesRegex(BlastError, "报告发布失败"):
                publish_blast_report(raw_source, mapping_source, output)
        self.assertEqual(output.read_text(encoding="utf-8"), "old report")
        self.assertEqual(mapping.read_text(encoding="utf-8"), "old mapping")
        self.assertEqual(list(self.directory.glob("*.backup-*")), [])

    def test_raw_publication_failure_removes_new_mapping_if_none_existed(self):
        output, mapping, raw_source, mapping_source = self.staged_reports(existing_mapping=False)
        real_replace = os.replace

        def replace(source, destination):
            if Path(destination) == output:
                raise PermissionError("raw report is locked")
            return real_replace(source, destination)

        with mock.patch("blast_support.os.replace", side_effect=replace):
            with self.assertRaises(BlastError):
                publish_blast_report(raw_source, mapping_source, output)
        self.assertEqual(output.read_text(encoding="utf-8"), "old report")
        self.assertFalse(mapping.exists())

    def test_mapping_publication_failure_keeps_previous_pair(self):
        output, mapping, raw_source, mapping_source = self.staged_reports()
        with mock.patch("blast_support.os.replace", side_effect=PermissionError("mapping is locked")):
            with self.assertRaises(BlastError):
                publish_blast_report(raw_source, mapping_source, output)
        self.assertEqual(output.read_text(encoding="utf-8"), "old report")
        self.assertEqual(mapping.read_text(encoding="utf-8"), "old mapping")

    def test_backup_failure_keeps_previous_pair(self):
        output, mapping, raw_source, mapping_source = self.staged_reports()
        with mock.patch("blast_support.shutil.copy2", side_effect=OSError("backup disk full")):
            with self.assertRaises(BlastError):
                publish_blast_report(raw_source, mapping_source, output)
        self.assertEqual(output.read_text(encoding="utf-8"), "old report")
        self.assertEqual(mapping.read_text(encoding="utf-8"), "old mapping")

    def test_rollback_failure_preserves_backup_and_reports_inconsistent_pair(self):
        output, mapping, raw_source, mapping_source = self.staged_reports()
        real_replace = os.replace

        def replace(source, destination):
            if Path(destination) == output or ".backup-" in Path(source).name:
                raise PermissionError("file is locked")
            return real_replace(source, destination)

        with mock.patch("blast_support.os.replace", side_effect=replace):
            with self.assertRaisesRegex(BlastError, "回滚失败.*旧映射备份保留在") as raised:
                publish_blast_report(raw_source, mapping_source, output)
        backups = list(self.directory.glob("*.backup-*"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_text(encoding="utf-8"), "old mapping")
        self.assertIn(str(backups[0]), str(raised.exception))
        self.assertEqual(output.read_text(encoding="utf-8"), "old report")
        self.assertEqual(mapping.read_text(encoding="utf-8"), "new mapping")

    def test_runner_drains_output_and_surfaces_nonzero_exit(self):
        result = run_blast([sys.executable, "-c", "import sys;sys.stdout.write('x'*100000);sys.stderr.write('warning')"])
        self.assertEqual(len(result.stdout), 100000)
        self.assertEqual(result.stderr, "warning")
        with self.assertRaisesRegex(BlastError, "退出码 3.*fixture error"):
            run_blast([sys.executable, "-c", "import sys;sys.stderr.write('fixture error');sys.exit(3)"])

    def test_command_and_output_paths_with_spaces_work_in_a_real_child(self):
        fake_program = self.directory / "fake blast.py"
        fake_program.write_text(
            "import pathlib,sys\n"
            "args=sys.argv\n"
            "assert args[args.index('-task')+1]=='blastn-short'\n"
            "query=pathlib.Path(args[args.index('-query')+1]).read_text()\n"
            "ids=[s[1:] for s in query.splitlines() if s.startswith('>')]\n"
            "out=[]\n"
            "for name in ids:\n"
            " out += ['# BLASTN fake', '# Query: '+name, '# 0 hits found']\n"
            "out += ['# BLAST processed '+str(len(ids))+' queries']\n"
            "pathlib.Path(args[args.index('-out')+1]).write_text('\\n'.join(out)+'\\n')\n",
            encoding="utf-8",
        )
        query = self.directory / "query with spaces.fasta"
        output = self.directory / "result with spaces.txt"
        mapping = write_query_fasta([(3, "ACGTACGTACGTACGTACGT")], query)
        command = [sys.executable, str(fake_program)] + build_blast_command("ignored", "db", query, output)[1:]
        run_blast(command)
        self.assertEqual(parse_blast_results(output, mapping), {"q000001": []})

    def test_cancellation_before_launch_does_not_create_process(self):
        event = threading.Event()
        event.set()
        with mock.patch("blast_support.subprocess.Popen") as popen:
            with self.assertRaises(BlastCancelled):
                run_blast([sys.executable, "-c", "pass"], event)
            popen.assert_not_called()

    def test_cancellation_reaps_real_running_process(self):
        event = threading.Event()
        children = []
        real_popen = subprocess.Popen

        def launch(*args, **kwargs):
            process = real_popen(*args, **kwargs)
            children.append(process)
            event.set()
            return process

        with mock.patch("blast_support.subprocess.Popen", side_effect=launch):
            with self.assertRaises(BlastCancelled):
                run_blast([sys.executable, "-c", "import time;time.sleep(30)"], event)
        self.assertEqual(len(children), 1)
        self.assertIsNotNone(children[0].poll())

    @unittest.skipIf(os.name == "nt", "POSIX signal-handler escalation check")
    def test_cancellation_kills_process_that_ignores_termination(self):
        event = threading.Event()
        ready = self.directory / "ready"
        children = []
        real_popen = subprocess.Popen

        def launch(*args, **kwargs):
            process = real_popen(*args, **kwargs)
            children.append(process)
            return process

        def cancel_after_ready():
            deadline = time.monotonic() + 3
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            event.set()

        watcher = threading.Thread(target=cancel_after_ready, daemon=True)
        watcher.start()
        code = (
            "import pathlib,signal,sys,time;"
            "signal.signal(signal.SIGTERM,signal.SIG_IGN);"
            "pathlib.Path(sys.argv[1]).write_text('ready');time.sleep(30)"
        )
        with mock.patch("blast_support.subprocess.Popen", side_effect=launch):
            with self.assertRaises(BlastCancelled):
                run_blast([sys.executable, "-c", code, str(ready)], event,
                          poll_interval=0.05, termination_timeout=0.1)
        watcher.join(timeout=4)
        self.assertTrue(ready.exists())
        self.assertEqual(children[0].returncode, -9)


if __name__ == "__main__":
    unittest.main()

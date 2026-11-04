"""Small real-BLAST regression test; skipped when BLAST+ is not installed."""

import importlib.util
from pathlib import Path
import queue
import shutil
import subprocess
import tempfile
import threading
import unittest

import pandas as pd

from blast_support import (
    build_blast_command,
    parse_blast_results,
    run_blast,
    summarize_hits,
    write_query_fasta,
)
from probe_io import read_probe_table, write_probe_table


QUERY = "ACGTCAGTGCATGACTAGCT"


def create_database(directory):
    """Build with makeblastdb, then relocate to a database path with spaces."""
    original = directory / "database"
    original.mkdir()
    reference = original / "reference.fasta"
    reference.write_text(
        f">subject_a\nGGCCGGCCTTCG{QUERY}CGGCTCCTCCGG\n"
        f">subject_b\nTGCCTCGGCCGT{QUERY}CGATGGCCTCCA\n",
        encoding="utf-8",
    )
    subprocess.run(
        [shutil.which("makeblastdb"), "-in", str(reference), "-dbtype", "nucl",
         "-parse_seqids", "-out", str(original / "tiny_database")],
        check=True, capture_output=True, text=True, timeout=30,
    )
    # BLAST interprets -db as a list even when subprocess passes one argument.
    # Moving the completed DB isolates that behavior from makeblastdb's parser.
    moved = directory / "database with spaces"
    original.rename(moved)
    return moved / "tiny_database"


@unittest.skipUnless(
    shutil.which("blastn") and shutil.which("makeblastdb"),
    "Real BLAST integration requires blastn and makeblastdb on PATH.",
)
class BlastIntegrationTests(unittest.TestCase):
    def test_twenty_base_queries_find_multiple_subjects_and_keep_zero_hits_separate(self):
        # These synthetic sequences contain no 7-base A/T homopolymer, so the
        # poly-A query has no seed on either strand, even with dust disabled.
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            database = create_database(directory)
            fasta = directory / "queries with spaces.fasta"
            report = directory / "results with spaces.tsv"
            query_rows = write_query_fasta(
                [(0, QUERY), (2, "A" * 20), (5, QUERY)], fasta
            )
            self.assertEqual(query_rows, {"q000001": 0, "q000002": 2, "q000003": 5})

            command = build_blast_command(
                shutil.which("blastn"), database, fasta, report
            )
            run_blast(command)
            results = parse_blast_results(report, query_rows)

            self.assertEqual(set(results), set(query_rows))
            self.assertEqual(results["q000002"], [])
            for query_id in ("q000001", "q000003"):
                with self.subTest(query_id=query_id):
                    exact = [hit for hit in results[query_id] if hit["full_length_exact"]]
                    subjects = {hit["subject"].removeprefix("lcl|") for hit in exact}
                    self.assertEqual(subjects, {"subject_a", "subject_b"})
                    for hit in exact:
                        self.assertEqual(hit["query_length"], 20)
                        self.assertEqual(hit["alignment_length"], 20)
                        self.assertEqual(hit["identity"], 100.0)
                        self.assertEqual(hit["query_coverage"], 100.0)
                    summary = summarize_hits(results[query_id])
                    self.assertEqual(summary["blast_hits_count"], 2)
                    self.assertEqual(summary["blast_full_match_count"], 2)
                    self.assertNotEqual(summary["specificity"], "高")

    def test_real_background_analysis_blast_and_merge_in_all_three_table_formats(self):
        # Construct only the worker state: no Tk window, widget mocks, fake
        # BLAST report, or replaced analysis/calculation code is involved.
        spec = importlib.util.spec_from_file_location(
            "fish_real_blast_pipeline", Path(__file__).resolve().parents[1] / "0C9.py"
        )
        app = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(app)
        ids = ["001", "duplicate", "duplicate", "sample A", "zero", "missing", "N", "short"]
        sequences = [QUERY] * 4 + ["A" * 20, "", "ACGTNACGTACGTACGTACG", "ACGT"]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            database = create_database(directory)
            tables = directory / "table files with spaces"
            tables.mkdir()
            for extension in (".csv", ".tsv", ".xlsx"):
                with self.subTest(extension=extension):
                    source = tables / ("input" + extension)
                    output = tables / ("analysis" + extension)
                    merged = tables / ("analysis_with_blast" + extension)
                    report = tables / ("blast" + extension + ".txt")
                    write_probe_table(pd.DataFrame({"id": ids, "sequence": sequences}), source)
                    ui = app.DNAProbeAnalyzerUI.__new__(app.DNAProbeAnalyzerUI)
                    ui.ui_update_queue = queue.Queue()
                    ui.pause_flag = threading.Event()
                    ui.cancel_flag = threading.Event()
                    ui.worker_threads = []
                    ui.analyzer = app.DNAProbeAnalyzer()
                    config = {
                        "input_file": str(source), "output_file": str(output),
                        "tm_method": "santalucia", "run_blast": True,
                        "blast_path": shutil.which("blastn"), "db_path": str(database),
                        "blast_output": str(report),
                    }
                    ui._start_worker(ui.run_analysis, config)
                    worker = ui.worker_threads[0]
                    worker.join(timeout=30)
                    if worker.is_alive():
                        ui.cancel_flag.set()
                        worker.join(timeout=5)
                        self.fail("The tiny-database analysis worker did not finish in 30 seconds")
                    messages = []
                    while not ui.ui_update_queue.empty():
                        messages.append(ui.ui_update_queue.get_nowait())
                    errors = [args for task, args in messages
                              if task == "show_message" and args[0] == "showerror"]
                    self.assertEqual(errors, [], messages)
                    self.assertIn(("update_status", ["分析完成"]), messages)
                    self.assertIn(("enable_buttons", []), messages)

                    analysis = read_probe_table(output)
                    result = read_probe_table(merged)
                    self.assertEqual(analysis["id"].tolist(), ids)
                    self.assertEqual(result["id"].tolist(), ids)
                    self.assertEqual(result.loc[:4, "blast_hits_count"].map(int).tolist(), [2, 2, 2, 2, 0])
                    self.assertEqual(result.loc[:4, "blast_full_match_count"].map(int).tolist(), [2, 2, 2, 2, 0])
                    self.assertEqual(result.loc[4, "specificity"], "无匹配（当前参数）")
                    self.assertEqual(result.loc[[5, 6, 7], "blast_query_id"].tolist(), ["", "", ""])
                    self.assertEqual(result.loc[[5, 6, 7], "blast_hits_count"].tolist(), ["", "", ""])
                    for position in (5, 6):
                        self.assertFalse(result.loc[position, "valid_sequence"])
                        self.assertIn("无效序列", result.loc[position, "specificity"])
                    self.assertIn("长度不足", result.loc[7, "specificity"])

                    sidecar = read_probe_table(str(report) + ".queries.csv")
                    self.assertEqual(sidecar["id"].tolist(), ids[:5])
                    self.assertEqual(sidecar["sequence"].tolist(), sequences[:5])
                    self.assertEqual(sidecar["source_row"].tolist(), [1, 2, 3, 4, 5])
                    query_ids = [f"q{number:06d}" for number in range(1, 6)]
                    self.assertEqual(sidecar["blast_query_id"].tolist(), query_ids)
                    self.assertEqual(result.loc[:4, "blast_query_id"].tolist(), query_ids)
                    parsed = parse_blast_results(report, query_ids)
                    self.assertEqual(parsed["q000005"], [])


if __name__ == "__main__":
    unittest.main()

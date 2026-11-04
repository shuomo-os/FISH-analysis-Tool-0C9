import tempfile
from pathlib import Path
import unittest
from unittest import mock

import pandas as pd

from probe_io import read_fasta_records, read_probe_table, write_probe_table


class ProbeTableTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def test_round_trip_preserves_identifiers_and_missing_sequences(self):
        expected = pd.DataFrame({
            "id": ["1", "001", "sample A", "001", "NA", "nan", ""],
            "sequence": ["ACGT", "", "AACCGG", "acgt", "", "nan", "ACGTN"],
            "valid_sequence": [True, False, True, True, False, False, False],
            "tm": [60.5, 0.0, 63.2, 60.5, 0.0, 0.0, 0.0],
        })
        for extension in (".csv", ".tsv", ".txt", ".xlsx", ".CSV"):
            with self.subTest(extension=extension):
                path = self.directory / f"probes{extension}"
                write_probe_table(expected, path)
                actual = read_probe_table(path, required_columns=("id", "sequence"))
                self.assertEqual(actual["id"].tolist(), expected["id"].tolist())
                self.assertEqual(actual["sequence"].tolist(), expected["sequence"].tolist())
                self.assertEqual(actual["valid_sequence"].tolist(), expected["valid_sequence"].tolist())
                self.assertAlmostEqual(actual["tm"].mean(), expected["tm"].mean())

    def test_csv_does_not_turn_missing_cells_into_nan_sequences(self):
        path = self.directory / "input.csv"
        path.write_text('id,sequence\n001,\nNA,NA\nnan,nan\n"sample, A",ACGT\n', encoding="utf-8")
        table = read_probe_table(path)
        self.assertEqual(table["id"].tolist(), ["001", "NA", "nan", "sample, A"])
        self.assertEqual(table["sequence"].tolist(), ["", "NA", "nan", "ACGT"])

    def test_tsv_is_selected_by_extension_even_when_csv_parser_would_succeed(self):
        path = self.directory / "input.tsv"
        path.write_text("id\tsequence\n0001\tACGT\n", encoding="utf-8")
        table = read_probe_table(path, required_columns=("sequence",))
        self.assertEqual(table.columns.tolist(), ["id", "sequence"])
        self.assertEqual(table.iloc[0]["id"], "0001")

    def test_default_allows_optional_id_and_arbitrary_reverse_complement_columns(self):
        path = self.directory / "input.csv"
        path.write_text("oligo\nACGT\n", encoding="utf-8")
        self.assertEqual(read_probe_table(path)["oligo"].tolist(), ["ACGT"])
        with self.assertRaisesRegex(ValueError, "sequence"):
            read_probe_table(path, required_columns=("sequence",))
        path.write_text("sequence\nACGT\n", encoding="utf-8")
        self.assertEqual(read_probe_table(path, required_columns="sequence")["sequence"].tolist(), ["ACGT"])

    def test_empty_file_is_rejected_but_header_only_table_is_preserved(self):
        path = self.directory / "empty.csv"
        path.write_text("", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "表头"):
            read_probe_table(path)
        path.write_text("id,sequence\n", encoding="utf-8")
        table = read_probe_table(path, required_columns=("sequence",))
        self.assertTrue(table.empty)
        self.assertEqual(table.columns.tolist(), ["id", "sequence"])

    def test_writer_creates_directory_and_unknown_extension_uses_csv(self):
        path = self.directory / "new" / "output.data"
        write_probe_table(pd.DataFrame({"id": ["001"], "sequence": ["ACGT"]}), path)
        self.assertEqual(read_probe_table(path).iloc[0]["id"], "001")

    def test_failed_save_does_not_truncate_existing_file(self):
        path = self.directory / "results.csv"
        path.write_text("previous results", encoding="utf-8")
        with mock.patch.object(pd.DataFrame, "to_csv", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                write_probe_table(pd.DataFrame({"id": ["001"]}), path)
        self.assertEqual(path.read_text(encoding="utf-8"), "previous results")
        self.assertEqual(list(self.directory.iterdir()), [path])


class FastaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.path = Path(self.temporary.name) / "targets.fasta"
        self.addCleanup(self.temporary.cleanup)

    def read(self, text):
        self.path.write_text(text, encoding="utf-8")
        return read_fasta_records(self.path)

    def test_multiple_targets_are_not_concatenated(self):
        self.assertEqual(
            self.read(">target_1 description\nATCG\n>target_2\nGCTA\n"),
            [("target_1", "ATCG"), ("target_2", "GCTA")],
        )

    def test_wrapped_sequences_bom_and_whitespace(self):
        self.assertEqual(
            self.read("\ufeff>rna\n a c g u \n\tga\tcu\n\n"),
            [("rna", "ACGUGACU")],
        )

    def test_invalid_sequence_characters_are_not_silently_removed(self):
        self.assertEqual(self.read(">target\nACGN-7?\n"), [("target", "ACGN-7?")])

    def test_duplicate_record_names_remain_separate(self):
        self.assertEqual(self.read(">same\nAT\n>same\nGC\n"), [("same", "AT"), ("same", "GC")])

    def test_malformed_or_empty_fasta_is_rejected(self):
        for text in ("", "\n", "ACGT\n", ">\nACGT\n", ">target\n", ">first\n>second\nACGT\n"):
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    self.read(text)


if __name__ == "__main__":
    unittest.main()

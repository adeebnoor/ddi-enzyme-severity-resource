#!/usr/bin/env python3
"""Meaningful parser/join checks for authorised-source reference reconstruction."""
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("rebuild_reference_membership.py")
spec = importlib.util.spec_from_file_location("reference_rebuild", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ReferenceParserTests(unittest.TestCase):
    def test_encoded_cid_offset_reversal_dedup(self):
        pairs, n = module.parse_reference("CID100001983,CID100004583\nCID100004583,CID100001983\nCID100001983,CID100004583\n", "trueDDI")
        self.assertEqual(n, 3)
        self.assertEqual(pairs, {(1983, 4583)})

    def test_encoded_cid_malformed_row_rejected(self):
        for text in ("CID100001983,missing\n", "CID100001983,CID100004583,extra\n", "1983,4583\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                module.parse_reference(text, "trueDDI")

    def test_ordinary_pubchem_cid_not_offset(self):
        with self.assertRaises(ValueError):
            module.parse_reference("CID1983,CID4583\n", "trueDDI")

    def test_drugbank_reversal_dedup_and_blank_lines(self):
        pairs, n = module.parse_reference("DB00945:DB00682\n\n DB00682 : DB00945 \n", "DrugBank")
        self.assertEqual(n, 2)
        self.assertEqual(pairs, {("DB00682", "DB00945")})

    def test_drugbank_bad_delimiter_or_accession(self):
        for text in ("DB00945,DB00682\n", "DB00945:00682\n", "DB00945:DB00682:DB00266\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                module.parse_reference(text, "DrugBank")

    def test_unmapped_drugbank_identifiers_not_testable(self):
        primary = [{"CID_A": "1983", "CID_B": "4583"}, {"CID_A": "1983", "CID_B": "9999"}]
        result = module.reference_membership(primary, {("DB00682", "DB00945")}, "DrugBank", {"1983": "DB00945", "4583": "DB00682"})
        self.assertEqual((result[0]["testable"], result[0]["reference_member"]), ("Yes", "Yes"))
        self.assertEqual((result[1]["testable"], result[1]["reference_member"]), ("No", "No"))

    def test_wrong_raw_snapshot_hash_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "wrong_snapshot.txt"
            source.write_text("CID100001983,CID100004583\n")
            result = subprocess.run([sys.executable, str(SCRIPT), "trueDDI", str(source)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Source hash differs", result.stderr)


if __name__ == "__main__":
    unittest.main()

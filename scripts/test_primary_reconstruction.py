#!/usr/bin/env python3
"""Small synthetic corner-case checks for deposited-snapshot aggregation."""
import tempfile
import unittest
from pathlib import Path
from rebuild_primary_snapshot import aggregate_records, assert_record_equivalence, verify_input


class PrimaryReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.bridge = [{"CID_pubchem": "2", "display_name": "Drug, two", "DDInterIDs": "DDInter2;DDInter22"},
                       {"CID_pubchem": "1", "display_name": "Drug one", "DDInterIDs": "DDInter1"}]
        self.rows = [{"CID_A": "2", "CID_B": "1", "enzyme_gene": "Z", "uniprot": "P2", "direction": "induction"},
                     {"CID_A": "1", "CID_B": "2", "enzyme_gene": "A (alias)", "uniprot": "P1", "direction": "transporter_inhibition"}]

    def test_reversal_sorting_unique_sets_and_bridge_values(self):
        row = aggregate_records(self.rows, self.bridge)[0]
        self.assertEqual((row["CID_A"], row["CID_B"], row["drug_B"], row["ddinter_ids_B"]),
                         ("1", "2", "Drug, two", "DDInter2;DDInter22"))
        self.assertEqual(row["enzymes"], "A (alias), Z")
        self.assertEqual(row["mechanism"], "CYP induction, transporter inhibition")

    def test_reversed_duplicate_rejected(self):
        with self.assertRaises(ValueError):
            aggregate_records(self.rows + [dict(self.rows[0], CID_A="1", CID_B="2")], self.bridge)

    def test_missing_or_ambiguous_bridge_rejected(self):
        for bridge in (self.bridge[:1], self.bridge + [self.bridge[0]],
                       [dict(self.bridge[0], display_name=""), self.bridge[1]]):
            with self.assertRaises(ValueError):
                aggregate_records(self.rows, bridge)

    def test_invalid_identifier_direction_and_self_pair_rejected(self):
        for replacement in ({"CID_A": "zero"}, {"CID_A": "0"}, {"CID_A": "1", "CID_B": "1"},
                            {"direction": "unknown"}, {"enzyme_gene": ""}):
            with self.assertRaises(ValueError):
                aggregate_records([dict(self.rows[0], **replacement)], self.bridge)

    def test_semantic_comparison_ignores_row_order_but_detects_changed_contents(self):
        extra_bridge = {"CID_pubchem": "3", "display_name": "Drug three", "DDInterIDs": "DDInter3"}
        extra_row = dict(self.rows[0], CID_A="2", CID_B="3")
        rebuilt = aggregate_records(self.rows + [extra_row], self.bridge + [extra_bridge])
        assert_record_equivalence(rebuilt, list(reversed(rebuilt)))
        with self.assertRaises(ValueError):
            assert_record_equivalence(rebuilt, [dict(rebuilt[0], enzymes="Other"), rebuilt[1]])
        with self.assertRaises(ValueError):
            assert_record_equivalence(rebuilt, rebuilt + rebuilt)

    def test_wrong_snapshot_hash_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "altered.csv"
            path.write_text("CID_A,CID_B\n1,2\n")
            with self.assertRaises(ValueError):
                verify_input(path, "attribution")


if __name__ == "__main__":
    unittest.main(verbosity=2)

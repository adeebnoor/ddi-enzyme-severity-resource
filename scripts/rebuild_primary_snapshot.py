#!/usr/bin/env python3
"""Rebuild primary records from the deposited attribution snapshot and bridge.

This downstream transformation verifies record content after canonical sorting.
It does not reconstruct selection from the historical 20,618-pair superset or
claim to recover the original primary CSV's row order. Outputs are local only.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT_SHA256 = {
    "attribution": "a43e37c6aac4074a823647f10c9fb476b022b4d34ffb7f9559a4bc31fa1ef6f6",
    "bridge": "98bdc6b180cf9491caaf19e353342c6703d0e0efddf1044732c0293cf505a8cc",
    "reference": "d879e618912303882a86e0900584b70ed1f76ae072172224ac2dfa5609686149",
}
FIELDS = ["drug_A", "drug_B", "CID_A", "CID_B", "ddinter_ids_A",
          "ddinter_ids_B", "enzymes", "mechanism"]
DIRECTIONS = {"inhibition": "CYP inhibition", "induction": "CYP induction",
              "transporter_inhibition": "transporter inhibition",
              "transporter_induction": "transporter induction"}


def read_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        result = list(reader)
    if not result or any(None in row or None in row.values() for row in result):
        raise ValueError(f"Empty or malformed CSV: {path.name}")
    return result


def pair_key(row):
    a, b = int(row["CID_A"]), int(row["CID_B"])
    if a <= 0 or b <= 0 or a == b:
        raise ValueError("Compound identifiers must be positive and distinct")
    return tuple(sorted((a, b)))


def aggregate_records(attributions, bridge):
    """Use CID keys, bridge display names and sets of preserved raw annotations."""
    lookup = {}
    for row in bridge:
        cid = int(row["CID_pubchem"])
        if cid <= 0 or cid in lookup:
            raise ValueError("Bridge compound identifiers must be positive and unique")
        lookup[cid] = row
    grouped, seen = {}, set()
    for row in attributions:
        key = pair_key(row)
        direction, label, accession = row["direction"], row["enzyme_gene"], row["uniprot"]
        if direction not in DIRECTIONS or not label.strip() or not accession.strip():
            raise ValueError("Unknown direction or empty source protein annotation")
        attribution = (key, label, accession, direction)
        if attribution in seen:
            raise ValueError("Duplicate unordered pair/protein/accession/direction row")
        seen.add(attribution)
        labels, mechanisms = grouped.setdefault(key, (set(), set()))
        labels.add(label)
        mechanisms.add(DIRECTIONS[direction])
    records = []
    for (a, b), (labels, mechanisms) in sorted(grouped.items()):
        if a not in lookup or b not in lookup:
            raise ValueError("Attribution compound has no identifier bridge row")
        ra, rb = lookup[a], lookup[b]
        if not ra["display_name"] or not rb["display_name"]:
            raise ValueError("Primary compounds require populated bridge display names")
        for row in (ra, rb):
            ids = row["DDInterIDs"].split(";")
            if not all(identifier.startswith("DDInter") for identifier in ids):
                raise ValueError("Primary compounds require DDInter identifiers")
        records.append(dict(zip(FIELDS, [ra["display_name"], rb["display_name"],
            str(a), str(b), ra["DDInterIDs"], rb["DDInterIDs"],
            ", ".join(sorted(labels)), ", ".join(sorted(mechanisms))])))
    return records


def assert_record_equivalence(rebuilt, reference):
    if any(set(row) != set(FIELDS) for row in reference):
        raise ValueError("Reference primary table has an unexpected schema")
    expected = {}
    for row in reference:
        key = pair_key(row)
        if key in expected:
            raise ValueError("Duplicate unordered reference pair")
        # Reference CID orientation is preserved; equivalence includes every value.
        expected[key] = row
    if {pair_key(row): row for row in rebuilt} != expected:
        raise ValueError("Rebuilt primary record contents differ from the reference")


def verify_input(path, kind):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != INPUT_SHA256[kind]:
        raise ValueError(f"{kind} snapshot SHA-256 differs from the deposited input")
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attribution", type=Path, default=ROOT / "data/enzyme_pair_attribution.csv")
    parser.add_argument("--bridge", type=Path, default=ROOT / "data/cid_severity_bridge.csv")
    parser.add_argument("--reference", type=Path, default=ROOT / "data/ddi_enzyme_database.csv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "data/derived")
    args = parser.parse_args()
    try:
        hashes = {kind: verify_input(getattr(args, kind), kind) for kind in INPUT_SHA256}
        attributions, bridge, reference = (read_rows(getattr(args, kind)) for kind in INPUT_SHA256)
        rebuilt = aggregate_records(attributions, bridge)
        assert_record_equivalence(rebuilt, reference)
        if len(rebuilt) != 1900 or len(attributions) != 3072:
            raise ValueError("Unexpected deposited snapshot row counts")
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Reconstruction rejected: {error}\n")
    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / "ddi_enzyme_database_rebuilt.csv"
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rebuilt)
    audit = {"input_sha256": hashes, "attribution_rows": len(attributions),
        "primary_pairs": len(rebuilt), "drug_count": len({cid for row in rebuilt for cid in pair_key(row)}),
        "bridge_rows": len(bridge),
        "bridge_joined_cids": sum(row["joined_to_ddinter"] == "yes" for row in bridge),
        "bridge_unjoined_cids": sum(row["joined_to_ddinter"] == "no" for row in bridge),
        "semantic_record_equivalence": True, "output_order": "ascending numerical unordered CID pair",
        "historical_csv_byte_equivalence_claimed": False,
        "historical_20618_superset_selection_reconstructed": False,
        "scope": "downstream transformation of deposited attribution and identifier snapshots",
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}
    (args.outdir / "primary_reconstruction_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print("Rebuilt 1,900 primary records from 3,072 attribution rows; canonical record contents match exactly.")
    print("Historical row order and selection from the 20,618-pair upstream superset are not reconstructed.")


if __name__ == "__main__":
    main()

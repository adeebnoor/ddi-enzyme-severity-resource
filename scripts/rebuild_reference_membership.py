#!/usr/bin/env python3
"""Rebuild historical reference membership from an authorised local source copy.

This is a membership audit, not a claim of independent clinical validation.
The archived author-labelled trueDDI source and DrugBank subset are not included
in the public record. Exact source hashes identify the historical input copies;
users must obtain authorised copies from the author under applicable terms.

Examples after placing this script in the repository scripts/ directory:
  python3 scripts/rebuild_reference_membership.py trueDDI /local/data-trueDDI-new.txt
  python3 scripts/rebuild_reference_membership.py DrugBank /local/HardProvenDDI.txt

Optional --grades data/derived/ddi_enzyme_database_with_severity.csv adds
per-grade aggregate counts after the reference DDInter fingerprint is verified.
Only aggregate CSV and an audit summary are written. Pair-level source content and flags are never exported.
Requires only the Python standard library.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = {
    "trueDDI": "04d70f5575199bd9095d0730c6171b58aecbfbbcac157445e14913191f486f02",
    "DrugBank": "8a36d3ccda8a54f8241f3db1cf83031972ee86090d9328901194b968ec86b057",
}
REFERENCE_GRADE_SHA = "aa59130abacc7ef80249c419bc693b9d417d91d2dfc489ff25adf17fb2fc7609"
EXPECTED = {"trueDDI": (1900, 229), "DrugBank": (1172, 73)}


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def pair(row):
    return tuple(sorted((int(row["CID_A"]), int(row["CID_B"]))))


def dump_csv(path, records):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def parse_reference(text, kind):
    reference = set()
    source_rows = 0
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        source_rows += 1
        if kind == "trueDDI":
            match = re.fullmatch(r"CID(\d+)\s*,\s*CID(\d+)\s*", line.strip())
            if not match:
                raise ValueError(f"Line {line_number}: expected two CID-prefixed identifiers.")
            encoded = tuple(map(int, match.groups()))
            if min(encoded) <= 100000000:
                raise ValueError(f"Line {line_number}: unexpected encoded CID range.")
            # The historical file uses CID100000000 + PubChem CID. This explicit
            # archived-input convention must not be applied to ordinary PubChem IDs.
            reference.add(tuple(sorted(n - 100000000 for n in encoded)))
        else:
            match = re.fullmatch(r"(DB\d+)\s*:\s*(DB\d+)\s*", line.strip())
            if not match:
                raise ValueError(f"Line {line_number}: expected two DrugBank accessions separated by colon.")
            reference.add(tuple(sorted(match.groups())))
    return reference, source_rows


def reference_membership(primary, reference, kind, cross=None):
    cross = cross or {}
    output = []
    for r in primary:
        a, b = pair(r)
        testable = kind == "trueDDI" or (str(a) in cross and str(b) in cross)
        key = (a, b) if kind == "trueDDI" else tuple(sorted((cross.get(str(a), ""), cross.get(str(b), ""))))
        output.append(dict(CID_A=a, CID_B=b, testable="Yes" if testable else "No",
                           reference_member="Yes" if testable and key in reference else "No"))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=sorted(SOURCE_SHA))
    parser.add_argument("source", type=Path)
    parser.add_argument("--primary", type=Path, default=ROOT / "data" / "ddi_enzyme_database.csv")
    parser.add_argument("--crosswalk", type=Path, default=ROOT / "data" / "cid_drugbank_crosswalk.csv")
    parser.add_argument("--grades", type=Path)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data" / "derived")
    args = parser.parse_args()
    digest = hashlib.sha256(args.source.read_bytes()).hexdigest()
    if digest != SOURCE_SHA[args.kind]:
        raise SystemExit(f"Source hash differs from the recovered {args.kind} archive. "
                         "This is not the reference input; do not treat its counts as reproduction.")
    reference, source_rows = parse_reference(args.source.read_text(encoding="utf-8-sig"), args.kind)
    primary = read(args.primary)
    if len({pair(r) for r in primary}) != len(primary):
        raise SystemExit("Primary table contains duplicate unordered pairs.")
    cross = {r["CID_pubchem"]: r["DrugBank_id"] for r in read(args.crosswalk)} if args.kind == "DrugBank" else {}
    output = reference_membership(primary, reference, args.kind, cross)
    testable = sum(r["testable"] == "Yes" for r in output)
    members = sum(r["reference_member"] == "Yes" for r in output)
    if (testable, members) != EXPECTED[args.kind]:
        raise SystemExit(f"Reference source does not reproduce expected mapped/matched counts: {testable}/{members}.")
    aggregate = [dict(grade="All", resource_pairs=len(primary), testable_pairs=testable, reference_pairs=members)]
    if args.grades:
        grades = {pair(r): r["severity"] for r in read(args.grades)}
        if set(grades) != {pair(r) for r in primary}:
            raise SystemExit("Grade table does not cover exactly the primary pair set.")
        canonical = "".join(f"{a},{b},{grades[(a, b)]}\n" for a, b in sorted(grades))
        if hashlib.sha256(canonical.encode()).hexdigest() != REFERENCE_GRADE_SHA:
            raise SystemExit("Joined grade fingerprint differs from the reference; run rebuild_severity.py.")
        for grade in ("Major", "Moderate", "Minor", "Unknown"):
            subset = [r for r in output if grades[(r["CID_A"], r["CID_B"])] == grade]
            aggregate.append(dict(grade=grade, resource_pairs=len(subset),
                                  testable_pairs=sum(r["testable"] == "Yes" for r in subset),
                                  reference_pairs=sum(r["reference_member"] == "Yes" for r in subset)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dump_csv(args.output_dir / f"{args.kind}_aggregate_local.csv", aggregate)
    metadata = dict(source_sha256=digest, source_rows=source_rows, source_unique_pairs=len(reference),
                    primary_pairs=len(primary), testable_pairs=testable, reference_pairs=members,
                    interpretation="Recovered historical reference membership; independence/clinical validation not established by this test.",
                    pair_outputs="Not exported; aggregate audit only.",
                    source_self_pairs=sum(a == b for a, b in reference))
    (args.output_dir / f"{args.kind}_reconstruction.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"{args.kind}: source {source_rows} rows / {len(reference)} unique pairs; "
          f"primary {len(primary)}; testable {testable}; matched {members}. Reference counts reproduced.")
    for r in aggregate:
        print(f"  {r['grade']}: {r['reference_pairs']} reference pairs / {r['testable_pairs']} testable / {r['resource_pairs']} resource pairs")


if __name__ == "__main__":
    main()

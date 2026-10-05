#!/usr/bin/env python3
"""Rebuild the DDInter severity column from your own DDInter 2.0 download.

Why this script exists
----------------------
DDInter 2.0 is distributed under its own CC BY-NC-SA 4.0 terms. Clinical grades
are excluded from this public record; obtain an authorised provider download
and follow the source terms when reconstructing or reusing them.
The public tables instead carry the DDInter identifiers of both drugs in every pair.
Download DDInter 2.0 yourself (https://ddinter2.scbdd.com/download/, free for
non-commercial use under its own terms) and run:

    python3 scripts/rebuild_severity.py path/to/ddinter_downloads/

The argument may be a directory or one or more CSV files. Every CSV must contain the
columns DDInterID_A, DDInterID_B and Level (other columns are ignored).

Rule (identical to the one used to build the resource)
------------------------------------------------------
For each pair, every combination of the two drugs' DDInter identifiers is looked up
as an unordered pair. If several DDInter rows match, the most severe grade is kept
under the ordering Minor < Unknown < Moderate < Major (see Methods). A pair with no
matching DDInter row is reported and written as "NotFound".

Outputs (written to data/derived/, which is not part of the public record)
--------------------------------------------------------------------------
  ddi_enzyme_database_with_severity.csv   primary table + `severity`
  enzyme_pair_severity.csv                long table + `severity`

The rebuilt column is checked against a SHA-256 fingerprint of the grades used in the
Data Descriptor, so you can confirm that your download reproduces them exactly.
Requires only the Python standard library.
"""
import csv
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = DATA / "derived"
ORDER = {"Minor": 0, "Unknown": 1, "Moderate": 2, "Major": 3}
FINGERPRINT = "aa59130abacc7ef80249c419bc693b9d417d91d2dfc489ff25adf17fb2fc7609"
EXPECTED = {"Major": 93, "Moderate": 419, "Minor": 50, "Unknown": 1338}


def rows(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def ddinter_files(args):
    files = []
    for a in args:
        p = Path(a)
        files += sorted(p.glob("*.csv")) if p.is_dir() else [p]
    if not files:
        sys.exit("No DDInter CSV files found.")
    return files


def load_ddinter(files):
    grade = {}
    for f in files:
        for r in rows(f):
            low = {k.strip().lower(): v for k, v in r.items() if k}
            try:
                a, b = low["ddinterid_a"].strip(), low["ddinterid_b"].strip()
                lvl = low["level"].strip()
            except KeyError:
                sys.exit(f"{f}: expected columns DDInterID_A, DDInterID_B, Level")
            if lvl not in ORDER:
                continue
            k = frozenset((a, b))
            if k not in grade or ORDER[lvl] > ORDER[grade[k]]:
                grade[k] = lvl
    return grade


def ids(cell):
    return [x.strip() for x in cell.split(";") if x.strip()]


def write(path, header, data):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(data)


def main(argv):
    if not argv:
        sys.exit(__doc__)
    grade = load_ddinter(ddinter_files(argv))
    primary = rows(DATA / "ddi_enzyme_database.csv")

    sev_by_cid, missing = {}, 0
    out = []
    for r in primary:
        found = [grade[frozenset((x, y))] for x in ids(r["ddinter_ids_A"])
                 for y in ids(r["ddinter_ids_B"]) if frozenset((x, y)) in grade]
        s = max(found, key=ORDER.get) if found else "NotFound"
        missing += s == "NotFound"
        sev_by_cid[tuple(sorted((int(r["CID_A"]), int(r["CID_B"]))))] = s
        out.append(list(r.values()) + [s])

    OUT.mkdir(exist_ok=True)
    write(OUT / "ddi_enzyme_database_with_severity.csv", list(primary[0]) + ["severity"], out)

    long_rows = rows(DATA / "enzyme_pair_attribution.csv")
    write(OUT / "enzyme_pair_severity.csv", list(long_rows[0]) + ["severity"],
          [list(r.values()) + [sev_by_cid[tuple(sorted((int(r["CID_A"]), int(r["CID_B"]))))]]
           for r in long_rows])

    canon = "".join(f"{a},{b},{s}\n" for (a, b), s in sorted(sev_by_cid.items()))
    counts = {g: sum(s == g for s in sev_by_cid.values()) for g in EXPECTED}
    ok = hashlib.sha256(canon.encode()).hexdigest() == FINGERPRINT
    print(f"pairs: {len(sev_by_cid)}  " + "  ".join(f"{g}: {n}" for g, n in counts.items())
          + f"  NotFound: {missing}")
    if ok:
        print("Severity column reproduces the Data Descriptor exactly (fingerprint match).")
        return 0
    print("WARNING: severity column differs from the Data Descriptor "
          f"(expected {EXPECTED}). Your DDInter download may be a different release.")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

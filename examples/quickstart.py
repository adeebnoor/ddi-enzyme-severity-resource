#!/usr/bin/env python3
"""
quickstart.py — six worked examples against the released data.

Usage:
    python3 examples/quickstart.py

Requires pandas (the resource itself needs nothing beyond the standard library;
pandas is used here only because it is what most reusers will reach for).

Each example prints a short, readable result so you can check the recipe does
what the comment says before adapting it.

Severity grades come from DDInter 2.0 (CC BY-NC-SA 4.0) and are not redistributed in
this record. Recipes 1-3 and 6 need them: first run
    python3 scripts/rebuild_severity.py <your DDInter 2.0 download>
which writes data/derived/. Recipes 4 and 5 run on the public files alone.
"""
import sys

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
DERIVED = DATA / "derived"
if not (DERIVED / "ddi_enzyme_database_with_severity.csv").exists():
    sys.exit("Run scripts/rebuild_severity.py <your DDInter 2.0 download> first "
             "(DDInter grades are not redistributed in this record).")


def banner(n, title):
    print(f"\n{'─' * 74}\n{n}. {title}\n{'─' * 74}")


# ---------------------------------------------------------------------------
banner(1, "Every Major-severity pair annotated with CYP2C9")

ddi = pd.read_csv(DERIVED / "ddi_enzyme_database_with_severity.csv")
cyp2c9_major = ddi[(ddi.severity == "Major")
                   & ddi.enzymes.str.contains("CYP2C9", na=False)]
print(f"{len(cyp2c9_major)} pairs\n")
print(cyp2c9_major[["drug_A", "drug_B", "enzymes", "mechanism"]]
      .head(8).to_string(index=False))


# ---------------------------------------------------------------------------
banner(2, "Severity occurrence profile of accession-backed genes")

pairs = pd.read_csv(DERIVED / "enzyme_pair_severity.csv")
audit = pd.read_csv(DATA / "protein_annotation_audit.csv")
resolved = audit[audit.annotation_resolution == "verified_unique_human_gene"]
pairs = pairs.merge(resolved[["source_uniprot", "official_primary_gene"]], left_on="uniprot", right_on="source_uniprot", how="left", validate="many_to_one")
profile = (pairs.dropna(subset=["official_primary_gene"]).groupby(["official_primary_gene", "severity"]).size()
           .unstack(fill_value=0)
           .reindex(columns=["Major", "Moderate", "Minor", "Unknown", "NotFound"],
                    fill_value=0))
profile["total"] = profile.sum(axis=1)
print(profile.sort_values("total", ascending=False).head(10).to_string())


# ---------------------------------------------------------------------------
banner(3, "Interactions for one drug, ranked by severity")

DRUG = "Tamoxifen"
order = {"Major": 0, "Moderate": 1, "Minor": 2, "Unknown": 3}
hits = ddi[(ddi.drug_A == DRUG) | (ddi.drug_B == DRUG)].copy()
hits["partner"] = hits.apply(
    lambda r: r.drug_B if r.drug_A == DRUG else r.drug_A, axis=1)
hits = hits.sort_values("severity", key=lambda s: s.map(order))
print(f"{len(hits)} interactions involving {DRUG}\n")
print(hits[["partner", "enzymes", "mechanism", "severity"]]
      .head(10).to_string(index=False))


# ---------------------------------------------------------------------------
banner(4, "Exact single-enzyme attribution subset")

# Source labels describe annotations, not independent external confirmation.
single = ddi[ddi.enzymes.str.split(",").map(len) == 1]
print(f"{len(single)} primary pairs have one source protein label")
print(single[["drug_A", "drug_B", "enzymes", "mechanism"]].head(8).to_string(index=False))


# ---------------------------------------------------------------------------
banner(5, "Joining to DrugBank via the crosswalk")

cross = pd.read_csv(DATA / "cid_drugbank_crosswalk.csv")
print(f"crosswalk: {len(cross)} PubChem CID -> DrugBank accession mappings\n")
cid_cols = [c for c in cross.columns if "cid" in c.lower()]
db_cols = [c for c in cross.columns if c not in cid_cols]
m = (ddi.merge(cross, left_on="CID_A", right_on=cid_cols[0], how="left")
        .merge(cross, left_on="CID_B", right_on=cid_cols[0], how="left",
               suffixes=("_A", "_B")))
both = m.dropna(subset=[f"{db_cols[0]}_A", f"{db_cols[0]}_B"])
print(f"{len(both)} of {len(ddi)} pairs have DrugBank accessions for both drugs.")
print("DrugBank interaction content itself is not redistributed (DrugBank terms);\n"
      "join these accessions to your own licensed DrugBank release.\n")
print(both[["drug_A", "drug_B", f"{db_cols[0]}_A", f"{db_cols[0]}_B"]]
      .head(6).to_string(index=False))


# ---------------------------------------------------------------------------
banner(6, "The right way to handle 'Unknown'")

n_unknown = (ddi.severity == "Unknown").sum()
print(f"{n_unknown:,} of {len(ddi):,} pairs ({100 * n_unknown / len(ddi):.0f}%) "
      f"are graded Unknown.\n")
print("Unknown means DDInter records no clinical grade — NOT that the pair is\n"
      "safe. For a supervised task, drop these rows or model them as missing;\n"
      "do not treat them as negatives.\n")
graded = ddi[ddi.severity.isin(["Major", "Moderate", "Minor"])]
print(f"Rows carrying observed grades (additional benchmark exclusions still required): {len(graded)}")
print(graded.severity.value_counts().to_string())

print("\nDone. See ../README.md for the full file inventory and the caveats "
      "that apply\nto each table, and ../figures/FIGURES.md for the figures.")

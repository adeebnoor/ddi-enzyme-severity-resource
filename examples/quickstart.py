#!/usr/bin/env python3
"""
quickstart.py — six worked examples against the released data.

Usage:
    python3 examples/quickstart.py

Requires pandas (the resource itself needs nothing beyond the standard library;
pandas is used here only because it is what most reusers will reach for).

Each example prints a short, readable result so you can check the recipe does
what the comment says before adapting it.
"""

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"


def banner(n, title):
    print(f"\n{'─' * 74}\n{n}. {title}\n{'─' * 74}")


# ---------------------------------------------------------------------------
banner(1, "Every Major-severity interaction mediated by CYP2C9")

ddi = pd.read_csv(DATA / "ddi_enzyme_database.csv")
cyp2c9_major = ddi[(ddi.severity == "Major")
                   & ddi.enzymes.str.contains("CYP2C9", na=False)]
print(f"{len(cyp2c9_major)} pairs\n")
print(cyp2c9_major[["drug_A", "drug_B", "enzymes", "mechanism"]]
      .head(8).to_string(index=False))


# ---------------------------------------------------------------------------
banner(2, "Severity profile of each enzyme (long format is easier to group)")

pairs = pd.read_csv(DATA / "enzyme_pair_severity.csv")
profile = (pairs.groupby(["enzyme_gene", "severity"]).size()
           .unstack(fill_value=0)
           .reindex(columns=["Major", "Moderate", "Minor", "Unknown"],
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
banner(4, "Only the externally corroborated pairs")

# in_trueDDI and in_trial are independent flags; requiring both is the
# strictest filter the resource supports.
strict = ddi[(ddi.in_trueDDI == "Yes") & (ddi.in_trial == "Yes")]
print(f"{len(strict)} pairs are in BOTH the trueDDI gold set and the "
      f"clinical-trial set\n")
print(strict.sort_values("trial_AEs", ascending=False)
      [["drug_A", "drug_B", "enzymes", "severity", "trial_AEs"]]
      .head(8).to_string(index=False))


# ---------------------------------------------------------------------------
banner(5, "Joining to DrugBank via the crosswalk")

cross = pd.read_csv(DATA / "cid_drugbank_crosswalk.csv")
flagged = pd.read_csv(DATA / "enzyme_pair_drugbank_flagged.csv")
print("The pair table already carries DrugBank accessions (DB_A, DB_B); use the "
      "crosswalk\nwhen you are starting from PubChem CIDs of your own.\n")
print(f"crosswalk: {len(cross)} CID → DrugBank mappings")
confirmed = flagged[flagged.drugbank_hard_proven == "Yes"]
print(f"pair rows independently confirmed in DrugBank: {len(confirmed)}\n")
print(confirmed[["drug_A", "drug_B", "DB_A", "DB_B", "enzyme_gene", "severity"]]
      .head(6).to_string(index=False))


# ---------------------------------------------------------------------------
banner(6, "The right way to handle 'Unknown'")

n_unknown = (ddi.severity == "Unknown").sum()
print(f"{n_unknown:,} of {len(ddi):,} pairs ({100 * n_unknown / len(ddi):.0f}%) "
      f"are graded Unknown.\n")
print("Unknown means DDInter records no clinical grade — NOT that the pair is\n"
      "safe. For a supervised task, drop these rows or model them as missing;\n"
      "do not treat them as negatives.\n")
graded = ddi[ddi.severity != "Unknown"]
print(f"Rows suitable for supervised use: {len(graded)}")
print(graded.severity.value_counts().to_string())

print("\nDone. See ../README.md for the full file inventory and the caveats "
      "that apply\nto each table, and ../figures/FIGURES.md for the figures.")

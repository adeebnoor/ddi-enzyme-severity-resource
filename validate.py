#!/usr/bin/env python3
"""
validate.py — re-derive every figure quoted in README.md and data/DATA_DICTIONARY.md
directly from the released CSV files.

Usage:
    python3 validate.py

Exits 0 if every documented figure matches the data, 1 otherwise.
Requires only the Python standard library.
"""

import csv
import re
import sys
from collections import Counter
from pathlib import Path

DATA = Path(__file__).parent / "data"
failures = []
checks = 0


def check(label, got, want):
    global checks
    checks += 1
    if got != want:
        failures.append(f"  {label}: documented {want!r}, computed {got!r}")


def rows(name):
    with open(DATA / name, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def summary(name):
    return {r["metric"]: r["value"] for r in rows(name)}


# ---------------------------------------------------------------- row counts
DOCUMENTED_ROWS = {
    "ddi_enzyme_database.csv": 1900,
    "enzyme_pair_severity.csv": 3072,
    "enzyme_pair_drugbank_flagged.csv": 1799,
    "enzyme_severity_stats.csv": 9,
    "enzyme_phenotype_enrichment.csv": 434,
    "drugbank_validation.csv": 73,
    "drugbank_validation_summary.csv": 8,
    "kegg_ddi_validation.csv": 1243,
    "kegg_enzyme_concordance.csv": 86,
    "kegg_validation_summary.csv": 13,
    "trial_validation.csv": 79,
    "trial_validation_summary.csv": 8,
    "liddi_coverage_comparison.csv": 7,
    "cid_severity_bridge.csv": 448,
    "cid_drugbank_crosswalk.csv": 255,
    "ddinter2_mechanism_annotations.csv": 11298,
    "mimic_coprescription_examples.csv": 13,
}

for fname, want in DOCUMENTED_ROWS.items():
    check(f"row count {fname}", len(rows(fname)), want)

# ------------------------------------------------- primary table composition
ddi = rows("ddi_enzyme_database.csv")

sev = Counter(r["severity"] for r in ddi)
check("severity Major", sev["Major"], 93)
check("severity Moderate", sev["Moderate"], 419)
check("severity Minor", sev["Minor"], 50)
check("severity Unknown", sev["Unknown"], 1338)
check("graded pairs", sum(v for k, v in sev.items() if k != "Unknown"), 562)

drugs = {d for r in ddi for d in (r["drug_A"], r["drug_B"])}
check("distinct drugs", len(drugs), 240)

# enzyme tokens: split on ; and , but keep parenthesised aliases attached
enzymes = set()
for r in ddi:
    for tok in re.split(r"[;,]", r["enzymes"]):
        tok = tok.strip()
        if tok and not tok.startswith("("):
            enzymes.add(tok)
check("distinct enzymes/transporters", len(enzymes), 35)

mech = Counter(m.strip() for r in ddi for m in r["mechanism"].split(","))
check("mechanism CYP inhibition", mech["CYP inhibition"], 1547)
check("mechanism transporter inhibition", mech["transporter inhibition"], 667)
check("mechanism CYP induction", mech["CYP induction"], 630)
check("mechanism transporter induction", mech["transporter induction"], 228)

check("in_trial = Yes", sum(r["in_trial"] == "Yes" for r in ddi), 79)
check("in_trueDDI = Yes", sum(r["in_trueDDI"] == "Yes" for r in ddi), 229)

# ------------------------------------------------------ per-enzyme statistics
stats = {r["enzyme"]: r for r in rows("enzyme_severity_stats.csv")}
DOCUMENTED_STATS = {
    "SLCO1B1 (OATP1B1)": ("12", "50.0", "5.24", "1.65", "16.64", "0.022"),
    "CYP2C9": ("94", "30.9", "2.86", "1.71", "4.76", "0.001"),
    "ABCC2 (MRP2)": ("18", "22.2", "1.44", "0.46", "4.48", "0.669"),
    "ABCB1 (P-gp)": ("149", "16.8", "1.08", "0.66", "1.76", "0.9"),
    "ABCC1": ("13", "15.4", "0.9", "0.2", "4.13", "1.0"),
    "CYP1A2": ("39", "10.3", "0.55", "0.19", "1.59", "0.559"),
    "CYP3A4": ("234", "10.7", "0.49", "0.3", "0.8", "0.017"),
    "CYP2C19": ("37", "8.1", "0.42", "0.13", "1.4", "0.316"),
    "CYP2D6": ("35", "5.7", "0.29", "0.07", "1.22", "0.22"),
}
for enz, (n, pct, orv, lo, hi, q) in DOCUMENTED_STATS.items():
    if enz not in stats:
        failures.append(f"  enzyme_severity_stats: missing row {enz!r}")
        continue
    s = stats[enz]
    for col, want in (("n_pairs", n), ("pct_major", pct), ("OR", orv),
                      ("CI_lo", lo), ("CI_hi", hi), ("q_fdr", q)):
        check(f"{enz}.{col}", float(s[col]), float(want))

fdr_sig = {e for e, s in stats.items() if float(s["q_fdr"]) < 0.05}
check("FDR-significant enzymes",
      fdr_sig, {"SLCO1B1 (OATP1B1)", "CYP2C9", "CYP3A4"})

# E_value must be populated exactly for the FDR-significant rows
check("E_value populated set",
      {e for e, s in stats.items() if s["E_value"].strip()}, fdr_sig)

# ------------------------------------------------- phenotype enrichment shape
pheno = rows("enzyme_phenotype_enrichment.csv")
check("phenotype distinct enzymes", len({r["enzyme"] for r in pheno}), 7)
check("phenotype distinct HPO terms", len({r["HPO"] for r in pheno}), 338)
check("phenotype all FDR-significant", all(float(r["q"]) < 0.05 for r in pheno), True)

# ------------------------------------------------------ validation summaries
db = summary("drugbank_validation_summary.csv")
check("DrugBank mapped pairs", db["enzyme DB-mapped pairs"], "1172")
check("DrugBank hard-proven confirmed", db["DrugBank hard-proven confirmed"], "73")
check("DrugBank confirmed Major", db["Major"], "11")
check("DrugBank confirmed Moderate", db["Moderate"], "37")
check("DrugBank confirmed Minor", db["Minor"], "2")
check("DrugBank confirmed Unknown", db["Unknown"], "23")
check("Micromedex mapped", db["Micromedex pairs mapped"], "72")
check("drugbank_validation.csv total == summary",
      len(rows("drugbank_validation.csv")),
      int(db["DrugBank hard-proven confirmed"]))
check("drugbank_validation.csv severity split",
      Counter(r["severity"] for r in rows("drugbank_validation.csv")),
      Counter({"Moderate": 37, "Unknown": 23, "Major": 11, "Minor": 2}))

# the DrugBank confirmation gradient must be reproducible from the pair table,
# not only asserted in the summary file
_pairs = {}
for r in rows("enzyme_pair_drugbank_flagged.csv"):
    k = (r["CID_A"], r["CID_B"])
    _pairs.setdefault(k, {"sev": r["severity"], "hp": False})
    if r["drugbank_hard_proven"] == "Yes":
        _pairs[k]["hp"] = True
check("unique DrugBank-mapped pairs", len(_pairs),
      int(db["enzyme DB-mapped pairs"]))
for grade, want_pct in (("Major", 25.0), ("Moderate", 15.0),
                        ("Minor", 7.7), ("Unknown", 2.7)):
    d = sum(v["sev"] == grade for v in _pairs.values())
    n = sum(v["sev"] == grade and v["hp"] for v in _pairs.values())
    check(f"DrugBank confirmation rate {grade}", round(100 * n / d, 1), want_pct)
    check(f"DrugBank confirmed count {grade}", n, int(db[grade]))
check("DrugBank confirmations sum to the summary total",
      sum(1 for v in _pairs.values() if v["hp"]),
      int(db["DrugBank hard-proven confirmed"]))

kg = summary("kegg_validation_summary.csv")
check("KEGG drugs mapped", float(kg["drugs_mapped_to_KEGG"]), 228.0)
check("KEGG DDI pairs", float(kg["KEGG_DDI_pairs_in_set"]), 3034.0)
check("KEGG testable pairs", float(kg["testable_pairs"]), 1243.0)
check("KEGG testable == kegg_ddi_validation rows",
      len(rows("kegg_ddi_validation.csv")), int(float(kg["testable_pairs"])))
check("KEGG Major confirm %", float(kg["Major_confirm_pct"]), 35.1)
check("KEGG Moderate confirm %", float(kg["Moderate_confirm_pct"]), 35.4)
check("KEGG Minor confirm %", float(kg["Minor_confirm_pct"]), 14.3)
check("KEGG Unknown confirm %", float(kg["Unknown_confirm_pct"]), 9.6)
check("KEGG graded-vs-unknown OR", float(kg["graded_vs_unknown_OR"]), 4.82)
check("KEGG enzyme concordance %", float(kg["enzyme_concordance_pct"]), 86.0)
check("KEGG concordance n == file rows",
      len(rows("kegg_enzyme_concordance.csv")),
      int(float(kg["enzyme_concordance_n"])))

conc = rows("kegg_enzyme_concordance.csv")
agree_pct = round(100 * sum(r["agree"] == "True" for r in conc) / len(conc), 1)
check("KEGG concordance % recomputed", agree_pct, 86.0)

tr = summary("trial_validation_summary.csv")
check("trial enzyme pairs total", float(tr["enzyme_pairs_total"]), 20618.0)
check("trial in clinical-trial set", float(tr["in_clinical_trial_DDI"]), 545.0)
check("trial pct", float(tr["pct_in_trials"]), 2.6)
check("trueDDI gold", float(tr["in_trueDDI_gold"]), 4233.0)
check("trueDDI pct", float(tr["pct_in_true"]), 20.5)
check("graded and in trials", float(tr["graded_and_in_trials"]), 79.0)
check("graded_and_in_trials == trial_validation.csv rows",
      len(rows("trial_validation.csv")), int(float(tr["graded_and_in_trials"])))
check("trial-validated Major", float(tr["trial_validated_Major"]), 7.0)
check("trial-validated Moderate", float(tr["trial_validated_Moderate"]), 29.0)

li = summary("liddi_coverage_comparison.csv")
check("LIDDI unique pairs", li["LIDDI_unique_pairs"], "4070")
check("GoldD3R mechanism pairs", li["GoldD3R_mechanism_pairs"], "12493")
check("shared drugs by CUI", li["shared_drugs_CUI"], "169")
check("shared interaction pairs", li["shared_interaction_pairs"], "0")
check("LIDDI drugs", li["LIDDI_drugs"], "265")
check("GoldD3R drugs", li["GoldD3R_drugs"], "1078")

# ------------------------------------------------- cross-file consistency
bridge = rows("cid_severity_bridge.csv")
check("bridge display_name unique per CID",
      len({r["CID_pubchem"] for r in bridge}), len(bridge))
# the populated display names are exactly the drugs that reach the primary table
check("bridge display_name set == primary-table drug set",
      {r["display_name"].strip() for r in bridge if r["display_name"].strip()},
      drugs)
check("bridge rows with a populated display_name",
      sum(1 for r in bridge if r["display_name"].strip()), 240)

pairs = rows("enzyme_pair_severity.csv")
check("pair table severity vocabulary",
      set(r["severity"] for r in pairs),
      {"Major", "Moderate", "Minor", "Unknown"})
check("pair table direction vocabulary",
      set(r["direction"] for r in pairs),
      {"inhibition", "induction", "transporter_inhibition", "transporter_induction"})
# the pair x enzyme direction counts must reconcile with the primary table's
# mechanism-class occurrences (bare inhibition/induction == CYP-mediated)
dirc = Counter(r["direction"] for r in pairs)
check("direction inhibition == mechanism CYP inhibition",
      dirc["inhibition"], mech["CYP inhibition"])
check("direction induction == mechanism CYP induction",
      dirc["induction"], mech["CYP induction"])
check("direction transporter_inhibition == mechanism transporter inhibition",
      dirc["transporter_inhibition"], mech["transporter inhibition"])
check("direction transporter_induction == mechanism transporter induction",
      dirc["transporter_induction"], mech["transporter induction"])
check("primary table severity vocabulary",
      set(sev), {"Major", "Moderate", "Minor", "Unknown"})

flagged = rows("enzyme_pair_drugbank_flagged.csv")
check("flagged table is a subset of pair table columns",
      {"CID_A", "drug_A", "CID_B", "drug_B", "enzyme_gene",
       "direction", "severity", "DB_A", "DB_B", "drugbank_hard_proven"}
      <= set(flagged[0]), True)

# every file non-empty and header-consistent
for f in sorted(DATA.glob("*.csv")):
    with open(f, encoding="utf-8", newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        widths = {len(r) for r in rdr if r}
    check(f"{f.name} ragged rows", widths - {len(header)}, set())

# ----------------------------------- documented data-quality facts
# These assert the caveats in README.md and DATA_DICTIONARY.md, so that a
# reader can confirm the limitations are real and still hold.

trial = rows("trial_validation.csv")
check("trial_validation severity breakdown",
      Counter(r["ddinter_severity"] for r in trial),
      Counter({"Unknown": 38, "Moderate": 29, "Major": 7, "Minor": 5}))
check("trial_validation graded subset", sum(
    1 for r in trial if r["ddinter_severity"] != "Unknown"), 41)
check("in_trial=Yes rows match trial_validation",
      sum(1 for r in ddi if r["in_trial"] == "Yes"), len(trial))

dd2 = rows("ddinter2_mechanism_annotations.csv")
check("ddinter2 layer is Major-only",
      set(r["severity"] for r in dd2), {"Major"})

kegg = rows("kegg_ddi_validation.csv")
check("kegg_flag vocabulary", set(r["kegg_flag"] for r in kegg), {"", "P"})
check("kegg_flag has no contraindication rows",
      sum(1 for r in kegg if r["kegg_flag"] == "CI"), 0)
check("kegg_flag precaution rows",
      sum(1 for r in kegg if r["kegg_flag"] == "P"), 204)

# SLCO1B1 appears under two labels with two different accessions
slco = {lab: {r["uniprot"] for r in pairs if r["enzyme_gene"] == lab}
        for lab in ("SLCO1B1", "SLCO1B1 (OATP1B1)")}
check("SLCO1B1 bare-label accession", slco["SLCO1B1"], {"Q9Y6L6"})
check("SLCO1B1 (OATP1B1) accession", slco["SLCO1B1 (OATP1B1)"], {"Q4U2R8"})
check("bare SLCO1B1 rows", sum(
    1 for r in pairs if r["enzyme_gene"] == "SLCO1B1"), 2)

# the long-format tables do NOT carry cleaned display names
pair_names = {n for r in pairs for n in (r["drug_A"], r["drug_B"])}
check("raw PubChem titles in enzyme_pair_severity",
      len(pair_names - drugs), 32)

# attributions that are a shared pharmacodynamic target rather than a
# metabolising enzyme or transporter
TARGET_ONLY = {"CYP19A1", "PTGS1", "PTGS2", "HPRT1", "XDH", "AOX1", "PGD",
               "SLC31A1"}


def enzyme_set(row):
    out = set()
    for part in row["enzymes"].split(";"):
        for sub in part.split(", "):
            sub = sub.strip()
            if sub:
                out.add(sub)
    return out


target_only = [r for r in ddi if enzyme_set(r) and enzyme_set(r) <= TARGET_ONLY]
check("pairs attributed only to a shared target", len(target_only), 26)
check("of those, graded Major",
      sum(1 for r in target_only if r["severity"] == "Major"), 2)

check("Micromedex pairs mapped", db["Micromedex pairs mapped"], "72")

# ------------------------------------------------------------------- report
print(f"validate.py — {checks} documented figures checked against the released data")
if failures:
    print(f"\nFAILED ({len(failures)}):")
    print("\n".join(failures))
    sys.exit(1)
print("All documented figures match the released data. OK")

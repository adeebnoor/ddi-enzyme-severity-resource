#!/usr/bin/env python3
"""
validate.py — re-derive every figure quoted in README.md and data/DATA_DICTIONARY.md
directly from the released files.

Usage:
    python3 validate.py

Two tiers of checks:

  * Public tier (always runs). Everything that can be re-derived from the files in
    this record, including the schema check that no specified restricted field is
    redistributed at pair level (no DDInter grade, no DrugBank or KEGG flag).

  * Severity tier (runs after `python3 scripts/rebuild_severity.py <DDInter download>`).
    DDInter 2.0 grades are CC BY-NC-SA 4.0 and are not redistributed here. Once you
    have rebuilt them from your own DDInter download into data/derived/, this tier
    re-derives every severity-dependent figure: the grade distribution, every
    per-enzyme count.

Exits 0 if every check passes, 1 otherwise. Requires only the Python standard library.
"""

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT / "data"
DERIVED = DATA / "derived"
TABLES = ROOT / "tables"
failures = []
checks = 0


def check(label, got, want):
    global checks
    checks += 1
    if got != want:
        failures.append(f"  {label}: documented {want!r}, computed {got!r}")


def rows(path):
    path = path if isinstance(path, Path) else DATA / path
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def summary(name):
    return {r["metric"]: r["value"] for r in rows(name)}


def tcsv(name):
    return rows(TABLES / name)


# ------------------------------------------------------------- statistics helpers
def gene_symbol(label):
    return re.split(r"[\s(]", label.strip(), maxsplit=1)[0].rstrip("?")


def fisher_two_sided(a, b, c, d):
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d
    lo, hi = max(0, c1 - r2), min(r1, c1)
    denom = math.comb(n, c1)
    prob = {x: math.comb(r1, x) * math.comb(r2, c1 - x) / denom for x in range(lo, hi + 1)}
    return min(1.0, sum(v for v in prob.values() if v <= prob[a] * (1 + 1e-7)))


def woolf(a, b, c, d):
    if 0 in (a, b, c, d):
        a, b, c, d = (x + 0.5 for x in (a, b, c, d))
    o = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    q, prev = [0.0] * m, 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q



# =============================================================== PUBLIC TIER
DOCUMENTED_ROWS = {
    "ddi_enzyme_database.csv": 1900,
    "enzyme_pair_attribution.csv": 3072,
    "cid_severity_bridge.csv": 448,
    "cid_drugbank_crosswalk.csv": 255,
    "enzyme_severity_stats.csv": 9,
    "enzyme_severity_stats_single_enzyme.csv": 9,
    "enzyme_phenotype_enrichment.csv": 434,
    "severity_summary.csv": 5,
    "drugbank_validation_summary.csv": 8,
    "kegg_validation_summary.csv": 13,
    "liddi_coverage_comparison.csv": 7,
}
for fname, want in DOCUMENTED_ROWS.items():
    check(f"row count {fname}", len(rows(fname)), want)
check("released CSV set", {f.name for f in DATA.glob("*.csv")}, set(DOCUMENTED_ROWS))

# ---- schema exclusion checks; not a guarantee of upstream source rights
FORBIDDEN = {"severity", "ddinter_severity", "our_severity", "level",
             "drugbank_hard_proven", "in_kegg_ddi", "kegg_flag", "kegg_class",
             "mechanism_text", "hadm_id", "in_trial", "trial_aes", "in_trueddi"}
for f in sorted(DATA.glob("*.csv")):
    with open(f, encoding="utf-8", newline="") as fh:
        header = {h.strip().lower() for h in next(csv.reader(fh))}
    check(f"{f.name}: no redistributed third-party column", header & FORBIDDEN, set())

# ---- primary table
ddi = rows("ddi_enzyme_database.csv")
drugs = {d for r in ddi for d in (r["drug_A"], r["drug_B"])}
check("distinct drugs", len(drugs), 240)
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
check("every pair carries DDInter identifiers for both drugs",
      all(r["ddinter_ids_A"].startswith("DDInter") and r["ddinter_ids_B"].startswith("DDInter")
          for r in ddi), True)

bridge = rows("cid_severity_bridge.csv")
by_cid = {r["CID_pubchem"]: r for r in bridge}
check("bridge one row per CID", len(by_cid), len(bridge))
check("bridge display_name set == primary-table drug set",
      {r["display_name"].strip() for r in bridge if r["display_name"].strip()}, drugs)
check("primary-table CIDs and DDInter IDs agree with the bridge",
      all(by_cid[r["CID_A"]]["display_name"] == r["drug_A"]
          and by_cid[r["CID_B"]]["display_name"] == r["drug_B"]
          and by_cid[r["CID_A"]]["DDInterIDs"] == r["ddinter_ids_A"]
          and by_cid[r["CID_B"]]["DDInterIDs"] == r["ddinter_ids_B"] for r in ddi), True)
check("CIDs mapping to more than one DDInter ID",
      sum(";" in r["DDInterIDs"] for r in bridge), 16)

# ---- long table
pairs = rows("enzyme_pair_attribution.csv")
dirc = Counter(r["direction"] for r in pairs)
check("direction vocabulary", set(dirc),
      {"inhibition", "induction", "transporter_inhibition", "transporter_induction"})
check("direction inhibition == mechanism CYP inhibition", dirc["inhibition"], mech["CYP inhibition"])
check("direction induction == mechanism CYP induction", dirc["induction"], mech["CYP induction"])
check("direction transporter_inhibition == mechanism transporter inhibition",
      dirc["transporter_inhibition"], mech["transporter inhibition"])
check("direction transporter_induction == mechanism transporter induction",
      dirc["transporter_induction"], mech["transporter induction"])
check("long table covers the 1,900 primary pairs",
      {frozenset((r["CID_A"], r["CID_B"])) for r in pairs},
      {frozenset((r["CID_A"], r["CID_B"])) for r in ddi})
slco = {lab: {r["uniprot"] for r in pairs if r["enzyme_gene"] == lab}
        for lab in ("SLCO1B1", "SLCO1B1 (OATP1B1)")}
check("SLCO1B1 bare-label accession", slco["SLCO1B1"], {"Q9Y6L6"})
check("SLCO1B1 (OATP1B1) accession", slco["SLCO1B1 (OATP1B1)"], {"Q4U2R8"})
check("raw PubChem titles in the long table",
      len({n for r in pairs for n in (r["drug_A"], r["drug_B"])} - drugs), 32)

TARGET_ONLY = {"CYP19A1", "PTGS1", "PTGS2", "HPRT1", "XDH", "AOX1", "PGD", "SLC31A1"}


def enzyme_set(row):
    return {s.strip() for part in row["enzymes"].split(";") for s in part.split(", ") if s.strip()}


target_only = [r for r in ddi if enzyme_set(r) and enzyme_set(r) <= TARGET_ONLY]
check("pairs attributed only to a shared target", len(target_only), 26)

# ---- aggregate severity counts
sev_sum = {r["grade"]: int(r["pairs"]) for r in rows("severity_summary.csv")}
check("severity summary", sev_sum,
      {"Major": 93, "Moderate": 419, "Minor": 50, "Unknown": 1338, "Total": 1900})
check("severity summary total == primary rows", sev_sum["Total"], len(ddi))


# ---- per-enzyme statistics: internally consistent from their own counts
def check_stats_file(fname, n_included):
    st = rows(fname)
    exact = []
    for s in st:
        a, n1 = int(s["n_major"]), int(s["n_pairs"])
        c, n2 = int(s["comparator_major"]), int(s["comparator_n"])
        e = s["gene_symbol"]
        check(f"{fname} {e} arms partition the included pairs", n1 + n2, n_included)
        o, lo, hi = woolf(a, n1 - a, c, n2 - c)
        check(f"{fname} {e} pct_major", s["pct_major"], f"{100 * a / n1:.1f}")
        check(f"{fname} {e} OR", s["OR"], f"{o:.2f}")
        check(f"{fname} {e} CI", (s["CI_lo"], s["CI_hi"]), (f"{lo:.2f}", f"{hi:.2f}"))
        p = fisher_two_sided(a, n1 - a, c, n2 - c)
        check(f"{fname} {e} p", s["p"], f"{p:.2e}")
        exact.append((s, o, p))
    for (s, o, _), q in zip(exact, bh([p for _, _, p in exact])):
        check(f"{fname} {s['gene_symbol']} q", s["q_fdr"], (f"{q:.2e}" if 0 < q < 0.0005 else f"{q:.3f}"))
    return st


stats = check_stats_file("enzyme_severity_stats.csv", 549)
check_stats_file("enzyme_severity_stats_single_enzyme.csv", 424)
check("FDR-significant enzymes",
      {s["gene_symbol"] for s in stats if float(s["q_fdr"]) < 0.05}, {"CYP2C9", "CYP3A4"})

# ---- phenotype layer
pheno = rows("enzyme_phenotype_enrichment.csv")
check("phenotype distinct enzymes", len({r["enzyme"] for r in pheno}), 7)
check("phenotype distinct HPO terms", len({r["HPO"] for r in pheno}), 338)
check("phenotype all FDR-significant", all(float(r["q"]) < 0.05 for r in pheno), True)

# ---- external-validation summaries (aggregate counts only) and their tables
db = summary("drugbank_validation_summary.csv")
for k, v in (("enzyme DB-mapped pairs", "1172"), ("DrugBank hard-proven confirmed", "73"),
             ("Major", "11"), ("Moderate", "37"), ("Minor", "2"), ("Unknown", "23"),
             ("Micromedex pairs mapped", "72"), ("Micromedex ∩ hard-proven", "35")):
    check(f"DrugBank summary {k}", db[k], v)
t6 = {r["grade"]: r for r in tcsv("SupplementaryTableS3_external_aggregate_counts.csv") if r["resource"] == "DrugBank"}
for g, pct in (("Major", "25.0"), ("Moderate", "15.0"), ("Minor", "7.7"), ("Unknown", "2.7")):
    n, k = int(t6[g]["testable_pairs"]), int(t6[g]["confirmed_pairs"])
    check(f"Table 6 confirmed {g} == summary", str(k), db[g])
    check(f"Table 6 % {g}", f"{100 * k / n:.1f}", pct)
check("Table 6 n sums to DrugBank-mapped pairs",
      sum(int(t6[g]["testable_pairs"]) for g in ("Major", "Moderate", "Minor", "Unknown")),
      int(db["enzyme DB-mapped pairs"]))

kg = summary("kegg_validation_summary.csv")
for k, v in (("drugs_mapped_to_KEGG", 228.0), ("KEGG_DDI_pairs_in_set", 3034.0),
             ("testable_pairs", 1243.0), ("Major_confirm_pct", 35.1),
             ("Moderate_confirm_pct", 35.4), ("Minor_confirm_pct", 14.3),
             ("Unknown_confirm_pct", 9.6), ("graded_vs_unknown_OR", 4.82),
             ("enzyme_concordance_pct", 86.0), ("enzyme_concordance_n", 86.0)):
    check(f"KEGG summary {k}", float(kg[k]), v)
t8 = {r["grade"]: r for r in tcsv("SupplementaryTableS3_external_aggregate_counts.csv") if r["resource"] == "KEGG"}
for g in ("Major", "Moderate", "Minor", "Unknown"):
    n, k = int(t8[g]["testable_pairs"]), int(t8[g]["confirmed_pairs"])
    check(f"Table 8 % {g} == summary", round(100 * k / n, 1), float(kg[f"{g}_confirm_pct"]))
check("Table 8 n sums to testable pairs",
      sum(int(t8[g]["testable_pairs"]) for g in ("Major", "Moderate", "Minor", "Unknown")),
      int(float(kg["testable_pairs"])))

li = summary("liddi_coverage_comparison.csv")
for k, v in (("LIDDI_unique_pairs", "4070"), ("GoldD3R_mechanism_pairs", "12493"),
             ("shared_drugs_CUI", "169"), ("shared_interaction_pairs", "0"),
             ("LIDDI_drugs", "265"), ("GoldD3R_drugs", "1078")):
    check(f"LIDDI {k}", li[k], v)

t4 = {r["enzyme"]: r for r in tcsv("SupplementaryTableS1_enzyme_statistics.csv")}
for st in stats:
    check(f"Supplementary S1 {st['enzyme']}", t4[st["enzyme"]], st)

# ---- source file
gold = [l for l in (ROOT / "sources" / "GoldD3R.txt").read_text().splitlines()
        if l.startswith("C")]
check("GoldD3R pairs", len(gold), 21897)
check("GoldD3R labelled pairs", sum("[" in l for l in gold), 12493)

for f in sorted(DATA.glob("*.csv")):
    with open(f, encoding="utf-8", newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        widths = {len(r) for r in rdr if r}
    check(f"{f.name} ragged rows", widths - {len(header)}, set())

# ---- structural release and provenance checks (not an upstream licence guarantee)
import hashlib
import json
key = lambda r: tuple(sorted((int(r["CID_A"]), int(r["CID_B"]))))
check("primary unordered pair uniqueness", len({key(r) for r in ddi}), len(ddi))
check("no self pairs", all(r["CID_A"] != r["CID_B"] for r in ddi), True)
check("long pair-gene-direction tuples unique", len({(key(r), r["enzyme_gene"], r["direction"]) for r in pairs}), len(pairs))
check("primary required columns", set(ddi[0]), {"drug_A","drug_B","CID_A","CID_B","ddinter_ids_A","ddinter_ids_B","enzymes","mechanism"})
for name in ("SupplementaryTableS1_enzyme_statistics.csv", "SupplementaryTableS2_single_enzyme_sensitivity.csv"):
    source = "enzyme_severity_stats.csv" if "S1_" in name else "enzyme_severity_stats_single_enzyme.csv"
    check(f"{name} identical to data statistics", tcsv(name), rows(source))
check("S4 identical to selected phenotype file", tcsv("SupplementaryTableS4_selected_phenotype.csv"), pheno)
source_manifest = rows(ROOT / "sources" / "sources.csv")
for entry in source_manifest:
    if entry["snapshot_path"] and entry["snapshot_status"] == "available_public":
        path = ROOT / entry["snapshot_path"]
        check(f"source manifest file exists {entry['source_id']}", path.is_file(), True)
        if path.is_file():
            check(f"source manifest hash {entry['source_id']}", hashlib.sha256(path.read_bytes()).hexdigest(), entry["sha256"])
    if entry["snapshot_status"] != "available_public":
        check(f"missing snapshot hash honest {entry['source_id']}", entry["sha256"], "unknown")
workflow = json.loads((ROOT / "workflow.yaml").read_text())
check("workflow original reconstruction status", workflow["original_source_reconstruction"], "incomplete_missing_upstream_snapshots")
for step in workflow["steps"]:
    if step.get("script"):
        path = ROOT / step["script"]
        check(f"workflow script exists {step['id']}", path.is_file(), True)
        check(f"workflow script hash {step['id']}", hashlib.sha256(path.read_bytes()).hexdigest(), step["script_sha256"])
html = (ROOT / "docs" / "index.html").read_text()
embedded = json.loads(re.search(r"const DATA=(\[.*?\]);", html).group(1))
check("browser primary record count", len(embedded), len(ddi))
check("browser embedded values match public primary table", [{k:v for k,v in r.items() if k != "severity"} for r in embedded], ddi)
check("browser has no embedded severity grades", {r["severity"] for r in embedded}, {""})
check("browser contains no withdrawn pair flags", any(x in html for x in ("in_trial", "trial_AEs", "in_trueDDI")), False)
check("browser no third-party scripts", bool(re.search(r"<script[^>]+src=", html)), False)

public_checks = checks

# ============================================================= SEVERITY TIER
severity_tier = (DERIVED / "enzyme_pair_severity.csv").exists()
if severity_tier:
    dsev = rows(DERIVED / "ddi_enzyme_database_with_severity.csv")
    sev = Counter(r["severity"] for r in dsev)
    check("rebuilt severity distribution", {g: sev[g] for g in ("Major", "Moderate", "Minor", "Unknown")},
          {g: sev_sum[g] for g in ("Major", "Moderate", "Minor", "Unknown")})
    check("target-only pairs graded Major",
          sum(1 for r in dsev if enzyme_set(r) and enzyme_set(r) <= TARGET_ONLY
              and r["severity"] == "Major"), 2)

    # every per-enzyme statistic re-derived from the rebuilt long table
    _sev, _inh = {}, {}
    for r in rows(DERIVED / "enzyme_pair_severity.csv"):
        key = tuple(sorted((r["CID_A"], r["CID_B"])))
        _sev[key] = r["severity"]
        if r["direction"] in ("inhibition", "transporter_inhibition"):
            _inh.setdefault(key, set()).add(gene_symbol(r["enzyme_gene"]))
    _incl = {k: g for k, g in _inh.items() if _sev[k] in {"Major", "Moderate", "Minor"}}
    check("graded inhibition-attributed pairs", len(_incl), 549)
    for fname, sub in (("enzyme_severity_stats.csv", _incl),
                       ("enzyme_severity_stats_single_enzyme.csv",
                        {k: g for k, g in _incl.items() if len(g) == 1})):
        released = {r["gene_symbol"]: r for r in rows(fname)}
        cnt = Counter(e for g in sub.values() for e in g)
        tested = {e for e, n in cnt.items() if n >= 10}
        check(f"{fname}: tested set (rebuilt)", set(released), tested)
        for e in tested & set(released):
            a = sum(_sev[k] == "Major" for k, g in sub.items() if e in g)
            n1 = sum(1 for g in sub.values() if e in g)
            c = sum(_sev[k] == "Major" for k, g in sub.items() if e not in g)
            s = released[e]
            check(f"{fname} {e} counts (rebuilt)",
                  (int(s["n_pairs"]), int(s["n_major"]), int(s["comparator_n"]), int(s["comparator_major"])),
                  (n1, a, len(sub) - n1, c))

# ------------------------------------------------------------------- report
print(f"validate.py — public tier: {public_checks} checks")
if severity_tier:
    print(f"validate.py — severity tier: {checks - public_checks} checks (DDInter grades rebuilt)")
else:
    print("validate.py — severity tier skipped: run scripts/rebuild_severity.py "
          "<your DDInter 2.0 download> to enable it")
if failures:
    print(f"\nFAILED ({len(failures)}):")
    print("\n".join(failures))
    sys.exit(1)
print("All checks pass. OK")

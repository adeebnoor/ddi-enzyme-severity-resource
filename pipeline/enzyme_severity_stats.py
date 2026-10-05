#!/usr/bin/env python3
"""Per-enzyme odds of Major severity, regenerated from the released long-format table.

Input : data/derived/enzyme_pair_severity.csv
        (created by scripts/rebuild_severity.py from your own DDInter 2.0 download;
        DDInter grades are not redistributed in this record)
Output: data/enzyme_severity_stats.csv               (tested enzymes, primary analysis)
        data/enzyme_severity_stats_single_enzyme.csv  (sensitivity: single-enzyme pairs only)

Analysis rule (stated in the Methods of the Data Descriptor):

  * Protein accessions are mapped to verified unique human genes using the
    frozen data/protein_annotation_audit.csv. Unresolved rows are excluded.
    This resolves protein identity, not the validity of a drug-protein relation.
  * Unit of analysis: the unique unordered drug pair.
  * Included pairs: severity in {Major, Moderate, Minor} and at least one
    inhibition-direction attribution (direction = inhibition or transporter_inhibition).
  * Tested enzymes: those with >= 10 included pairs.
  * Exposed arm for enzyme E: included pairs with an inhibition attribution to E.
    Comparator arm: all other included pairs (no inhibition attribution to E).
    The two arms are disjoint; a pair inhibiting several enzymes is in the exposed
    arm of each of them and in the comparator arm of every other enzyme.
  * Outcome: Major vs (Moderate or Minor).
  * Two-sided Fisher exact test; odds ratio with log-method (Woolf) 95% CI
    (Haldane-Anscombe +0.5 only if a cell is zero); Benjamini-Hochberg q across
    the tested enzymes.

Requires: Python 3.10+, scipy.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
INHIBITION = {"inhibition", "transporter_inhibition"}
GRADED = {"Major", "Moderate", "Minor"}
MIN_PAIRS = 10
# Common aliases are display aids; grouping uses accession-backed primary genes.
DISPLAY = {"SLC22A6": "SLC22A6 (OAT1)", "SLC22A8": "SLC22A8 (OAT3)",
           "SLCO1B1": "SLCO1B1 (OATP1B1)", "ABCB1": "ABCB1 (P-gp)",
           "ABCC2": "ABCC2 (MRP2)", "ABCC4": "ABCC4 (MRP4)", "ABCG2": "ABCG2 (BCRP)"}


def accession_genes():
    with (ROOT / "data" / "protein_annotation_audit.csv").open(encoding="utf-8", newline="") as fh:
        audit = list(csv.DictReader(fh))
    mapping = {}
    for row in audit:
        if row["annotation_resolution"] != "verified_unique_human_gene":
            continue
        if row["organism"] != "Homo sapiens (Human)" or not row["official_primary_gene"] or len(row["official_primary_gene"].split()) != 1:
            raise ValueError("Invalid verified gene-resolution record.")
        accession, gene = row["source_uniprot"], row["official_primary_gene"]
        if accession in mapping and mapping[accession] != gene:
            raise ValueError("Conflicting accession-to-gene resolution.")
        mapping[accession] = gene
    return mapping


def odds_ratio(a, b, c, d):
    if 0 in (a, b, c, d):
        a, b, c, d = (x + 0.5 for x in (a, b, c, d))
    o = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def load():
    sev, genes = {}, defaultdict(set)
    mapping = accession_genes()
    excluded_rows = 0
    src = ROOT / "data" / "derived" / "enzyme_pair_severity.csv"
    if not src.exists():
        raise SystemExit("Run scripts/rebuild_severity.py <your DDInter 2.0 download> first.")
    with open(src, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            pair = tuple(sorted((r["CID_A"], r["CID_B"])))
            sev[pair] = r["severity"]
            if r["direction"] in INHIBITION:
                gene = mapping.get(r["uniprot"])
                if gene:
                    genes[pair].add(gene)
                else:
                    excluded_rows += 1
    included = {p: g for p, g in genes.items() if sev[p] in GRADED}
    print(f"Unresolved inhibition-attribution rows excluded: {excluded_rows}")
    return sev, included


def analyse(sev, included):
    count = defaultdict(int)
    for g in included.values():
        for e in g:
            count[e] += 1
    tested = sorted(e for e, n in count.items() if n >= MIN_PAIRS)
    rows = []
    for e in tested:
        exp = [p for p, g in included.items() if e in g]
        com = [p for p, g in included.items() if e not in g]
        a = sum(sev[p] == "Major" for p in exp); b = len(exp) - a
        c = sum(sev[p] == "Major" for p in com); d = len(com) - c
        o, lo, hi = odds_ratio(a, b, c, d)
        p = fisher_exact([[a, b], [c, d]])[1]
        rows.append(dict(enzyme=DISPLAY.get(e, e), gene_symbol=e, n_pairs=len(exp),
                         n_major=a, pct_major=round(100 * a / len(exp), 1),
                         comparator_n=len(com), comparator_major=c,
                         OR=o, CI_lo=lo, CI_hi=hi, p=p))
    for r, q in zip(rows, bh([r["p"] for r in rows])):
        r["q_fdr"] = q
        r["direction_effect"] = "enriched_Major" if r["OR"] > 1 else "depleted_Major"
    return sorted(rows, key=lambda r: -r["OR"])


def write(rows, name):
    cols = ["enzyme", "gene_symbol", "n_pairs", "n_major", "pct_major", "comparator_n",
            "comparator_major", "OR", "CI_lo", "CI_hi", "p", "q_fdr",
            "direction_effect"]
    fmt = {"OR": "{:.2f}", "CI_lo": "{:.2f}", "CI_hi": "{:.2f}", "p": "{:.2e}",
           "q_fdr": "{:.3f}"}
    with open(ROOT / "data" / name, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([(f"{r[c]:.2e}" if c == "q_fdr" and 0 < r[c] < 0.0005 else fmt[c].format(r[c])) if c in fmt and r[c] is not None
                        else ("" if r[c] is None else r[c]) for c in cols])


if __name__ == "__main__":
    sev, included = load()
    primary = analyse(sev, included)
    write(primary, "enzyme_severity_stats.csv")
    single = {p: g for p, g in included.items() if len(g) == 1}
    write(analyse(sev, single), "enzyme_severity_stats_single_enzyme.csv")
    print(f"included pairs: {len(included)}; single-enzyme pairs: {len(single)}")
    for r in primary:
        print(f"{r['enzyme']:20s} n={r['n_pairs']:3d} Major={r['n_major']:3d} "
              f"OR={r['OR']:.2f} ({r['CI_lo']:.2f}-{r['CI_hi']:.2f}) "
              f"P={r['p']:.2e} q={r['q_fdr']:.3g}")

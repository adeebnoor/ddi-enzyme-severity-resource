# Figures

Every figure here is generated from the CSV files in `data/` by
[`../scripts/make_figures.py`](../scripts/make_figures.py). Nothing is drawn by
hand and no value is typed in by hand, so a figure cannot drift away from the
data it describes. Regenerate them at any time with:

```bash
python3 scripts/make_figures.py
```

Each figure is written as **PNG** (200 dpi, for the README and for print) and
**SVG** (vector, for resizing or for a manuscript). The underlying numbers are
always available as a table — the CSV the figure was built from is named below.

## Design notes

- The palette is colour-vision-deficiency safe and was verified with a
  contrast/ΔE validator rather than by eye. Blue and red are the two poles of a
  diverging scale (depleted / enriched); grey is the de-emphasis colour.
- **Grey means "not significant".** Results that do not survive FDR correction
  are deliberately recessive so they cannot be misread as trends.
- Severity uses a single-hue ordinal ramp (light Minor → dark Major). `Unknown`
  is a neutral grey because it is the *absence* of a grade, not a level of the
  scale.
- Text never wears a data colour; identity comes from the coloured mark beside
  the label, so the figures survive greyscale printing.

---

## Figure 1 — Which enzymes carry the severe interactions

![Per-enzyme odds of Major severity](fig1_enzyme_forest.png)

Odds ratio, with 95% confidence interval, that a pair attributed to a given
enzyme carries a Major severity grade, against all other pairs. Log scale;
the vertical rule is OR = 1 (no association).

Three results survive Benjamini–Hochberg correction: **SLCO1B1 (OATP1B1)** and
**CYP2C9** are enriched for Major-severity interactions, and **CYP3A4** — the
enzyme most often assumed to dominate clinical DDI risk — is *depleted*. The
other six enzymes are null results and are shown in grey.

Read the SLCO1B1 result with care: it rests on 12 pairs, and the confidence
interval runs from 1.65 to 16.64.

Source: `data/enzyme_severity_stats.csv`

---

## Figure 2 — Assigned severity predicts independent confirmation

![Confirmation rate by severity grade in DrugBank and KEGG](fig2_severity_confirmation.png)

The proportion of pairs at each severity grade that are independently confirmed
in two external gold standards, neither of which was used to build the resource.

The DrugBank gradient is monotone across grades (25.0% → 15.0% → 7.7% → 2.7%,
Cochran–Armitage trend *p* = 1.8 × 10⁻¹⁸). In KEGG, Major and Moderate are
effectively indistinguishable (35.1% vs 35.4%), but graded pairs still confirm
far more often than ungraded ones (graded-vs-`Unknown` OR 4.82).

This is the resource's central external check: severity labels inherited
through harmonisation carry information that an independent curator agrees with.

Sources: `data/enzyme_pair_drugbank_flagged.csv` (rates recomputed from the pair
table, not copied from a summary), `data/kegg_validation_summary.csv`

---

## Figure 3 — What the 1,900 enzyme-annotated pairs are graded as

![Severity composition of the primary table](fig3_severity_composition.png)

Composition of the primary table by severity grade. **70% of pairs are
`Unknown`** — DDInter records no clinical grade for them. That is a statement
about coverage, not about safety, and treating `Unknown` as a negative class
will bias any downstream model.

Source: `data/ddi_enzyme_database.csv`

---

## Figure 4 — How the interactions are mediated

![Mechanism class occurrences](fig4_mechanism_classes.png)

Occurrences of each mechanism class across the primary table. Inhibition
dominates induction, and CYP-mediated mechanisms dominate transporter-mediated
ones. A pair may carry more than one class, so occurrences exceed the 1,900
pairs.

These four counts reconcile exactly with the `direction` column of
`data/enzyme_pair_severity.csv` — `validate.py` asserts it.

Source: `data/ddi_enzyme_database.csv`

---

## Figure 5 — Strongest enzyme–adverse-event signatures

![Top enzyme-phenotype enrichments](fig5_phenotype_enrichment.png)

The 12 strongest of 434 FDR-significant enzyme × adverse-event (HPO term)
associations, by log₂ fold enrichment over the background model.

**These are descriptive signatures mediated by drug class, not causal
enzyme-attributable risk.** An enzyme that handles a particular drug class will
inherit that class's adverse-event profile; the association says nothing about
the enzyme causing the event. The figure is included because the signatures are
useful for hypothesis generation and because hiding them would be worse than
labelling them clearly.

Source: `data/enzyme_phenotype_enrichment.csv`

---

## Figure 6 — The interactive browser

![Screenshot of the interactive browser](fig6_browser_screenshot.png)

`docs/index.html` opened in a browser. All 1,900 records are searchable and
filterable by drug, enzyme/transporter, severity grade and trial-observation
status, with CSV export. The file is entirely self-contained: no external
scripts, no stylesheets, no network calls, works offline.

Source: `docs/index.html`

---

## Graphical abstract

![Graphical abstract](graphical_abstract.png)

Summary of the harmonisation pipeline and headline findings, as submitted with
the manuscript.

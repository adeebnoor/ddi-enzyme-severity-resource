# Tables

Machine-readable versions of Tables 1–8 of the Data Descriptor, as UTF-8 CSV, plus a
formatted workbook containing all eight.

| File | Content |
|---|---|
| `Table1_harmonisation_flow.csv` | Flow of the GoldD3R–DDInter harmonisation and definition of the analysis set |
| `Table2_mechanism_severity_contingency.csv` | Primary-mechanism by severity contingency table (N = 963) |
| `Table3_hypothesis_odds_ratios.csv` | Hypothesis-specific odds ratios (N = 963), with BH q-values and E-values |
| `Table4_enzyme_odds_ratios.csv` | Enzyme-resolved odds of Major severity among inhibition pairs |
| `Table5_phenotype_enrichments.csv` | Representative adverse-event phenotype enrichments per enzyme |
| `Table6_drugbank_validation_by_severity.csv` | External validation against the DrugBank hard-proven set |
| `Table7_trueddi_validation_by_severity.csv` | External validation against the trueDDI reference set |
| `Table8_kegg_validation_by_severity.csv` | External validation against the curated KEGG DDI set |
| `SupplementaryTables.xlsx` | All eight tables, one per sheet, with captions |

## Which of these are recomputable from this deposit

**Tables 4, 6, 7 and 8 are fully recomputable** from the CSV files in `data/`:

- Table 4 — `enzyme_severity_stats.csv`, and the `n pairs` and `% Major` columns reproduce
  exactly from `enzyme_pair_severity.csv` (unique non-`Unknown` pairs, inhibition direction).
- Table 6 — recomputable from `enzyme_pair_drugbank_flagged.csv`; `validate.py` asserts the
  25.0 / 15.0 / 7.7 / 2.7% gradient directly from that file.
- Table 7 — recomputable from the `in_trueDDI` column of `ddi_enzyme_database.csv`.
- Table 8 — recomputable from `kegg_ddi_validation.csv`.

**Tables 1, 2, 3 and 5 are reported values.** They derive from the GoldD3R↔DDInter overlap
(1,604 matched pairs; 1,064 with a mechanism label; N = 963) and from the full
enzyme–phenotype observation matrix. Those intermediate tables are upstream of this deposit
and are not redistributed here, so these four cannot be recomputed from `data/` alone. They
are included in machine-readable form so that the published values can at least be read,
cited and cross-checked without transcribing them from the PDF.

## Provenance note

These files were regenerated from the author's submission package. In that package all six
original table files carried the name of the **preceding** table — `Table1_harmonisation_flow.csv`
contained the DrugBank validation, `Table2_mechanism_severity_contingency.csv` contained the
harmonisation flow, and so on through a one-position rotation, in both the CSVs and the six
sheets of the workbook. The filenames and sheet names here have been corrected to match their
contents; no value was altered. Tables 7 and 8 are new, computed from the released data.

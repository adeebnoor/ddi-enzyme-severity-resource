# Final figures

Regenerate with `python3 scripts/make_figures.py` (Matplotlib). PNG, SVG and PDF files are produced; no original predictive-model ROC plot is retained.

| Stem | Manuscript use | Exact basis |
| --- | --- | --- |
| fig1_provenance | Figure 1 | Diagram explicitly separates GoldD3R context from unresolved enzyme inputs and local grade reconstruction |
| fig2_severity_composition | Figure 2 | data/severity_summary.csv; Unknown is missing grading, not safety |
| fig3_enzyme_forest | Figure 3 | data/enzyme_severity_stats.csv; accession-backed human genes; three of nine q values <0.05 (CYP2C9, CYP3A4, SLC22A6/OAT1); true SLCO1B1 below eligibility threshold |
| fig4_external_concordance | Figure 4 | SupplementaryTableS3_external_aggregate_counts.csv; three resource-specific denominators; descriptive recovered DrugBank/trueDDI and reported KEGG memberships |
| supp_fig1_mechanism_classes | Supplementary Figure 1 | Mechanism occurrences in primary CSV; multiple classes per pair |

The earlier graphical abstract and browser screenshot have been removed because they carried retired claims or schemas. The provenance figure represents missing inputs with dashed outlines; it does not imply a verified GoldD3R-to-enzyme derivation. Figure regeneration reproduces deposited aggregate summaries; exact recovered DrugBank/trueDDI raw-source membership counts can also be rebuilt locally. The selected phenotype figure is withdrawn; no HPO layer remains.

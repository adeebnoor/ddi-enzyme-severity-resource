# Final figures

Regenerate with `python3 scripts/make_figures.py` (Matplotlib). PNG, SVG and PDF files are produced; no original predictive-model ROC plot is retained.

| Stem | Manuscript use | Exact basis |
| --- | --- | --- |
| fig1_provenance | Figure 1 | Diagram explicitly separates GoldD3R context from unresolved enzyme inputs and local grade reconstruction |
| fig2_severity_composition | Figure 2 | data/severity_summary.csv; Unknown is missing grading, not safety |
| fig3_enzyme_forest | Figure 3 | data/enzyme_severity_stats.csv; two of nine q values <0.05; OATP1B1 q=0.055 |
| fig4_external_concordance | Figure 4 | SupplementaryTableS3_external_aggregate_counts.csv; resource-specific denominators; historical descriptive counts |
| supp_fig1_phenotype_selected | Supplementary Figure 1 | Twelve highest log2FE among 434 selected associations; full correction universe absent |
| supp_fig2_mechanism_classes | Supplementary Figure 2 | Mechanism occurrences in primary CSV; multiple classes per pair |

The earlier graphical abstract and browser screenshot have been removed because they carried retired claims or schemas. The provenance figure represents missing inputs with dashed outlines; it does not imply a verified GoldD3R-to-enzyme derivation. Figure regeneration reproduces the deposited summaries, not unavailable original source calculations.

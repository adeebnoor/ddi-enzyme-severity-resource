# An enzyme-resolved dataset linking drug interaction mechanisms to clinical severity grades

Version **2.1.0**, review snapshot, 5 October 2026. Repository: <https://github.com/adeebnoor/ddi-enzyme-severity-resource>.

The released author-derived export contains **1,900 unique unordered drug pairs**, **240 drugs**, **35 source protein labels**, and **3,072 pair–protein–direction rows**. DDInter 2.0 grades are rebuilt locally from a user's own download; public files contain identifiers and aggregate statistics. This is a research resource, not a prescribing tool.

**Reproducibility boundary:** the public files, their structural consistency, local identifier-based severity reconstruction, and the exploratory per-enzyme calculations are verifiable. The original enzyme-detail input files and their transformations are absent. The relationship of that branch to the separately deposited GoldD3R export is not established. This release therefore does **not** claim complete reconstruction from original source resources. Exact available file hashes and missing metadata are recorded in [sources/sources.csv](sources/sources.csv) and [workflow.yaml](workflow.yaml).

**Archive status:** DOI **10.5281/zenodo.23169653** is reserved for an **unpublished Zenodo draft**. It is not a published, resolving dataset citation. The current accessible record is this repository; citation metadata is in [CITATION.cff](CITATION.cff).

![Provenance and reproducibility boundaries](figures/fig1_provenance.png)

## Release inventory

All CSV files are UTF-8 with one header row. Column definitions and caveats are in [data/DATA_DICTIONARY.md](data/DATA_DICTIONARY.md).

| File | Rows | Role |
| --- | ---: | --- |
| data/ddi_enzyme_database.csv | 1,900 | Primary pairs, identifiers, source protein labels and mechanism classes |
| data/enzyme_pair_attribution.csv | 3,072 | Long-format pair × protein label × direction |
| data/cid_severity_bridge.csv | 448 | Candidate PubChem–DDInter identifier bridge; 240 display names populated |
| data/cid_drugbank_crosswalk.csv | 255 | PubChem–DrugBank identifier-only bridge |
| data/severity_summary.csv | 5 | Aggregate expected grade counts, including total |
| data/enzyme_severity_stats.csv | 9 | Exploratory per-enzyme statistics |
| data/enzyme_severity_stats_single_enzyme.csv | 9 | Single-enzyme sensitivity statistics |
| data/enzyme_phenotype_enrichment.csv | 434 | Previously selected phenotype associations; full tested universe absent |
| data/drugbank_validation_summary.csv | 8 | Historical reported aggregate comparison |
| data/kegg_validation_summary.csv | 13 | Historical reported aggregate comparison |
| data/liddi_coverage_comparison.csv | 7 | Historical reported coverage comparison |

[Tables 1–2 and supplementary tables S1–S4](tables/README.md) provide the release inventory, source boundaries, exploratory statistics, reported external counts and selected phenotype layer. The workbook contains the same numeric records. Prior hypothesis tables are retired. Prior clinical-trial and trueDDI flags are withdrawn: the supplied candidate sources could not establish their provenance or reproduce all flags.

## Local use and validation

```bash
python3 validate.py
sha256sum -c CHECKSUMS.sha256

# Acquire your own DDInter 2.0 CSVs under its source terms:
python3 scripts/rebuild_severity.py path/to/ddinter_downloads/
python3 validate.py
python3 pipeline/enzyme_severity_stats.py
python3 examples/quickstart.py
```

`validate.py` and `rebuild_severity.py` use the Python standard library. Statistics require SciPy; figures require Matplotlib; the workbook requires openpyxl; examples require pandas. The environment used for this review is recorded in `requirements-review.txt`. Locally rebuilt grades go to ignored `data/derived/` and are excluded from public releases.

The join treats DDInter drug identifiers as unordered pairs, considers all combinations of multiple identifiers, and keeps the maximum matched grade under **Minor < Unknown < Moderate < Major**. Missing source pairs become `NotFound`. The expected joined output is **93 Major, 419 Moderate, 50 Minor and 1,338 Unknown**. Its SHA-256 fingerprint is `aa59130abacc7ef80249c419bc693b9d417d91d2dfc489ff25adf17fb2fc7609`. A match verifies this output column; it does not establish the original full DDInter download version or reproduce missing enzyme-detail inputs.

Open `docs/index.html` locally to search drugs, filter exact source protein labels, sort, and export filtered CSVs. The HTML is self-contained. You can load DDInter CSVs in the browser; no files or grades are uploaded. The browser checks the same joined-grade fingerprint where the browser cryptography API is available. `scripts/make_browser.py` regenerates its embedded records from the primary public CSV.

## Descriptive coverage and exploratory reuse

![Severity coverage](figures/fig2_severity_composition.png)

`Unknown` is missing clinical grading in the matched source, not evidence of absence or safety. `NotFound` is distinct: the local download has no matching DDInter row. Neither belongs in a negative supervised-learning class.

![Exploratory per-enzyme comparison](figures/fig3_enzyme_forest.png)

The calculation uses unique unordered pairs with a Major, Moderate or Minor grade and at least one inhibition-direction attribution (549 pairs). For each gene symbol, its exposed arm contains pairs with inhibition attribution to that symbol; its disjoint comparator contains all other included pairs. Multi-enzyme pairs enter each relevant enzyme's exposed arm. Nine symbols have at least 10 included pairs. Two-sided Fisher tests, Woolf 95% intervals with +0.5 only for zero cells, and Benjamini–Hochberg adjustment across nine tests are implemented in `pipeline/enzyme_severity_stats.py`. The sensitivity set contains 424 pairs with only one inhibition-attributed symbol.

CYP2C9 (OR 2.78, q=0.002) and CYP3A4 (OR 0.44, q=0.005) reach the predefined FDR threshold in this overlap; **SLCO1B1/OATP1B1 does not** (OR 3.92, q=0.055). The other seven rows do not meet the threshold. These comparisons illustrate reuse and are not causal enzyme-specific clinical risk estimates. E-values have been removed.

![Reported aggregate external concordance](figures/fig4_external_concordance.png)

DrugBank and KEGG aggregate counts are retained as **historical descriptive comparisons**. Their original snapshots, exact acquisition dates and extraction workflows are unavailable. Agreement and arithmetic can be checked against the deposited aggregate exports; source-level independent reanalysis cannot. The figure shows a distinct denominator for each resource and grade. In KEGG, Major (35.1%) and Moderate (35.4%) have similar reported rates. No claim of fully independent new validation is made.

The 434 phenotype associations are **selected-only exploratory output**, not technical validation. They cover seven source protein labels and 338 HPO identifiers. Original ontology snapshot, full observation matrix, non-significant tests and the correction universe are unavailable. The published q values can be read but not independently recomputed from this selected subset. The layer is unsuitable for estimating causal enzyme-attributable adverse-event risk.

## Identifier and attribution caveats

- Join long-format rows on CIDs, not drug names: 32 long-table names differ from cleaned display names.
- Sixteen candidate CIDs carry multiple DDInter identifiers. Preserve the explicit max-grade join rule.
- SLCO1B1 appears as a bare symbol and with the OATP1B1 alias under two accessions (`Q9Y6L6` and `Q4U2R8`). These raw values are retained; the statistics merge by bare source symbol. Gene-label stripping is a pragmatic grouping rule, not a validated HGNC remapping.
- `SLCO1B3? (OATP)` retains a source uncertainty marker. Normalising its label does not resolve that uncertainty.
- Some labels denote shared pharmacodynamic targets: 26 primary pairs contain only labels in the documented target set. Prodrug–metabolite pairs and combination-product or brand names also occur. Benchmark users should define exclusions explicitly.

## Rights and source terms

MIT applies to author-written code. CC BY 4.0 applies to the author's original documentation, figures and contribution to derived tables **to the extent the author owns those rights**. It does not override third-party terms or certify missing upstream permissions. The exact rights of the unresolved original enzyme-detail snapshot are not established by this package. Consult [LICENSE](LICENSE) and the source manifest before redistribution or commercial reuse.

DDInter clinical grades are excluded from the public CSVs and browser; use your own authorised download under DDInter's terms. Pair-level DrugBank, KEGG and MIMIC content are excluded. Private historical exports, raw trueDDI material and local derived grades must not be included in GitHub or Zenodo release files.

Questions and corrections: [repository issues](https://github.com/adeebnoor/ddi-enzyme-severity-resource/issues).

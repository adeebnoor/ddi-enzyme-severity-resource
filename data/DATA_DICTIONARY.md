# Data dictionary

Enzyme-resolved DDI research resource, review snapshot 2.2.0. Author code is MIT; author contributions are licensed only within rights the author owns. Missing upstream source permissions remain unresolved.

**Third-party content.** DDInter 2.0 severity grades (CC BY-NC-SA 4.0) are not
redistributed: the public tables carry DDInter identifiers, and
`scripts/rebuild_severity.py` adds the grades from the user's own DDInter
download (written to `data/derived/`, outside the public record). Pair-level
DrugBank and KEGG content and MIMIC-IV content are not redistributed; only
aggregate counts are.

Every file is UTF-8 CSV with a single header row and comma separators; text
fields containing commas are double-quoted. Row counts exclude the header and
were verified against the released files. **Column names below are the literal
header strings** — they are reproduced exactly as they appear in each file.

Primary-table drug names are cleaned display names; the long table retains source titles. For machine linkage use
PubChem CIDs (`CID_A`, `CID_B`, `CID_pubchem`) or DDInter identifiers rather than name matching;
`cid_severity_bridge.csv` is the authoritative name source.

Missing values are empty fields (not `NA`, not `NaN`). Per-enzyme numeric display fields are rounded; the released counts and deterministic script permit recalculation. E-values were removed.

---

## Primary table

### `ddi_enzyme_database.csv` — 1,900 rows

The primary resource: one row per drug pair.

| Column | Type | Description |
| --- | --- | --- |
| `drug_A` | string | First interacting drug (cleaned display name) |
| `drug_B` | string | Second interacting drug (cleaned display name) |
| `CID_A` | integer | PubChem Compound ID of drug A |
| `CID_B` | integer | PubChem Compound ID of drug B |
| `ddinter_ids_A` | string | DDInter 2.0 identifier(s) of drug A; semicolon-separated where the compound maps to several (16 CIDs do) |
| `ddinter_ids_B` | string | DDInter 2.0 identifier(s) of drug B |
| `enzymes` | string | Source-attributed protein label(s); may include shared targets rather than a proven pharmacokinetic mechanism, comma-separated gene symbols; transporter symbols carry a parenthesised common alias, e.g. `ABCB1 (P-gp)` |
| `mechanism` | string | Mechanism class(es), comma-separated. Historical controlled values: `CYP inhibition`, `CYP induction`, `transporter inhibition`, `transporter induction`. The CYP-labelled classes include general Enzyme and shared-target relations under the recovered D3 rules; they are not exclusively CYP pharmacokinetic mechanisms |

Composition: 240 distinct drugs; 35 distinct source protein labels and 35 accession-backed current human genes;
mechanism-class occurrences CYP inhibition 1,547, transporter inhibition 667,
CYP induction 630, transporter induction 228. Clinical-trial and reference pair flags are excluded. Correct recovered trueDDI and DrugBank raw snapshots reproduce historical memberships; aggregate-only results are in reference_membership_aggregate.csv.

**Severity.** `scripts/rebuild_severity.py` writes
`data/derived/ddi_enzyme_database_with_severity.csv`, this table plus `severity`:
the DDInter grade `Major`, `Moderate`, `Minor` or `Unknown` (`NotFound` if the
user's download lacks the pair). Each pair takes the most severe grade over all
combinations of its two drugs' DDInter identifiers, under the ordering
Minor < Unknown < Moderate < Major. With DDInter 2.0 as used in the Data
Descriptor the result is 93 `Major`, 419 `Moderate`, 50 `Minor`, 1,338 `Unknown`
(aggregate counts in `severity_summary.csv`), and the script confirms an exact
match with a SHA-256 fingerprint. `Unknown` means DDInter records no grade —
**not** that the interaction is absent or safe.

---

## Pair × enzyme tables

### `enzyme_pair_attribution.csv` — 3,072 rows

Long format: one row per pair × enzyme × direction.

| Column | Type | Description |
| --- | --- | --- |
| `CID_A` | integer | PubChem Compound ID of drug A |
| `drug_A` | string | Drug A name. **Not the cleaned display name** — 32 of the 240 names in this file are raw PubChem systematic titles (e.g. `Azane;cyclobutane-1,1-dicarboxylic acid;platinum` for carboplatin, `(E,Z)-Tamoxifen`, `Methotrexate, (A+-)-`). Join on `CID_A`/`CID_B`, never on names |
| `CID_B` | integer | PubChem Compound ID of drug B |
| `drug_B` | string | Drug B name — same caveat as `drug_A` |
| `enzyme_gene` | string | Single source protein label; not a validated causal attribution or independently verified HGNC mapping |
| `uniprot` | string | Preserved source UniProt accession. Current official audit resolves Q4U2R8 to SLC22A6, Q8TCC7 to SLC22A8, and Q9Y6L6 to SLCO1B1; two raw gene labels conflict with their accessions. Group using protein_annotation_audit.csv, not raw label stripping |
| `direction` | enum | Mechanism class for this pair × enzyme row: `inhibition`, `induction`, `transporter_inhibition`, `transporter_induction`. The bare terms denote historical general Enzyme-rule attributions, including metabolism, substrate and has_target relations; the `transporter_`-prefixed terms denote historical transporter-rule attributions |

Occurrences of `direction` are 1,547 `inhibition`, 667
`transporter_inhibition`, 630 `induction` and 228 `transporter_induction`, which
reconcile exactly with the mechanism-class occurrences in
`ddi_enzyme_database.csv`. `rebuild_severity.py` writes the same table with a
`severity` column to `data/derived/enzyme_pair_severity.csv`; that file is the
input to `pipeline/enzyme_severity_stats.py`.

---

## Exploratory usage layer: statistics

### `enzyme_severity_stats.csv` — 9 rows

Per-gene association with `Major` severity, using verified accession-backed human identities for nine genes with at least **10 Major/Moderate/Minor inhibition-direction pairs**. The frozen official audit resolves all 35 accessions; unresolved annotations would be excluded (zero here). True SLCO1B1 has two graded pairs and is ineligible; corrected SLC22A6/OAT1 has 12 graded pairs, including six Major.

| Column | Type | Description |
| --- | --- | --- |
| `enzyme` | string | Enzyme / transporter gene symbol |
| `n_pairs` | integer | Pairs attributed to this enzyme |
| `n_major` | integer | Of those, pairs graded `Major` |
| `pct_major` | float | `n_major` / `n_pairs` × 100 |
| `gene_symbol` | string | Normalised gene symbol used for aggregation |
| `comparator_n` | integer | Graded inhibition-attributed pairs with no inhibition attribution to this enzyme (comparator arm) |
| `comparator_major` | integer | Of those, graded `Major` |
| `OR` | float | Odds ratio for `Major` severity, this enzyme's pairs vs the comparator arm |
| `CI_lo` | float | Lower bound, 95% confidence interval |
| `CI_hi` | float | Upper bound, 95% confidence interval |
| `p` | float | Uncorrected two-sided *p*-value |
| `q_fdr` | float | Benjamini–Hochberg FDR-adjusted *q*-value |
| `direction_effect` | enum | `enriched_Major` or `depleted_Major` |

`CYP2C9` (OR 2.78, q=0.002), `CYP3A4` (OR 0.44, q=0.005) and `SLC22A6 (OAT1)` (OR 5.24, 95% CI 1.65–16.64, q=0.022) reach FDR<0.05 in this exploratory snapshot. The other six rows do not. Grouping uses verified accession-backed human genes; exposed and comparator arms are disjoint. True SLCO1B1 has two graded pairs and fails the minimum-10 threshold. The prior merged raw SLCO-labelled biological estimate is withdrawn because it conflated SLCO1B1 and SLC22A6. A non-significant row's `direction_effect` describes its point estimate only, not a finding. Protein identity resolution does not validate the underlying drug-pair relation.

## Descriptive historical reference comparisons

### `drugbank_validation_summary.csv` — 8 rows (`metric`, `value`)

Exact author-held input snapshot recovered (16,316 unique pairs, dated correspondence 16 January 2017), reproducing 1,172 enzyme pairs
DrugBank-mapped; 73 hard-proven confirmed (11 `Major`, 37 `Moderate`, 2
`Minor`, 23 `Unknown`); 72 Micromedex pairs mapped, 35 in the intersection with
the hard-proven set. Confirmation rate by grade: `Major` 25.0%, `Moderate`
15.0%, `Minor` 7.7%, `Unknown` 2.7%. The Micromedex counts are historical reported aggregates; their licensed original input is not reconstructed.

### `kegg_validation_summary.csv` — 13 rows (`metric`, `value`)

Reported counts from an unavailable original source snapshot: 228 drugs mapped to KEGG; 3,034 KEGG DDI pairs; 1,243 testable. Confirmation by
grade: `Major` 35.1%, `Moderate` 35.4%, `Minor` 14.3%, `Unknown` 9.6%;
Cochran–Armitage *z* = 10.49, *p* = 9.6 × 10⁻²⁶; graded-vs-`Unknown` OR 4.82
(*p* = 8.0 × 10⁻²³); enzyme concordance 86.0% (*n* = 86).

### `liddi_coverage_comparison.csv` — 7 rows (`metric`, `value`)

Coverage overlap with the independent LIDDI corpus: 4,070 LIDDI unique pairs vs
12,493 GoldD3R mechanism pairs; 265 LIDDI drugs vs 1,078 GoldD3R drugs; 169
drugs shared by CUI; **0 shared interaction pairs**, so all 4,070 LIDDI pairs
are novel relative to GoldD3R. The original matching intermediate is unavailable; this overlap statement cannot be independently source-recomputed from the deposited files.

---

## Identifier bridges and linkage layers

### `cid_severity_bridge.csv` — 448 rows

| Column | Type | Description |
| --- | --- | --- |
| `CID_pubchem` | integer | PubChem Compound ID |
| `drug_name` | string | PubChem compound title (verbatim) |
| `DDInterIDs` | string | Matched DDInter drug ID(s) |
| `joined_to_ddinter` | enum | `yes` / `no` |
| `display_name` | string | Cleaned display name — **authoritative** cleaned name used in the primary table; long-table names may differ |

### `cid_drugbank_crosswalk.csv` — 255 rows

| Column | Type | Description |
| --- | --- | --- |
| `CID_pubchem` | integer | PubChem Compound ID |
| `DrugBank_id` | string | DrugBank accession, e.g. `DB00339` |

### `severity_summary.csv` — 5 rows (`grade`, `pairs`)

Aggregate number of pairs per DDInter grade (Major 93, Moderate 419, Minor 50,
Unknown 1,338, Total 1,900). Counts only; the pair-level grades are rebuilt
locally.

---

## Known naming variation

Raw source labels are preserved and must not be grouped by stripping aliases alone:

- `SLCO1B1 (OATP1B1)` uses Q4U2R8, officially SLC22A6/OAT1. Bare `SLCO1B1` uses Q9Y6L6 and resolves to SLCO1B1; they are different proteins.
- `SLCO1B3? (OATP)` uses Q8TCC7, officially SLC22A8/OAT3. The uncertainty marker does not repair the misidentification.
- The frozen official crosswalk resolves all 35 accessions to current reviewed human genes; 33 raw symbols agree and two conflict. This resolves identifier identity only, leaving historical pair-level attributions subject to the source boundary.

---

## Provenance and reproducibility status

The exact available snapshots, SHA-256 values and explicit missing original
metadata are in `sources/sources.csv`. Runnable release processing steps and exact
script hashes are in `workflow.yaml`. Database publications identify a resource,
not the exact historical input version. The deposited GoldD3R export is contextual;
no documented transformation connects it to the primary enzyme-detail branch.
Original D3 Java/query logic is recovered: enzyme match combines an inhibits/induces relation with metabolism, substrate or has_target relation, protein type Enzyme, querying both drug orders. Transporter match combines transported-by with inhibits/induces and protein type Transporter. Two complete historical candidate graphs are recovered and assessed. The first matches 2256/3072 deposited CID–UniProt–direction keys, with 816 missing and 1675 extra within deposited pairs; the second yields zero core keys. Neither reconstructs exact current inputs or the original superset selection.

The correct author-held trueDDI reference and DrugBank subset are recovered with exact hashes and dated archival correspondence. Their membership aggregates reproduce, but provider releases, actual export creation dates and clinical independence are not established by this audit. KEGG original snapshot remains unavailable. The incomplete phenotype layer is withdrawn entirely.

### `reference_membership_aggregate.csv` — 10 rows

Columns: `reference` (trueDDI or DrugBank), `grade` (All, Major, Moderate, Minor, Unknown), `resource_pairs`, `testable_pairs`, `reference_pairs` (integer counts), and `interpretation` (descriptive membership scope). trueDDI All = 229/1,900; DrugBank All = 73/1,172 testable out of 1,900 resource pairs. Each resource's per-grade denominators sum to its All row. Exact authorized source files can regenerate aggregates with scripts/rebuild_reference_membership.py; pair flags and raw input content are not exported.

## What the enzyme attribution does and does not mean

The `enzymes` / `enzyme_gene` columns record historical protein attribution. Recovered rules admit metabolism, substrate and shared-target relations; the historical CYP-labelled classes are not proof of a exclusively CYP-mediated pharmacokinetic mechanism. The separate current accession audit resolves protein identity, not the accuracy of this attribution.

Twenty-six primary pairs contain only labels in an analyst-defined target/atypical annotation set (`CYP19A1`, `PTGS1`, `PTGS2`, `HPRT1`, `XDH`, `AOX1`, `PGD`, `SLC31A1`); two are graded Major. This is a reproducible screening count, not a validated classification that every such pair acts only through shared pharmacodynamic targets. Prodrug/metabolite and combination-product labels also occur. Benchmark users should review these categories and define exclusions for their task explicitly.

### `protein_annotation_audit.csv` — 35 rows

Separate curation record; original pair and attribution snapshots are unchanged. `source_label` and `source_uniprot` join the raw long table. `canonical_uniprot` and `canonical_accession_status` report the official accession (all current primary accessions); `official_primary_gene`, `organism` and `entry_reviewed` identify the current primary gene and species (all reviewed Homo sapiens). `source_symbol_matches_official` is No only for Q4U2R8 and Q8TCC7. `long_rows` counts original attribution rows. `annotation_resolution` is `verified_unique_human_gene` for all 35; unresolved or ambiguous entries must be excluded from gene analyses. `retrieved_date`, `uniprot_release` and `source_url` freeze acquisition metadata (2026-10-05; 2026_03). `interpretation` explicitly separates protein identity from validation of the drug-pair relation.

`sources/uniprot_accession_gene_snapshot.tsv` contains the exact official 35-entry response, with HTTP release headers preserved beside it. UniProt Consortium annotations are CC BY 4.0 under their source terms. Browser-generated exports add `official_genes`, a semicolon-separated aggregate of accession-backed identities per unordered primary pair; raw `enzymes` remains separate.

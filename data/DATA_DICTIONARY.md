# Data dictionary

Enzyme-resolved DDI research resource, review snapshot 2.1.0. Author code is MIT; author contributions are licensed only within rights the author owns. Missing upstream source permissions remain unresolved.

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
| `mechanism` | string | Mechanism class(es), comma-separated. Controlled values: `CYP inhibition`, `CYP induction`, `transporter inhibition`, `transporter induction` |

Composition: 240 distinct drugs; 35 distinct source protein labels (34 after bare-symbol grouping);
mechanism-class occurrences CYP inhibition 1,547, transporter inhibition 667,
CYP induction 630, transporter induction 228. Clinical-trial and trueDDI flags were withdrawn because their exact provenance could not be established.

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
| `uniprot` | string | UniProt accession for that protein. **Not a reliable normalisation key** — `SLCO1B1 (OATP1B1)` carries `Q4U2R8` while bare `SLCO1B1` carries `Q9Y6L6`, although both denote the same transporter |
| `direction` | enum | Mechanism class for this pair × enzyme row: `inhibition`, `induction`, `transporter_inhibition`, `transporter_induction`. The bare terms denote CYP-mediated mechanisms; the `transporter_`-prefixed terms denote transporter-mediated ones |

Occurrences of `direction` are 1,547 `inhibition`, 667
`transporter_inhibition`, 630 `induction` and 228 `transporter_induction`, which
reconcile exactly with the mechanism-class occurrences in
`ddi_enzyme_database.csv`. `rebuild_severity.py` writes the same table with a
`severity` column to `data/derived/enzyme_pair_severity.csv`; that file is the
input to `pipeline/enzyme_severity_stats.py`.

---

## Exploratory usage layer: statistics

### `enzyme_severity_stats.csv` — 9 rows

Per-enzyme association with `Major` severity, for the nine enzymes carrying at
least **10 non-`Unknown`, inhibition-direction pairs**. (Fourteen enzymes have
≥10 pairs overall; the inhibition and non-`Unknown` restrictions are what reduce
the set to nine.)

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

Only `CYP2C9` (OR 2.78, *q* = 0.002) and `CYP3A4` (OR 0.44, *q* = 0.005)
reach FDR significance; `SLCO1B1 (OATP1B1)` is nominal only (OR 3.92,
*p* = 0.018, *q* = 0.055). Gene labels are normalised to the bare symbol before
aggregation, and the comparator is all other graded inhibition-attributed pairs
(disjoint arms); see `pipeline/enzyme_severity_stats.py`. The other seven rows do not reach the FDR threshold and are released for completeness; the
`direction_effect` label on a non-significant row describes the point estimate
only and must not be read as a finding.

### `enzyme_phenotype_enrichment.csv` — 434 rows

FDR-significant enzyme × adverse-event associations, covering 7 enzymes and
338 distinct HPO terms.

| Column | Type | Description |
| --- | --- | --- |
| `enzyme` | string | Enzyme / transporter gene symbol |
| `HPO` | string | Human Phenotype Ontology term identifier, e.g. `HP:0012410` |
| `AE_name` | string | HPO term label |
| `observed` | integer | Observed enzyme–phenotype observation count |
| `expected` | float | Expected count under the background model |
| `log2FE` | float | log₂ fold enrichment of `observed` over `expected` |
| `p` | float | Uncorrected *p*-value |
| `q` | float | FDR-adjusted *q*-value |

**Selected-only exploratory layer.** These previously reported signatures are not causal enzyme-attributable risk estimates. The original ontology snapshot, full observation matrix, all non-significant tests and correction universe are absent. Raw p and adjusted q are reported values, not independently recomputable from this selected-only export.

---

## Historical reported aggregate comparisons

### `drugbank_validation_summary.csv` — 8 rows (`metric`, `value`)

Reported counts from an unavailable original source snapshot: 1,172 enzyme pairs
DrugBank-mapped; 73 hard-proven confirmed (11 `Major`, 37 `Moderate`, 2
`Minor`, 23 `Unknown`); 72 Micromedex pairs mapped, 35 in the intersection with
the hard-proven set. Confirmation rate by grade: `Major` 25.0%, `Moderate`
15.0%, `Minor` 7.7%, `Unknown` 2.7%; Cochran–Armitage trend *p* = 1.8 × 10⁻¹⁸.

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

Enzyme labels retain source-level variation and should be normalised before
aggregating:

- `SLCO1B1` and `SLCO1B1 (OATP1B1)` both occur and denote the same transporter,
  **under two different UniProt accessions** (`Q9Y6L6`, the reviewed SwissProt
  entry, and `Q4U2R8`). The statistics merge them
  (*n* = 14); the raw labels are retained in this file.
- `SLCO1B3? (OATP)` retains an upstream uncertainty marker (`?`) indicating a
  tentative attribution in the source annotation.
- Transporter symbols carry parenthesised aliases (`ABCB1 (P-gp)`,
  `ABCG2 (BCRP)`, `ABCC2 (MRP2)`, `ABCC4 (MRP4)`, `SLCO1A2 (OATP1A2)`); CYP,
  UGT, SLC and other symbols do not.

---

## Provenance and reproducibility status

The exact available snapshots, SHA-256 values and explicit missing original
metadata are in `sources/sources.csv`. Runnable release processing steps and exact
script hashes are in `workflow.yaml`. Database publications identify a resource,
not the exact historical input version. The deposited GoldD3R export is contextual;
no documented transformation connects it to the primary enzyme-detail branch.
Original enzyme-detail files, clinical-grading download snapshot and external
comparison snapshots are absent. Structural validation and local grade rebuilding
do not constitute complete original-source reconstruction.

## What the enzyme attribution does and does not mean

The `enzymes` / `enzyme_gene` columns record the protein to which the source
annotation attributes the interaction. For the great majority of rows that is a
metabolising enzyme or a membrane transporter, which is the intended reading.
A small number of rows are different in kind and should be handled explicitly:

- **Shared pharmacodynamic target, not a pharmacokinetic mechanism.** 26 of the
  1,900 pairs (1.4%) are attributed *only* to a gene of this type — `PTGS2`
  (COX-2) for celecoxib + ketorolac, `CYP19A1` for anastrozole + letrozole
  (both are aromatase inhibitors, so `CYP19A1` is their shared target rather
  than the enzyme metabolising either one), `HPRT1`, `XDH`, `AOX1`, `PGD`,
  `PTGS1`, `SLC31A1`. Two of the 26 are graded `Major`.
- **Prodrug–metabolite pairs.** 5-fluorouracil + capecitabine and
  azathioprine + 6-mercaptopurine appear as interacting pairs, but in each case
  one member is a prodrug of the other. These co-occurrences are genuine in the
  source data but are not drug–drug interactions in the mechanistic sense.


Benchmark builders should exclude both categories; the row counts above make
that a cheap filter.

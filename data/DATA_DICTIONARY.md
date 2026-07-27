# Data dictionary

Enzyme-resolved DDI mechanism–severity resource. Released under CC BY 4.0.

Every file is UTF-8 CSV with a single header row and comma separators; text
fields containing commas are double-quoted. Row counts exclude the header and
were verified against the released files. **Column names below are the literal
header strings** — they are reproduced exactly as they appear in each file.

Drug names throughout are cleaned display names. For machine linkage use
`CID_pubchem` (PubChem) or DrugBank accessions rather than name matching;
`cid_severity_bridge.csv` is the authoritative name source.

Missing values are empty fields (not `NA`, not `NaN`). In
`enzyme_severity_stats.csv` the `E_value` column is deliberately empty for
enzymes that do not survive FDR correction.

---

## Primary table

### `ddi_enzyme_database.csv` — 1,900 rows

The primary resource: one row per drug pair.

| Column | Type | Description |
| --- | --- | --- |
| `drug_A` | string | First interacting drug (cleaned display name) |
| `drug_B` | string | Second interacting drug (cleaned display name) |
| `enzymes` | string | Responsible CYP enzyme(s) / transporter(s), comma-separated gene symbols; transporter symbols carry a parenthesised common alias, e.g. `ABCB1 (P-gp)` |
| `mechanism` | string | Mechanism class(es), comma-separated. Controlled values: `CYP inhibition`, `CYP induction`, `transporter inhibition`, `transporter induction` |
| `severity` | enum | DDInter clinical grade: `Major`, `Moderate`, `Minor`, `Unknown`. `Unknown` means DDInter records no grade — **not** that the interaction is absent or safe |
| `in_trial` | enum | `Yes` / `No` — pair appears in the clinical-trial DDI set |
| `trial_AEs` | integer | Number of trial adverse-event reports for the pair (0 when `in_trial` = `No`) |
| `in_trueDDI` | enum | `Yes` / `No` — pair appears in the trueDDI gold set |

Composition: 240 distinct drugs; 35 distinct enzymes/transporters; severity
93 `Major`, 419 `Moderate`, 50 `Minor`, 1,338 `Unknown`; mechanism-class
occurrences CYP inhibition 1,547, transporter inhibition 667, CYP induction
630, transporter induction 228; 79 rows `in_trial` = `Yes`, 229 rows
`in_trueDDI` = `Yes`.

---

## Pair × enzyme tables

### `enzyme_pair_severity.csv` — 3,072 rows

Long format: one row per pair × enzyme × direction.

| Column | Type | Description |
| --- | --- | --- |
| `CID_A` | integer | PubChem Compound ID of drug A |
| `drug_A` | string | Drug A name. **Not the cleaned display name** — 32 of the 240 names in this file are raw PubChem systematic titles (e.g. `Azane;cyclobutane-1,1-dicarboxylic acid;platinum` for carboplatin, `(E,Z)-Tamoxifen`, `Methotrexate, (A+-)-`). Join on `CID_A`/`CID_B`, never on names |
| `CID_B` | integer | PubChem Compound ID of drug B |
| `drug_B` | string | Drug B name — same caveat as `drug_A` |
| `enzyme_gene` | string | Single responsible enzyme / transporter gene symbol |
| `uniprot` | string | UniProt accession for that protein. **Not a reliable normalisation key** — `SLCO1B1 (OATP1B1)` carries `Q4U2R8` while bare `SLCO1B1` carries `Q9Y6L6`, although both denote the same transporter |
| `direction` | enum | Mechanism class for this pair × enzyme row: `inhibition`, `induction`, `transporter_inhibition`, `transporter_induction`. The bare terms denote CYP-mediated mechanisms; the `transporter_`-prefixed terms denote transporter-mediated ones |
| `severity` | enum | `Major` / `Moderate` / `Minor` / `Unknown` |

Occurrences of `direction` are 1,547 `inhibition`, 667
`transporter_inhibition`, 630 `induction` and 228 `transporter_induction`, which
reconcile exactly with the mechanism-class occurrences in
`ddi_enzyme_database.csv`.

### `enzyme_pair_drugbank_flagged.csv` — 1,799 rows

The pair × enzyme table restricted to DrugBank-mappable rows, with the
hard-proven flag attached.

| Column | Type | Description |
| --- | --- | --- |
| `CID_A`, `drug_A`, `CID_B`, `drug_B` | | As above |
| `enzyme_gene` | string | Responsible enzyme / transporter |
| `direction` | enum | As in `enzyme_pair_severity.csv`: `inhibition`, `induction`, `transporter_inhibition`, `transporter_induction` |
| `severity` | enum | `Major` / `Moderate` / `Minor` / `Unknown` |
| `DB_A` | string | DrugBank accession of drug A |
| `DB_B` | string | DrugBank accession of drug B |
| `drugbank_hard_proven` | enum | `Yes` / `No` — pair present in the DrugBank hard-proven set |

---

## Statistics

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
| `OR` | float | Odds ratio for `Major` severity, this enzyme vs all others |
| `CI_lo` | float | Lower bound, 95% confidence interval |
| `CI_hi` | float | Upper bound, 95% confidence interval |
| `p` | float | Uncorrected two-sided *p*-value |
| `q_fdr` | float | Benjamini–Hochberg FDR-adjusted *q*-value |
| `E_value` | float | E-value for unmeasured confounding; **empty when `q_fdr` ≥ 0.05** |
| `direction_effect` | enum | `enriched_Major` or `depleted_Major` |

Only `SLCO1B1 (OATP1B1)` (OR 5.24, *q* = 0.022), `CYP2C9` (OR 2.86,
*q* = 0.001) and `CYP3A4` (OR 0.49, *q* = 0.017) reach FDR significance. The
other six rows are null results and are released for completeness; the
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

**Caveat.** These are drug-class-**mediated** descriptive signatures, not causal
enzyme-attributable risk estimates. The file contains the FDR-significant
associations only, not the full underlying observation set.

---

## External validation

### `drugbank_validation.csv` — 73 rows

Pairs independently confirmed in the DrugBank hard-proven set.

| Column | Type | Description |
| --- | --- | --- |
| `drug_A`, `drug_B` | string | Interacting drugs |
| `DrugBank_A`, `DrugBank_B` | string | DrugBank accessions |
| `enzymes` | string | Responsible enzyme(s) / transporter(s) |
| `severity` | enum | Severity grade assigned by this resource |
| `drugbank_hard_proven` | enum | `Yes` for every row in this file |

### `drugbank_validation_summary.csv` — 8 rows (`metric`, `value`)

Counts behind the severity → confirmation gradient: 1,172 enzyme pairs
DrugBank-mapped; 73 hard-proven confirmed (11 `Major`, 37 `Moderate`, 2
`Minor`, 23 `Unknown`); 72 Micromedex pairs mapped, 35 in the intersection with
the hard-proven set. Confirmation rate by grade: `Major` 25.0%, `Moderate`
15.0%, `Minor` 7.7%, `Unknown` 2.7%; Cochran–Armitage trend *p* = 1.8 × 10⁻¹⁸.

### `kegg_ddi_validation.csv` — 1,243 rows

Testable pairs against the manually curated KEGG DDI set.

| Column | Type | Description |
| --- | --- | --- |
| `drug_A`, `drug_B` | string | Interacting drugs, both mapped to KEGG DRUG |
| `our_severity` | enum | Severity assigned by this resource |
| `in_kegg_ddi` | boolean | `True` if the pair is in the curated KEGG DDI set |
| `kegg_flag` | enum | KEGG severity flag: `CI` = contraindication, `P` = precaution; empty if not in KEGG. **All 204 flagged rows are `P`** — no contraindication-level pair survives the mapping |
| `kegg_class` | string | KEGG mechanism annotation, e.g. `CYP inhibition: CYP3A4`; empty if not in KEGG |

### `kegg_enzyme_concordance.csv` — 86 rows

| Column | Type | Description |
| --- | --- | --- |
| `drug` | string | Drug carrying both our and KEGG enzyme annotations |
| `our_enzymes` | string | Enzyme(s) attributed by this pipeline (list literal) |
| `kegg_enzymes` | string | Enzyme(s) from KEGG METABOLISM annotation (list literal) |
| `agree` | boolean | `True` if at least one enzyme is shared |

### `kegg_validation_summary.csv` — 13 rows (`metric`, `value`)

228 drugs mapped to KEGG; 3,034 KEGG DDI pairs; 1,243 testable. Confirmation by
grade: `Major` 35.1%, `Moderate` 35.4%, `Minor` 14.3%, `Unknown` 9.6%;
Cochran–Armitage *z* = 10.49, *p* = 9.6 × 10⁻²⁶; graded-vs-`Unknown` OR 4.82
(*p* = 8.0 × 10⁻²³); enzyme concordance 86.0% (*n* = 86).

### `trial_validation.csv` — 79 rows

**All** enzyme-annotated pairs observed in clinical-trial adverse-event reports,
not only the graded ones: 7 `Major`, 29 `Moderate`, 5 `Minor` and **38
`Unknown`**. The graded subset is 41 pairs.

| Column | Type | Description |
| --- | --- | --- |
| `drug_A`, `drug_B` | string | Interacting drugs (display names) |
| `drug_A_pubchem_title`, `drug_B_pubchem_title` | string | Verbatim PubChem compound titles |
| `enzymes` | string | Responsible enzyme(s) / transporter(s) |
| `ddinter_severity` | enum | DDInter clinical grade |
| `n_trial_AEs` | integer | Number of trial adverse-event reports for the pair |

### `trial_validation_summary.csv` — 8 rows (`metric`, `value`)

20,618 enzyme pairs total; 545 (2.6%) in the clinical-trial DDI set; 4,233
(20.5%) in the trueDDI gold set.

**Two caveats.** (1) The metric labelled `graded_and_in_trials = 79` counts
*every* trial-observed pair in this deposit, not only graded ones; 41 of the 79
carry a grade (7 `Major`, 29 `Moderate`, 5 `Minor`) and 38 are `Unknown`.
(2) The 2.6% and 20.5% rates are computed over a 20,618-pair enzyme-attribution
superset that is **not part of this deposit**; within the deposit the
corresponding rates are 79/1,900 = 4.2% and 229/1,900 = 12.1%. The two sets of
numbers are not comparable.

### `liddi_coverage_comparison.csv` — 7 rows (`metric`, `value`)

Coverage overlap with the independent LIDDI corpus: 4,070 LIDDI unique pairs vs
12,493 GoldD3R mechanism pairs; 265 LIDDI drugs vs 1,078 GoldD3R drugs; 169
drugs shared by CUI; **0 shared interaction pairs**, so all 4,070 LIDDI pairs
are novel relative to GoldD3R. LIDDI is complementary, not redundant.

---

## Identifier bridges and linkage layers

### `cid_severity_bridge.csv` — 448 rows

| Column | Type | Description |
| --- | --- | --- |
| `CID_pubchem` | integer | PubChem Compound ID |
| `drug_name` | string | PubChem compound title (verbatim) |
| `DDInterIDs` | string | Matched DDInter drug ID(s) |
| `joined_to_ddinter` | enum | `yes` / `no` |
| `display_name` | string | Cleaned display name — **authoritative** name used across all other files |

### `cid_drugbank_crosswalk.csv` — 255 rows

| Column | Type | Description |
| --- | --- | --- |
| `CID_pubchem` | integer | PubChem Compound ID |
| `DrugBank_id` | string | DrugBank accession, e.g. `DB00339` |

### `ddinter2_mechanism_annotations.csv` — 11,298 rows

DDInter 2.0 mechanism layer, retrieved from the DDInter 2.0 interaction server
(<https://ddinter2.scbdd.com>).

| Column | Type | Description |
| --- | --- | --- |
| `drug_a_name`, `drug_b_name` | string | Interacting drugs |
| `drugbankID_a`, `drugbankID_b` | string | DrugBank accessions |
| `internalID_a`, `internalID_b` | string | DDInter drug identifiers |
| `severity` | enum | DDInter severity grade |
| `metabolism` | 0/1 | Pharmacokinetic mechanism flag — metabolism |
| `synergistic_effect` | 0/1 | Pharmacodynamic mechanism flag — synergism |
| `antagonistic_effect` | 0/1 | Pharmacodynamic mechanism flag — antagonism |
| `absorption` | 0/1 | Pharmacokinetic mechanism flag — absorption |
| `distribution` | 0/1 | Pharmacokinetic mechanism flag — distribution |
| `excretion` | 0/1 | Pharmacokinetic mechanism flag — excretion |
| `others` | 0/1 | Mechanism flag — other / unclassified |
| `mechanism_text` | string | DDInter free-text mechanism description |

**Caveat.** Every one of the 11,298 rows carries `severity = Major`. This is a
Major-only slice, **not a complete export**. It is provided as a documentation and
identifier-linkage layer only and carries no analytical claim; do not compute
severity distributions or mechanism frequencies from it.

### `mimic_coprescription_examples.csv` — 13 rows

| Column | Type | Description |
| --- | --- | --- |
| `hadm_id` | integer | MIMIC-IV demo hospital admission identifier |
| `drug_A`, `drug_B` | string | Co-prescribed interacting drugs |
| `enzymes` | string | Responsible enzyme(s) / transporter(s) |
| `severity` | enum | Severity grade assigned by this resource |

**Caveat.** Drawn from the 100-patient MIMIC-IV demo cohort and illustrative
only; these rows support no prevalence, rate or risk estimate.

---

## Known naming variation

Enzyme labels retain source-level variation and should be normalised before
aggregating:

- `SLCO1B1` and `SLCO1B1 (OATP1B1)` both occur and denote the same transporter,
  **under two different UniProt accessions** (`Q9Y6L6`, the reviewed SwissProt
  entry, and `Q4U2R8`). Only `SLCO1B1 (OATP1B1)` was carried into
  `enzyme_severity_stats.csv`; merging the two changes that row from *n* = 12,
  50.0% Major, OR 5.24 (1.65–16.64) to *n* = 14, 42.9% Major,
  OR 3.92 (1.33–11.57), *p* = 0.018.
- `SLCO1B3? (OATP)` retains an upstream uncertainty marker (`?`) indicating a
  tentative attribution in the source annotation.
- Transporter symbols carry parenthesised aliases (`ABCB1 (P-gp)`,
  `ABCG2 (BCRP)`, `ABCC2 (MRP2)`, `ABCC4 (MRP4)`, `SLCO1A2 (OATP1A2)`); CYP,
  UGT, SLC and other symbols do not.

---

## Provenance of source datasets

| Source | Role | Identifier |
| --- | --- | --- |
| GoldD3R (D3-derived) | Mechanistic DDI corpus (UMLS CUI-indexed) | — |
| DDInter 2.0 | Clinical severity grades; mechanism layer | doi:10.1093/nar/gkae726 |
| PubChem | Compound identifiers and titles | doi:10.1093/nar/gkac956 |
| UniProt | Protein accessions | doi:10.1093/nar/gkac1052 |
| Human Phenotype Ontology | Adverse-event vocabulary | doi:10.1093/nar/gkaa1043 |
| DrugBank 5.0 | Independent validation (hard-proven set) | doi:10.1093/nar/gkx1037 |
| KEGG DRUG | Independent validation (curated DDI set) | Kanehisa *et al.* |
| MIMIC-IV demo | Illustrative co-prescription | doi:10.13026/dp1f-ex47 |
| LIDDI | Coverage comparison | doi:10.1007/978-3-319-25010-6_18 |
| Micromedex (Merative) | Commercial cross-check — counts only, not redistributed | — |


---

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
  5-fluorouracil + capecitabine additionally carries the largest trial
  adverse-event count in the deposit (465), which is the artefact you would
  expect from oncology co-coding rather than an interaction signal.

Benchmark builders should exclude both categories; the row counts above make
that a cheap filter.

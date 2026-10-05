# Review snapshot 2.1.0 — 5 October 2026

The release is narrowed to 11 public CSVs with explicit source boundaries. Clinical-trial/trueDDI primary flags and exports, hypothesis tables, unsupported prediction/ROC material, E-values, stale graphical abstract and browser screenshot are retired. All nine primary and sensitivity statistics remain available; SLCO1B1/OATP1B1 primary q rounds to 0.055. Tiny sensitivity q values use scientific notation rather than misleading zero display.

Sources/sources.csv identifies exact available public snapshots and unknown original metadata. Workflow.yaml identifies runnable steps and exact script SHA-256 values. The original enzyme-detail input, its transformations and permissions remain unresolved; GoldD3R is a separate contextual branch. DrugBank/KEGG are reported historical aggregate comparisons. Phenotypes are selected-only exploratory results.

Verification of this review snapshot: **276 public checks**, **23 additional privately reconstructed grade checks (299 total)**, reference joined-grade fingerprint match, byte-identical primary and sensitivity CSVs after SciPy regeneration, all six usage recipes completed, and **15 independent browser JavaScript checks** for search, exact-label filters, sorting, CSV round-trip, local grade import and fingerprint. The original supplied package also reproduced its 255/283 checks in a separate private baseline. No restricted source or locally reconstructed clinical grades were copied into the public record.

DOI 10.5281/zenodo.23169653 is reserved for an unpublished draft, not a published archival citation. Author code MIT; original author contributions CC BY only within owned rights; no blanket source rights certification.

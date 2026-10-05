# Source manifest

`sources.csv` identifies exact available snapshots with SHA-256, and explicitly records unavailable input metadata as `unknown`. Receipt of the author's package on 5 October 2026 is not substituted for an original retrieval date. A database publication or current download page is not an identifier for an exact historical snapshot.

`GoldD3R.txt` is the exact author-supplied contextual export: 21,897 pair lines, 12,493 labelled lines, 1,599 CUIs overall and 1,078 in labelled lines. Ignore the two log trailers. D3 original article: https://doi.org/10.1093/jamia/ocw128; corrigendum: https://doi.org/10.1093/jamia/ocz061. The export date, original acquisition metadata and source permissions are not independently established. No documented transformation connects this file to the 1,900-pair enzyme-detail branch.

The released primary and long-format tables are exact author-derived snapshots. Their hashes permit identity and integrity checks, not end-to-end reconstruction. Original enzyme files, full clinical-grading download and external comparison snapshots are missing. Clinical-trial and trueDDI layers have been excluded from the public release rather than presented as resolved validation.

See `workflow.yaml` for runnable processing steps, exact code hashes, outputs and boundaries. Rebuild clinical grades only from your own authorised DDInter download. Private test slices and historical graded files are excluded from the release.

#!/usr/bin/env python3
"""Recreate the separate protein identity audit from the frozen official snapshot.

This is an offline annotation-resolution step. Raw pair-level annotations remain
unchanged. Resolution of protein identity does not verify a drug-pair relation.
The frozen UniProt Consortium response is release 2026_03, retrieved 2026-10-05,
and used under the provider's CC BY 4.0 terms. Its exact hash prevents a current
or incomplete lookup from silently receiving these archived acquisition dates.
"""
import csv
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT_SHA = "ce9213b58819e4286207950e5eef9b94323363bd372fa4e77896a46c124d2277"


def main():
    snapshot = ROOT / "sources" / "uniprot_accession_gene_snapshot.tsv"
    if hashlib.sha256(snapshot.read_bytes()).hexdigest() != SNAPSHOT_SHA:
        raise SystemExit("Official protein snapshot hash differs from the frozen audit input.")
    with snapshot.open(encoding="utf-8", newline="") as fh:
        official_rows = list(csv.DictReader(fh, delimiter="\t"))
    official = {r["Entry"]:r for r in official_rows}
    if len(official) != len(official_rows):
        raise SystemExit("Official snapshot has ambiguous repeated accessions.")
    with (ROOT / "data" / "enzyme_pair_attribution.csv").open(encoding="utf-8", newline="") as fh:
        raw = list(csv.DictReader(fh))
    counts = Counter((r["enzyme_gene"],r["uniprot"]) for r in raw)
    output = []
    for (label, accession), n in sorted(counts.items()):
        entry = official.get(accession, {})
        gene = entry.get("Gene Names (primary)", "")
        resolved = entry.get("Organism") == "Homo sapiens (Human)" and gene and len(gene.split()) == 1
        output.append(dict(source_label=label, source_uniprot=accession,
            canonical_uniprot=entry.get("Entry", ""),
            canonical_accession_status="current_primary_accession" if accession == entry.get("Entry") else "unresolved",
            official_primary_gene=gene, organism=entry.get("Organism", ""), entry_reviewed=entry.get("Reviewed", ""),
            source_symbol_matches_official="Yes" if label.split(" ")[0].rstrip("?") == gene else "No",
            long_rows=n, annotation_resolution="verified_unique_human_gene" if resolved else "unresolved_excluded",
            retrieved_date="2026-10-05", uniprot_release="2026_03",
            source_url="https://rest.uniprot.org/uniprotkb/" + accession + ".tsv",
            interpretation="Protein identifier resolution only; drug-pair relation is not independently validated"))
    with (ROOT / "data" / "protein_annotation_audit.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"Protein audit: {len(output)} source label/accession records; "
          f"{sum(r['annotation_resolution'] == 'verified_unique_human_gene' for r in output)} resolved; "
          f"{sum(r['source_symbol_matches_official'] == 'No' for r in output)} source-symbol conflicts.")


if __name__ == "__main__":
    main()

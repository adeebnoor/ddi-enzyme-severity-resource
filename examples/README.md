# Worked examples

```bash
python3 examples/quickstart.py
```

Requires **pandas**. The resource itself needs nothing beyond the Python
standard library — `validate.py` and `scripts/make_figures.py` both run without
pandas — but the recipes here use it because it is what most reusers reach for.

[`quickstart.py`](quickstart.py) runs six recipes and prints a short, readable
result for each, so you can confirm a recipe does what its comment says before
adapting it.

| # | Recipe | What it demonstrates |
|---|---|---|
| 1 | Every Major-severity interaction mediated by CYP2C9 | Filtering the primary table on severity and enzyme. `enzymes` holds comma-separated symbols, so match with `str.contains`, not equality. |
| 2 | Severity profile of each enzyme | Why the long-format `enzyme_pair_severity.csv` beats the primary table for aggregation — one row per pair × enzyme means `groupby` just works. |
| 3 | All interactions for one drug, ranked by severity | Drugs appear in either `drug_A` or `drug_B`; this normalises to a "partner" column and sorts by clinical grade rather than alphabetically. |
| 4 | Only externally corroborated pairs | Combining the `in_trueDDI` and `in_trial` flags — the strictest filter the resource supports. |
| 5 | Joining to DrugBank | The pair table already carries accessions (`DB_A`, `DB_B`); the crosswalk is for when you are starting from your own PubChem CIDs. |
| 6 | The right way to handle `Unknown` | The most important caveat in the resource, made concrete. `Unknown` means *no grade recorded*, not *no interaction* — drop it or model it as missing, never as a negative. |

## Things worth knowing before you adapt these

- **Match enzymes with `str.contains`, not `==`.** The `enzymes` column holds
  comma-separated symbols, so `enzymes == "CYP3A4"` silently misses every pair
  where CYP3A4 acts alongside another enzyme.
- **Normalise enzyme symbols first if you aggregate.** `SLCO1B1` and
  `SLCO1B1 (OATP1B1)` are the same transporter under two labels.
- **Join on identifiers, not names.** Use `CID_pubchem` or DrugBank accessions.
  A few display names are PubChem brand titles (`Arthrotec`, `Anzemet`,
  `Respules`) and `Arthrotec` is a combination product.
- **`562` is the number of rows you can supervise on**, not 1,900. The rest are
  `Unknown`.

Full caveats are in the [repository README](../README.md#scope-caveats-and-known-limitations);
column definitions are in [`../data/DATA_DICTIONARY.md`](../data/DATA_DICTIONARY.md).

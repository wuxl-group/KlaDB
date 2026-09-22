# Motif data

This directory contains the public, record-level inputs and derived tables for the current KlaDB motif analysis.

| File | Contents |
| --- | --- |
| `KlaDB_motif_21aa.csv.gz` | 61,390 standard-amino-acid 21-aa windows centred on Kla (`K` at position 11). Fields: `site_id`, `species`, `uniprot_id`, `position`, and `window_21`. Rows remain at the site-record level; identical windows from different site records are not collapsed. |
| `KlaDB_motif_processing_summary.csv` | Per-species processing summary used for Supplementary Table S1: 64,072 site records, 61,390 retained 21-aa windows, and 2,682 excluded records. |
| `KlaDB_motif_enrichment_taxonomic.csv.gz` | Background-corrected enriched residues for all species and the three taxonomic groups, used for Supplementary Table S2. |
| `KlaDB_motif_enrichment_species.csv.gz` | Background-corrected enriched residues for each species, used for Supplementary Table S3. |

Records requiring terminal padding, or lacking a standard 21-aa sequence with the expected central lysine, were excluded from motif analysis. The enriched-residue tables retain non-central residues with Benjamini–Hochberg-adjusted `P < 0.05` and log2 enrichment `> 0`. They report the positive/background window totals and odds ratio for each residue-position test; `Displayed in figures` uses the more stringent log2-enrichment threshold `> 0.3`.

Reference proteome FASTA files are third-party UniProt source data and are not redistributed. Download them using the pinned proteome identifiers in `Scripts/motif_analysis/download_background_proteomes.py`, then run `Scripts/motif_analysis/export_enrichment_tables.py` to regenerate the two background-corrected tables. The original motif table can be regenerated from the current database export using `prepare_motif_windows.py`.

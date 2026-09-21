# Motif data

`KlaDB_motif_21aa.csv.gz` contains 61,390 standard-amino-acid, 21-aa windows centered on Kla (`K` at position 11). The windows were generated from the 64,072 current database site records. Records requiring terminal padding or having an invalid center/sequence were excluded from motif analysis.

Fields are `site_id`, `species`, `uniprot_id`, `position`, and `window_21`. Sequence rows are kept at the site-record level; identical windows from distinct site records are not collapsed.

The source folder supplied for revision contained 61,533 windows from an earlier extraction. The repository file was regenerated from the current database and is therefore the authoritative motif input for this release.

Run the scripts in `../Scripts/motif_analysis/` to regenerate the motif table, download reference proteomes, draw sequence logos, or perform background-corrected enrichment analysis.

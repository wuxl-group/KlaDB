# Database

`KlaDB_full_curated_sites.csv.gz` contains 64,072 curated Kla site records from the current KlaDB database. Each row represents one database site record and retains a stable `site_id`. When more than one study supports a record, PMIDs and dataset codes are joined with semicolons.

`KlaDB_studies.csv` contains publication- and dataset-level metadata for 28 studies, including organism, sample context, enrichment, mass spectrometry, search software, thresholds, and reported/retained counts.

The website's complete site-evidence download can contain more than 64,072 rows because one site may be supported by multiple study datasets. This repository's site file is deliberately one row per KlaDB site record.

See `data_dictionary.csv` for field definitions. Missing source metadata are reported as `Not reported` where applicable.

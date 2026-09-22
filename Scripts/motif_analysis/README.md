# Motif analysis scripts

From the repository root:

```bash
python Scripts/motif_analysis/prepare_motif_windows.py
python Scripts/motif_analysis/plot_sequence_logos.py
python Scripts/motif_analysis/download_background_proteomes.py --output-dir background_proteomes
python Scripts/motif_analysis/background_enrichment.py
python Scripts/motif_analysis/export_enrichment_tables.py
```

The public motif table is generated only from standard 21-aa windows centered on Kla. `export_enrichment_tables.py` writes the selected taxonomic and species-level background-corrected tables in `Motif/`; these are the machine-readable sources of Supplementary Tables S2 and S3. The UniProt downloads are intentionally not committed because they are third-party source data and can be regenerated with the pinned proteome identifiers in the script.

Dependencies are listed in `requirements.txt`.

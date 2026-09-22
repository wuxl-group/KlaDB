"""Create the public 21-aa motif table from the KlaDB site export."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sites", type=Path, default=Path("Database/KlaDB_full_curated_sites.csv.gz"))
    parser.add_argument("--output", type=Path, default=Path("Motif/KlaDB_motif_21aa.csv.gz"))
    parser.add_argument("--summary", type=Path, default=Path("Motif/KlaDB_motif_processing_summary.csv"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    kept = excluded = 0
    raw_by_species = defaultdict(int)
    kept_by_species = defaultdict(int)
    fields = ["site_id", "species", "uniprot_id", "position", "window_21"]
    with open_text(args.sites, "rt") as source, open_text(args.output, "wt") as destination:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for row in reader:
            raw_by_species[row["species"]] += 1
            window = row["window_51"].strip().upper()
            motif = window[15:36] if len(window) == 51 else ""
            if len(motif) != 21 or motif[10] != "K" or any(aa not in STANDARD_AA for aa in motif):
                excluded += 1
                continue
            writer.writerow({**{field: row[field] for field in fields[:-1]}, "window_21": motif})
            kept_by_species[row["species"]] += 1
            kept += 1
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", encoding="utf-8", newline="") as summary_file:
        writer = csv.DictWriter(summary_file, fieldnames=[
            "Species", "Raw sequences", "Retained 21-aa peptides", "Removed sequences", "Retention rate (%)"
        ])
        writer.writeheader()
        for species in sorted(raw_by_species):
            raw = raw_by_species[species]
            retained = kept_by_species[species]
            writer.writerow({
                "Species": species,
                "Raw sequences": raw,
                "Retained 21-aa peptides": retained,
                "Removed sequences": raw - retained,
                "Retention rate (%)": round(retained / raw * 100, 2),
            })
        writer.writerow({
            "Species": "Total",
            "Raw sequences": sum(raw_by_species.values()),
            "Retained 21-aa peptides": kept,
            "Removed sequences": excluded,
            "Retention rate (%)": round(kept / (kept + excluded) * 100, 2),
        })
    print(f"wrote={kept} excluded={excluded} output={args.output} summary={args.summary}")


if __name__ == "__main__":
    main()

"""Create the public 21-aa motif table from the KlaDB site export."""

import argparse
import csv
import gzip
from pathlib import Path


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sites", type=Path, default=Path("Database/KlaDB_full_curated_sites.csv.gz"))
    parser.add_argument("--output", type=Path, default=Path("Motif/KlaDB_motif_21aa.csv.gz"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    kept = excluded = 0
    fields = ["site_id", "species", "uniprot_id", "position", "window_21"]
    with open_text(args.sites, "rt") as source, open_text(args.output, "wt") as destination:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for row in reader:
            window = row["window_51"].strip().upper()
            motif = window[15:36] if len(window) == 51 else ""
            if len(motif) != 21 or motif[10] != "K" or any(aa not in STANDARD_AA for aa in motif):
                excluded += 1
                continue
            writer.writerow({**{field: row[field] for field in fields[:-1]}, "window_21": motif})
            kept += 1
    print(f"wrote={kept} excluded={excluded} output={args.output}")


if __name__ == "__main__":
    main()

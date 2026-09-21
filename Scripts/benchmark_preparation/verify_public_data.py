"""Verify public KlaDB database, motif, and benchmark files."""

import csv
import gzip
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def main():
    sites = list(rows(ROOT / "Database" / "KlaDB_full_curated_sites.csv.gz"))
    studies = list(rows(ROOT / "Database" / "KlaDB_studies.csv"))
    motifs = list(rows(ROOT / "Motif" / "KlaDB_motif_21aa.csv.gz"))
    assert len(sites) == 64072
    assert len(studies) == 28
    assert len({row["site_id"] for row in sites}) == len(sites)
    assert all(len(row["window_21"]) == 21 and row["window_21"][10] == "K" for row in motifs)

    sample_ids = set()
    label_counts = {"positive": 0, "putative_negative": 0}
    for label, dirname in [("positive", "positive_sets"), ("putative_negative", "putative_negative_sets")]:
        for path in (ROOT / "Benchmark" / dirname).glob("*.csv.gz"):
            for row in rows(path):
                assert row["sample_id"] not in sample_ids
                sample_ids.add(row["sample_id"])
                label_counts[label] += 1
    scores = list(rows(ROOT / "Benchmark" / "model_scores.csv.gz"))
    assert {row["sample_id"] for row in scores} == sample_ids
    assert len(sample_ids) == 108175
    print({"sites": len(sites), "studies": len(studies), "motifs": len(motifs), **label_counts, "benchmark": len(sample_ids)})


if __name__ == "__main__":
    main()

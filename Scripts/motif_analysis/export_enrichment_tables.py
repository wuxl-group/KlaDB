"""Export the public background-corrected motif-enrichment tables.

Reference proteomes are downloaded separately with
``download_background_proteomes.py`` and are intentionally not versioned in
this repository.  This script writes the selected enriched residues used in
Supplementary Tables S2 and S3 to ``Motif/``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from background_enrichment import (
    GROUP_MAP,
    count_21aa_sequences,
    count_background_proteome,
    enrichment_analysis,
    load_positive_sequences,
)


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
OUTPUT_DIR = PROJECT_DIR / "Motif"


def display_name(identifier: str) -> str:
    return identifier.replace("_", " ")


def selected_rows(frame: pd.DataFrame, dataset: str, dataset_id: str, level: str) -> pd.DataFrame:
    result = frame[
        (frame["position"] != 0)
        & (frame["qvalue"] < 0.05)
        & (frame["log2_enrichment"] > 0)
    ].copy()
    ranked = frame[frame["position"] != 0].copy()
    ranked["rank"] = ranked["log2_enrichment"].rank(method="first", ascending=False).astype(int)
    result = result.merge(ranked[["position", "aa", "rank"]], on=["position", "aa"], how="left")
    denominator = (result["positive_n"] - result["positive_count"]) * result["background_count"]
    result["odds_ratio"] = np.where(
        denominator != 0,
        result["positive_count"] * (result["background_n"] - result["background_count"]) / denominator,
        np.nan,
    )
    result["displayed"] = (result["qvalue"] < 0.05) & (result["log2_enrichment"] > 0.3)
    return pd.DataFrame({
        "Dataset": dataset,
        "Dataset_ID": dataset_id,
        "Level": level,
        "Rank by log2 enrichment": result["rank"].astype(int),
        "Position": result["position"].astype(int),
        "Amino acid": result["aa"],
        "Motif feature": result["aa"] + result["position"].map(lambda value: f"{value:+d}").str.replace("+", "", regex=False),
        "Positive count": result["positive_count"].astype(int),
        "Positive total": result["positive_n"].astype(int),
        "Positive frequency": result["positive_freq"],
        "Background count": result["background_count"].astype(int),
        "Background total": result["background_n"].astype(int),
        "Background frequency": result["background_freq"],
        "Odds ratio": result["odds_ratio"],
        "Log2 enrichment": result["log2_enrichment"],
        "P value": result["pvalue"],
        "Adjusted P value (BH)": result["qvalue"],
        "Displayed in figures": result["displayed"],
    })


def main() -> None:
    positive = load_positive_sequences()
    species_counts: dict[str, np.ndarray] = {}
    background_counts: dict[str, np.ndarray] = {}
    positive_totals: dict[str, int] = {}
    background_totals: dict[str, int] = {}
    species_tables: list[pd.DataFrame] = []

    for identifier, sequences in positive.items():
        pos_counts, pos_total = count_21aa_sequences(sequences)
        bg_counts, bg_total, _ = count_background_proteome(identifier)
        if bg_counts is None or bg_total == 0:
            raise RuntimeError(f"No usable reference-proteome background for {identifier}.")
        species_counts[identifier] = pos_counts
        background_counts[identifier] = bg_counts
        positive_totals[identifier] = pos_total
        background_totals[identifier] = bg_total
        stats = enrichment_analysis(pos_counts, pos_total, bg_counts, bg_total, identifier)
        species_tables.append(selected_rows(stats, display_name(identifier), identifier, "species"))

    group_tables: list[pd.DataFrame] = []
    all_identifiers = sorted(species_counts)
    for label, members, dataset_id, level in [
        ("All species", all_identifiers, "All_species", "overall"),
        ("Mammals", [key for key in all_identifiers if GROUP_MAP[key] == "Mammals"], "Mammals", "taxonomic_group"),
        ("Plants", [key for key in all_identifiers if GROUP_MAP[key] == "Plants"], "Plants", "taxonomic_group"),
        ("Insects", [key for key in all_identifiers if GROUP_MAP[key] == "Insects"], "Insects", "taxonomic_group"),
    ]:
        pos_counts = sum((species_counts[key] for key in members), np.zeros((21, 20), dtype=int))
        bg_counts = sum((background_counts[key] for key in members), np.zeros((21, 20), dtype=int))
        pos_total = sum(positive_totals[key] for key in members)
        bg_total = sum(background_totals[key] for key in members)
        stats = enrichment_analysis(pos_counts, pos_total, bg_counts, bg_total, dataset_id)
        group_tables.append(selected_rows(stats, label, dataset_id, level))

    group_output = pd.concat(group_tables, ignore_index=True)
    species_output = pd.concat(species_tables, ignore_index=True)
    for frame, path in [
        (group_output, OUTPUT_DIR / "KlaDB_motif_enrichment_taxonomic.csv.gz"),
        (species_output, OUTPUT_DIR / "KlaDB_motif_enrichment_species.csv.gz"),
    ]:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, index=False)
        print(f"Wrote {len(frame)} rows: {path}")


if __name__ == "__main__":
    main()

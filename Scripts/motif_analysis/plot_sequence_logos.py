"""Draw species, taxonomic-group, and pooled Kla sequence logos."""

import argparse
import gzip
from pathlib import Path

import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


AAS = list("ACDEFGHIKLMNPQRSTVWY")
POSITIONS = list(range(-10, 11))
GROUPS = {
    "Cavia porcellus": "Mammals", "Homo sapiens": "Mammals",
    "Mus musculus": "Mammals", "Rattus norvegicus": "Mammals",
    "Sus scrofa": "Mammals", "Glycine max": "Plants",
    "Oryza sativa": "Plants", "Triticum aestivum": "Plants",
    "Frankliniella occidentalis": "Insects",
}


def ppm(sequences):
    counts = pd.DataFrame(0.0, index=POSITIONS, columns=AAS)
    for sequence in sequences:
        for index, aa in enumerate(sequence):
            counts.loc[POSITIONS[index], aa] += 1
    return counts.div(counts.sum(axis=1), axis=0).fillna(0)


def draw(sequences, title, output):
    fig, axis = plt.subplots(figsize=(10, 3))
    logo = logomaker.Logo(ppm(sequences), ax=axis)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    axis.set(title=f"{title} (n={len(sequences):,})", xlabel="Position relative to lactylated lysine", ylabel="Frequency", xlim=(-10.5, 10.5))
    axis.set_xticks(POSITIONS)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("Motif/KlaDB_motif_21aa.csv.gz"))
    parser.add_argument("--output-dir", type=Path, default=Path("generated_figures/motif"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input)
    for species, frame in data.groupby("species", sort=True):
        draw(frame["window_21"].tolist(), species, args.output_dir / f"{species.replace(' ', '_')}_logo.png")
    data["group"] = data["species"].map(GROUPS)
    for group, frame in data.groupby("group", sort=True):
        draw(frame["window_21"].tolist(), group, args.output_dir / f"{group}_logo.png")
    draw(data["window_21"].tolist(), "All species", args.output_dir / "All_species_logo.png")


if __name__ == "__main__":
    main()

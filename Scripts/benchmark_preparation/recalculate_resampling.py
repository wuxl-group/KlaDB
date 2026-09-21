"""Recalculate deterministic 30-seed benchmark metrics from public files."""

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (accuracy_score, average_precision_score,
                             brier_score_loss, confusion_matrix, f1_score,
                             matthews_corrcoef, precision_score, roc_auc_score)


MODELS = {
    "DeepKla": "deepkla_score",
    "Auto-Kla": "autokla_score",
    "HybridKla": "hybridkla_score",
    "PCBert-Kla": "pcbert_kla_score",
}


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "little")


def load_samples(benchmark_dir):
    frames = []
    for label, dirname in [(1, "positive_sets"), (0, "putative_negative_sets")]:
        for path in sorted((benchmark_dir / dirname).glob("*.csv.gz")):
            frame = pd.read_csv(path, low_memory=False)
            frame["label"] = label
            frames.append(frame)
    base = pd.concat(frames, ignore_index=True)
    scores = pd.read_csv(benchmark_dir / "model_scores.csv.gz", low_memory=False)
    data = base.merge(scores, on="sample_id", validate="one_to_one")
    if len(data) != len(base) or data["sample_id"].duplicated().any():
        raise ValueError("Sample/score join failed")
    return data


def select_ids(data, group, ratio, seed):
    subset = data if group == "All species" else data[data["species"].eq(group)]
    positive = subset[subset["label"].eq(1)]["sample_id"]
    negative = subset[subset["label"].eq(0)]["sample_id"]
    if ratio == "all":
        return set(pd.concat([positive, negative]))
    count = min(len(negative), int(ratio) * len(positive))
    selected = negative.sample(n=count, replace=False, random_state=stable_seed(group, ratio, seed))
    return set(pd.concat([positive, selected]))


def metrics(frame, score_column):
    y = frame["label"].to_numpy(int)
    score = frame[score_column].to_numpy(float)
    pred = (score >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "AUC": roc_auc_score(y, score), "AUPRC": average_precision_score(y, score),
        "MCC": matthews_corrcoef(y, pred), "Sn": tp / (tp + fn),
        "Sp": tn / (tn + fp), "Precision": precision_score(y, pred, zero_division=0),
        "F1": f1_score(y, pred, zero_division=0), "ACC": accuracy_score(y, pred),
        "Brier": brier_score_loss(y, score), "TN": tn, "FP": fp, "FN": fn, "TP": tp,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", type=Path, default=Path("Benchmark"))
    parser.add_argument("--output-dir", type=Path, default=Path("generated_analysis"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = load_samples(args.benchmark_dir)
    groups = ["All species", *sorted(data["species"].unique())]
    rows = []
    for group in groups:
        for ratio in (1, 5, 10, "all"):
            for seed in range(30):
                ids = select_ids(data, group, ratio, seed)
                selected = data[data["sample_id"].isin(ids)]
                for model, score_column in MODELS.items():
                    rows.append({
                        "scope": "common_strict", "group": group, "ratio": str(ratio),
                        "seed": seed, "model": model, "threshold": 0.5,
                        "samples": len(selected), "positive": int(selected["label"].sum()),
                        "putative_negative": int(selected["label"].eq(0).sum()),
                        **metrics(selected, score_column),
                    })
    result = pd.DataFrame(rows)
    result.to_csv(args.output_dir / "resampling_30_seeds.csv", index=False)

    one = result[result["ratio"].eq("1")]
    metric_names = ["AUC", "AUPRC", "MCC", "Sn", "Sp", "Precision", "F1", "ACC", "Brier"]
    summary = []
    for (group, model), frame in one.groupby(["group", "model"], sort=True):
        row = {
            "scope": "common_strict", "negative_positive_ratio": "1:1",
            "group": group, "model": model, "resamples": len(frame),
            "positive_per_resample": int(frame["positive"].iloc[0]),
            "putative_negative_per_resample": int(frame["putative_negative"].iloc[0]),
            "samples_per_resample": int(frame["samples"].iloc[0]), "threshold": 0.5,
        }
        for metric in metric_names:
            row[metric] = frame[metric].mean()
            row[f"{metric}_low"] = frame[metric].quantile(0.025)
            row[f"{metric}_high"] = frame[metric].quantile(0.975)
        summary.append(row)
    pd.DataFrame(summary).to_csv(args.output_dir / "benchmark_all_species.csv", index=False, float_format="%.10f")
    print(f"Wrote {len(result)} resampling rows and {len(summary)} summary rows to {args.output_dir}")


if __name__ == "__main__":
    main()

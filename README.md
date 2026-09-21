# KlaDB public data and reproducibility resources

This repository accompanies KlaDB, a manually curated multi-species database of experimentally reported lysine lactylation (Kla) sites.

## Contents

- `Database/`: 64,072 curated site records and metadata for 28 source studies.
- `Motif/`: 21-aa Kla-centered windows suitable for motif analysis.
- `Benchmark/`: the external leakage-screened benchmark, model scores, 30-seed 1:1 results, sensitivity summaries, and pinned references to the four evaluated model repositories.
- `Scripts/`: motif-analysis and benchmark-preparation scripts.

## Clone the repository

The benchmark data and summary results can be viewed without downloading the model repositories. To also retrieve the four referenced model repositories, clone KlaDB with its Git submodules:

```bash
git clone --recurse-submodules https://github.com/wuxl-group/KlaDB.git
```

If KlaDB was cloned previously without submodules, initialize them with:

```bash
git submodule update --init --recursive
```

## Referenced prediction-model repositories

The original model repositories are included as Git submodules under `Benchmark/kla_prediction_models/`. The parent KlaDB repository pins a specific commit of each upstream repository, allowing the referenced source state to be identified without copying or rewriting the upstream history.

| Model | Submodule path | Upstream repository | Use in the KlaDB benchmark |
| --- | --- | --- | --- |
| DeepKla | `Benchmark/kla_prediction_models/DeepKla/` | <https://github.com/linDing-group/DeepKla> | Evaluated from the released repository and available model artifacts. |
| Auto-Kla | `Benchmark/kla_prediction_models/Auto-Kla/` | <https://github.com/tubic/Auto-Kla> | Evaluated from the released repository and available model artifacts. |
| HybridKla | `Benchmark/kla_prediction_models/HybridKla/` | <https://github.com/kongxianzw/HybridKla> | The published workflow was reconstructed locally; the submodule records the upstream source repository used as the reference. |
| PCBert-Kla | `Benchmark/kla_prediction_models/PCBert-Kla/` | <https://github.com/ZhangHongqi215/PCBert-Kla> | The model was trained locally according to the published workflow; the submodule records the upstream source repository used as the reference. |

The exact benchmark settings, local reconstruction/training status, environments, thresholds, input hashes, and available weight information are recorded in `Benchmark/model_settings.json` and `Benchmark/model_environment.json`. The submodules alone should not be interpreted as guaranteeing one-command reproduction of every model environment or pretrained checkpoint.

## Access

- KlaDB website: <https://kladb.wuxl-group.com>
- Online data downloads: <https://kladb.wuxl-group.com/download>
- Study provenance: <https://kladb.wuxl-group.com/studies>

## Licenses

Data files are released under CC BY 4.0; analysis code is released under the MIT License. See `DATA_LICENSE.md` and `LICENSE`.

The four model submodules are third-party projects. Their authorship, copyright, citation requirements, and licenses remain those specified in the respective upstream repositories. The KlaDB MIT License does not relicense third-party submodule content.

## Integrity

`SHA256SUMS.txt` lists checksums for the regular KlaDB release files. Git records the exact commit of each submodule. No ZIP bundle is required: files can be downloaded directly from the repository.

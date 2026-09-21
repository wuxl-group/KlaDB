# Benchmark

The primary comparison uses the `common_strict` external leakage-screened benchmark. Samples homologous to any model training sequence were excluded at 70% identity and at least 80% bidirectional coverage. Positives are experimentally reported Kla sites; negatives are putative/unlabeled lysines not reported as Kla sites.

For the main comparison, every positive was paired with one putative/unlabeled negative in 30 deterministic 1:1 resamples. `benchmark_all_species.csv` reports the 30-resample means and empirical 95% intervals. Threshold-dependent metrics use a fixed probability threshold of 0.5.

## Files

- `positive_sets/` and `putative_negative_sets/`: six-species common-strict sample pools.
- `model_scores.csv.gz`: one probability from each of the four evaluated models for every common-strict sample.
- `resampling_30_seeds.csv`: per-seed results for 1:1, 1:5, 1:10, and all-candidate analyses.
- `benchmark_all_species.csv`: 1:1 30-seed summary used by the website and Figure 4.
- `sensitivity_summary.csv`: identity-threshold and terminal-window sensitivity analysis.
- `filtering_audit.csv`: sample counts through the leakage-control workflow.
- `model_settings.json`: training-sequence provenance and checksums.
- `model_environment.json`: model source, local reconstruction status, thresholds, frameworks, CUDA/GPU context, and input/weight checksums where available.

The 1:5, 1:10, and all-candidate results are sensitivity analyses; the balanced 1:1 results are the primary comparison because the published models were generally developed and evaluated with approximately balanced classes.

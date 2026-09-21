# Benchmark preparation and verification

`build_strict_benchmarks.py` documents the all-species, all-label CD-HIT workflow and the 70%, 80%, and 90% identity filters. Rebuilding from raw inputs additionally requires the four published models' training sets and Docker for the pinned CD-HIT image.

`prepare_model_inputs.py` converts a strict dataset into model-specific input formats.

From the repository root, verify the published tables and regenerate all 30-seed metrics with:

```bash
python Scripts/benchmark_preparation/verify_public_data.py
python Scripts/benchmark_preparation/recalculate_resampling.py
```

Regenerated files are written to `generated_analysis/` and do not overwrite the distributed results.

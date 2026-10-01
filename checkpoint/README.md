# Frozen assets

The inference exporter creates tensor-only weights plus a SHA256 manifest.
It reads trusted local experiment checkpoints, not arbitrary downloaded pickle files.

```bash
pixi run python scripts/export_assets.py /path/to/trusted/step2_experiment
```

The manifest identifies all six selected policy branches per route, codec/decoder,
DiMA frozen reference, normalization, length distribution, calibration and official
reference FASTA. Inference uses weights_only=True and strict state-dict loading.

The approved frozen inference weights/reference bundle is distributed as a
[release asset](https://github.com/nicetone9/amp-challenge-submission/releases/tag/frozen-pool-seed42-20260930),
not as Git blobs. `inference-checkpoint-seed42.tar.gz` has SHA256
`10006178df027959ff2c25a619b9e2e4440d3083a737090a056f23d20bee8e12`.
It contains only manifest-listed generator inputs plus the reference-only bundle.
It does not contain full training/test data or external predictor weights.

For explicit model resampling, download and verify the archive, then extract its
contents under `checkpoint/` in a fresh clone. The entry verifies both manifests
before sampling and uses tensor-only strict loading. No random-weight fallback is
provided. Public asset availability is not a successful model-resampling certificate.

The default fixed-pool filtering entry needs no model checkpoint or CUDA GPU.

External oracle assets belong in external/ANIA and external/HemoPI2, with
external/manifest.json containing a sha256 mapping of exact relative files.
No oracle downloads execute automatically. Do not substitute versions silently.

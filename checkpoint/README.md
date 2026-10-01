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
not as Git blobs. Use `inference-checkpoint-seed42-with-notices.tar.gz`, SHA256
`72f9327ee1f5ed9dd113d7901b026d241eb3e09da53030f6c78a4f2c2af171a4`.
It preserves full DiMA, ESM and AMP Challenge license notices. Model tensors and
input manifests match the first archive exactly. The earlier archive remains for
provenance; use the notices-inclusive archive for new downloads.
It contains manifest-listed generator inputs, the reference-only bundle and notices.
It does not contain full training/test data or external predictor weights.
The DiMA decoder includes parameters derived from the pretrained ESM head;
the full ESM encoder is not bundled. Original-code MIT does not relicense third-party
weights/data or resolve incomplete source-data redistribution clearance.

For explicit model resampling, download and verify the archive, then extract its
contents under `checkpoint/` in a fresh clone. The entry verifies both manifests
before sampling and uses tensor-only strict loading. No random-weight fallback is
provided. Public asset availability is not a successful model-resampling certificate.

The default fixed-pool filtering entry needs no model checkpoint or CUDA GPU.

External oracle assets belong in external/ANIA and external/HemoPI2, with
external/manifest.json containing a sha256 mapping of exact relative files.
No oracle downloads execute automatically. Do not substitute versions silently.

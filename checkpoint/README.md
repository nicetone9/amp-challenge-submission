# Frozen assets

The inference exporter creates tensor-only weights plus a SHA256 manifest.
It reads trusted local experiment checkpoints, not arbitrary downloaded pickle files.

```bash
pixi run python scripts/export_assets.py /path/to/trusted/step2_experiment
```

The manifest identifies all six selected policy branches per route, codec/decoder,
DiMA frozen reference, normalization, length distribution, calibration and official
reference FASTA. Inference uses weights_only=True and strict state-dict loading.

Assets exist in the Lane staging checkout but are ignored by Git until weight/data
redistribution terms and the delivery mechanism are approved. A new clone therefore
fails clearly rather than substituting random weights or replaying old FASTA files.

External oracle assets belong in external/ANIA and external/HemoPI2, with
external/manifest.json containing a sha256 mapping of exact relative files.
No oracle downloads execute automatically. Do not substitute versions silently.

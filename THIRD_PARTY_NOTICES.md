# Third-party notices

- DiMA, https://github.com/MeshchaninovViacheslav/DiMA,
  commit 18f4a2e67988efe1cc5593fcaec1ece8f09fdac0.
  MIT, copyright 2024 Petr Grinberg; full notice in third_party/DiMA_LICENSE.
  Checked against the upstream LICENCE at that exact commit on 2026-09-30.
  Five upstream source files are in src/amp_submission/vendor. The score estimator's
  import is made package-relative. No structural components are included.
- AMP Challenge template, https://github.com/szczurek-lab/amp-challenge-2027.
  scripts/verify_submission.py retained unchanged, BSD-3-Clause notice retained in
  third_party/AMP_CHALLENGE_LICENSE. Its SHA256 is
  3f2eb1bd61200abfccf07d86e9c226d57f3d12abcf25715af1d90f41531942cf.
- ESM-2, https://huggingface.co/facebook/esm2_t33_650M_UR50D and
  https://github.com/facebookresearch/esm, MIT, copyright Meta Platforms, Inc.
  and affiliates; full notice in third_party/ESM_LICENSE. The full frozen encoder
  is not bundled. The trained DiMA decoder was initialized from the upstream
  ESM LM head and tied embedding/output table; those derived parameters are
  included in the optional inference checkpoint.
- ProtT5 weights are external and not bundled.
- ANIA and HemoPI2 are external dependencies, not vendored or relicensed here.
  HemoPI2 source is GPL-3.0; redistribution/combination and all weight terms need review.
- MishaLaskin/vqvae and mahdip72/vq_encoder_decoder were method references;
  their repository contents are not copied into this package.

The repository MIT license applies only to original code and documentation, not
third-party data, predictions, weights or upstream source under its own notice.

Built wheels include this notice and both full upstream license files under
amp_submission/third_party. Original-code MIT terms do not replace those notices.
No license is asserted here for unreleased nonpublic training data; see docs/DATA.md.

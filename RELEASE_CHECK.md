# Release check — paper-refreshment-designs-replication

Date: 2026-09-29T04:38:29Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-refreshment-designs-replication/tar.gz/main`.

- ref: `main`; commit: `3ae4e29a1ee382697fd07be729e43eb13c85fcb5`
- archive SHA-256: `6f91f7aead1ffe4f1ba4ab5b7da18ab5a7fee5e51e92ea76b17431ecaa7acb1c`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 79; listed in FILE_MANIFEST.txt: 79; missing from snapshot: 0; not listed in manifest: 0
- content scan of the published snapshot: 0 file(s) matched
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (491s); log and sessionInfo kept; comparison below
- third-party reproduction: none; this record is the author's own re-execution.

## Environment of the clean run

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Sequoia 15.7.3

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] lpSolve_5.6.23

loaded via a namespace (and not attached):
[1] compiler_4.6.0

Python side:
Python 3.14.7
```

## Regenerated vs shipped outputs

```
file                                   max |diff|          cells/tokens  status
check_design_rank_output.txt           0.000e+00           -             identical
check_funnel_examples_output.txt       0.000e+00           -             identical
check_multiwave_completion_output.txt  1.133e-14           550           within tolerance
sim1_results.csv                       2.992e-17           221           within tolerance
sim2_plim.csv                          0.000e+00           -             identical
sim2_results.csv                       9.714e-17           68            within tolerance
sim3_diagnostics.csv                   5.502e-18           143           within tolerance
sim3_estimators.csv                    9.992e-16           976           within tolerance
sim3_pairs.csv                         3.036e-17           77            within tolerance
sim3_partD.csv                         9.758e-19           36            within tolerance
sim3_partD_bootstrap.csv               0.000e+00           -             identical
sim3_summary.txt                       0.000e+00           -             identical
sim3_variance_check.csv                9.975e-18           208           within tolerance

files compared: 13; identical: 5; within tolerance: 8; FAILED: 0; largest numeric difference: 1.133e-14; tol: 1e-08
comparison passed
```

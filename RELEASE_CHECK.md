# Release check — paper-refreshment-designs-replication

Date: 2026-09-25T05:33:16Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-refreshment-designs-replication/tar.gz/d516967f9065be73186a67438241e3a8c28db8e2`.

- ref: `d516967f9065be73186a67438241e3a8c28db8e2`; commit: `d516967f9065be73186a67438241e3a8c28db8e2`
- archive SHA-256: `7c2c71dbb087324e738540e1a34f43d357c9f0d31a932aed9c39523cd5e6a453`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 72; listed in FILE_MANIFEST.txt: 72; missing from snapshot: 0; not listed in manifest: 0
- content scan of the published snapshot: 0 file(s) matched
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (449s); log and sessionInfo kept; comparison below
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

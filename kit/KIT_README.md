# The JLPS pipeline behind Section 10 — data, access, and what each number comes from

Section 10 of this paper reports an application to the Japanese Life Course Panel Surveys (JLPS).
**No individual record, and no file derived from individual records, is in this archive or in its
Git history.** What is here is the code that produces the section's numbers, a synthetic input with
the same structure as the real one, and the expected outputs of running the code on that synthetic
input. Everything the paper quotes is an aggregate with cells below ten observations suppressed.

## 1. The data, and how a third party obtains them

The file analysed is the **integrated release (waves 1–19) distributed to participants of the panel
survey project**, not a public release obtained by application, and it therefore carries no archive
study number. A reader with a licence of their own obtains comparable waves by applying to the
**Social Science Japan Data Archive (SSJDA, University of Tokyo)**. The study numbers recorded in the
author's companion archive for the Japanese-language paper are `PY130` for JLPS-Y and `PM130` for
JLPS-M, with `PY130_add2` for the 2011 additional sample and `PY130_refresh` for the 2019 refresh
sample. **Confirm the current study numbers and wave coverage with SSJDA before relying on them**:
they are reproduced here from the author's own earlier archive, not verified against the catalogue
for this paper, and the releases differ from the integrated file, so respondent counts and cleaning
need not match. Byte-identical reproduction of Section 10 requires the same integrated release; from
an SSJDA public release this is a re-analysis with the same design and definitions, not a replication.

Restricted data are never placed inside the project folder. The pipeline reads them from a directory
outside it, named by the environment variable `P1_DATA_DIR`, and `00_config.R` aborts if that
directory resolves inside the project tree.

## 2. What a licensed user runs

```sh
export P1_DATA_DIR=<the directory holding the integrated .dta>
export P1_PROJECT_DIR=<this archive's kit/ directory>
Rscript R/04_build_analysis.R          # derives the wave-5 risk set and the item universe
Rscript R/15_panelcond_designs.R       # the five designs, the diagnostics, the funnel ratios
Rscript R/15b_item_flags.R             # routing and hygiene flags for the item-nonresponse family
Rscript R/11_negcontrol_w13.R          # the 2019-episode negative-control battery (Appendix table)
Rscript R/11b_negcontrol_w13_boot.R --B 2000   # its pooled value with a person-level bootstrap SE
Rscript check_manuscript_values_T2.R <results dir>   # reprints every number quoted in Section 10
```

Requires R (≥ 4.3) with `haven`, `data.table`, and `panelcond`
(`remotes::install_github("sokubo/panelcond@v0.1.4")`). The reported design run used 0.1.2; every
standard error and diagnostic quoted in Section 10 comes from the joint person bootstrap of
`15_panelcond_designs.R`, which recomputes every arm and every group share in each replicate, so the
correction of the analytic entry-wave and diagnostic variances in 0.1.4 (they now include the
estimated shares) does not change any quoted number. The design run takes about 13 minutes with 500
bootstrap replications over two dose definitions; `11b` takes a few minutes with 2,000 replications.

`11_negcontrol_w13.R` and `11b_negcontrol_w13_boot.R` read an item map,
`results/03_negcontrol_sets_v1.csv`, listing for each of the 23 battery items its variable name in the
2007, 2011 and 2019 entry questionnaires. The map is built from the provider's variable labels and is
therefore **not** shipped here (like `varmap_w1-19.csv`); the 23 items are listed in the paper's
appendix table, and the author supplies the map privately to licensed users on request.

Without a licence, the same code runs on the synthetic input, self-contained:

```sh
export LC_ALL=C.UTF-8                    # the item labels are Japanese; see the archive README
export P1_DATA_DIR=/tmp/p1synth P1_PROJECT_DIR=$PWD P1_RESULTS_DIR=/tmp/p1synth_out
Rscript R/15_make_synthetic.R            # regenerates synthetic/*.dta byte-identically
Rscript R/15_make_synthetic_varmap.R     # the w1 <-> w5 variable map for the synthetic items
export P1_VARMAP=$P1_DATA_DIR/varmap_synthetic.csv
Rscript R/04_build_analysis.R
Rscript R/15_panelcond_designs.R --B 50
Rscript R/15b_item_flags.R
diff -r /tmp/p1synth_out synthetic_out    # should differ only in 15_env.txt (paths and timestamps)
```

Run the commands from this `kit/` directory; each script finds `00_config.R` and its other companion
files in its own folder (`R/`). Both `P1_DATA_DIR` and `P1_PROJECT_DIR` must be set: the configuration
has no default location and stops with a message if either is missing. `15_make_synthetic.R` fixes the
date stamp that `haven::write_dta()` writes into the Stata header (otherwise the current time), so the
regenerated input has the SHA-256 given in the archive README. The expected outputs in `synthetic_out/`
were produced with panelcond 0.1.4. With 0.1.2 or 0.1.3 the analytic columns differ (`se_analytic` of
the entry-wave arms in `15_designs_items.csv`; `T1_stat`, `T1_p`, `T2_stat`, `T2_p` in
`15_tests_items.csv`) because those versions omitted the share terms; every bootstrap column and every
summary file is identical, and no quoted number uses an analytic column.

Ten seconds end to end. It exercises every code path the real run uses, including the entry-wave
correction (36 of 42 synthetic items have an entry-wave counterpart). The synthetic run is a
**pipeline check only**: its numbers are not the paper's and must never be reported as such — its
analysis arms are 2,682 and 963 where the paper's are 2,797 and 963.

The real run additionally reads the provider's variable list, `varmap_w1-19.csv`, which maps each
item to its column in every wave. That file is the data provider's documentation and is **not**
shipped here; a licensed user already has it with the data, and points `P1_VARMAP` at it.
`R/15_make_synthetic_varmap.R` builds the equivalent for the synthetic input, so the public pipeline
needs nothing from outside this archive.

## 3. Where each number in Section 10 comes from

`check_manuscript_values_T2.R` performs every comparison below and prints the result; on 23 September
2026 it reproduced all 48 quantities from the aggregate outputs of the real run (34 from scripts 04, 11
and 15; 14 from the person-level bootstrap of script 11b, which run once its outputs exist).

| Quantity in Section 10 | Output file | Column |
|---|---|---|
| 2,797 continuing survivors; 963 entrants; 574 and 540 survivors | `15_arms.csv` | `n_old_S`, `n_new_total`, `n_new_sm`, `n_new_sm1` |
| 500 bootstrap replications | `15_arms.csv` | `B` |
| 268 items with the same coding at entry | `15_arms.csv` | `n_ec_items` |
| 470 substantive items with variation | `15_detection_counts.csv` | `n_items`, naive row |
| flagged 28 / 19 / 13 over all items; 22 and 21 on 265 | `15_detection_counts.csv` | `n_affected` by `estimator` |
| flagged 19 / 15 / 9 / 22 / 21 on the common set of 265 | `15_designs_items.csv` | `class3 == "affected"`, restricted to items with an `ec` row |
| four items flagged by all five designs | `15_designs_items.csv` | intersection over `estimator` (DQ26, DQ39, DQ44_4A, DQ55_Q) |
| 13.1% and 8.2% diagnostic rejections | `15_diagnostics_summary.csv` | `share_T1_p05`, `share_T2_p05` |
| 19% for item nonresponse | `15_diagnostics_summary.csv` | `share_T2_p05`, family `B_itemnonresp` (.186) |
| extreme-category fall of .21–.25; midpoint rise of .19–.21 | `15_designs_items.csv` | `d_std` for `pdq_ext_share`, `pdq_mid_share` |
| employment, 4.9 points, losing significance when standardised | `15_designs_items.csv` | `DQ02`, `estimate` and `q` for `ec` and `ec_adj` |
| seven of 405 items with a mass ratio above one; 1.96 largest | `15_designs_items.csv` | `funnel_ratio` |
| 23 negative-control items; median \|d\| .033; 1 rejection; 9 equivalences; pooled −.035 (independence SE .0095) | `11_negcontrol_w13_summary.csv` | all columns |
| person-level bootstrap SE .0186 (2.0 times the independence SE; design effect 3.84); 95% interval [−.072, .001]; −.028 (SE .0185) without the rejected item | `11b_negcontrol_w13_pooled.csv` | `se_person_bootstrap`, `design_effect`, `ci95_lo`, `ci95_hi`; the row with the BH-rejected item dropped |
| mean correlation .116 between item contrasts across replications | `11b_negcontrol_w13_corr.csv` | `mean_offdiag_corr` |
| implied bias .017–.052 over Γ ∈ [.5, 1.5]; .036–.108 (.072Γ) at the far end of the interval | `11b_loading_grid.csv` | `correction`, `ci95_lo` (sign reversed) |
| the 23 item contrasts of the appendix table | `11_negcontrol_w13.csv` | `d`, `se`, `q` |
| n = 2,638 and n = 619 at the 2019 episode | `02_response_by_wave_cn.csv` | wave 13 by cohort |

## 4. The aggregate outputs themselves

They are **not** shipped here. The paper's Section 10 reports them in summary, but releasing the
item-level files is a finer-grained disclosure than the paper makes, and the disclosure conditions of
the survey and of the data provider have not been confirmed for that. Suppression of cells below ten
observations does not by itself authorise release. A licensed user regenerates them with the sequence
in §2; the author can supply them privately on request.

## 5. The freeze record

Because the real-data run cannot be repeated by a third party, the author records it: the identity of
the input (file name, byte size, SHA-256, dimensions — **never its content**), the SHA-256 of every
script executed and every output produced, the environment, and a comparison of every number
transcribed into the manuscript with the outputs. `15_env.txt` already carries the input identity and
the environment. That record is supplied privately on request; it is the author's own
frozen re-execution and is **not** third-party reproduction, and the paper says so.

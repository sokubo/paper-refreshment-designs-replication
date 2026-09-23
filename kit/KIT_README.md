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
**Social Science Japan Data Archive (SSJDA, University of Tokyo)**. The current releases are the study
numbers `PY160` (JLPS-Y) and `PM160` (JLPS-M): the JLPS user page of the Center for Social Research and
Data Archives (https://csrda.iss.u-tokyo.ac.jp/socialresearch/user/, checked 23 September 2026)
recommends that users apply for "the latest data set (PY160/PM160)", with corrections reflected
through wave 16. The numbers `PY130`/`PM130` (with `PY130_add2` for the 2011 additional sample and
`PY130_refresh` for the 2019 refresh sample) recorded in the author's earlier archive refer to older
releases. **Confirm the study numbers and wave coverage with SSJDA before applying.** Every SSJDA
release differs from the integrated file analysed here, so respondent counts and cleaning need not
match. Byte-identical reproduction of Section 10 requires the same integrated release; from
an SSJDA public release this is a re-analysis with the same design and definitions, not a replication.

Restricted data are never placed inside the project folder. The pipeline reads them from a directory
outside it, named by the environment variable `P1_DATA_DIR`, and `00_config.R` aborts if that
directory resolves inside the project tree.

## 2. What a licensed user runs

```sh
export P1_DATA_DIR=<the directory holding the integrated .dta>
export P1_PROJECT_DIR=<this archive's kit/ directory>
export P1_INPUT_DTA=<file name of the integrated .dta inside $P1_DATA_DIR/raw>   # one named input
Rscript R/15_test_mass_ratio.R         # unit test of the mass-domination helper (no data read)
Rscript R/04_build_analysis.R          # derives the wave-5 risk set and the item universe
Rscript R/15_panelcond_designs.R       # the five designs, the diagnostics, the mass-domination ratios
Rscript R/15b_item_flags.R             # routing and hygiene flags for the item-nonresponse family
Rscript R/11_negcontrol_w13.R          # the 2019-episode negative-control battery (Appendix table)
Rscript R/11b_negcontrol_w13_boot.R --B 2000   # its pooled value with a person-level bootstrap SE
Rscript check_manuscript_values_T2.R <results dir> --selftest   # every number in Section 10, then the checker's self-test
```

**One named input.** Every script reads the file named by `P1_INPUT_DTA`; if it is unset, the single
master `.dta` in `$P1_DATA_DIR/raw` is used and the run stops if there is more than one candidate (no
naming heuristic chooses between versions). Every script prints the file's name, size, MD5 and SHA-256
and writes them to its run record (`15_env.txt`, `11_env.txt`, `11b_env.txt`; the `11b` record of the frozen run was
written before the size line was added and carries the name, MD5 and SHA-256). `04` saves the input's
fingerprint, the row positions and cohorts of the wave-5 risk set, and the respondent identifiers, and
`15` stops unless it reads the same file and reconstructs the same persons in the same order. The
checker requires the three run records to name one and the same file with the SHA-256 given in §5.

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

The synthetic example exercises the scripts of the 2011 episode (`04`, `15`, `15b`). It does **not**
exercise the 2019 battery scripts (`11`, `11b`), which need the confidential item map described above,
and it does not contain a category reported by stayers but by no fresh entrant; that case, a
zero/zero codebook category, a rare category and an ordinary one are tested by
`R/15_test_mass_ratio.R` on constructed data.

Run the commands from this `kit/` directory; each script finds `00_config.R` and its other companion
files in its own folder (`R/`). Both `P1_DATA_DIR` and `P1_PROJECT_DIR` must be set: the configuration
has no default location and stops with a message if either is missing. `15_make_synthetic.R` fixes the
date stamp that `haven::write_dta()` writes into the Stata header (otherwise the current time), so the
regenerated input has the SHA-256 given in the archive README. The expected outputs in `synthetic_out/`
were produced with panelcond 0.1.4. With 0.1.2 or 0.1.3 the analytic columns differ (`se_analytic` of
the entry-wave arms in `15_designs_items.csv`; `T1_stat`, `T1_p`, `T2_stat`, `T2_p` in
`15_tests_items.csv`) because those versions omitted the share terms; every bootstrap column and every
summary file is identical, and no quoted number uses an analytic column.

About half a minute end to end. It exercises every code path of the 2011-episode scripts that the
real run uses, including the entry-wave correction (36 of 42 synthetic items have an entry-wave
counterpart) and the exclusion of nominal codes from the counts (the synthetic marital-status item
`DQ43`, as in the real run). The synthetic run is a
**pipeline check only**: its numbers are not the paper's and must never be reported as such — its
analysis arms are 2,682 and 963 where the paper's are 2,797 and 963.

The real run additionally reads the provider's variable list, `varmap_w1-19.csv`, which maps each
item to its column in every wave. That file is the data provider's documentation and is **not**
shipped here; a licensed user already has it with the data, and points `P1_VARMAP` at it.
`R/15_make_synthetic_varmap.R` builds the equivalent for the synthetic input, so the public pipeline
needs nothing from outside this archive.

## 3. Where each number in Section 10 comes from

`check_manuscript_values_T2.R` performs every comparison below and prints the result. It requires all
sixteen aggregate files and run records listed at its top (a missing file stops it), requires every key
to identify exactly one row, and first checks that the run records of scripts 15, 11 and 11b name the
same input file with the SHA-256 given in §5. On 24 September 2026 (JST) it reproduced all 89 quantities from the aggregate outputs of the frozen rerun, and its self-test detected all nine corruptions. One statement of Section 10 is a description rather than a number and is not asserted: that most positive/zero items are asked only of a subgroup (read from the item labels).

| Quantity in Section 10 | Output file | Column |
|---|---|---|
| 2,797 continuing survivors; 963 entrants; 574 and 540 survivors | `15_arms.csv` | `n_old_S`, `n_new_total`, `n_new_sm`, `n_new_sm1` |
| 500 bootstrap replications | `15_arms.csv` | `B` |
| 268 items with the same coding at entry | `15_arms.csv` | `n_ec_items` |
| 478 items in the inventory; 76 outside the counts: the 44 of `R/15_item_exclude.csv`, 31 further items with 04's routing flag (three flagged items are also on the list), and one follow-up of a flagged question (DQ46Y) | `15_designs_items.csv` | `exclude_flag`, `nr_routing_flag`, `routing_followup`, `count_exclude` (= any of the three) |
| 402 items for which a mean contrast is meaningful; 401 for the matched designs (one item without variation in the matched arms); 227 with an entry wave | `15_detection_counts.csv` | `n_items` by `estimator` |
| flagged 22 / 16 / 9 / 20 / 19 | `15_detection_counts.csv` | `n_affected` by `estimator` |
| flagged 15 / 12 / 6 / 20 / 19 on the common set of 227 | `15_designs_items.csv` | `class3 == "affected"`, restricted to items with an `ec` row that enter the counts |
| two items flagged by all five designs | `15_designs_items.csv` | intersection over `estimator` (DQ26, DQ44_4A) |
| diagnostic rejections: 30 of 227 (13.2%) and 36 of 401 (9.0%) | `15_diagnostics_summary.csv` | `n_T1_reject`/`n_T1_tests`, `n_T2_reject`/`n_T2_tests`, `share_T1_p05`, `share_T2_p05`, family `A_substantive` |
| 66 of 396 (16.7%) for item nonresponse; 50 of the 66 from the grids DQ58C, DQ04(3), DQ09 and DQ08D | `15_diagnostics_summary.csv`; `15_tests_items.csv` | `n_T2_reject`/`n_T2_tests`, `share_T2_p05`, family `B_itemnonresp`; `T2_p_boot < .05` by variable |
| extreme-category fall of .21–.25; midpoint rise of .19–.21 | `15_designs_items.csv` | `d_std` for `pdq_ext_share`, `pdq_mid_share` |
| employment: continuing respondents more often employed by 4.9 points, flagged by the entry-wave correction (q = .02) and survival matching (q = .05), not standardised (q = .15) or naive (q = .38) | `15_designs_items.csv`; `04_item_meta.csv` | `DQ02`, `estimate`, `q` and `class3` for `naive`, `sm`, `ec`, `ec_adj` (BH families: the items that enter the counts, per estimator). `DQ02` keeps its raw codes, 1 = working and 2 = not working, so the estimate −.049 is a lower not-working share among continuing respondents; the checker confirms the coding from the reach rate of the follow-up `DQ02_1`, which only those not working are asked |
| mass-domination diagnostic: items with two to nine observed categories; items and categories with a positive/zero category; finite exceedances; exceedances in categories with at least ten fresh entrants; exceedances only in sparser categories; the two maxima | `15_mass_diagnostic_summary.csv` | `n_items_eligible`, `n_zero_denom_items`, `n_zero_denom_categories`, `n_finite_gt1`, `n_supported_gt1`, `n_sparse_only_gt1`, `max_finite_ratio`, `max_supported_ratio` (exact dose) |
| the flagged items themselves: the largest finite ratio (spouse's occupation, `dq44_2l`, 2.256); party identification (`DQ30`, 1.955, sparse categories only); the two supported exceedances (`DQ45A`, 1.192; `DQ08B_4`, 1.020); the two positive/zero items with a routing flag | `15_mass_diagnostic_flags.csv` | one row per item: `funnel_p` (retention for the item's own population: survivors who answered, as a share of the cohort, divided by the share of fresh entrants who reached the question), `funnel_ratio`, `funnel_ratio_supported`, `funnel_zero_denom`, `funnel_sparse_gt1` |
| 23 negative-control items; median \|d\| .033; 1 rejection; 9 equivalences; pooled −.035 (independence SE .0095) | `11_negcontrol_w13_summary.csv` | all columns |
| person-level bootstrap SE .0186 (2.0 times the independence SE; design effect 3.84); 95% interval [−.072, .001]; −.028 (SE .0185) without the rejected item | `11b_negcontrol_w13_pooled.csv` | `se_person_bootstrap`, `design_effect`, `ci95_lo`, `ci95_hi`; the row with the BH-rejected item dropped |
| mean correlation .116 between item contrasts across replications | `11b_negcontrol_w13_corr.csv` | `mean_offdiag_corr` |
| implied bias .017–.052 over Γ ∈ [.5, 1.5]; .036–.108 (.072Γ) at the far end of the interval | `11b_loading_grid.csv` | `correction`, `ci95_lo` (sign reversed) |
| the 23 item contrasts of the appendix table | `11_negcontrol_w13.csv` | `d`, `se`, `q` |
| n = 2,638 and n = 619 at the 2019 episode | `11b_negcontrol_w13_pooled.csv` | `n_2007_t13`, `n_2011_t9` |
| the input file behind all of them | `15_env.txt`, `11_env.txt`, `11b_env.txt` | `input`, `input sha256` (must agree with §5) |

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
transcribed into the manuscript with the outputs. The run records `15_env.txt`, `11_env.txt` and
`11b_env.txt` carry the input identity and the environment. The input behind Section 10 is the integrated
file `ZQ115AQ212BQ116CQ111DQ211EQ115FQ108GQ112HQ112IQ107JQ109KQ105LQ106MQ107NQ304OQ103PQ203QQ201RQ102.dta` (115,702,847 bytes; SHA-256 `d7fea333353e78bacde801045ae433ee3f24ae8d6b3e5af48d5720fec55306c5`); the checker refuses aggregate outputs
whose run records name any other file. The run behind the earlier version 0.5 of the paper read two files: the design script read
this one, and the preprocessing and battery scripts read the preceding release `…RQ101.dta` of the same
integrated file (a first-file-in-the-directory rule, now removed). A local, value-free comparison of the two
releases (same persons in the same order; every variable compared, names matched without regard to letter
case) found no difference in any value, including the 80 sampling-design variables whose names differ only
in letter case (among them the cohort identifier `CN`). The releases differ only in that letter case, in the
value labels (not the values or variable labels) of 13 variables, none of which this pipeline reads (three
and two variables of waves 2 and 6, whose only use here is the response markers, and eight of wave 19), and in
four variables added at wave 19. That record is supplied privately on request; it is the author's own
frozen re-execution and is **not** third-party reproduction, and the paper says so.

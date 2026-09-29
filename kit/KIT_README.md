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
Rscript R/15_test_style_items.R        # unit test of the response-style helpers (no data read)
# R/15_style_items.csv lists the rating-scale items of the response-style composites (read by 04 and 15; see §3b)
# R/15_item_scale.csv declares the scale, special codes, recodes, indicators, derived items and entry-wave class of
#   every wave-5 item; R/15_entry_overrides.csv and R/15_entry_subgroups.csv the 2007 counterparts and subgroups (§3c)
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
export P1_ITEM_SCALE=$P1_DATA_DIR/item_scale_synthetic.csv       # the synthetic items' rows of the item specification table (§3c)
export P1_STYLE_ITEMS=$PWD/synthetic/style_items_synthetic.csv   # the synthetic items' entry in the style-item list (§3b)
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

About a minute end to end. It exercises the principal computational paths of the 2011-episode
workflow, including the entry-wave correction (78 of the 89 universe columns have an entry-wave
counterpart, 77 of them in the main entry-wave analysis) and, since v1.0, every rule of the item specification
table (§3c): the merged categories of `DQ15`, the 0/1 indicators of the nominal questions `DQ30`, `DQ43`, `DQ57B`,
`DQ57C`, `DQ03_1` and `DQ03_2`, the merged party category of `DQ30` (code 7 names a different party at the two waves
and is folded into "other" at both), the routing codes of `DQ03_1` set to not applicable and its working students
returned to the questionnaire's category 3, the don't-know indicator `DQ42_G` that blanks its siblings, the filter of
`DQ03_5ES`, the scale-external code of `DQ17A`, the four clock times and the relationship duration built from their
components at both waves, the 2007 subgroups of the class-C items `DQ03_1`, `DQ03_2`, `DQ03_13` and `dq54_1_months`,
the entry-wave recodes of `DQ57A`, `DQ57D` and `DQ03_2` (the 2007 after-code 10 = agriculture folded into "other"), and
one class-D item (`DQ103`) whose entry-wave estimators stay outside the main counts. It does not contain every support pattern: the positive/zero case, the
missing-mass allocation and the other edge cases of the mass diagnostic are checked by `R/15_test_mass_ratio.R`,
and the synthetic input has no routing-flagged item. The synthetic run is a
**pipeline check only**: its numbers are not the paper's and must never be reported as such — its
analysis arms are 2,682 and 963 where the paper's are 2,797 and 963.

The real run additionally reads the provider's variable list, `varmap_w1-19.csv`, which maps each
item to its column in every wave. That file is the data provider's documentation and is **not**
shipped here; a licensed user already has it with the data, and points `P1_VARMAP` at it.
`R/15_make_synthetic_varmap.R` builds the equivalent for the synthetic input, together with the synthetic rows of the
item specification table, so the public pipeline needs nothing from outside this archive.

## 3. Where each number in Section 10 comes from

*The table below describes the run behind manuscript v1.0 (the licensed run of 29 September 2026 with the item
specification table of §3c). The run behind versions 0.5 to 0.10 (24 September 2026) used the inferred universe of 478
items; its counts are no longer quoted.*

`check_manuscript_values_T2.R` performs every comparison below and prints the result. It requires all
twenty-two aggregate files and run records listed at its top (a missing file stops it), requires every key
to identify exactly one row, first checks that the run records of scripts 15, 11 and 11b name the
same input file with the SHA-256 given in §5, and checks that the run record of 15 carries the MD5 of the item
specification table shipped in `R/`. On 29 September 2026 (JST) it reproduced all 147 quantities from the aggregate outputs of the licensed run, and its self-test detected all thirteen corruptions.

| Quantity in Section 10 | Output file | Column |
|---|---|---|
| 2,797 continuing survivors; 963 entrants; 574 and 540 survivors | `15_arms.csv` | `n_old_S`, `n_new_total`, `n_new_sm`, `n_new_sm1` |
| 500 bootstrap replications | `15_arms.csv` | `B` |
| the item universe: 540 wave-5 variables in the table, 53 excluded, 15 nominal questions, 6 derived items; 523 columns = 284 binary + 102 ordered + 41 continuous items + 90 indicators + 6 derived items; no variable outside the table | `R/15_item_scale.csv` (its MD5 in `15_env.txt`); `04_item_meta.csv`; `04_items_unlisted.csv` | `scale`, `construct`; `in_universe`, `spec_scale`, `construct`, `qgroup`; the unlisted file is empty |
| 320 columns with an entry-wave counterpart; 269 in the main entry-wave analysis (classes A, B, C); 16 of the 320 routing-flagged and 35 class D; in the counts 181 (A + B) + 87 (C) = 268, the 35 class-D columns outside them; admitting class D raises the two entry-wave counts by two each (29 → 31, 23 → 25); the 2007 subgroups of class C | `15_arms.csv`; `15_ec_status.csv`; `15_ec_counts_by_class.csv`; `15_ec_scope_sensitivity.csv` | `n_ec_items`, `n_ec_items_main`; `ec_ok`, `ec_class`, `ec_main`, `ec_subgroup`; `n_items` by `cls`; `n_affected` by `scope` |
| 33 columns outside every count: 31 flagged as potentially incomparable because of differential item nonresponse or routing (04's flag: item-nonresponse gap above 15 points), one follow-up of a flagged question (DQ46Y), and one sensitivity variant (`dq57d_hr_s12`, `in_counts = FALSE`). The flag is a screening rule; the questionnaire filters have not been verified | `15_designs_items.csv` | `exclude_flag` (= not `in_counts`), `nr_routing_flag`, `routing_followup`, `count_exclude` (= any of the three, or, for `ec` and `ec_adj`, outside the main entry-wave analysis) |
| the same counts with the flag recomputed at gap thresholds of 10 and 20 points (identical flagged set and counts) and with no routing exclusion (522 columns; 31 / 21 / 14 / 29 / 23; the same three items flagged by all five designs) | `15_routing_sensitivity.csv` | one row per threshold: `n_routing`, `n_followup`, `n_items_*`, `aff_*`, `common_*`, `n_all5`, `all5_items`; `n_disagree_with_04` (= 0: the gap recomputed in 15 at .15 reproduces 04's flag) |
| 490 columns for which a mean contrast is meaningful; 489 for the matched designs (one item without variation in the matched arms); 268 in the entry-wave counts | `15_detection_counts.csv` | `n_items` by `estimator` |
| flagged 30 / 21 / 13 / 29 / 23 | `15_detection_counts.csv` | `n_affected` by `estimator` |
| flagged 18 / 16 / 9 / 29 / 23 on the common set of 268 | `15_designs_items.csv` | `class3 == "affected"`, restricted to columns with an `ec` row that enter the counts |
| three items flagged by all five designs (subjective social position, owner-occupancy of a detached house, anxiety about married life as a reason for remaining single) | `15_designs_items.csv` | intersection over `estimator` (DQ26, DQ39__1, DQ55_Q) |
| diagnostic rejections: 39 of 268 (14.6%) and 43 of 489 (8.8%) | `15_diagnostics_summary.csv` | `n_T1_reject`/`n_T1_tests`, `n_T2_reject`/`n_T2_tests`, `share_T1_p05`, `share_T2_p05`, family `A_substantive` |
| 73 of 409 (17.8%) for item nonresponse; 50 of the 73 from the grids DQ58C, DQ04(3), DQ09 and DQ08D | `15_diagnostics_summary.csv`; `15_tests_items.csv` | `n_T2_reject`/`n_T2_tests`, `share_T2_p05`, family `B_itemnonresp`; `T2_p_boot < .05` by variable |
| response-style composites: the item set (46 rating-scale items on the list; verified codes; the common battery of 42 items with the same codes at entry, 25 with a neutral midpoint) and the 40 items that the earlier rule admits on the present universe (26 frequency scales, occupational rank, smoking, drinking, education and other classifications; no nominal code, those being indicators now) | `15_style_items_audit.csv` | `set`, `codes_verified`, `entry_wave_ok`, `in_main`, `in_main_midpoint`, `in_v07_rule`, `reason`; the codes and their labels as read from the file (§3b) |
| extreme-category and midpoint margins over the five designs (common battery) | `15_designs_items.csv` | `d_std` for `pdq_ext_share`, `pdq_mid_share`, family `P_style` |
| the same margins under other item sets: all rating items (same-wave designs), rating and frequency scales, the bipolar five-point scales only, and the rule used up to v0.7 (every item with four to seven consecutively numbered codes) as run on the present universe and on the common battery | `15_style_sensitivity.csv` | `battery`, `indicator`, `estimator`, `n_items_w5`, `n_items_w1`, `d_std`; the `rating_common (main)` rows repeat the item-file values; the `v07_asrun` rows are the earlier rule applied to the v1.0 universe (a fall of .21 to .26 and a rise of .20 to .22; on the v0.9 universe the rule gave .21 to .25 and .19 to .21) |
| employment: continuing respondents more often employed by 4.9 points, flagged by the entry-wave correction (q = .03) and survival matching (q = .05), not standardised (q = .14) or naive (q = .38) | `15_designs_items.csv`; `04_item_meta.csv` | `DQ02`, `estimate`, `q` and `class3` for `naive`, `sm`, `ec`, `ec_adj` (BH families: the columns that enter the counts, per estimator). `DQ02` keeps its raw codes, 1 = working and 2 = not working, so the estimate −.049 is a lower not-working share among continuing respondents; the checker confirms the coding from the reach rate of the follow-up `DQ02_2`, which only those not working are asked |
| mass-domination diagnostic: items with two to nine observed categories; items and categories with a positive/zero category; finite exceedances; exceedances in categories with at least ten fresh entrants; exceedances only in sparser categories; the two maxima | `15_mass_diagnostic_summary.csv` | `n_items_eligible`, `n_zero_denom_items`, `n_zero_denom_categories`, `n_finite_gt1`, `n_supported_gt1`, `n_sparse_only_gt1`, `max_finite_ratio`, `max_supported_ratio` (exact dose) |
| fresh item nonresponse: all 18 flags (10 positive/zero columns, one category each, 9 of them asked only of a subgroup; 8 finite exceedances) are reconcilable by the fresh missing mass; median missing mass .018 over the 410 eligible columns; largest needed mass .0065; reach rates differ by 9.1 points at most (the owner-occupied-housing follow-ups DQ39_A–E), the unmarried block (DQ50 and its follow-ups) next | `15_mass_diagnostic_summary.csv`; `15_mass_diagnostic_flags.csv`; `15_tests_items.csv` | `n_flagged_not_reconcilable`, `n_finite_gt1_not_reconcilable`, `n_zero_denom_not_reconcilable`, `n_supported_gt1_not_reconcilable`, `median_fresh_missing_mass` |
| the flagged columns themselves: the two largest finite ratios, both rare events exceeding one only in sparse categories (mother died in the past year, `DQ09_D`, 1.369; expects to have taken over the family business in ten years, `DQ56_C`, 1.270); the two supported exceedances (`DQ45A`, 1.192; `DQ08B_4`, 1.020); the six positive/zero columns on how a respondent with a partner met that person (`DQ54_2*`) and the two with a routing flag (`DQ49_2P`, `DQ49_2Z`); per column, the fresh missing mass `r_M`, the mass needed to cover the stayers under the identity map, and whether the map remains feasible; the reach and answer shares of both arms | `15_mass_diagnostic_flags.csv` | one row per item: `funnel_p` (retention for the item's own population: survivors who answered, as a share of the cohort, divided by the share of fresh entrants who reached the question), `funnel_ratio`, `funnel_ratio_supported`, `funnel_zero_denom`, `funnel_sparse_gt1` |
| 23 negative-control items; median \|d\| .033; 1 rejection; 9 equivalences; pooled −.035 (independence SE .0095) | `11_negcontrol_w13_summary.csv` | all columns |
| person-level bootstrap SE .0186 (2.0 times the independence SE; design effect 3.84); 95% interval [−.072, .001]; −.028 (SE .0185) without the rejected item | `11b_negcontrol_w13_pooled.csv` | `se_person_bootstrap`, `design_effect`, `ci95_lo`, `ci95_hi`; the row with the BH-rejected item dropped |
| mean correlation .116 between item contrasts across replications | `11b_negcontrol_w13_corr.csv` | `mean_offdiag_corr` |
| implied bias .017–.052 over Γ ∈ [.5, 1.5]; .036–.108 (.072Γ) at the far end of the interval | `11b_loading_grid.csv` | `correction`, `ci95_lo` (sign reversed) |
| the 23 item contrasts of the appendix table | `11_negcontrol_w13.csv` | `d`, `se`, `q` |
| n = 2,638 and n = 619 at the 2019 episode | `11b_negcontrol_w13_pooled.csv` | `n_2007_t13`, `n_2011_t9` |
| the input file behind all of them | `15_env.txt`, `11_env.txt`, `11b_env.txt` | `input`, `input sha256` (must agree with §5) |

### 3a. The mass diagnostic and fresh item nonresponse

The mass ratio of Section 10 is a complete-case quantity. `15_panelcond_designs.R` computes, per item, `p` as the
survivors who gave a substantive answer, per member of the cohort, divided by the share of fresh entrants who reached
the item (their reach rate stands in for eligibility; an item nonresponse recorded as missing rather than with a
no-answer code cannot be told from routing, see the comment in `04_build_analysis.R`); `Q` as the distribution of
the stayers' substantive answers; and `F2` as the distribution of the fresh entrants who gave one. Reading the ratio
as the population inequality of Corollary 4 requires, beyond refreshment validity, that fresh answerers be
representative of the eligible subgroup. The codes are treated as in `04_build_analysis.R`: a not-applicable code marks a respondent who did not reach the
item; DK, refusal and no-answer codes are reached but non-substantive; every other value is a substantive answer.
Survivors who reached the item without a substantive answer are folded into `p` and treated like attriters. Since the
fresh missing mass (`funnel_missing_mass` = fresh entrants who reached the item without a substantive answer, as a
share of those who reached it) could sit anywhere, `mass_ratio_item()` also evaluates the
identity map with that mass allocated freely: `funnel_needed_mass` = Σ_a max{l_a − r_a, 0} with `l_a` = p·Q(a) and
`r_a` = (1 − r_M)·F2(a); `funnel_identity_feasible` is `needed ≤ missing`. The reach and answer shares of both arms
(`reach_new`, `reach_old_S`, `answer_new_given_reach`, `answer_old_S_given_reach`) are written for every flagged item
and for every item in `15_tests_items.csv`. These are population compatibility calculations on complete-case
shares, not calibrated tests; the paper reports the ratios as descriptive diagnostics.

### 3b. The response-style composites and their item set

The two person-level composites of Section 10 are the share of a respondent's substantive answers that fall in an
extreme category and the share that fall in the neutral middle category. Since 24 September 2026 they are computed
over a **fixed list of rating-scale items**, `R/15_style_items.csv` (`set = rating`: agreement, satisfaction
and evaluation scales with verbal anchors; `set = frequency`: frequency scales, used only in the sensitivity file).
The list was written before the licensed run and revised on the value labels that the run's audit file printed
(the neutral middle of the two standard-of-living items; attention to politics and advice from co-workers moved to
or added to the frequency set; the frequency of meeting a partner removed, its first code being a status), and the
run was repeated with the revised list; the list's MD5 in `15_env.txt` identifies the version behind Section 10. The list declares each item's number of substantive codes `k` and whether its middle code is a
neutral category. `R/15_style_items.R` verifies every listed item against the value labels in the file: the
substantive codes (all labels minus the don't-know, refusal, no-answer and not-applicable codes, and the
`NAP_OVERRIDE` codes of `00_config.R`) must be exactly `1..k`, the observed values must lie in `1..k` (at the entry wave
too, for the common battery), and a declared midpoint must be an odd `k` whose middle label names a neutral category
("どちらともいえない", "ふつう", "変わらない", ...). An item that fails the code check enters no composite; one that fails only the
midpoint check enters the extreme-category composite only. The list is read strictly (a malformed line stops the run, and so does a `k` that is not a whole number, such as 5.9,
which is not read as 5; every such stop names the item),
and its MD5 is recorded by `04` in the derived file and by `15` in `15_env.txt`; `15` stops if `04` used another list. `04_build_analysis.R`
computes the composites over all verified rating items (the `rating_all` set); `15_panelcond_designs.R` recomputes
them, checks that it reproduces `04`, and uses as the **main composites** the *common battery*: the verified rating
items whose entry-wave counterpart has the same codes, so that the entry-wave and the comparison-wave composites,
and all five designs, use one item set. `15_style_items_audit.csv` records, for every listed item, its codes and
labels as read from the file, whether it was verified, its entry-wave counterpart, and which composites it enters;
it also lists the items that the earlier rule admitted but that are not on the list. The unit test
`R/15_test_style_items.R` (36 checks, no data read) covers the reading of the list, including the malformed lists
that must stop the run, the verification of codes and midpoint labels, the range check of observed values, the
person-level shares, and the exclusion reason written to the audit file when zero, one or several items have
entry-wave values outside `1..k`. (Up to version 0.8, two or more such items stopped `15` with an error while that
reason was written; the licensed run and the synthetic input have no such item, so no output changed.)

Up to manuscript v0.7 the composites were built from every universe item with four to seven consecutively numbered
substantive codes, whatever its meaning, which admitted marital status, occupational rank, education and other
classifications. That rule is retained only as the `v07_asrun` rows of `15_style_sensitivity.csv` (the w5 composite
over all such items paired with the w1 composite over those with an entry-wave counterpart, as it was run; since v1.0
the rule is applied to the universe of §3c, on which the nominal codes are indicators with two codes and no longer
fall under it) and the `v07_common` rows (the same rule on the common battery), next to `rating_all` (same-wave designs), `rating_common`
(= the main rows), `agree_common` (the five-point agreement items only, which are asked of every respondent;
the job-characteristics items and the job and marriage satisfaction items reach only respondents with a job or a
spouse), `bipolar_common` (the five-point agreement and satisfaction scales) and `ratingfreq_all` /
`ratingfreq_common` (frequency scales added). Each row gives the number of items in the
composite at each wave, the estimate, its bootstrap SE (the same joint person bootstrap as every other quantity)
and the standardised margin `d_std`.

### 3c. The item specification table and the entry-wave classes (v1.0)

Up to version 0.9, `04_build_analysis.R` inferred each item's scale from the number of its value labels (two:
binary; three to nine: ordinal; otherwise continuous) and recognised special codes from label strings only. An audit
of the 2011 questionnaire (`review_round06/SCALE_AUDIT_2011_ja.md`) found nominal items averaged as if ordered, codes
outside the scale entering the means ("public sector", "other", "no parent at the time", "not fixed"), don't-know
options that had become items of their own, and clock and duration variables whose distributed recodes differ between
waves. Since version 1.0 the treatment of every wave-5 item is declared in **`R/15_item_scale.csv`** and read by `04`
and `15` (`R/15_item_scale.R`): its scale (`binary`, `ordinal`, `continuous`, `nominal`, `exclude`), the substantive
codes it must show (`codes_expected`, checked against the file; a mismatch stops the run), the special codes beyond the
label rule (`nap_add`, `dk_add`, `nr_add`, `ref_add`), a recode applied at both waves (`recode`; e.g. the two "currently
zero" categories of smoking and drinking merged), the construction of 0/1 indicators for nominal single-choice
questions (`construct = indicators`, one column `<var>__<code>` per category; the question's item-nonresponse and
don't-know indicators are attached to its first category), derived items built from components by one rule at both
waves (`clock:` AM/PM, hour and minute to decimal hours — return-home times earlier than the leave-home time and
bed times before noon are moved past midnight, 12 pm is read as midnight; `months:` years and months to months), a
filter (`filter`; e.g. the annual-salary amount kept for those who chose that pay form, the health-insurance
categories blanked for those who answered "don't know"), and whether the item enters the counts (`in_counts`). A
D-series variable not listed in the table is left out and named in `04_items_unlisted.csv`; `04_item_scale_resolved.csv`
records the codes, labels and special codes of every universe item as resolved, for review and for freezing into
`codes_expected`. The table's MD5 is recorded by `04` in the derived file and by `15` in `15_env.txt`, and `15` stops if
`04` used another table.

The same table carries the **entry-wave comparability class** of each item (`ec_class`, from
`review_round06/EC_ITEMS_REVIEW_ja.md`, which compares the 2007 and 2011 questionnaires item by item): A (identical
wording, options, codes and respondents), B (minor differences of layout or context), C (identical wording but a
different filter or reference, comparable within the same subgroup), D (different reference period, option set,
format or wording). The **main entry-wave analysis counts classes A, B and C** (`ec_main`); the entry-wave estimators
of class-D items and of items without a class are computed and written but stay outside the counts (`class3 =
"excluded (entry-wave class D)"`), and `15_ec_scope_sensitivity.csv` gives the counts with and without them
(`15_ec_counts_by_class.csv` by class). For class-C items the entry-wave term is computed within the 2007 subgroup
named in `ec_subgroup` and defined in **`R/15_entry_subgroups.csv`** (workers, employees, the married, respondents
with a partner, parents at 2007), so that the restriction B4 of the paper is a restriction on that subgroup.
**`R/15_entry_overrides.csv`** names the 2007 counterpart where the provider's map has none or the codes differ
(a recode of the 2007 codes, or the 2007 components of a derived item), and the same recodes and special codes of
the table are applied to the 2007 variable before its label set is compared with the 2011 one; `15_ec_status.csv`
lists every universe item with its counterpart, class and subgroup, and `15_ec_unavailable.csv` the reason where no
entry-wave term exists. The former `R/15_item_exclude.csv` is superseded by `in_counts`.

The label sets of the 2007 and 2011 sides were compared code by code on 29 September 2026 against the provider's
value-label file (metadata only) and both questionnaires, for all 254 entry-wave pairs of the licensed run. Three
substantive differences were found and are handled in the tables: party identification (`DQ30`) has a different party
at code 7 in the two years (2011 みんなの党, 2007 新党日本), so the table folds code 7 into 8 = other at both waves
(`recode = 7=8`; 8 indicators); the 2007 occupation precodes (`JC_2`, `ZQ55_2`) carry an after-code 10 = agriculture
that the 2011 list does not have, folded into 8 = other by `15_entry_overrides.csv` (`map:10=8`); and work status
(`DQ03_1`, `JC_1`) carries the provider's routing codes 10 and 11 (not working; not answer options, and already
represented by `DQ02` and `DQ02_1`) and 12 (students working non-regularly, split out of option 3), which the table
sets to not applicable and returns to 3 respectively, so that the indicators are the questionnaire's nine options over
current workers at both waves. All other differences are wording only (e.g. 家族従業者 / 家族従事者, ○選択 / 選択).

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
frozen re-execution and is **not** third-party reproduction, and the paper says so. The run behind Section 10 of
versions 0.5 to 0.10 used the versions of `R/15_panelcond_designs.R` and `R/15_style_items.R` in tag `paper-v0.8`. Version 0.9 changes them
only in how an exclusion reason is appended to the audit table and in how the list is read: a `k` that is not a
whole number is rejected, a reading problem stops the run after the reader has returned rather than inside it, and
the other checks of the list stop with a message that names the item (§3b); the licensed run has no item with entry-wave values outside `1..k` and its list has
whole-number `k`, so neither change alters any output of that run, and the outputs on the synthetic input are
byte-identical to those of version 0.8. **Version 1.0 (kit of 28--29 September 2026) changes the item universe itself**
(§3c: the item specification table, the nominal indicators, the derived items, the filters and the entry-wave
classes). The licensed run was repeated on 29 September 2026 (JST) with these scripts (`04_build_analysis.R` v6.1,
`15_panelcond_designs.R` v1.0, `15b_item_flags.R`; R 4.6.0 on macOS, `panelcond` 0.1.2, `data.table` 1.18.4, `haven`
2.5.5; 500 bootstrap replications, seed 20260915; about 15 minutes) on the same input file, with the item
specification table of MD5 `5b056b5a425c4d471cda9e0f8e8b8f13` and the style-item list of MD5
`19c5597ee76037603652da37de5ce037`, both recorded in `15_env.txt`, and its aggregate outputs are the run behind
Section 10 of manuscript v1.0 and tag `paper-v1.0`; `check_manuscript_values_T2.R` reproduced all 147 quoted
quantities from them. The 2019-episode battery (scripts 11 and 11b) was not repeated: it does not read the item
specification table, and its run records name the same input. The frozen run of 24 September 2026 (478 items, the run
behind versions 0.5 to 0.10) is superseded; its outputs are kept by the author but are no longer quoted. The synthetic
input was regenerated with the items that exercise the new rules (its SHA-256 in the archive README changes
accordingly; regenerated again on 29 September 2026 with `DQ03_1`/`JC_1`, `DQ03_2`/`JC_2` and the nine-party coding
of `DQ30`/`ZQ30`).

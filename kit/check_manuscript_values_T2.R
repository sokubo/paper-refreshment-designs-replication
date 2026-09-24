#!/usr/bin/env Rscript
# ============================================================
# T2 — every number quoted in Section 10 and in the negative-control appendix table, reproduced
# mechanically from the aggregate outputs of the JLPS pipeline (the counterpart of the simulation
# checker sims/check_manuscript_values.R). Inputs are aggregates only (cells below ten suppressed);
# no individual record is read.
#   Rscript check_manuscript_values_T2.R [<results directory>]            check the quoted values
#   Rscript check_manuscript_values_T2.R [<results directory>] --selftest check, then check the checker
#     default results directory: ../../P1_jlps_diagnosis/results
# Every file in REQUIRED must be present; a missing file stops the run (nothing is skipped). A
# comparison fails if a value is missing, non-finite, of the wrong length, or taken from a key that
# does not identify exactly one row. The three run records (15_env.txt, 11_env.txt, 11b_env.txt) must
# name the same input file with the same SHA-256, equal to the one declared in INPUT below.
# ============================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
SELFTEST <- "--selftest" %in% args; args <- setdiff(args, "--selftest")
RES0 <- if (length(args)) args[1] else file.path("..", "..", "P1_jlps_diagnosis", "results")

## the licensed input all quoted numbers come from (identity only, never content)
INPUT <- list(file   = "ZQ115AQ212BQ116CQ111DQ211EQ115FQ108GQ112HQ112IQ107JQ109KQ105LQ106MQ107NQ304OQ103PQ203QQ201RQ102.dta",
              sha256 = "d7fea333353e78bacde801045ae433ee3f24ae8d6b3e5af48d5720fec55306c5")

REQUIRED <- c("15_arms.csv", "15_detection_counts.csv", "15_diagnostics_summary.csv", "15_designs_items.csv",
              "15_mass_diagnostic_summary.csv", "15_mass_diagnostic_flags.csv", "15_tests_items.csv", "15_routing_sensitivity.csv", "15_env.txt", "04_item_meta.csv",
              "11_negcontrol_w13_summary.csv", "11_negcontrol_w13.csv", "11_env.txt",
              "11b_negcontrol_w13_pooled.csv", "11b_negcontrol_w13_corr.csv", "11b_loading_grid.csv", "11b_env.txt")

## the comparison rule: NULL when the value reproduces the quoted one, otherwise the reason it does not
compare <- function(got, want, tol = 0) {
  if (!length(want) || anyNA(want)) return("no expected value declared")
  if (length(got) != length(want)) return(sprintf("%d value(s), expected %d", length(got), length(want)))
  if (is.character(want)) return(if (anyNA(got) || !identical(as.character(got), want)) "differs" else NULL)
  g <- suppressWarnings(as.numeric(got))
  if (any(!is.finite(g))) return("missing or non-finite value")
  if (any(abs(g - as.numeric(want)) > tol)) "differs" else NULL
}

run_checks <- function(RES, verbose = TRUE) {
  miss <- REQUIRED[!file.exists(file.path(RES, REQUIRED))]
  if (length(miss)) stop("missing aggregate outputs in ", RES, ": ", paste(miss, collapse = ", "), call. = FALSE)
  rd <- function(f) fread(file.path(RES, f), encoding = "UTF-8")
  ok <- 0L; bad <- 0L; failed <- character(0)
  say <- function(...) if (verbose) cat(...)
  record <- function(what, got, want, why) {
    if (is.null(why)) ok <<- ok + 1L else { bad <<- bad + 1L; failed <<- c(failed, paste(what, why, sep = " :: ")) }
    say(sprintf("%-58s %-22s %-22s %s\n", what, paste(format(got), collapse = " "),
                paste(format(want), collapse = " "), if (is.null(why)) "ok" else paste("MISMATCH:", why)))
  }
  chk <- function(what, got, want, tol = 0) record(what, got, want, compare(got, want, tol))
  ## exactly one row for a key; otherwise the failure is recorded and a one-row NA table returned
  one <- function(dt, what) {
    if (nrow(dt) == 1L) return(dt)
    record(paste("key identifies one row:", what), nrow(dt), 1L, "key does not identify exactly one row")
    dt[NA_integer_]
  }
  uniq <- function(dt, keys, what) {
    d <- dt[, .N, by = keys][N > 1L, .N]
    record(paste("unique keys:", what), d, 0L, if (d == 0L) NULL else "duplicated key rows")
  }
  env_input <- function(f) {
    x <- readLines(file.path(RES, f), warn = FALSE)
    g <- function(tag) { v <- sub(paste0("^", tag, ":\\s*"), "", grep(paste0("^", tag, ":"), x, value = TRUE)); if (length(v) == 1L) v else NA_character_ }
    c(file = g("input"), sha256 = g("input sha256"))
  }

  say(sprintf("%-58s %-22s %-22s %s\n", "quantity (Section 10)", "from outputs", "quoted in paper", ""))
  say(strrep("-", 112), "\n")

  ## --- the input: one file, the same for every script -------------------------------------------
  e15 <- env_input("15_env.txt"); e11 <- env_input("11_env.txt"); e11b <- env_input("11b_env.txt")
  shas <- c(e15[["sha256"]], e11[["sha256"]], e11b[["sha256"]])
  chk("15, 11 and 11b record one and the same input SHA-256", as.numeric(!anyNA(shas) && length(unique(shas)) == 1L), 1)
  chk("input file named by 15 / 11 / 11b", c(e15[["file"]], e11[["file"]], e11b[["file"]]), rep(INPUT$file, 3))
  chk("input SHA-256 recorded by 15 / 11 / 11b", c(e15[["sha256"]], e11[["sha256"]], e11b[["sha256"]]), rep(INPUT$sha256, 3))

  arms <- one(rd("15_arms.csv")[dose_def == "exact"], "15_arms exact")
  dc   <- rd("15_detection_counts.csv")[dose_def == "exact" & family == "A_substantive"]
  dg   <- rd("15_diagnostics_summary.csv")[dose_def == "exact"]
  it   <- rd("15_designs_items.csv")[dose_def == "exact"]
  if (!"count_exclude" %in% names(it)) stop("15_designs_items.csv has no count_exclude column: rerun R/15_panelcond_designs.R with this kit", call. = FALSE)
  ms   <- one(rd("15_mass_diagnostic_summary.csv")[dose_def == "exact"], "15_mass_diagnostic_summary exact")
  nc   <- one(rd("11_negcontrol_w13_summary.csv"), "11_negcontrol_w13_summary")
  uniq(dc, "estimator", "15_detection_counts (exact, A_substantive)")
  uniq(dg, "family", "15_diagnostics_summary (exact)")
  uniq(it, c("family", "var", "estimator"), "15_designs_items (exact)")

  ## --- arms and bootstrap -----------------------------------------------------
  chk("continuing survivors of five waves", arms$n_old_S, 2797)
  chk("fresh entrants at the 2011 refreshment", arms$n_new_total, 963)
  chk("entrants surviving four further waves", arms$n_new_sm, 574)
  chk("entrants surviving five further waves", arms$n_new_sm1, 540)
  chk("bootstrap replications", arms$B, 500)
  chk("items with the same coding at entry", arms$n_ec_items, 268)

  ## --- waterfall over all items with variation --------------------------------
  g <- function(e, col) one(dc[estimator == e], paste("detection counts", e))[[col]]
  chk("substantive items in the counts (naive arm)", g("naive", "n_items"), COUNTS$n_naive)
  chk("flagged: naive / survival matching / symmetric matching / entry-wave / standardised",
      c(g("naive", "n_affected"), g("sm", "n_affected"), g("ssm", "n_affected"), g("ec", "n_affected"), g("ec_adj", "n_affected")), COUNTS$flagged)
  chk("items in the counts: survival / symmetric matching", c(g("sm", "n_items"), g("ssm", "n_items")), COUNTS$n_matched)
  chk("items with an available entry wave", g("ec", "n_items"), 227)

  ## --- restricted to the common set of 227 ------------------------------------
  A <- it[family == "A_substantive"]
  if (!"routing_followup" %in% names(A)) stop("15_designs_items.csv has no routing_followup column: rerun R/15_panelcond_designs.R with this kit", call. = FALSE)
  u <- unique(A[, .(var, exclude_flag, count_exclude, routing_followup)])
  chk("items in the inventory", nrow(u), 478)
  chk("outside the counts: listed codes / routing only", c(sum(u$exclude_flag), sum(u$count_exclude & !u$exclude_flag)), COUNTS$excluded)
  rtf <- unique(A[, .(var, nr_routing_flag, exclude_flag)])
  chk("routing-flagged items that are also on the exclusion list", sum(rtf$nr_routing_flag %in% TRUE & rtf$exclude_flag), 3)
  chk("follow-up of a routing-flagged question (outside the counts)", paste(sort(u[routing_followup & !exclude_flag, var]), collapse = ", "), COUNTS$followup)
  ## the detection counts cover items with variation that are outside neither 15_item_exclude.csv (nominal codes,
  ## date and clock-time components, duplicate recodes) nor the routing flag (count_exclude); the common set is
  ## defined the same way
  A <- A[class3 != "n/a (no variation)" & !count_exclude]
  common <- A[estimator == "ec", unique(var)]
  chk("items in the detection counts with an entry-wave counterpart", length(common), 227)
  cnt <- sapply(c("naive", "sm", "ssm", "ec", "ec_adj"),
                function(e) A[estimator == e & var %in% common & class3 == "affected", uniqueN(var)])
  chk("common set: naive / sm / ssm / ec / ec_adj", cnt, COUNTS$common)
  five <- Reduce(intersect, lapply(c("naive", "sm", "ssm", "ec", "ec_adj"),
                 function(e) A[estimator == e & class3 == "affected", unique(var)]))
  chk("items flagged by all five designs", length(five), COUNTS$n_all5)
  chk("the items flagged by all five designs", paste(sort(five), collapse = ", "), COUNTS$all5_items)
  say("      (the items: ", paste(sort(five), collapse = ", "), ")\n", sep = "")

  ## --- diagnostics ------------------------------------------------------------
  dA <- one(dg[family == "A_substantive"], "diagnostics A"); dB <- one(dg[family == "B_itemnonresp"], "diagnostics B")
  chk("T_NS: rejections / items with the statistic (entry wave needed)", c(dA$n_T1_reject, dA$n_T1_tests), DIAG$ns)
  chk("T_SD: rejections / items with the statistic", c(dA$n_T2_reject, dA$n_T2_tests), DIAG$sd)
  chk("T_SD on item nonresponse: rejections / indicators", c(dB$n_T2_reject, dB$n_T2_tests), DIAG$sd_nr)
  chk("T_NS rejections, substantive items (13.2%)", one(dg[family == "A_substantive"], "diagnostics A")$share_T1_p05, 0.132, 5e-4)
  chk("T_SD rejections, substantive items (9.0%)", one(dg[family == "A_substantive"], "diagnostics A")$share_T2_p05, 0.090, 5e-4)
  ## item nonresponse: the rejections that come from four question grids answered or skipped as a block
  tsd <- rd("15_tests_items.csv")[dose_def == "exact" & family == "B_itemnonresp" & count_exclude == FALSE]
  uniq(tsd, "var", "15_tests_items (exact, item nonresponse)")
  rj <- tsd[T2_p_boot < .05, toupper(var)]
  chk("item-nonresponse rejections from the grids DQ58C, DQ04(3), DQ09, DQ08D / all", c(sum(grepl("^(DQ58C_|DQ04_3|DQ09_|DQ08D_)", rj)), length(rj)), DIAG$grid4)
  chk("T_SD rejections, item nonresponse (share)", one(dg[family == "B_itemnonresp"], "diagnostics B")$share_T2_p05, DIAG$share_nr, 5e-4)

  ## --- response-style margins -------------------------------------------------
  st <- it[family == "P_style"]
  ext <- st[grepl("ext_share", var), d_std]; mid <- st[grepl("mid_share", var), d_std]
  chk("extreme-category use: five designs", length(ext), 5)
  chk("midpoint use: five designs", length(mid), 5)
  chk("extreme-category use, range over the five designs", round(range(ext), 2), c(-0.25, -0.21))
  chk("midpoint use, range over the five designs",        round(range(mid), 2), c(0.19, 0.21))

  ## --- the employment item ----------------------------------------------------
  emp <- A[var == "DQ02"]
  chk("employment, entry-wave corrected (percentage points)",
      round(100 * one(emp[estimator == "ec"], "DQ02 ec")$estimate, 1), -4.9, 0.05)
  chk("employment, unstandardised: q (.02) below .10",
      round(one(emp[estimator == "ec"], "DQ02 ec")$q, 2), 0.02, 5e-3)
  chk("employment, standardised: q above .10",
      as.numeric(one(emp[estimator == "ec_adj"], "DQ02 ec_adj")$q > 0.10), 1)
  chk("employment, standardised: q (.15)", round(one(emp[estimator == "ec_adj"], "DQ02 ec_adj")$q, 2), 0.15, 5e-3)
  chk("employment, naive: q (.38)", round(one(emp[estimator == "naive"], "DQ02 naive")$q, 2), 0.38, 5e-3)
  ## direction: DQ02 keeps its raw codes (1 = working, 2 = not working); its follow-up DQ02_1 is reached only by those
  ## not working, and the naive contrast over the wave-5 risk set (dose "any") equals the difference in that reach rate,
  ## so a negative estimate means that continuing respondents are MORE often employed
  im <- rd("04_item_meta.csv"); f1 <- one(im[toupper(var) == "DQ02_1"], "04_item_meta DQ02_1")
  nv_any <- one(rd("15_designs_items.csv")[dose_def == "any" & family == "A_substantive" & var == "DQ02" & estimator == "naive"], "DQ02 naive, dose any")$estimate
  chk("DQ02 coding: naive contrast (any dose) = difference in the not-working share", round(nv_any - (f1$base_rate_T1 - f1$base_rate_T0), 3), 0, 1.5e-3)
  chk("employment: the entry-wave-corrected contrast is negative (more often employed)", as.numeric(one(emp[estimator == "ec"], "DQ02 ec")$estimate < 0), 1)
  chk("employment, survival matching: q (.05), flagged", c(round(one(emp[estimator == "sm"], "DQ02 sm")$q, 2), as.numeric(one(emp[estimator == "sm"], "DQ02 sm")$class3 == "affected")), c(0.05, 1), 5e-3)

  ## --- mass-domination diagnostic (exact dose, substantive items with two to nine observed categories) ----
  chk("items with two to nine observed categories",          ms$n_items_eligible,      MASS$eligible)
  chk("items with a positive/zero category (infinite ratio)", ms$n_zero_denom_items,    MASS$zero_denom_items)
  chk("positive/zero categories in total",                   ms$n_zero_denom_categories, MASS$zero_denom_categories)
  chk("items whose finite ratio exceeds one",                 ms$n_finite_gt1,          MASS$finite_gt1)
  chk("largest finite ratio",                                 ms$max_finite_ratio,      MASS$max_finite, 5e-4)
  chk("items exceeding one in a supported category (>= 10)",  ms$n_supported_gt1,       MASS$supported_gt1)
  chk("items exceeding one only in sparse categories",        ms$n_sparse_only_gt1,     MASS$sparse_only_gt1)
  chk("largest ratio in a supported category",                ms$max_supported_ratio,   MASS$max_supported, 5e-4)
  fl <- rd("15_mass_diagnostic_flags.csv")[dose_def == "exact"]
  uniq(fl, "var", "15_mass_diagnostic_flags (exact)")
  chk("item with the largest finite ratio", one(fl[funnel_ratio == max(funnel_ratio, na.rm = TRUE)], "largest finite ratio")$var, MASS$max_finite_item)
  chk("items exceeding one in a supported category", sort(fl[funnel_ratio_supported > 1, var]), MASS$supported_items)
  pid <- one(fl[var == "DQ30"], "party identification (DQ30) among the flags")
  chk("party identification: finite ratio", pid$funnel_ratio, MASS$party_id, 5e-4)
  chk("siblings for help in finding work (DQ08B_4): supported ratio", one(fl[var == "DQ08B_4"], "DQ08B_4")$funnel_ratio_supported, 1.020, 5e-4)
  chk("spouse prepares meals (DQ45A): supported ratio", one(fl[var == "DQ45A"], "DQ45A")$funnel_ratio_supported, 1.192, 5e-4)
  chk("party identification exceeds one only in sparse categories", as.numeric(pid$funnel_ratio_supported <= 1), 1)
  ## positive/zero items in the marriage-history block flagged for differential item nonresponse or routing
  rt <- unique(it[family == "A_substantive", .(var, nr_routing_flag)])
  chk("positive/zero items on how a respondent with a fiancé(e) or partner met that person (DQ54_2*)", sum(grepl("^DQ54_2", toupper(fl[funnel_zero_denom > 0, var]))), 6)
  chk("positive/zero items with the routing flag", sort(fl[funnel_zero_denom > 0 & var %in% rt[nr_routing_flag %in% TRUE, var], var]), MASS$zero_denom_routing)
  ## fresh item nonresponse: flags that survive the allocation of the fresh missing mass (identity map infeasible even then)
  if (!"funnel_identity_feasible" %in% names(fl)) stop("15_mass_diagnostic_flags.csv has no funnel_identity_feasible column: rerun R/15_panelcond_designs.R with this kit", call. = FALSE)
  chk("flags not reconcilable by the fresh missing mass: all / finite > 1 / positive-zero / supported > 1",
      c(ms$n_flagged_not_reconcilable, ms$n_finite_gt1_not_reconcilable, ms$n_zero_denom_not_reconcilable, ms$n_supported_gt1_not_reconcilable), MASS$not_reconcilable)
  chk("the same, recounted from the flag list", c(sum((fl$funnel_ratio > 1 | fl$funnel_zero_denom > 0) & fl$funnel_identity_feasible %in% FALSE, na.rm = TRUE),
      sum(fl$funnel_ratio > 1 & fl$funnel_identity_feasible %in% FALSE, na.rm = TRUE)), MASS$not_reconcilable[1:2])
  chk("median fresh missing mass among eligible items (1.8%)", ms$median_fresh_missing_mass, MASS$median_missing, 5e-4)
  chk("the one flag not reconcilable is the minute of bedtime (DQ57DZ), a positive/zero item with no fresh missing mass",
      c(fl[funnel_identity_feasible %in% FALSE, var], as.character(fl[var == "DQ57DZ", funnel_zero_denom] > 0), as.character(fl[var == "DQ57DZ", funnel_missing_mass])), c("DQ57DZ", "TRUE", "0"))
  chk("largest mass needed to cover the stayers among the flagged items (.0074)", max(fl$funnel_needed_mass, na.rm = TRUE), 0.0074, 5e-5)
  tst <- rd("15_tests_items.csv")[dose_def == "exact" & family == "A_substantive" & !is.na(funnel_ncat)]
  gap <- tst$reach_new - tst$reach_old_S
  chk("reach rates of the two arms differ by less than ten points over the 419 items (max gap in points, one decimal)", round(100 * max(abs(gap)), 1), 9.1, 5e-2)
  chk("largest reach gaps: owner-occupied housing follow-ups (DQ39_*) and the unmarried block (DQ50)",
      c(all(grepl("^DQ39_", tst$var[abs(gap) > 0.08])), any(grepl("^DQ50$", tst$var[gap > 0.07]))), c(TRUE, TRUE))
  chk("items named in the text: reconcilable (spouse's occupation, party identification, meals, siblings)",
      as.numeric(fl[match(c("dq44_2l", "DQ30", "DQ45A", "DQ08B_4"), var), funnel_identity_feasible]), MASS$named_feasible)
  ## routing-threshold sensitivity (the flag recomputed at 10 and 20 points and with no routing exclusion)
  rs <- rd("15_routing_sensitivity.csv")[dose_def == "exact"]
  uniq(rs, "threshold", "15_routing_sensitivity (exact)")
  r15 <- one(rs[abs(threshold - 0.15) < 1e-9], "routing sensitivity, threshold .15")
  chk("routing gap recomputed in 15 at .15 reproduces 04's flag (disagreements)", r15$n_disagree_with_04, 0)
  chk("threshold .15 reproduces the main counts: naive / sm / ssm / ec / ec_adj", c(r15$aff_naive, r15$aff_sm, r15$aff_ssm, r15$aff_ec, r15$aff_ec_adj), COUNTS$flagged)
  chk("threshold .15 reproduces the routing exclusions and the follow-up", c(r15$n_routing, r15$n_followup), c(COUNTS$excluded[2] - 1, 1))
  sens <- function(th) { r <- one(if (is.finite(th)) rs[abs(threshold - th) < 1e-9] else rs[!is.finite(threshold)], paste("routing sensitivity", th))
    c(r$n_routing + r$n_followup, r$aff_naive, r$aff_sm, r$aff_ssm, r$aff_ec, r$aff_ec_adj, r$n_all5) }
  chk("threshold .10: excluded / flagged x 5 / all five", sens(0.10), SENS$t10)
  chk("threshold .20: excluded / flagged x 5 / all five", sens(0.20), SENS$t20)
  chk("no routing exclusion: excluded / flagged x 5 / all five", sens(Inf), SENS$none)
  rnone <- one(rs[!is.finite(threshold)], "routing sensitivity, none")
  chk("no routing exclusion: items in the counts (434) and the two items flagged by all five designs", c(rnone$n_items_naive, rnone$all5_items), c(434, COUNTS$all5_items_semicolon))

  ## --- 2019 episode, negative-control battery ---------------------------------
  chk("negative-control items at the 2019 episode", nc$n_items, 23)
  chk("median absolute contrast", nc$median_abs_d, 0.033, 5e-4)
  chk("rejections at BH q < .10", nc$n_reject_q10, 1)
  chk("TOST equivalences", nc$n_equiv, 9)
  chk("inverse-variance pooled contrast", round(nc$pooled_d_invvar, 3), -0.035, 5e-4)
  chk("its standard error under independence", round(nc$pooled_se, 4), 0.0095, 5e-5)

  ## --- appendix table: the 23 item contrasts (rounded to three decimals as printed) --------------------
  it23 <- rd("11_negcontrol_w13.csv")
  uniq(it23, c("stem", "contrast"), "11_negcontrol_w13 (item x contrast)")
  printed_d  <- c(-.033, -.032, .015, .040, -.074, -.092, -.103, .009, -.033, -.051, -.014, -.012, -.015, .025, .023, .022, -.020, -.176, -.069, .078, -.106, -.096, -.069)
  printed_se <- c(.045, .047, .045, .045, .045, .045, .046, .049, .045, .045, .045, .045, .045, .045, .045, .045, .045, .045, .046, .049, .045, .045, .046)
  chk("appendix table: 23 rows", nrow(it23), 23)
  chk("appendix table: d (rounded)", round(it23$d, 3), printed_d, 5e-4)
  chk("appendix table: SE (rounded)", round(it23$se, 3), printed_se, 5e-4)
  chk("appendix table: the one item with q < .10 is row 18 (telephone)", which(it23$q < .10), 18L)

  ## --- person-level bootstrap of the pooled value (11b) --------------------------------------------------
  p11b <- rd("11b_negcontrol_w13_pooled.csv")
  uniq(p11b, "estimator", "11b_negcontrol_w13_pooled")
  pb  <- one(p11b[estimator == "inverse-variance pooled (as in script 11)"], "11b pooled (as in 11)")
  pbx <- one(p11b[estimator == "inverse-variance pooled, BH-rejected item(s) dropped"], "11b pooled, rejected item dropped")
  chk("11b reproduces the pooled value of 11", round(pb$estimate, 4), round(nc$pooled_d_invvar, 4), 5e-5)
  chk("11b reproduces the independence SE of 11", round(pb$se_independence, 4), round(nc$pooled_se, 4), 5e-5)
  chk("2019 episode: 2007 cohort at tenure 13, 2011 cohort at tenure 9", c(pb$n_2007_t13, pb$n_2011_t9), c(2638, 619))
  chk("replications of the person bootstrap", pb$B, 2000)
  chk("person-bootstrap SE of the pooled value", round(pb$se_person_bootstrap, 4), 0.0186, 5e-5)
  chk("its ratio to the independence SE (2.0)", round(pb$se_person_bootstrap / pb$se_independence, 1), 2.0, 5e-2)
  chk("design effect (variance ratio)", round(pb$design_effect, 2), 3.84, 5e-3)
  chk("95% percentile interval, lower end", round(pb$ci95_lo, 3), -0.072, 5e-4)
  chk("95% percentile interval, upper end", round(pb$ci95_hi, 3), 0.001, 5e-4)
  chk("pooled value without the rejected item", round(pbx$estimate, 3), -0.028, 5e-4)
  chk("its bootstrap SE", round(pbx$se_person_bootstrap, 4), 0.0185, 5e-5)
  chk("mean correlation between item contrasts",
      round(one(rd("11b_negcontrol_w13_corr.csv"), "11b corr")$mean_offdiag_corr, 3), 0.116, 5e-4)
  lg <- rd("11b_loading_grid.csv"); uniq(lg, "Gamma", "11b_loading_grid")
  gam <- function(G, col) one(lg[abs(Gamma - G) < 1e-9], paste("loading grid, Gamma =", G))[[col]]
  chk("implied bias, point estimate, Gamma = .5 and 1.5", round(-c(gam(.5, "correction"), gam(1.5, "correction")), 3), c(.017, .052), 5e-4)
  chk("implied bias at the interval's far end, Gamma = .5 and 1.5", round(-c(gam(.5, "ci95_lo"), gam(1.5, "ci95_lo")), 3), c(.036, .108), 5e-4)
  chk("far end of the interval as a multiple of Gamma (.072)", round(-gam(1, "ci95_lo"), 3), .072, 5e-4)

  say(strrep("-", 112), "\n")
  say(sprintf("checked %d quantities: %d reproduced, %d mismatched\n", ok + bad, ok, bad))
  list(ok = ok, bad = bad, failed = unique(failed))
}

## values quoted in Section 10 that are set from the frozen rerun (counts, diagnostic denominators, mass diagnostic)
COUNTS <- list(n_naive = 402, n_matched = c(401, 401), flagged = c(22, 16, 9, 20, 19), excluded = c(44, 32), followup = "DQ46Y",
               common = c(15, 12, 6, 20, 19), n_all5 = 2, all5_items = "DQ26, DQ44_4A", all5_items_semicolon = "DQ26; DQ44_4A")
DIAG <- list(ns = c(30, 227), sd = c(36, 401), sd_nr = c(66, 396), share_nr = 0.167, grid4 = c(50, 66))   # rejections and denominators (frozen rerun)
SENS <- list(t10 = c(32, 22, 16, 9, 20, 19, 2), t20 = c(32, 22, 16, 9, 20, 19, 2), none = c(0, 23, 16, 10, 20, 19, 2))   # routing-threshold sensitivity (frozen rerun)
MASS <- list(max_finite_item = "dq44_2l", supported_items = c("DQ08B_4", "DQ45A"), eligible = 419, zero_denom_items = 14, zero_denom_categories = 16,
             not_reconcilable = c(1, 0, 1, 0), median_missing = 0.0177, named_feasible = c(1, 1, 1, 1),
             finite_gt1 = 13, max_finite = 2.256, supported_gt1 = 2, sparse_only_gt1 = 11, max_supported = 1.192,
             party_id = 1.955, zero_denom_routing = c("DQ49_2P", "DQ49_2Z"))

res <- run_checks(RES0)

if (SELFTEST) {
  cat("\n== self-test: the checker must fail on each corruption below ==\n")
  st_ok <- 0L; st_bad <- 0L
  expect_fail <- function(what, mutate) {
    tmp <- tempfile("chk_"); dir.create(tmp); file.copy(file.path(RES0, REQUIRED), tmp)
    mutate(tmp)
    r <- tryCatch(run_checks(tmp, verbose = FALSE), error = function(e) NULL)
    ## detected if the run stops, or if some check fails (or fails for another reason) that did not on the intact outputs
    failed <- is.null(r) || length(setdiff(r$failed, res$failed)) > 0L
    if (failed) st_ok <<- st_ok + 1L else st_bad <<- st_bad + 1L
    cat(sprintf("  %-72s %s\n", what, if (failed) "detected" else "NOT DETECTED"))
    unlink(tmp, recursive = TRUE)
  }
  csv_edit <- function(f, fn) function(d) { x <- fread(file.path(d, f), encoding = "UTF-8"); fwrite(fn(x), file.path(d, f)) }
  expect_fail("a required file is missing (11b_loading_grid.csv)", function(d) file.remove(file.path(d, "11b_loading_grid.csv")))
  expect_fail("a required file is missing (15_mass_diagnostic_summary.csv)", function(d) file.remove(file.path(d, "15_mass_diagnostic_summary.csv")))
  expect_fail("one appendix-table row deleted (11_negcontrol_w13.csv)", csv_edit("11_negcontrol_w13.csv", function(x) x[-5]))
  expect_fail("a detection-count row duplicated (15_detection_counts.csv)", csv_edit("15_detection_counts.csv", function(x) rbind(x, x[1])))
  expect_fail("the exact-dose arms row removed (15_arms.csv)", csv_edit("15_arms.csv", function(x) x[dose_def != "exact"]))
  expect_fail("a quoted value set to missing (15_arms.csv, n_old_S)", csv_edit("15_arms.csv", function(x) x[dose_def == "exact", n_old_S := NA]))
  expect_fail("the style margins of one design removed (15_designs_items.csv)",
              csv_edit("15_designs_items.csv", function(x) x[!(family == "P_style" & estimator == "ssm")]))
  expect_fail("another input named by 11_env.txt", function(d) {
    x <- readLines(file.path(d, "11_env.txt")); x <- sub("^input sha256: .*", paste("input sha256:", strrep("0", 64)), x)
    writeLines(x, file.path(d, "11_env.txt")) })
  ## the comparison rule itself: empty, missing, non-finite or wrong-length values never pass
  z <- c(is.null(compare(numeric(0), numeric(0))), is.null(compare(numeric(0), 1)), is.null(compare(NA, 1)),
         is.null(compare(1, NA)), is.null(compare(c(1, 2), 1)), is.null(compare(Inf, 1)), is.null(compare(NA_character_, "a")))
  hit <- !any(z) && is.null(compare(c(1, 2), c(1, 2))) && is.null(compare(1.0004, 1, 5e-4))
  if (hit) st_ok <- st_ok + 1L else st_bad <- st_bad + 1L
  cat(sprintf("  %-72s %s\n", "empty, missing, non-finite or wrong-length values never pass", if (hit) "detected" else "NOT DETECTED"))
  cat(sprintf("self-test: %d of %d corruptions detected\n", st_ok, st_ok + st_bad))
  if (st_bad > 0L) quit(status = 2)
}
if (res$bad > 0L) quit(status = 1)

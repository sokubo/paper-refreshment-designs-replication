#!/usr/bin/env Rscript
# ============================================================
# T2 — Section 10 と付録表(負対照バッテリー)の引用値を、キットの集計出力から機械的に再現する。
#   prose と出力が乖離しないようにするための照合スクリプト(T1 の check_manuscript_values.R の対応物)。
#   入力はすべて**集計値**(N<10 抑制済み)。個票は読まない。
#   実行: Rscript check_manuscript_values_T2.R [<results ディレクトリ>]
#     既定の results ディレクトリ: ../../P1_jlps_diagnosis/results
#   必要ファイル: 15_arms.csv, 15_detection_counts.csv, 15_diagnostics_summary.csv,
#                 15_designs_items.csv, 11_negcontrol_w13_summary.csv
# ============================================================
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
RES <- if (length(args)) args[1] else file.path("..", "..", "P1_jlps_diagnosis", "results")
need <- c("15_arms.csv","15_detection_counts.csv","15_diagnostics_summary.csv",
          "15_designs_items.csv","11_negcontrol_w13_summary.csv")
miss <- need[!file.exists(file.path(RES, need))]
if (length(miss)) stop("missing aggregate outputs in ", RES, ": ", paste(miss, collapse=", "))

arms <- fread(file.path(RES,"15_arms.csv"))[dose_def=="exact"]
dc   <- fread(file.path(RES,"15_detection_counts.csv"))[dose_def=="exact" & family=="A_substantive"]
dg   <- fread(file.path(RES,"15_diagnostics_summary.csv"))[dose_def=="exact"]
it   <- fread(file.path(RES,"15_designs_items.csv"))[dose_def=="exact"]
nc   <- fread(file.path(RES,"11_negcontrol_w13_summary.csv"))

ok <- 0L; bad <- 0L
chk <- function(what, got, want, tol = 0) {
  hit <- if (is.character(want)) identical(as.character(got), want) else
         all(abs(as.numeric(got) - as.numeric(want)) <= tol)
  if (isTRUE(hit)) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("%-58s %-22s %-22s %s\n", what,
              paste(format(got), collapse=" "), paste(format(want), collapse=" "),
              if (isTRUE(hit)) "ok" else "MISMATCH"))
}
cat(sprintf("%-58s %-22s %-22s %s\n", "quantity (Section 10)", "from outputs", "quoted in paper", ""))
cat(strrep("-", 112), "\n")

## --- arms and bootstrap -----------------------------------------------------
chk("continuing survivors of five waves", arms$n_old_S, 2797)
chk("fresh entrants at the 2011 refreshment", arms$n_new_total, 963)
chk("entrants surviving four further waves", arms$n_new_sm, 574)
chk("entrants surviving five further waves", arms$n_new_sm1, 540)
chk("bootstrap replications", arms$B, 500)
chk("items with the same coding at entry", arms$n_ec_items, 268)

## --- waterfall over all items with variation --------------------------------
g <- function(e, col) dc[estimator == e][[col]]
chk("substantive items with variation (naive arm)", g("naive","n_items"), 470)
chk("flagged: naive", g("naive","n_affected"), 28)
chk("flagged: survival matching", g("sm","n_affected"), 19)
chk("flagged: symmetric survival matching", g("ssm","n_affected"), 13)
chk("items with an available entry wave", g("ec","n_items"), 265)
chk("flagged: entry-wave correction", g("ec","n_affected"), 22)
chk("flagged: entry-wave correction, standardised", g("ec_adj","n_affected"), 21)

## --- restricted to the common set of 265 ------------------------------------
A <- it[family == "A_substantive"]
common <- A[estimator == "ec", unique(var)]
cnt <- sapply(c("naive","sm","ssm","ec","ec_adj"),
              function(e) A[estimator == e & var %in% common & class3 == "affected", uniqueN(var)])
chk("common set of 265: naive / sm / ssm / ec / ec_adj", cnt, c(19,15,9,22,21))
five <- Reduce(intersect, lapply(c("naive","sm","ssm","ec","ec_adj"),
               function(e) A[estimator == e & class3 == "affected", unique(var)]))
chk("items flagged by all five designs", length(five), 4)
cat("      (the four items: ", paste(sort(five), collapse=", "), ")\n", sep="")

## --- diagnostics ------------------------------------------------------------
chk("T_NS rejections, substantive items", dg[family=="A_substantive", share_T1_p05], 0.131, 5e-4)
chk("T_SD rejections, substantive items", dg[family=="A_substantive", share_T2_p05], 0.082, 5e-4)
chk("T_SD rejections, item nonresponse (paper says 19%)",
    round(dg[family=="B_itemnonresp", share_T2_p05], 2), 0.19, 5e-3)

## --- response-style margins -------------------------------------------------
st <- it[family == "P_style"]
ext <- st[grepl("ext_share", var), range(d_std)]; mid <- st[grepl("mid_share", var), range(d_std)]
chk("extreme-category use, range over the five designs", round(ext,2), c(-0.25,-0.21))
chk("midpoint use, range over the five designs",        round(mid,2), c(0.19,0.21))
cat(sprintf("      extreme %.4f..%.4f   midpoint %.4f..%.4f  (paper: a fall of .21 to .25, a rise of .19 to .21)\n",
            ext[1], ext[2], mid[1], mid[2]))

## --- the employment item ----------------------------------------------------
emp <- A[var == "DQ02"]
chk("employment, entry-wave corrected (percentage points)",
    round(100 * emp[estimator=="ec", estimate], 1), -4.9, 0.05)
chk("employment, standardised: q above .10",
    as.numeric(emp[estimator=="ec_adj", q] > 0.10), 1)

## --- mass-domination flags --------------------------------------------------
fr <- unique(A[!is.na(funnel_ratio), .(var, funnel_ratio)])
chk("items with a mass ratio above one", fr[funnel_ratio > 1, .N], 7)
chk("largest mass ratio (party identification)", round(fr[, max(funnel_ratio)], 2), 1.96, 0.005)
cat("      ratios above one: ", paste(sprintf("%s %.2f", fr[funnel_ratio>1][order(-funnel_ratio), var],
                                              fr[funnel_ratio>1][order(-funnel_ratio), funnel_ratio]),
                                      collapse="; "), "\n", sep="")

## --- 2019 episode, negative-control battery ---------------------------------
chk("negative-control items at the 2019 episode", nc$n_items, 23)
chk("median absolute contrast", nc$median_abs_d, 0.033, 5e-4)
chk("rejections at BH q < .10", nc$n_reject_q10, 1)
chk("TOST equivalences", nc$n_equiv, 9)
chk("inverse-variance pooled contrast", round(nc$pooled_d_invvar, 3), -0.035, 5e-4)
chk("its standard error under independence", round(nc$pooled_se, 4), 0.0095, 5e-5)

## --- appendix table: the 23 item contrasts (rounded to three decimals as printed) --------------------
if (file.exists(file.path(RES, "11_negcontrol_w13.csv"))) {
  it23 <- fread(file.path(RES, "11_negcontrol_w13.csv"), encoding = "UTF-8")
  printed_d  <- c(-.033, -.032, .015, .040, -.074, -.092, -.103, .009, -.033, -.051, -.014, -.012, -.015, .025, .023, .022, -.020, -.176, -.069, .078, -.106, -.096, -.069)
  printed_se <- c(.045, .047, .045, .045, .045, .045, .046, .049, .045, .045, .045, .045, .045, .045, .045, .045, .045, .045, .046, .049, .045, .045, .046)
  chk("appendix table: 23 rows", nrow(it23), 23)
  chk("appendix table: d (rounded)", round(it23$d, 3), printed_d, 5e-4)
  chk("appendix table: SE (rounded)", round(it23$se, 3), printed_se, 5e-4)
  chk("appendix table: the one item with q < .10 is row 18 (telephone)", which(it23$q < .10), 18L)
}
## --- person-level bootstrap of the pooled value (11b) --------------------------------------------------
f11b <- file.path(RES, "11b_negcontrol_w13_pooled.csv")
if (file.exists(f11b)) {
  pb <- fread(f11b)[1]
  chk("11b reproduces the pooled value of 11", round(pb$estimate, 4), round(nc$pooled_d_invvar, 4), 5e-5)
  chk("11b reproduces the independence SE of 11", round(pb$se_independence, 4), round(nc$pooled_se, 4), 5e-5)
  pbx <- fread(f11b)[estimator == "inverse-variance pooled, BH-rejected item(s) dropped"]
  chk("replications of the person bootstrap", pb$B, 2000)
  chk("person-bootstrap SE of the pooled value", round(pb$se_person_bootstrap, 4), 0.0186, 5e-5)
  chk("its ratio to the independence SE (2.0)", round(pb$se_person_bootstrap / pb$se_independence, 1), 2.0, 5e-2)
  chk("design effect (variance ratio)", round(pb$design_effect, 2), 3.84, 5e-3)
  chk("95% percentile interval, lower end", round(pb$ci95_lo, 3), -0.072, 5e-4)
  chk("95% percentile interval, upper end", round(pb$ci95_hi, 3), 0.001, 5e-4)
  chk("pooled value without the rejected item", round(pbx$estimate, 3), -0.028, 5e-4)
  chk("its bootstrap SE", round(pbx$se_person_bootstrap, 4), 0.0185, 5e-5)
  fc <- file.path(RES, "11b_negcontrol_w13_corr.csv")
  if (file.exists(fc)) chk("mean correlation between item contrasts", round(fread(fc)$mean_offdiag_corr, 3), 0.116, 5e-4)
  fg <- file.path(RES, "11b_loading_grid.csv")
  if (file.exists(fg)) {
    lg <- fread(fg)
    chk("implied bias, point estimate, Gamma = .5 and 1.5", round(-lg[Gamma %in% c(.5, 1.5)]$correction, 3), c(.017, .052), 5e-4)
    chk("implied bias at the interval's far end, Gamma = .5 and 1.5", round(-lg[Gamma %in% c(.5, 1.5)]$ci95_lo, 3), c(.036, .108), 5e-4)
    chk("far end of the interval as a multiple of Gamma (.072)", round(-lg[Gamma == 1]$ci95_lo, 3), .072, 5e-4)
  }
} else cat("      11b_negcontrol_w13_pooled.csv not found: run R/11b_negcontrol_w13_boot.R first\n")

cat(strrep("-", 112), "\n")
cat(sprintf("checked %d quantities: %d reproduced, %d mismatched\n", ok + bad, ok, bad))
if (bad > 0) quit(status = 1)

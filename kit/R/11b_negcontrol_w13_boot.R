# ============================================================
# P1 JLPS kit — 11b: w13 incumbent negative-control battery, covariance-aware pooling (2026-09-23)
#
# Purpose.  Script 11 pools the 23 item contrasts d_k (2007 cohort at tenure 13 vs 2011 cohort at tenure 9,
#   common birth-year support, 3-year birth band x sex strata) by inverse-variance weights and reports
#   SE = sqrt(1 / sum(1/se_k^2)).  That SE treats the 23 items as independent, but all items are measured on the
#   same respondents.  This script recomputes the same item contrasts and pooled values and adds a
#   PERSON-LEVEL BOOTSTRAP (persons resampled with replacement within each of the two cohorts; every item
#   contrast, every stratum weight and the pooling weights recomputed in each replicate), which keeps the
#   covariance between items.  It also reports the equal-weight mean, the median, the pooled value without the
#   item(s) rejected at BH q < .10 in script 11, the bootstrap correlation of the item contrasts, and a
#   loading-sensitivity grid Gamma x pooled value.
#
# Definitions are identical to script 11 (same item map, same w13-participation rule, same cohort coding,
#   same support, same strata, same >= 10 per stratum rule).  Point estimates therefore reproduce 11's
#   11_negcontrol_w13_summary.csv (checked below and written to 11b_check_vs_11.csv).
#
# Inputs : RAW_DIR master .dta (outside the project folder), results/03_negcontrol_sets_v1.csv (item map).
# Outputs (aggregates only; no cell below MIN_CELL): results/11b_negcontrol_w13_pooled.csv,
#          results/11b_negcontrol_w13_corr.csv, results/11b_loading_grid.csv, results/11b_check_vs_11.csv,
#          results/11b_env.txt
# Run    : cd <P1 root> && Rscript analysis/R/11b_negcontrol_w13_boot.R [--B 2000] [--seed 20260923]
#          (about 2-5 minutes with B = 2000)
# ============================================================

.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
suppressMessages({library(haven); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
B_BOOT <- as.integer(getarg("--B", "2000"))
SEED11B <- as.integer(getarg("--seed", "20260923"))
MIN_ANSWERS_W13 <- 3
GAMMA_GRID <- c(0, 0.5, 0.75, 1, 1.25, 1.5, 2)

main <- function() {
  t0 <- Sys.time()
  nc <- fread(file.path(RESULTS_DIR, "03_negcontrol_sets_v1.csv"))
  nc <- nc[complete3 == 1]
  f_all <- list.files(RAW_DIR, pattern = "\\.dta$", full.names = TRUE)
  f <- f_all[!grepl("online", basename(f_all), ignore.case = TRUE)][1]
  d <- as.data.table(read_dta(f))
  gc_ <- function(nm) { hit <- names(d)[toupper(trimws(names(d))) == toupper(nm)]; if (!length(hit)) return(NULL); d[[hit[1]]] }
  num <- function(x) if (is.null(x)) NULL else suppressWarnings(as.numeric(zap_labels(x)))
  cn <- as.integer(zap_labels(gc_("CN")))
  SPEC_PAT <- "わからない|分からない|ＤＫ|DK|答えたくない|回答したくない|無回答|不明|非該当"
  clean <- function(nm) {
    x <- gc_(nm); if (is.null(x)) return(NULL)
    vl <- attr(x, "labels", exact = TRUE); xv <- num(x)
    if (!is.null(vl)) { spec <- unname(vl[grepl(SPEC_PAT, names(vl))]); xv[xv %in% spec] <- NA }
    xv
  }
  lvars <- grep("^[Ll][Qq][0-9]", names(d), value = TRUE)
  lmat <- as.matrix(d[, lapply(.SD, function(x) suppressWarnings(as.numeric(zap_labels(x)))), .SDcols = lvars])
  in13 <- rowSums(!is.na(lmat)) >= MIN_ANSWERS_W13
  grp <- rep(NA_character_, nrow(d))
  grp[cn == 1 & in13] <- "c2007_t13"; grp[cn == 2 & in13] <- "c2011_t9"
  yb <- num(gc_("ybirth")); sx <- num(gc_("sex"))
  if (is.null(yb) || is.null(sx)) stop("sex / ybirth not found")
  lo <- max(min(yb[grp == "c2007_t13"], na.rm = TRUE), min(yb[grp == "c2011_t9"], na.rm = TRUE))
  hi <- min(max(yb[grp == "c2007_t13"], na.rm = TRUE), max(yb[grp == "c2011_t9"], na.rm = TRUE))
  sup <- !is.na(yb) & yb >= lo & yb <= hi
  strat <- as.integer(interaction(3 * floor(yb / 3), sx, drop = TRUE))
  I1 <- which(grp == "c2007_t13" & sup); I0 <- which(grp == "c2011_t9" & sup)
  cat(sprintf("diag: support [%s, %s]; persons in the two incumbent cohorts on the support: %d and %d\n", lo, hi, length(I1), length(I0)))

  ## item vectors: 2007 cohort measured in ZQ (w1), 2011 cohort in DQ (w5), as in script 11
  items <- list(); labs <- character(0)
  for (i in seq_len(nrow(nc))) {
    yZ <- clean(nc$w1_Z[i]); yD <- clean(nc$w5_D[i])
    if (is.null(yZ) || is.null(yD)) next
    y <- rep(NA_real_, nrow(d)); y[grp == "c2007_t13" & !is.na(grp)] <- yZ[grp == "c2007_t13" & !is.na(grp)]
    y[grp == "c2011_t9" & !is.na(grp)] <- yD[grp == "c2011_t9" & !is.na(grp)]
    items[[length(items) + 1]] <- y; labs <- c(labs, nc$stem[i])
  }
  K <- length(items); cat("diag: items =", K, "\n")

  ## stratified standardised difference on resampled index sets (identical rule to script 11's sdiff)
  sdiff_idx <- function(y, i1, i0) {
    a_all <- y[i1]; b_all <- y[i0]; sa <- strat[i1]; sb <- strat[i0]
    ka <- !is.na(a_all); kb <- !is.na(b_all); a_all <- a_all[ka]; sa <- sa[ka]; b_all <- b_all[kb]; sb <- sb[kb]
    est <- 0; wsum <- 0
    for (st in sort(unique(c(sa, sb)))) {
      a <- a_all[sa == st]; b <- b_all[sb == st]
      if (length(a) < 10 || length(b) < 10) next
      sp <- sqrt((var(a) * (length(a) - 1) + var(b) * (length(b) - 1)) / (length(a) + length(b) - 2))
      if (!is.finite(sp) || sp == 0) next
      vv <- 1 / length(a) + 1 / length(b)
      est <- est + ((mean(a) - mean(b)) / sp) / vv; wsum <- wsum + 1 / vv
    }
    if (wsum == 0) return(c(NA_real_, NA_real_))
    c(est / wsum, sqrt(1 / wsum))
  }
  pooled <- function(i1, i0, drop = integer(0)) {
    r <- vapply(seq_len(K), function(k) sdiff_idx(items[[k]], i1, i0), numeric(2))
    dd <- r[1, ]; se <- r[2, ]; keep <- setdiff(which(!is.na(dd)), drop)
    w <- 1 / se[keep]^2
    c(invvar = sum(w * dd[keep]) / sum(w), se_indep = sqrt(1 / sum(w)),
      equal = mean(dd[keep]), se_equal_indep = sqrt(sum(se[keep]^2)) / length(keep), median = median(dd[keep]), dd)
  }
  full <- pooled(I1, I0)
  dfull <- full[-(1:5)]; sefull <- vapply(seq_len(K), function(k) sdiff_idx(items[[k]], I1, I0)[2], 0)
  p <- 2 * pnorm(-abs(dfull / sefull)); q <- p.adjust(p, "BH"); rej <- which(q < 0.10)
  full_x <- pooled(I1, I0, drop = rej)

  ## person bootstrap within cohort
  set.seed(SEED11B)
  bm <- matrix(NA_real_, B_BOOT, 5 + K); bx <- matrix(NA_real_, B_BOOT, 3)
  for (b in seq_len(B_BOOT)) {
    i1 <- sample(I1, length(I1), replace = TRUE); i0 <- sample(I0, length(I0), replace = TRUE)
    r <- pooled(i1, i0); bm[b, ] <- r
    rx <- pooled(i1, i0, drop = rej); bx[b, ] <- rx[c("invvar", "equal", "median")]
    if (b %% 250 == 0) cat(sprintf("   replicate %d / %d (%.1f min)\n", b, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  sdb <- function(v) sd(v, na.rm = TRUE); ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
  out <- data.table(
    estimator = c("inverse-variance pooled (as in script 11)", "equal-weight mean", "median of item contrasts",
                  "inverse-variance pooled, BH-rejected item(s) dropped", "equal-weight mean, BH-rejected item(s) dropped",
                  "median, BH-rejected item(s) dropped"),
    n_items = c(K, K, K, K - length(rej), K - length(rej), K - length(rej)),
    estimate = round(c(full["invvar"], full["equal"], full["median"], full_x["invvar"], full_x["equal"], full_x["median"]), 5),
    se_independence = round(c(full["se_indep"], full["se_equal_indep"], NA, full_x["se_indep"], full_x["se_equal_indep"], NA), 5),
    se_person_bootstrap = round(c(sdb(bm[, 1]), sdb(bm[, 3]), sdb(bm[, 5]), sdb(bx[, 1]), sdb(bx[, 2]), sdb(bx[, 3])), 5),
    ci95_lo = round(c(ci(bm[, 1])[1], ci(bm[, 3])[1], ci(bm[, 5])[1], ci(bx[, 1])[1], ci(bx[, 2])[1], ci(bx[, 3])[1]), 5),
    ci95_hi = round(c(ci(bm[, 1])[2], ci(bm[, 3])[2], ci(bm[, 5])[2], ci(bx[, 1])[2], ci(bx[, 2])[2], ci(bx[, 3])[2]), 5),
    B = B_BOOT, seed = SEED11B, n_2007_t13 = length(I1), n_2011_t9 = length(I0))
  out[, design_effect := round((se_person_bootstrap / se_independence)^2, 3)]
  write_aggregate(out, "11b_negcontrol_w13_pooled.csv", exempt = c("n_items", "B", "seed", "n_2007_t13", "n_2011_t9"))
  C <- suppressWarnings(cor(bm[, -(1:5)], use = "pairwise.complete.obs")); off <- C[upper.tri(C)]
  corr <- data.table(n_items = K, mean_offdiag_corr = round(mean(off, na.rm = TRUE), 4), median_offdiag_corr = round(median(off, na.rm = TRUE), 4),
                     min_offdiag_corr = round(min(off, na.rm = TRUE), 4), max_offdiag_corr = round(max(off, na.rm = TRUE), 4),
                     share_positive = round(mean(off > 0, na.rm = TRUE), 3), B = B_BOOT)
  write_aggregate(corr, "11b_negcontrol_w13_corr.csv", exempt = c("n_items", "B"))
  grid <- data.table(Gamma = GAMMA_GRID, correction = round(GAMMA_GRID * full["invvar"], 5),
                     ci95_lo = round(GAMMA_GRID * ci(bm[, 1])[1], 5), ci95_hi = round(GAMMA_GRID * ci(bm[, 1])[2], 5),
                     note = "correction = Gamma x pooled incumbent negative-control contrast (SD units); subtract from the target increment")
  write_aggregate(grid, "11b_loading_grid.csv", exempt = c("Gamma"))
  s11 <- file.path(RESULTS_DIR, "11_negcontrol_w13_summary.csv")
  chk <- data.table(quantity = c("pooled_d_invvar", "pooled_se_independence", "n_items", "n_reject_q10"),
                    value_11b = c(round(full["invvar"], 4), round(full["se_indep"], 4), K, length(rej)))
  if (file.exists(s11)) { s <- fread(s11); chk[, value_11 := c(s$pooled_d_invvar[1], s$pooled_se[1], s$n_items[1], s$n_reject_q10[1])]
    chk[, agree := abs(value_11b - value_11) < 1e-4] }
  write_aggregate(chk, "11b_check_vs_11.csv", exempt = c("value_11b", "value_11"))
  writeLines(c(paste("time:", format(Sys.time())), R.version.string, paste("data.table", packageVersion("data.table")),
               paste("haven", packageVersion("haven")), paste("input:", basename(f)), paste("B:", B_BOOT, "seed:", SEED11B),
               paste("rejected item index (BH q < .10):", paste(rej, collapse = ","))), file.path(RESULTS_DIR, "11b_env.txt"))
  print(out); print(corr); print(grid); print(chk)
  cat(sprintf("\n== 11b done (%.1f min). Please share results/11b_*.csv and 11b_env.txt (aggregates only). ==\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

run_guarded("11b_negcontrol_w13_boot", main)

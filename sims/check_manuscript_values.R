# ============================================================
# check_manuscript_values.R -- every simulation and design-rank value quoted in the manuscript, asserted against
# a named cell of an output file.  Run from sims/ after the simulations and checks:
#     Rscript check_manuscript_values.R [path/to/main.qmd]
# Each assertion (i) selects exactly one row of an output file by its key columns, (ii) formats the value as it
# is printed in the manuscript, and (iii) requires that exact string (with its units stated in the label) to
# appear in main.qmd.  Range statements are asserted from the minimum and maximum of the rows they cover.
# Exits with status 1 on any failure; --selftest confirms that a perturbed value and a missing row make it fail.
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
selftest <- "--selftest" %in% args
qmd_path <- if (length(setdiff(args, "--selftest"))) setdiff(args, "--selftest")[1] else file.path("..", "manuscript", "main.qmd")
QMD <- paste(readLines(qmd_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
n_ok <- 0; n_fail <- 0
ok <- function(cond, label) {
  if (isTRUE(cond)) { n_ok <<- n_ok + 1 } else { n_fail <<- n_fail + 1; cat("FAIL:", label, "\n") }
}
has <- function(txt, label) ok(grepl(txt, QMD, fixed = TRUE, useBytes = TRUE), sprintf("%s -- expected text: %s", label, txt))
one <- function(df, ...) {                      # exactly one row matching all key = value pairs
  keys <- list(...); sel <- rep(TRUE, nrow(df))
  for (k in names(keys)) sel <- sel & (as.character(df[[k]]) == as.character(keys[[k]]))
  if (sum(sel) != 1) stop(sprintf("key selects %d rows: %s", sum(sel), paste(names(keys), keys, collapse = ", ")))
  df[sel, , drop = FALSE]
}
d3 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.3f", x))            # .123 / -.123
d2 <- function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))
m3 <- function(x) sub("^-", intToUtf8(0x2212), d3(x))                                 # table minus sign
tx3 <- function(x) sprintf("$%s$", d3(x))

run_checks <- function(s1, s2, p2, e3, g3, v3, pr3, D12, D3, drk) {
  ## ---------------- Simulation 1 (sim1_results.csv; outcome units, h in SD of Y2*) ----------------
  nb <- function(b, p) s1$naive_bias[s1$beta == b & s1$p_target == p]
  ok(max(abs(nb(0, .6)), abs(nb(0, .8))) <= .004, "Sim1: naive bias at beta = 0 at most .004")
  has(sprintf("$%s$--$%s$ at $\\beta = .3$ and $%s$--$%s$ at $\\beta = .6$ when $p = .8$", d3(min(nb(.3, .8))), d3(max(nb(.3, .8))), d3(min(nb(.6, .8))), d3(max(nb(.6, .8)))), "Sim1: naive bias, p = .8")
  has(sprintf("and $%s$--$%s$ and $%s$--$%s$ when $p = .6$", d3(min(nb(.3, .6))), d3(max(nb(.3, .6))), d3(min(nb(.6, .6))), d3(max(nb(.6, .6)))), "Sim1: naive bias, p = .6")
  pw <- function(b, p) unique(round(s1$pop_width[s1$beta == b & s1$p_target == p], 10))
  for (b in c(0, .3, .6)) for (p in c(.6, .8)) ok(length(pw(b, p)) == 1, sprintf("Sim1: population width unique in (beta %s, p %s)", b, p))
  has(sprintf("widths $0$, $%s$ and $%s$ at $p = .8$ and $0$, $%s$ and $%s$ at $p = .6$", d3(pw(.3, .8)), d3(pw(.6, .8)), d3(pw(.3, .6)), d3(pw(.6, .6))), "Sim1: population identified-set widths")
  ok(pw(0, .6) == 0 && pw(0, .8) == 0, "Sim1: beta = 0 population set is a point")
  mw <- function(p) range(s1$mean_width[s1$p_target == p]); rw <- function(p) range(s1$relaxed_width[s1$p_target == p])
  has(sprintf("mean widths $%s$--$%s$ at $p = .8$ and $%s$--$%s$ at $p = .6$", d2(mw(.8)[1]), d2(mw(.8)[2]), d2(mw(.6)[1]), d2(mw(.6)[2])), "Sim1: relaxed plug-in widths")
  has(sprintf("($%s$--$%s$ at $p = .8$, $%s$--$%s$ at $p = .6$)", d2(rw(.8)[1]), d2(rw(.8)[2]), d2(rw(.6)[1]), d2(rw(.6)[2])), "Sim1: population relaxed widths")
  cv <- function(p) range(s1$coverage[s1$p_target == p])
  ok(cv(.6)[2] == 1, "Sim1: max coverage at p = .6 is 1")
  has(sprintf("is $%s$--$1$ at $p = .6$ and $%s$--$%s$ at $p = .8$", d3(cv(.6)[1]), d3(cv(.8)[1]), d3(cv(.8)[2])), "Sim1: coverage ranges")
  has(sprintf("width $%s\\sigma$, against mean plug-in widths of $%s$--$%s$", d2(unique(s1$relaxed_width[s1$beta == 0 & s1$p_target == .8])),
              d2(min(s1$mean_width[s1$beta == 0 & s1$p_target == .8])), d2(max(s1$mean_width[s1$beta == 0 & s1$p_target == .8]))), "Remark 2: relaxed width vs plug-in")
  has(sprintf("widths between $0$ and $%s\\sigma$ at $p = .8$", d2(pw(.6, .8))), "Section 5: location-set widths in Simulation 1")

  ## ---------------- Simulation 2 (sim2_results.csv, sim2_plim.csv; population-SD units) ----------------
  m <- one(s2, cell = "matched"); gcell <- one(s2, cell = "generational"); dr <- one(s2, cell = "drift3")
  pinc <- one(p2, contrast = "increment", generational = FALSE)
  ok(abs(pinc$target - 0.15 * (exp(-2) - exp(-3)) / sqrt(2)) < 1e-12, "Sim2: increment target in population-SD units")
  has(sprintf("$(\\tau(13) - \\tau(9))/\\sqrt2 = .15(e^{-2} - e^{-3})/\\sqrt2 = %s$", sub("^0", "", sprintf("%.4f", pinc$target))), "Sim2: increment target")
  has(sprintf("$(\\tau(13) - \\tau(1))/\\sqrt2 = %s$", sub("^0", "", sprintf("%.4f", one(p2, contrast = "level", generational = FALSE)$target))), "Sim2: level target")
  has(sprintf("biased by $%s$ (population limit $%s$)", d3(m$bias_inc_naive), d3(pinc$plim_naive - pinc$target)), "Sim2: naive increment")
  has(sprintf("has bias $%s$ (Monte Carlo standard error $%s$; population limit zero) and RMSE $%s$", d3(m$bias_inc_trueG), d3(m$mcse_inc_trueG), d3(m$rmse_inc_trueG)), "Sim2: true-loading correction")
  ok(abs(pinc$asy_bias_trueG) < 1e-12, "Sim2: true loading is asymptotically unbiased (raw-scale transport)")
  has(sprintf("the unit loading leaves $%s$ (limit $%s$)", d3(m$bias_inc_unitG), d3(pinc$asy_bias_unitG)), "Sim2: unit loading")
  has(sprintf("has bias $%s$ (standard error $%s$) in the age-matched design and $%s$ in the generational design", d3(m$bias_lvl_trueG), d3(m$mcse_lvl_trueG), d3(gcell$bias_lvl_trueG)), "Sim2: level contrasts")
  has(sprintf("$-g/\\sqrt2 = %s$", d3(-0.35 / sqrt(2))), "Sim2: generational offset")
  has(sprintf("stays unbiased there ($%s$)", if (abs(gcell$bias_inc_trueG) < .0005) "-.000" else d3(gcell$bias_inc_trueG)), "Sim2: generational increment")
  qr <- range(s2$Q_reject_rate[s2$cell != "drift3"])
  has(sprintf("rejects at $%s$--$%s$ at nominal $.05$ without drift and at $%s$ with three drifted items", d3(qr[1]), d3(qr[2]), d2(dr$Q_reject_rate)), "Sim2: Q test")
  has(sprintf("in $%d\\%%$ of replications, against $%d\\%%$", round(100 * dr$item_flag_rate), round(100 * m$item_flag_rate)), "Sim2: item screen")
  ok(max(s2$false_flag_rate) < .01, "Sim2: per-item false-flag rate below .01")
  has(sprintf("(bias $%s$); screening before pooling ($%s$) and median pooling ($%s$)", d3(dr$bias_inc_trueG), d3(dr$bias_inc_screen), d3(dr$bias_inc_median)), "Sim2: drift contamination")

  ## ---------------- Simulation 3 (sim3_*.csv; outcome units) ----------------
  E <- function(rn, es, g = 0.3) one(e3, regime = rn, gamma = g, estimator = es)
  G <- function(rn, g = 0.3) one(g3, regime = rn, gamma = g)
  labs <- c(MCAR = "MCAR", MAR_X = "MAR on $X$", MNAR_trait = "MNAR trait", MNAR_nonstat = "MNAR non-stationary",
            MNAR_state = "MNAR state-dependent", MNAR_both = "MNAR, both failures")
  bold <- list(MNAR_nonstat = c("sm", "ssm", "TNS"), MNAR_state = c("sm", "ec", "TSD"), MNAR_both = c("ssm", "ec", "TNS", "TSD"))
  fb <- function(x) { s <- m3(x); if (x > 0 && s != ".000") s <- paste0("+", s); s }
  for (rn in names(labs)) {
    cell <- function(es) { r <- E(rn, es); v <- if (es == "ipw") fb(r$bias_own) else sprintf("%s (%s)", fb(r$bias_own), d2(r$cover95))
                           if (es %in% bold[[rn]]) paste0("**", v, "**") else v }
    rj <- function(st) { v <- d3(G(rn)[[paste0("rej_", st)]]); if (st %in% bold[[rn]]) paste0("**", v, "**") else v }
    row <- sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", labs[rn], cell("naive"), cell("sm"), cell("ssm"), cell("ec"), cell("ipw"), rj("TNS"), rj("TSD"))
    has(row, sprintf("Sim3 table row %s", rn))
  }
  has(sprintf("at most %s and of the rejection rates at most %s", sprintf("%.3f", ceiling(max(e3$mcse_bias) * 1000) / 1000), sprintf("%.3f", ceiling(max(g3$mcse_rej_TNS, g3$mcse_rej_TSD) * 1000) / 1000)), "Sim3: Monte Carlo SE bounds")
  allowed <- rbind(E("MCAR", "sm"), E("MCAR", "ssm"), E("MCAR", "ec"), E("MAR_X", "sm"), E("MAR_X", "ssm"), E("MAR_X", "ec"),
                   E("MNAR_trait", "sm"), E("MNAR_trait", "ssm"), E("MNAR_trait", "ec"), E("MNAR_nonstat", "ec"), E("MNAR_state", "ssm"))
  ok(max(abs(allowed$bias_own)) < .003, "Sim3: |bias| < .003 where the restriction holds")
  ok(d2(min(allowed$cover95)) == ".94" && d2(max(allowed$cover95)) == ".96", "Sim3: coverage .94 to .96 where the restriction holds")
  has("biases below $0.003$ in absolute value (outcome units) and 95% coverage between $0.94$ and $0.96$", "Sim3: text of allowed regimes")
  sp <- function(rn, es) { r <- E(rn, es); s <- sub("^-", "-", sprintf("%.3f", r$bias_own)); if (r$bias_own > 0) s <- paste0("+", s); sprintf("$%s$, %s$%s$", s, "", sprintf("%.2f", r$cover95)) }
  has(sprintf("under non-stationarity ($%s$, coverage $%s$) and under state dependence ($+%s$, $%s$)", sprintf("%.3f", E("MNAR_nonstat", "sm")$bias_own), sprintf("%.2f", E("MNAR_nonstat", "sm")$cover95), sprintf("%.3f", E("MNAR_state", "sm")$bias_own), sprintf("%.2f", E("MNAR_state", "sm")$cover95)), "Sim3/Sec8: SM")
  has(sprintf("under state dependence ($+%s$, $%s$) but not under non-stationarity ($%s$, $%s$)", sprintf("%.3f", E("MNAR_state", "ec")$bias_own), sprintf("%.2f", E("MNAR_state", "ec")$cover95), sprintf("%.3f", E("MNAR_nonstat", "ec")$bias_own), sprintf("%.2f", E("MNAR_nonstat", "ec")$cover95)), "Sim3/Sec8: EC")
  has(sprintf("under non-stationarity ($%s$, $%s$) but not under state dependence ($+%s$, $%s$)", sprintf("%.3f", E("MNAR_nonstat", "ssm")$bias_own), sprintf("%.2f", E("MNAR_nonstat", "ssm")$cover95), sprintf("%.3f", E("MNAR_state", "ssm")$bias_own), sprintf("%.2f", E("MNAR_state", "ssm")$cover95)), "Sim3/Sec8: SSM")
  sz <- unlist(lapply(c("MCAR", "MAR_X", "MNAR_trait"), function(rn) c(G(rn)$rej_TNS, G(rn)$rej_TSD)))
  has(sprintf("reject at rates between $%s$ and $%s$", sprintf("%.3f", min(sz)), sprintf("%.3f", max(sz))), "Sec8: diagnostic size")
  has(sprintf("rejects at $%s$ under non-stationarity and at $%s$ under state dependence", sprintf("%.2f", G("MNAR_nonstat")$rej_TNS), sprintf("%.3f", G("MNAR_state")$rej_TNS)), "Sec8: T_NS")
  has(sprintf("at $%s$ under state dependence and at $%s$ under non-stationarity", sprintf("%.2f", G("MNAR_state")$rej_TSD), sprintf("%.2f", G("MNAR_nonstat")$rej_TSD)), "Sec8: T_SD")
  ipw <- sapply(c("MNAR_trait", "MNAR_nonstat", "MNAR_state", "MNAR_both"), function(rn) E(rn, "ipw")$bias_own)
  has(sprintf("biased by $%s$ to $%s$ under every regime with selection on the unobserved trait", sprintf("%.2f", max(ipw)), sprintf("%.2f", min(ipw))), "Sec8: IPW")
  tr <- E("MNAR_trait", "ec"); tr0 <- E("MNAR_trait", "ec", 0)
  has(sprintf("the survivor target averages $%s$ against a population shift of $%s$", d3(tr$target_surv), d3(tr$target_pop)), "Sim3: targets")
  ok(tr$target_surv < tr$target_pop, "Sim3: survivors carry smaller shifts")
  has(sprintf("unbiased for the survivor effect $h^{S}$ ($%s$) and biased for the population shift ($%s$, or %d%% of $\\tau(5)$)", d3(tr$bias_own), d3(tr$bias_pop), round(100 * abs(tr$bias_pop) / 0.30 / (1 - exp(-1)))), "Sim3: EC targets")
  has(sprintf("the population bias is $%s$", d3(tr0$bias_pop)), "Sim3: gamma = 0 population bias")
  pt <- one(pr3, regime = "MNAR_trait")
  ok(pt$identical_Y0_entry == TRUE || pt$identical_Y0_entry == "TRUE", "Sim3: common random numbers across gamma blocks")
  has(sprintf("the paired difference of the population biases, $%s$,", d3(pt$ec_bias_pop_g030 - pt$ec_bias_pop_g0)), "Sim3: paired difference")
  ok(pt$mcse_diff_bias_pop < 1e-4, "Sim3: paired difference MC error below .0001")
  has(sprintf("is $%s$ with heterogeneous shifts and $%s$ when the composition term", d3(pt$rej_TSD_g030), d3(pt$rej_TSD_g0)), "Sim3: T_SD by gamma")
  has(sprintf("(paired difference $%s$, Monte Carlo standard error $%s$)", d3(pt$diff_rej), d3(pt$mcse_diff_rej)), "Sim3: T_SD paired difference")
  vt <- one(v3, regime = "MNAR_trait", gamma = 0.3)
  has(sprintf("(for example $%s$ and $%s$ under MNAR on the trait)", d3(vt$cover_ec_if), d3(vt$cover_ec_within)), "Sim3: coverage with and without share terms")
  D1 <- D12[grepl("^D1", D12$part), ]; stopifnot(nrow(D1) == 1)
  has(sprintf("averages $%s$, against an empirical variance of $%s$ and a theoretical value of $%s$, and gives coverage $%s$ and a $\\mathcal{T}_{\\mathrm{NS}}$ rejection rate of $%s$",
              sub("^0", "", sprintf("%.5f", D1$mean_var_if_ec)), sub("^0", "", sprintf("%.5f", D1$var_emp_ec)), sub("^0", "", sprintf("%.5f", D1$var_theory_ec)), d3(D1$cover_ec_if), d3(D1$rej_TNS_if)), "Sim3 D1: influence-function results")
  has(sprintf("averages $%s$ and gives $%s$ and $%s$", sub("^0", "", sprintf("%.5f", D1$mean_var_within_ec)), d3(D1$cover_ec_within), d3(D1$rej_TNS_within)), "Sim3 D1: within-group results")
  rel <- c(abs(D3$mean_se_if_ec / D3$mean_se_boot_ec - 1), abs(D3$mean_se_if_TNS / D3$mean_se_boot_TNS - 1))
  ok(max(rel) < .01, "Sim3 D3: bootstrap and influence-function SEs within 1%")
  ok(nrow(D3) == 3, "Sim3 D3: three blocks")

  ## ---------------- design-rank checks (check_design_rank_output.txt) ----------------
  num <- function(pat) as.integer(sub(pat, "\\1", grep(pat, drk, value = TRUE)[1]))
  ok(any(grepl("all checks passed", drk)), "design rank: all checks passed")
  has(sprintf("on all %s staggered trapezoids", format(num("^ *designs: ([0-9]+) .*$"), big.mark = ",")), "design rank: grid size")
  has(sprintf("On %d of these designs", num("^ *tenure-chain condition holds but nu > d: ([0-9]+) designs.*$")), "design rank: tenure-chain counterexamples")
  has(sprintf("up to eight periods of follow-up, %d attain", num("^ *K >= 4, w < Delta_2 - d, but nu = d: ([0-9]+) designs.*$")), "design rank: K >= 4 designs")
  invisible(NULL)
}

rd <- function(f) read.csv(f, stringsAsFactors = FALSE)
inputs <- list(s1 = rd("sim1_results.csv"), s2 = rd("sim2_results.csv"), p2 = rd("sim2_plim.csv"),
               e3 = rd("output/sim3_estimators.csv"), g3 = rd("output/sim3_diagnostics.csv"), v3 = rd("output/sim3_variance_check.csv"),
               pr3 = rd("output/sim3_pairs.csv"), D12 = rd("output/sim3_partD.csv"), D3 = rd("output/sim3_partD_bootstrap.csv"),
               drk = readLines("check_design_rank_output.txt"))
do.call(run_checks, inputs)
cat(sprintf("%d assertions passed, %d failed\n", n_ok, n_fail))
if (selftest) {
  res <- c()
  for (case in c("perturb", "droprow")) {
    x <- inputs; n_fail <- 0; n_ok <- 0
    if (case == "perturb") x$s2$bias_inc_naive[x$s2$cell == "matched"] <- x$s2$bias_inc_naive[x$s2$cell == "matched"] + 0.01
    if (case == "droprow") x$e3 <- x$e3[-1, ]
    r <- tryCatch({ sink(tempfile()); do.call(run_checks, x); sink(); n_fail > 0 }, error = function(e) { sink(); TRUE })
    res <- c(res, r); cat(sprintf("selftest %-8s: %s\n", case, if (r) "fails as it should" else "DID NOT FAIL"))
  }
  if (!all(res)) quit(status = 1)
}
if (n_fail > 0) quit(status = 1)

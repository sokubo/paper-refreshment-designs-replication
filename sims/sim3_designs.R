## sim3_designs.R -- Simulation 3 of the refreshment-designs paper (Sections 7-8 and the Appendix)
##
## Estimators of the survivors' mean conditioning shift at comparison wave t (continuing cohort observed with
## k = 4 prior interviews; fresh cohort entering at t):
##   naive : continuing respondents at t                 vs. all fresh first-wave respondents
##   sm    : survival matching                           vs. fresh respondents who go on to respond at t+1..t+k
##   ssm   : symmetric survival matching                 continuing arm must also respond at t+1; fresh arm must
##           respond at t+1..t+k+1, so that both arms pass k+1 response decisions, one of them taken at t
##   ec    : entry-wave correction (the observable functional of Das, Toepoel and van Soest 2011, Assumption 3):
##           continuing respondents at t minus the cohort's own entry-wave selection differential, vs. all fresh
##   ipw   : inverse-probability weighting on the observed entry covariate X (a population-effect estimator)
## Diagnostics: T_NS = SM - EC and T_SD = SSM - SM, two-sided at the 5% level.
##
## Data-generating process (designed; no survey microdata):  Y_w(0) = 0.5 X + 0.5 U + eps_w, X ~ Bernoulli(.5),
## U ~ N(0,1), eps_w ~ N(0,1) i.i.d. over waves.  Conditioning shift at tenure s = k + 1:
## tau_i = 0.30 (1 - exp(-k/4)) (1 + gamma U).  The response decision taken at wave w (governing response at w+1)
## has dropout logit  a0 + aX X + b_w U + aE eps_w  with a0 = -2.2.  b_w = aU at every wave of the fresh cohort and
## at the continuing cohort's waves t, t+1, ...; b_w = aU + shift_old at the continuing cohort's waves before t
## (so that shift_old > 0 makes trait selection STRONGER in the continuing cohort's selection period).
##
## Variance estimation.  Every estimator is a difference of subset means within two independent cohorts.  For a
## subset mean Ybar_B = sum B_i Y_i / sum B_i the influence function is B_i (Y_i - mu_B) / P(B); a linear
## combination of subset means within a cohort has the corresponding linear combination of influence functions,
## and its variance is estimated by sum(IF_i^2) / n^2.  Because the influence functions include the estimated
## group shares, the variance of the entry-wave correction contains the share term
##     p (1 - p) (mu_{E,S} - mu_{E,NS})^2 / n_old,
## which the within-group formula used in the previous version of this script omitted; the diagnostics contain
## the analogous share terms.  Both formulas are reported (columns *_within) so that the change can be audited.
##
## Parts
##   A  six attrition regimes x two effect-heterogeneity settings (gamma = 0.30, 0), R replications each.  The two
##      gamma blocks of a regime use the SAME seed (common random numbers): the data are identical except for the
##      conditioning shift, so paired differences between the blocks isolate the effect of heterogeneity.
##   D1 high-separation null experiment: U ~ Bernoulli(.5), stable latent response U + N(0, .1^2) at entry and at t,
##      survival S = U in both cohorts, no conditioning, 500 entrants per cohort, R_D1 replications.
##   D2 exact counterexample: E = S, Y = 1 among survivors, constant fresh outcome; EC - const = p_hat, whose
##      variance is p(1 - p)/n; the within-group formula returns zero.
##   D3 analytic (influence-function) versus person-bootstrap standard errors (persons resampled within cohort,
##      every arm and share recomputed), in D1 and in the MNAR-trait and MNAR-state regimes of Part A.
##
## Usage:  Rscript sim3_designs.R [R] [seed]      (defaults R = 1000, seed = 20260903)
## Output: output/sim3_estimators.csv, sim3_diagnostics.csv, sim3_variance_check.csv, sim3_pairs.csv,
##         sim3_partD.csv, sim3_summary.txt

args <- commandArgs(trailingOnly = TRUE)
R    <- if (length(args) >= 1) as.integer(args[1]) else 1000L
seed <- if (length(args) >= 2) as.integer(args[2]) else 20260903L
R_D1 <- 5000L; R_BOOT <- 200L; B_BOOT <- 200L
dir.create("output", showWarnings = FALSE)
t_start <- Sys.time()

n_old <- 4800; n_new <- 960
tau_inf <- 0.30; kappa <- 4
tau_k <- function(k) tau_inf * (1 - exp(-k / kappa))
a0 <- -2.2
SD_Y0 <- sqrt(0.25 * 0.25 + 0.25 + 1)        # Var(0.5 X + 0.5 U + eps) = 0.0625 + 0.25 + 1 -> SD 1.146

regimes <- list(
  MCAR         = list(aX = 0.0, aU = 0.0, aE = 0.0, shift_old = 0.0),
  MAR_X        = list(aX = 0.8, aU = 0.0, aE = 0.0, shift_old = 0.0),
  MNAR_trait   = list(aX = 0.4, aU = 0.8, aE = 0.0, shift_old = 0.0),
  MNAR_nonstat = list(aX = 0.4, aU = 0.8, aE = 0.0, shift_old = 0.6),
  MNAR_state   = list(aX = 0.4, aU = 0.4, aE = 0.8, shift_old = 0.0),
  MNAR_both    = list(aX = 0.4, aU = 0.8, aE = 0.8, shift_old = 0.6)
)
GAMMAS <- c(0.30, 0.00)

## ---------- data generator ------------------------------------------------------------------------------------
simulate_cohort <- function(n, k_prior, reg, period_shift, k_future, gam) {
  X <- rbinom(n, 1, 0.5); U <- rnorm(n)
  W <- k_prior + 1 + k_future
  eps <- matrix(rnorm(n * W), n, W)
  Y0  <- 0.5 * X + 0.5 * U + eps
  t   <- k_prior + 1
  tau_i <- tau_k(k_prior) * (1 + gam * U)
  Y_t <- Y0[, t] + tau_i
  lp  <- a0 + reg$aX * X + reg$aE * eps
  if (k_prior > 0) lp[, seq_len(k_prior)] <- lp[, seq_len(k_prior), drop = FALSE] + (reg$aU + period_shift) * U
  lp[, t:W] <- lp[, t:W, drop = FALSE] + reg$aU * U
  keep <- matrix(rbinom(n * W, 1, 1 - plogis(lp)), n, W)   # keep[, w] = responds at w + 1 given the decision at w
  S_prior  <- if (k_prior > 0) as.integer(rowSums(keep[, seq_len(k_prior), drop = FALSE]) == k_prior) else rep(1L, n)
  S_next   <- keep[, t]
  S_future <- if (k_future > 0) as.integer(rowSums(keep[, t:(t + k_future - 1), drop = FALSE]) == k_future) else rep(1L, n)
  S_future_m <- if (k_future > 1) as.integer(rowSums(keep[, t:(t + k_future - 2), drop = FALSE]) == k_future - 1) else S_future
  data.frame(X, U, Y0_entry = Y0[, 1], Y_t, tau_i, S_prior, S_next, S_future, S_future_m)
}

## ---------- influence functions --------------------------------------------------------------------------------
if_sub <- function(y, B) {                     # IF of the mean of y over the subset B (B logical, length n)
  out <- numeric(length(B)); pB <- mean(B)
  if (!any(B)) return(out + NA_real_)
  yb <- y[B]; out[B] <- (yb - mean(yb)) / pB; out
}
v_if <- function(IF) sum(IF^2) / length(IF)^2

ipw_mean <- function(y, s, X) {
  pS <- tapply(s, X, mean)[as.character(X)]
  w <- 1 / pmax(pS, 0.02)
  sum(w[s == 1] * y[s == 1]) / sum(w[s == 1])
}
vw <- function(v) if (length(v) >= 2) var(v) / length(v) else NA_real_

estimators <- function(old, new) {
  io <- old$S_prior == 1; is <- io & old$S_next == 1
  nm <- new$S_future_m == 1; m1 <- new$S_future == 1
  all_o <- rep(TRUE, nrow(old)); all_n <- rep(TRUE, nrow(new))
  yo <- old$Y_t[io]; yo_s <- old$Y_t[is]; yn <- new$Y_t; yn_m <- yn[nm]; yn_m1 <- yn[m1]
  naive <- mean(yo) - mean(yn)
  sm    <- mean(yo) - mean(yn_m)
  ssm   <- mean(yo_s) - mean(yn_m1)
  ipw   <- ipw_mean(old$Y_t, old$S_prior, old$X) - mean(yn)
  sel_entry <- mean(old$Y0_entry[io]) - mean(old$Y0_entry)
  ec    <- (mean(yo) - sel_entry) - mean(yn)
  ## influence functions, continuing (o) and fresh (n) cohorts
  IF_yS <- if_sub(old$Y_t, io); IF_ES <- if_sub(old$Y0_entry, io); IF_E <- if_sub(old$Y0_entry, all_o)
  IF_ySS <- if_sub(old$Y_t, is)
  IF_n <- if_sub(yn, all_n); IF_nm <- if_sub(yn, nm); IF_nm1 <- if_sub(yn, m1)
  v_naive <- v_if(IF_yS) + v_if(IF_n)
  v_sm    <- v_if(IF_yS) + v_if(IF_nm)
  v_ssm   <- v_if(IF_ySS) + v_if(IF_nm1)
  v_ec_o  <- v_if(IF_yS - IF_ES + IF_E); v_n <- v_if(IF_n)
  v_ec    <- v_ec_o + v_n
  v_TNS   <- v_if(IF_ES - IF_E) + v_if(IF_nm - IF_n)            # T_NS = (Ebar_S - Ebar) - (Ybar_m - Ybar)
  v_TSD   <- v_if(IF_ySS - IF_yS) + v_if(IF_nm1 - IF_nm)        # T_SD = (Ybar_SS1 - Ybar_S) - (Ybar_m1 - Ybar_m)
  ## within-group formulas of the previous version (no share terms), kept for the audit columns
  es <- old$Y0_entry[io]; ens <- old$Y0_entry[!io]; p <- mean(io)
  v_ec_within <- (var(yo) - 2 * (1 - p) * cov(yo, es) + (1 - p)^2 * var(es)) / sum(io) + (1 - p)^2 * vw(ens) + vw(yn)
  q_m <- mean(nm)
  v_TNS_within <- (1 - p)^2 * (vw(es) + vw(ens)) + (1 - q_m)^2 * (vw(yn[!nm]) + vw(yn_m))
  r_s <- mean(old$S_next[io] == 1); yo_ns <- old$Y_t[io & old$S_next == 0]; q1 <- mean(m1[nm])
  v_TSD_within <- (1 - r_s)^2 * (vw(yo_s) + vw(yo_ns)) + (1 - q1)^2 * (vw(yn[nm & !m1]) + vw(yn_m1))
  c(naive = naive, sm = sm, ssm = ssm, ipw = ipw, ec = ec, TNS = sm - ec, TSD = ssm - sm,
    target_pop = mean(old$tau_i), target_surv = mean(old$tau_i[io]), target_surv_s = mean(old$tau_i[is]),
    se_naive = sqrt(v_naive), se_sm = sqrt(v_sm), se_ssm = sqrt(v_ssm), se_ec = sqrt(v_ec),
    se_TNS = sqrt(v_TNS), se_TSD = sqrt(v_TSD),
    se_ec_within = sqrt(v_ec_within), se_TNS_within = sqrt(v_TNS_within), se_TSD_within = sqrt(v_TSD_within),
    p_surv = p, delta_E = mean(es) - mean(ens))
}
mcse_p <- function(p, R) sqrt(p * (1 - p) / R)

boot_se <- function(old, new, B) {
  bm <- t(replicate(B, {
    e <- estimators(old[sample.int(nrow(old), replace = TRUE), , drop = FALSE], new[sample.int(nrow(new), replace = TRUE), , drop = FALSE])
    e[c("ec", "TNS", "TSD")]
  }))
  apply(bm, 2, sd)
}

## ---------- Part A -------------------------------------------------------------------------------------------
k <- 4
est <- c("naive", "sm", "ssm", "ipw", "ec")
tgt <- c(naive = "target_surv", sm = "target_surv", ssm = "target_surv_s", ipw = "target_pop", ec = "target_surv")
seA <- c(naive = "se_naive", sm = "se_sm", ssm = "se_ssm", ipw = NA, ec = "se_ec")
runs <- list()
for (j in seq_along(regimes)) for (g in GAMMAS) {
  rn <- names(regimes)[j]; reg <- regimes[[rn]]
  set.seed(seed + 1000L * j)                                  # same seed for both gamma blocks: common random numbers
  o <- t(replicate(R, {
    old <- simulate_cohort(n_old, k, reg, reg$shift_old, k_future = 1, gam = g)
    new <- simulate_cohort(n_new, 0, reg, 0, k_future = k + 1, gam = g)
    estimators(old, new)
  }))
  runs[[paste(rn, g, sep = "@")]] <- list(rn = rn, g = g, o = o)
}
summA <- do.call(rbind, lapply(runs, function(x) {
  o <- x$o
  cov <- sapply(est, function(e) if (is.na(seA[e])) NA else mean(abs(o[, e] - o[, tgt[e]]) <= 1.96 * o[, seA[e]]))
  data.frame(regime = x$rn, gamma = x$g, estimator = est, target = unname(tgt[est]),
             bias_own  = sapply(est, function(e) mean(o[, e] - o[, tgt[e]])),
             mcse_bias = sapply(est, function(e) sd(o[, e] - o[, tgt[e]]) / sqrt(R)),
             bias_pop  = sapply(est, function(e) mean(o[, e] - o[, "target_pop"])),
             rmse_own  = sapply(est, function(e) sqrt(mean((o[, e] - o[, tgt[e]])^2))),
             sd_est    = sapply(est, function(e) sd(o[, e])),
             mean_se   = sapply(est, function(e) if (is.na(seA[e])) NA else mean(o[, seA[e]])),
             cover95   = cov, mcse_cover = mcse_p(cov, R),
             target_pop = mean(o[, "target_pop"]), target_surv = mean(o[, "target_surv"]), target_surv_s = mean(o[, "target_surv_s"]),
             R = R, row.names = NULL)
}))
summC <- do.call(rbind, lapply(runs, function(x) {
  o <- x$o
  data.frame(regime = x$rn, gamma = x$g,
             rej_TNS = mean(abs(o[, "TNS"] / o[, "se_TNS"]) > 1.96), rej_TSD = mean(abs(o[, "TSD"] / o[, "se_TSD"]) > 1.96),
             rej_TNS_within = mean(abs(o[, "TNS"] / o[, "se_TNS_within"]) > 1.96), rej_TSD_within = mean(abs(o[, "TSD"] / o[, "se_TSD_within"]) > 1.96),
             mean_TNS = mean(o[, "TNS"]), mean_TSD = mean(o[, "TSD"]), R = R, row.names = NULL)
}))
summC$mcse_rej_TNS <- mcse_p(summC$rej_TNS, R); summC$mcse_rej_TSD <- mcse_p(summC$rej_TSD, R)
vchk <- do.call(rbind, lapply(runs, function(x) {
  o <- x$o
  data.frame(regime = x$rn, gamma = x$g,
             var_emp_ec = var(o[, "ec"]), var_if_ec = mean(o[, "se_ec"]^2), var_within_ec = mean(o[, "se_ec_within"]^2),
             var_emp_TNS = var(o[, "TNS"]), var_if_TNS = mean(o[, "se_TNS"]^2), var_within_TNS = mean(o[, "se_TNS_within"]^2),
             var_emp_TSD = var(o[, "TSD"]), var_if_TSD = mean(o[, "se_TSD"]^2), var_within_TSD = mean(o[, "se_TSD_within"]^2),
             cover_ec_if = mean(abs(o[, "ec"] - o[, "target_surv"]) <= 1.96 * o[, "se_ec"]),
             cover_ec_within = mean(abs(o[, "ec"] - o[, "target_surv"]) <= 1.96 * o[, "se_ec_within"]),
             mean_p_surv = mean(o[, "p_surv"]), mean_delta_E = mean(o[, "delta_E"]), R = R, row.names = NULL)
}))
## paired (common-random-number) comparison of the two gamma blocks: MNAR-trait T_SD rejection and EC population bias
pairs <- do.call(rbind, lapply(names(regimes), function(rn) {
  a <- runs[[paste(rn, 0.3, sep = "@")]]$o; b <- runs[[paste(rn, 0, sep = "@")]]$o
  ra <- abs(a[, "TSD"] / a[, "se_TSD"]) > 1.96; rb <- abs(b[, "TSD"] / b[, "se_TSD"]) > 1.96
  pa <- a[, "ec"] - a[, "target_pop"]; pb <- b[, "ec"] - b[, "target_pop"]
  data.frame(regime = rn, rej_TSD_g030 = mean(ra), rej_TSD_g0 = mean(rb), diff_rej = mean(ra) - mean(rb),
             mcse_diff_rej = sd(ra - rb) / sqrt(R), discordant_pairs = sum(ra != rb),
             ec_bias_pop_g030 = mean(pa), ec_bias_pop_g0 = mean(pb), mcse_diff_bias_pop = sd(pa - pb) / sqrt(R),
             identical_Y0_entry = isTRUE(all.equal(a[, "delta_E"], b[, "delta_E"])), R = R)
}))
write.csv(summA, "output/sim3_estimators.csv", row.names = FALSE)
write.csv(summC, "output/sim3_diagnostics.csv", row.names = FALSE)
write.csv(vchk, "output/sim3_variance_check.csv", row.names = FALSE)
write.csv(pairs, "output/sim3_pairs.csv", row.names = FALSE)
cat(sprintf("Part A done (%.1f min)\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

## ---------- Part D -------------------------------------------------------------------------------------------
make_D1 <- function(n) {
  uo <- rbinom(n, 1, .5); un <- rbinom(n, 1, .5)
  eo <- uo + rnorm(n, sd = .1); en <- un + rnorm(n, sd = .1)
  old <- data.frame(X = 0L, U = uo, Y0_entry = eo, Y_t = eo, tau_i = 0, S_prior = uo, S_next = 1L, S_future = 1L, S_future_m = 1L)
  new <- data.frame(X = 0L, U = un, Y0_entry = en, Y_t = en, tau_i = 0, S_prior = 1L, S_next = 1L, S_future = un, S_future_m = un)
  list(old = old, new = new)
}
set.seed(seed + 50000L)
d1 <- t(replicate(R_D1, { dd <- make_D1(500); estimators(dd$old, dd$new) }))
D1 <- data.frame(part = "D1 high-separation null (n = 500 per cohort)", R = R_D1,
                 var_theory_ec = 2 * (.25 + .1^2) / 500, var_emp_ec = var(d1[, "ec"]),
                 mean_var_if_ec = mean(d1[, "se_ec"]^2), mean_var_within_ec = mean(d1[, "se_ec_within"]^2),
                 cover_ec_if = mean(abs(d1[, "ec"]) <= 1.96 * d1[, "se_ec"]), cover_ec_within = mean(abs(d1[, "ec"]) <= 1.96 * d1[, "se_ec_within"]),
                 rej_TNS_if = mean(abs(d1[, "TNS"] / d1[, "se_TNS"]) > 1.96), rej_TNS_within = mean(abs(d1[, "TNS"] / d1[, "se_TNS_within"]) > 1.96))
D1$mcse_cover_ec_if <- mcse_p(D1$cover_ec_if, R_D1); D1$mcse_rej_TNS_if <- mcse_p(D1$rej_TNS_if, R_D1)

## D2: exact counterexample, one realisation and the Monte Carlo variance of p_hat
set.seed(seed + 60000L)
n2 <- 1000
d2 <- t(replicate(2000, {
  S <- rbinom(n2, 1, .5)
  old <- data.frame(X = 0L, U = 0, Y0_entry = S, Y_t = ifelse(S == 1, 1, NA), tau_i = 0, S_prior = S, S_next = 1L, S_future = 1L, S_future_m = 1L)
  new <- data.frame(X = 0L, U = 0, Y0_entry = 0, Y_t = rep(0, n2), tau_i = 0, S_prior = 1L, S_next = 1L, S_future = 1L, S_future_m = 1L)
  e <- estimators(old, new); c(ec = e[["ec"]], p = mean(S), v_if = e[["se_ec"]]^2, v_within = e[["se_ec_within"]]^2)
}))
stopifnot(max(abs(d2[, "ec"] - d2[, "p"])) < 1e-12, max(abs(d2[, "v_within"])) < 1e-12,
          max(abs(d2[, "v_if"] - d2[, "p"] * (1 - d2[, "p"]) / n2)) < 1e-12)
D2 <- data.frame(part = "D2 exact counterexample (E = S, Y = 1 among survivors, constant fresh; n = 1000)", R = 2000,
                 var_theory_ec = .25 / n2, var_emp_ec = var(d2[, "ec"]), mean_var_if_ec = mean(d2[, "v_if"]), mean_var_within_ec = 0,
                 cover_ec_if = mean(abs(d2[, "ec"] - .5) <= 1.96 * sqrt(d2[, "v_if"])), cover_ec_within = mean(abs(d2[, "ec"] - .5) <= 0),
                 rej_TNS_if = NA, rej_TNS_within = NA, mcse_cover_ec_if = NA, mcse_rej_TNS_if = NA)

## D3: influence-function versus person-bootstrap standard errors
D3 <- list()
boot_block <- function(label, gen, Rb, B) {
  z <- t(replicate(Rb, {
    dd <- gen(); e <- estimators(dd$old, dd$new); sb <- boot_se(dd$old, dd$new, B)
    c(ec = e[["ec"]], TNS = e[["TNS"]], TSD = e[["TSD"]], tgt = e[["target_surv"]],
      se_ec = e[["se_ec"]], se_ec_within = e[["se_ec_within"]], se_ec_boot = sb[["ec"]],
      se_TNS = e[["se_TNS"]], se_TNS_within = e[["se_TNS_within"]], se_TNS_boot = sb[["TNS"]])
  }))
  data.frame(block = label, R = Rb, B = B,
             sd_emp_ec = sd(z[, "ec"]), mean_se_if_ec = mean(z[, "se_ec"]), mean_se_boot_ec = mean(z[, "se_ec_boot"]), mean_se_within_ec = mean(z[, "se_ec_within"]),
             cover_ec_if = mean(abs(z[, "ec"] - z[, "tgt"]) <= 1.96 * z[, "se_ec"]), cover_ec_boot = mean(abs(z[, "ec"] - z[, "tgt"]) <= 1.96 * z[, "se_ec_boot"]),
             sd_emp_TNS = sd(z[, "TNS"]), mean_se_if_TNS = mean(z[, "se_TNS"]), mean_se_boot_TNS = mean(z[, "se_TNS_boot"]), mean_se_within_TNS = mean(z[, "se_TNS_within"]))
}
set.seed(seed + 70000L)
D3[[1]] <- boot_block("D1 high-separation null", function() make_D1(500), R_BOOT, B_BOOT)
for (rn in c("MNAR_trait", "MNAR_state")) {
  reg <- regimes[[rn]]
  D3[[length(D3) + 1]] <- boot_block(paste0("Part A ", rn, ", gamma = 0.3"), function() list(
    old = simulate_cohort(n_old, k, reg, reg$shift_old, k_future = 1, gam = 0.3),
    new = simulate_cohort(n_new, 0, reg, 0, k_future = k + 1, gam = 0.3)), R_BOOT, B_BOOT)
}
D3 <- do.call(rbind, D3)
partD <- list(D12 = rbind(D1, D2), D3 = D3)
write.csv(partD$D12, "output/sim3_partD.csv", row.names = FALSE)
write.csv(partD$D3, "output/sim3_partD_bootstrap.csv", row.names = FALSE)

sink("output/sim3_summary.txt")
cat("Simulation 3 (R =", R, ", seed =", seed, "; R_D1 =", R_D1, "; bootstrap R =", R_BOOT, "x B =", B_BOOT, ")\n")
cat("tau(4) =", round(tau_k(4), 4), "outcome units (untreated SD", round(SD_Y0, 3), "-> tau(4) =", round(tau_k(4) / SD_Y0, 3), "SD)\n\n")
cat("Part A: estimators (bias against own target; coverage with influence-function standard errors)\n")
print(summA, digits = 3, row.names = FALSE)
cat("\nPart A: diagnostics (rejection at 5%; influence-function denominators, and within-group denominators of the previous version)\n")
print(summC, digits = 3, row.names = FALSE)
cat("\nPart A: variance check\n"); print(vchk, digits = 3, row.names = FALSE)
cat("\nPart A: common-random-number pairs of the two gamma blocks\n"); print(pairs, digits = 3, row.names = FALSE)
cat("\nPart D1-D2\n"); print(partD$D12, digits = 4, row.names = FALSE)
cat("\nPart D3: influence-function vs person-bootstrap standard errors\n"); print(partD$D3, digits = 4, row.names = FALSE)
sink()
cat("run time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n")
cat("done\n")

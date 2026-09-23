# ============================================================
# Simulation 2 of the refreshment-designs paper: Theorem 4 (two refreshments, negative-control correction)
# and Corollary 7 (generational refreshment).  Version 3.
#
# Data-generating mechanism.  Three cohorts enter at waves 1, 5, 13 (n = 3,000 each) and are compared at wave 13
# (tenures 13, 9, 1).  Person trait u ~ N(0,1); per-wave retention plogis(1.6 + 0.6 u), so cohort 1 has passed 12
# response decisions at wave 13, cohort 5 eight and cohort 13 none.  Target item y = g_e + tau(s) + u + e,
# e ~ N(0,1), tau(s) = 0.15 (1 - exp(-(s-1)/4)).  Negative-control battery of K = 23 items measured once at entry,
# x_k = 0.6 u + eta_k, eta_k ~ N(0,1).  Generational cell: cohort 13 carries g = 0.35 on the target item (not on the
# battery).  Contaminated cell: three battery items of cohort 1 carry recall drift +0.15.  300 replications.
#
# Estimands (population-SD units; the target item's population SD is sqrt(2)):
#   incumbent increment  (tau(13) - tau(9)) / sqrt(2) = 0.15 (e^-2 - e^-3) / sqrt(2) = 0.00907376
#   level contrast       (tau(13) - tau(1)) / sqrt(2) = 0.15 (1 - e^-3) / sqrt(2)    = 0.10078530
#
# Methods.  All corrections are computed on the RAW scale and divided by the fixed population SD sqrt(2)
# afterwards, so that the selection transport (A2) holds exactly on the scale on which it is applied:
#   naive          raw contrast of survivors' means / sqrt(2)
#   true loading   raw contrast - (1/0.6) x mean_k(raw negative-control contrast); the raw loading 1/0.6 is the
#                  ratio of the target's to the battery's loading on u, i.e. exact (A2)
#   unit loading   raw contrast - (sqrt(2)/sqrt(1.36)) x mean_k(raw negative-control contrast): a unit loading on the
#                  population-standardised scale (population SDs sqrt(2) and sqrt(1.36)), which under-corrects by
#                  the factor 0.6 sqrt(2)/sqrt(1.36) = 0.7276
#   screened       true loading, with items flagged by the per-item outlier screen removed before pooling
#   median-pooled  true loading, with the median instead of the mean of the item contrasts
# The previous version (v2) standardised every contrast by the pair-specific pooled SD of the two survivor groups
# and applied a loading fixed at the population-SD ratio 1.374369.  Because selection changes the target's and the
# battery's SDs differently, that loading is not the one for the pair-specific denominators; the population limits
# printed at the end reproduce its asymptotic residual biases (+0.006318 for the increment, +0.015550 for the level).
#
# Battery homogeneity: Q statistic on the standardised item contrasts (chi-square with K - 1 df); per-item outlier
# flag |z| > 2.5 with z centred at the median item contrast (a screen; the corrections above pool by the mean
# unless stated).
# Run: Rscript sim2_multicohort.R     (writes sim2_results.csv and sim2_plim.csv next to this script)
# ============================================================
set.seed(20260825)
R_REPS <- 300; N_C <- 3000; K <- 23; LAM <- 0.6
tau <- function(s) 0.15 * (1 - exp(-(s - 1) / 4))
SD_Y <- sqrt(2); SD_X <- sqrt(LAM^2 + 1)
G_TRUE <- 1 / LAM                       # raw-scale loading (exact A2)
G_UNIT <- SD_Y / SD_X                   # unit loading on the population-standardised scale, expressed on the raw scale
G_V2   <- (1 / sqrt(2)) / (LAM / sqrt(LAM^2 + 1))   # the v2 loading for pair-standardised contrasts (1.374369)
DINC_TRUE <- (tau(13) - tau(9)) / SD_Y; LVL_TRUE <- (tau(13) - tau(1)) / SD_Y
stopifnot(abs(DINC_TRUE - 0.00907376) < 5e-9, abs(LVL_TRUE - 0.1007853) < 5e-8)

one_rep <- function(generational = FALSE, nc_drift = FALSE) {
  res <- list()
  for (cohort in c(1, 5, 13)) {
    u <- rnorm(N_C)
    stay <- rep(TRUE, N_C)
    for (w in seq_len(13 - cohort)) stay <- stay & (runif(N_C) < plogis(1.6 + 0.6 * u))
    g <- if (generational && cohort == 13) 0.35 else 0
    x <- matrix(LAM * u, N_C, K) + matrix(rnorm(N_C * K), N_C, K)
    if (nc_drift && cohort == 1) x[, 1:3] <- x[, 1:3] + 0.15          # recall drift in three items
    res[[as.character(cohort)]] <- list(u = u, stay = stay, g = g, x = x)
  }
  yobs <- function(cc, s) { r <- res[[as.character(cc)]]; (r$g + tau(s) + r$u + rnorm(N_C))[r$stay] }
  xobs <- function(cc) { r <- res[[as.character(cc)]]; r$x[r$stay, , drop = FALSE] }
  std <- function(a, b) (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
  yC <- yobs(1, 13); yD <- yobs(5, 9); yE <- yobs(13, 1)            # same draw order as v2
  xC <- xobs(1); xD <- xobs(5); xE <- xobs(13)
  ## raw contrasts
  c_inc <- mean(yC) - mean(yD); c_lvl <- mean(yC) - mean(yE)
  n_inc <- colMeans(xC) - colMeans(xD); n_lvl <- colMeans(xC) - colMeans(xE)
  ## standardised item contrasts for the homogeneity test and the screen
  nk <- vapply(1:K, function(k) std(xC[, k], xD[, k]), 0)
  se_k <- sqrt(1 / nrow(xC) + 1 / nrow(xD))
  Q <- sum(((nk - mean(nk)) / se_k)^2)
  zk <- (nk - median(nk)) / se_k
  keep <- abs(zk) <= 2.5
  ## v2 method (pair-standardised contrasts, loading 1.374369), for the audit columns
  C_inc_std <- std(yC, yD); C_lvl_std <- std(yC, yE)
  nk_lvl <- vapply(1:K, function(k) std(xC[, k], xE[, k]), 0)
  c(inc_naive  = c_inc / SD_Y,
    inc_true   = (c_inc - G_TRUE * mean(n_inc)) / SD_Y,
    inc_unit   = (c_inc - G_UNIT * mean(n_inc)) / SD_Y,
    inc_screen = (c_inc - G_TRUE * mean(n_inc[keep])) / SD_Y,
    inc_median = (c_inc - G_TRUE * median(n_inc)) / SD_Y,
    lvl_true   = (c_lvl - G_TRUE * mean(n_lvl)) / SD_Y,
    inc_v2     = C_inc_std - G_V2 * mean(nk), lvl_v2 = C_lvl_std - G_V2 * mean(nk_lvl),
    Q = Q, reject = as.numeric(Q > qchisq(.95, K - 1)),
    item_flag = as.numeric(any(abs(zk[1:3]) > 2.5)), false_flag = mean(abs(zk[4:K]) > 2.5), n_screened = sum(!keep))
}

cells <- list(
  matched      = list(generational = FALSE, nc_drift = FALSE),
  generational = list(generational = TRUE,  nc_drift = FALSE),
  drift3       = list(generational = FALSE, nc_drift = TRUE))
res <- do.call(rbind, lapply(names(cells), function(nm) {
  cfg <- cells[[nm]]
  m <- replicate(R_REPS, one_rep(cfg$generational, cfg$nc_drift))
  b <- function(v, truth) mean(m[v, ]) - truth
  mc <- function(v) sd(m[v, ]) / sqrt(R_REPS)
  data.frame(cell = nm,
    bias_inc_naive  = b("inc_naive", DINC_TRUE),
    bias_inc_unitG  = b("inc_unit", DINC_TRUE),
    bias_inc_trueG  = b("inc_true", DINC_TRUE), mcse_inc_trueG = mc("inc_true"),
    rmse_inc_trueG  = sqrt(mean((m["inc_true", ] - DINC_TRUE)^2)),
    bias_inc_screen = b("inc_screen", DINC_TRUE), bias_inc_median = b("inc_median", DINC_TRUE),
    bias_lvl_trueG  = b("lvl_true", LVL_TRUE), mcse_lvl_trueG = mc("lvl_true"),
    bias_inc_v2 = b("inc_v2", DINC_TRUE), bias_lvl_v2 = b("lvl_v2", LVL_TRUE),
    Q_reject_rate  = mean(m["reject", ]),
    item_flag_rate = mean(m["item_flag", ]), false_flag_rate = mean(m["false_flag", ]), mean_n_screened = mean(m["n_screened", ]),
    R = R_REPS)
}))

## population limits (deterministic): selection moments of u after k decisions
f <- function(k, j) integrate(function(u) u^j * dnorm(u) * plogis(1.6 + .6 * u)^k, -Inf, Inf, rel.tol = 1e-11)$value
mom <- sapply(c(12, 8, 0), function(k) { pr <- f(k, 0); mu <- f(k, 1) / pr; c(p = pr, mu = mu, V = f(k, 2) / pr - mu^2) })
colnames(mom) <- c("cohort_1", "cohort_5", "cohort_13")
plim_row <- function(pair, gen) {                         # pair 2 = increment (cohort 1 vs 5), 3 = level (cohort 1 vs 13)
  dmu <- mom["mu", 1] - mom["mu", pair]
  ttrue <- if (pair == 2) tau(13) - tau(9) else tau(13) - tau(1)
  gshift <- if (pair == 3 && gen) -0.35 else 0
  raw <- ttrue + gshift + dmu; nraw <- LAM * dmu
  vy <- 1 + mean(mom["V", c(1, pair)]); vx <- LAM^2 * mean(mom["V", c(1, pair)]) + 1
  data.frame(contrast = if (pair == 2) "increment" else "level", generational = gen, target = ttrue / SD_Y,
             plim_naive = raw / SD_Y, plim_trueG = (raw - G_TRUE * nraw) / SD_Y, plim_unitG = (raw - G_UNIT * nraw) / SD_Y,
             plim_v2 = raw / sqrt(vy) - G_V2 * nraw / sqrt(vx), gamma_pair_specific = sqrt(vx) / (LAM * sqrt(vy)))
}
plim <- rbind(plim_row(2, FALSE), plim_row(3, FALSE), plim_row(3, TRUE))
plim$asy_bias_trueG <- plim$plim_trueG - plim$target; plim$asy_bias_unitG <- plim$plim_unitG - plim$target
plim$asy_bias_v2 <- plim$plim_v2 - plim$target; plim$asy_bias_naive <- plim$plim_naive - plim$target
## assertions: exact transport on the raw scale; the v2 residuals quoted in the header; Monte Carlo agreement
stopifnot(abs(plim$asy_bias_trueG[1:2]) < 1e-12,
          abs(plim$asy_bias_v2[1] - 0.006318) < 5e-6, abs(plim$asy_bias_v2[2] - 0.015550) < 5e-6,
          abs(plim$gamma_pair_specific[1] - 1.423883) < 5e-6, abs(plim$gamma_pair_specific[2] - 1.400333) < 5e-6,
          abs(res$bias_inc_trueG[1:2]) < 4 * res$mcse_inc_trueG[1:2], abs(res$bias_lvl_trueG[1]) < 4 * res$mcse_lvl_trueG[1])

cat(sprintf("true (population-SD units): tau(13)-tau(9) = %.8f, tau(13)-tau(1) = %.8f\n", DINC_TRUE, LVL_TRUE))
print(res, digits = 4)
cat("\npopulation limits\n"); print(mom, digits = 10); print(plim, digits = 7)
write.csv(res, "sim2_results.csv", row.names = FALSE)
write.csv(plim, "sim2_plim.csv", row.names = FALSE)

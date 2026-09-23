# ============================================================
# Simulation 1 of the refreshment-designs paper (Theorem 1, Corollary 2, Theorem 3): the plug-in funnel.  Version 2.
#  DGM: (Y1, Y2*) bivariate normal, means 0, variances 1, correlation .5; outcome selection pi = plogis(a + beta Y2*),
#       beta in {0, .3, .6}, a tuned to retention p in {.6, .8}; stayers report Y2 = Y2* + h, h in {0, .2};
#       n = 4,000 plus a refreshment sample of 2,000 reporting Y2* directly; 200 replications per cell.
#  Estimand: the location shift h.
#  Method (a RELAXED plug-in procedure, not the population identified set):
#       Gaussian kernel densities (bw.nrd0), evaluated exactly on a 512-point grid; candidate set {h' : sup_y p_hat q_hat(y + h') / f2_hat(y) <= 1 + kappa} with the supremum
#       taken only over the trimmed support {f2_hat >= 0.10 max f2_hat} and slack kappa = 0.10; candidate grid
#       h' in [-1, 1.5] step .01 (v1 used [-0.5, 0.7], which capped the upper end in two cells); coverage counts an
#       endpoint equal to the truth as covering it (tolerance 1e-9).
#       Naive comparator: quantile alignment ignoring attrition, median over quantiles .10(.05).90.
#  Performance: coverage of the true h by the relaxed plug-in set, mean width and endpoints, naive bias.
#  Population benchmarks (deterministic; no sampling): the population identified set of Theorem 3 for this DGM,
#       h + [0, Delta*] (a point when beta = 0), and the population analogue of the relaxed procedure (the same
#       trim and slack applied to the true densities, without smoothing), so that the plug-in widths can be compared
#       with what identification alone would give.
# Run: Rscript sim1_bounds.R      (no confidential data; writes sim1_results.csv next to this script)
# ============================================================
set.seed(20260825)
R_REPS <- 200; N <- 4000; N_R <- 2000; RHO <- 0.5; KAPPA <- 0.10
HGRID <- (-100:150) / 100          # k/100 exactly as the literals 0 and 0.2 are represented (v1 used seq(), whose
                                   # grid point near 0.2 was 0.2 + 1e-16 and so excluded the truth when it was an endpoint)

## Gaussian kernel density with the bw.nrd0 bandwidth, evaluated EXACTLY on 512 points of [-4, 4.7]
## (v1 used stats::density, whose FFT binning changed in R 4.4.0 -- see its `old.coords` argument -- so that
## R 4.3 and R >= 4.4 gave plug-in sets differing in a few replications; the exact sum is version-independent)
kde_grid <- function(x, from = -4, to = 4.7, n = 512) {
  g <- seq(from, to, length.out = n); bw <- bw.nrd0(x)
  list(x = g, y = vapply(g, function(z) mean(dnorm((z - x) / bw)) / bw, 0))
}

run_cell <- function(beta, h_true, p_target) {
  tune_a <- function() {
    fn <- function(a) { y2 <- rnorm(200000); mean(plogis(a + beta * y2)) - p_target }
    uniroot(fn, c(-6, 6))$root
  }
  a <- tune_a()
  out <- replicate(R_REPS, {
    y1 <- rnorm(N); y2 <- RHO * y1 + sqrt(1 - RHO^2) * rnorm(N)
    S <- rbinom(N, 1, plogis(a + beta * y2)) == 1
    rep_y2 <- y2[S] + h_true                      # stayers' conditioned reports
    ref <- rnorm(N_R); ref <- RHO * rnorm(N_R) + sqrt(1 - RHO^2) * rnorm(N_R) # fresh Y2*
    p_hat <- mean(S)
    dq <- kde_grid(rep_y2); df2 <- kde_grid(ref)
    f2_at <- approxfun(df2$x, df2$y, yleft = 0, yright = 0)
    q_at  <- approxfun(dq$x, dq$y, yleft = 0, yright = 0)
    ygrid <- df2$x[df2$y >= 0.10 * max(df2$y)]    # trimmed support of f2
    feas <- vapply(HGRID, function(hp) {
      r <- p_hat * q_at(ygrid + hp) / f2_at(ygrid)
      max(r, na.rm = TRUE) <= 1 + KAPPA
    }, TRUE)
    if (!any(feas)) return(c(NA, NA, NA, NA))
    hL <- min(HGRID[feas]); hU <- max(HGRID[feas])
    # 素朴: 分位整列(選択を無視した点推定)
    qs <- seq(.1, .9, .05)
    h_naive <- median(quantile(rep_y2, qs) - quantile(ref, qs))
    c(hL, hU, h_naive, p_hat)
  })
  ok <- !is.na(out[1, ])
  data.frame(beta = beta, h_true = h_true, p_target = p_target,
             n_ok = sum(ok),
             coverage = mean(out[1, ok] <= h_true + 1e-9 & h_true - 1e-9 <= out[2, ok]),
             coverage_v1_rule = mean(out[1, ok] <= h_true & h_true <= out[2, ok]),
             mean_width = mean(out[2, ok] - out[1, ok]),
             mean_hL = mean(out[1, ok]), mean_hU = mean(out[2, ok]),
             naive_bias = mean(out[3, ok]) - h_true,
             theory_width_note = beta * 1)      # βσ²(σ=1)
}

grid <- expand.grid(beta = c(0, 0.3, 0.6), h_true = c(0, 0.2), p_target = c(0.6, 0.8))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  run_cell(grid$beta[i], grid$h_true[i], grid$p_target[i])))
## ---------- population benchmarks (no random numbers are drawn below) ----------
a_exact <- function(beta, p) uniroot(function(a) integrate(function(y) plogis(a + beta * y) * dnorm(y), -Inf, Inf, rel.tol = 1e-10)$value - p,
                                     c(-10, 10), tol = 1e-12)$root
## log of p q(y + h') / f2(y) with eps = h' - h:  log plogis(a + beta (y + eps)) - eps y - eps^2 / 2
logratio_max <- function(a, beta, eps, ylim = c(-60, 60)) {
  f <- function(y) plogis(a + beta * (y + eps), log.p = TRUE) - eps * y - eps^2 / 2
  if (beta == 0 && eps != 0) return(Inf)                       # unbounded in one tail
  if (eps < 0 || eps > beta) return(Inf)                         # unbounded in the upper / lower tail
  optimize(f, ylim, maximum = TRUE, tol = 1e-12)$objective
}
pop_set <- function(beta, p) {                                   # identified set of eps = h' - h (Theorem 3: an interval)
  if (beta == 0) return(c(0, 0))
  a <- a_exact(beta, p)
  up <- uniroot(function(e) logratio_max(a, beta, e) - 0, c(1e-9, beta - 1e-9), tol = 1e-12)$root
  if (logratio_max(a, beta, beta) <= 0) up <- beta
  c(0, up)
}
relaxed_set <- function(beta, p, kappa = KAPPA, Y = sqrt(2 * log(10))) {   # trim |y| <= Y (10% of the mode), slack kappa
  a <- a_exact(beta, p); yy <- seq(-Y, Y, length.out = 4001); eg <- seq(-1.5, 2, by = 0.0005)
  ok <- vapply(eg, function(e) max(plogis(a + beta * (yy + e), log.p = TRUE) - e * yy - e^2 / 2) <= log(1 + kappa), TRUE)
  range(eg[ok])
}
pb <- t(mapply(function(beta, p) c(pop_set(beta, p), relaxed_set(beta, p)), grid$beta, grid$p_target))
res$pop_hL <- res$h_true + pb[, 1]; res$pop_hU <- res$h_true + pb[, 2]; res$pop_width <- pb[, 2] - pb[, 1]
res$relaxed_hL <- res$h_true + pb[, 3]; res$relaxed_hU <- res$h_true + pb[, 4]; res$relaxed_width <- pb[, 4] - pb[, 3]
stopifnot(all(res$pop_width <= res$theory_width_note + 1e-9))   # Delta* <= beta sigma^2 (Corollary 2 worked example)

print(res, digits = 3)
write.csv(res, "sim1_results.csv", row.names = FALSE)

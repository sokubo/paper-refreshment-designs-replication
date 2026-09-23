# ============================================================
# Deterministic checks of the worked examples and counterexamples of Sections 3-6 of the refreshment-designs
# paper.  Version 3 (adds block 12; version 2 added blocks 5-11).  No randomness; grids in log space; stops on any failed assertion.
#   (1) f2 = N(0,1), Q = N(-1,1), p = .5: identified set {-1} (Theorem 3's example of a degenerate interval).
#   (2) Gaussian-logistic example of Corollary 2: sigma^2 = beta = 1, alpha = -1 (retention .30): Delta* = 1 attained.
#   (3) Corollary 3(i): outcome-independent attrition, f2 = Q = N(0,1): {0} at every p < 1 (the constraint fails only
#       at |y| > log(1/p)/h', so a finite grid shows it only for h' beyond log(1/p)/40; cf. Remark 2 on trimming).
#   (4) Corollary 3(ii)-(iii): Laplace and logistic R(eps) = e^|eps| (half-width log(1/p)); Cauchy closed form.
#   (5) Central binding (Theorem 3, Corollary 2): f2 = N(0,1), Q = N(0, s^2), p = s < 1.  The implied selection
#       exp{-(s^-2 - 1) y^2 / 2} decays in both tails and the stayers' density has lighter tails, yet the ratio at
#       every h' != 0 is exp{h'^2 / (2(1 - s^2))} > 1: the identified set is {0}.
#   (6) Tail class does not determine the set: f(x) propto (1 + x^2)^-1 exp{2 e^{-x^2} cos(20 pi x)} (Cauchy tails)
#       and the same modulation of the logistic density (exponential tails).  At p = .5, shift .1 is feasible and .05
#       is not: the identified set is disconnected.  For the modulated logistic, R(eps) >= e^|eps| (outer bound of
#       Corollary 3(ii)) but R(eps) != e^|eps|.
#   (7) Corollary 3(iii)/Corollary 2(iii): if log f2 is L-Lipschitz, every |eps| <= L^-1 log(1/sup pi_bar) is feasible.
#   (8) Corollary 4: the image of a feasible map need only CONTAIN the reported categories; a non-identity map is
#       point identified when a category goes unreported.
#   (9) Corollary 5: under stationary selection a general map can be excluded although it satisfies both episodes'
#       domination constraints (c'(y) = 2y); in the location family the restriction does not bind.
#  (10) Proposition 1: a covariate that predicts retention only tightens the set (Laplace strata, p(x) = .8, .2).
#  (11) Theorem 4(a) with a bounded common loading: episode constraints must be intersected jointly in Gamma.
#  (12) Theorem 4(a) under marginal versus joint target-battery information: with a battery that predicts retention,
#       the loading set computed from the margins is an outer bound for the joint data (binary X, Laplace Y*).
# ============================================================
ys <- seq(-40, 40, by = 0.002)
lsup <- function(v) max(v)                                         # log sup over the grid
ok <- function(cond, msg) { cat(sprintf("  %s  %s\n", if (cond) "ok  " else "FAIL", msg)); if (!cond) stop(msg) }

cat("(1) f2 = N(0,1), Q = N(-1,1), p = .5\n")
r1 <- sapply(c(-1.2, -1.05, -1, -0.95, -0.8), function(hp) log(0.5) + lsup(dnorm(ys + hp, -1, 1, log = TRUE) - dnorm(ys, log = TRUE)))
ok(all(r1[-3] > 0) && abs(r1[3] - log(0.5)) < 1e-12, "only h' = -1 is feasible on the grid: identified set {-1}")

cat("(2) Gaussian-logistic, alpha = -1, beta = 1, sigma^2 = 1\n")
alpha <- -1; beta <- 1
p_ret <- integrate(function(y) plogis(alpha + beta * y) * dnorm(y), -Inf, Inf)$value
lr <- function(D) plogis(alpha + beta * (ys + D), log.p = TRUE) + dnorm(ys + D, log = TRUE) - dnorm(ys, log = TRUE)
ok(abs(p_ret - 0.3) < 0.01, sprintf("retention %.4f", p_ret))
ok(lsup(lr(1)) <= 1e-12 && lsup(lr(1.05)) > 0 && lsup(lr(-0.02)) > 0, "Delta = 1 = beta sigma^2 feasible; 1.05 and -0.02 not (lower- and upper-tail violations visible on |y| <= 40)")
ok(abs(exp(optimize(function(y) plogis(alpha + beta * (y + 1), log.p = TRUE) + dnorm(y + 1, log = TRUE) - dnorm(y, log = TRUE), c(-20, 20), maximum = TRUE)$objective) - exp(-1/2)) < 1e-6,
   "interior supremum at Delta = 1 equals e^{-1/2}")

cat("(3) outcome-independent attrition, f2 = Q = N(0,1)\n")
for (p in c(0.3, 0.5, 0.9)) {
  h <- 2 * log(1 / p) / 40
  ok(log(p) + lsup(dnorm(ys + h, log = TRUE) - dnorm(ys, log = TRUE)) > 0 && log(p) + lsup(dnorm(ys - h, log = TRUE) - dnorm(ys, log = TRUE)) > 0,
     sprintf("p = %.1f: h' = +-%.3f infeasible on |y| <= 40 (the ratio is exp{-h' y - h'^2/2}, unbounded)", p, h))
}

cat("(4) Laplace, logistic and Cauchy\n")
dlap <- function(y) -abs(y) - log(2)
hw <- function(logf, p, grid) { f <- sapply(grid, function(h) log(p) + lsup(logf(ys + h) - logf(ys)) <= 1e-9); max(grid[f]) }
for (p in c(0.3, 0.5, 0.8)) {
  hL <- hw(dlap, p, seq(0, 3, by = 0.01)); hG <- hw(function(y) dlogis(y, log = TRUE), p, seq(0, 3, by = 0.01))
  hC <- hw(function(y) dcauchy(y, log = TRUE), p, seq(0, 5, by = 0.01))
  cat(sprintf("   p = %.1f: Laplace %.2f, logistic %.2f (theory log(1/p) = %.2f) | Cauchy %.2f (theory (1-p)/sqrt(p) = %.2f)\n", p, hL, hG, log(1 / p), hC, (1 - p) / sqrt(p)))
  ok(abs(hL - log(1 / p)) < 0.011 && abs(hG - log(1 / p)) < 0.011 && abs(hC - (1 - p) / sqrt(p)) < 0.02, sprintf("p = %.1f half-widths agree with the closed forms", p))
}
for (e in c(0.3, 1, 2)) ok(abs(lsup(dlogis(ys + e, log = TRUE) - dlogis(ys, log = TRUE)) - e) < 1e-3, sprintf("logistic: log R(%.1f) = %.1f (1-Lipschitz log density with tail slopes -+1)", e, e))

cat("(5) central binding: f2 = N(0,1), Q = N(0, s^2), p = s\n")
for (s in c(0.5, 0.8)) {
  pibar <- exp(-(s^-2 - 1) * ys^2 / 2)
  ok(max(pibar) <= 1 && abs(s * dnorm(0, 0, s) / dnorm(0) - 1) < 1e-12, sprintf("s = %.1f: implied selection exp{-(s^-2-1)y^2/2} is admissible, sup = 1 at y = 0", s))
  for (hp in c(0.01, 0.1, 1)) {
    num <- log(s) + lsup(dnorm(ys + hp, 0, s, log = TRUE) - dnorm(ys, log = TRUE)); an <- hp^2 / (2 * (1 - s^2))
    ok(abs(num - an) < 1e-5 && an > 0, sprintf("s = %.1f, h' = %.2f: log sup ratio %.6f = h'^2/(2(1-s^2)) > 0: infeasible", s, hp, num))
  }
}
cat("   -> identified set {0} although selection decays in both tails and Q has strictly lighter tails than F2\n")

cat("(6) modulated Cauchy and logistic densities (p = .5)\n")
mod <- function(x) 2 * exp(-x^2) * cos(20 * pi * x)
lg_c <- function(x) -log1p(x^2) + mod(x); lg_l <- function(x) dlogis(x, log = TRUE) + mod(x)
## analytic bound for eps = .1 (one period of the modulation): the Cauchy and logistic log densities are 1-Lipschitz, and
## m(x + .1) - m(x) = 2 cos(20 pi x) (e^{-(x+.1)^2} - e^{-x^2}) is at most 2 x .1 x sqrt(2/e) in absolute value
bound <- 0.1 * (1 + 2 * sqrt(2 / exp(1)))
for (nm in c("Cauchy", "logistic")) {
  lg <- if (nm == "Cauchy") lg_c else lg_l
  R01 <- lsup(lg(ys + 0.1) - lg(ys)); R005 <- lsup(lg(ys + 0.05) - lg(ys)); w <- lg(0.1) - lg(0.05)
  cat(sprintf("   %s tails: log R(.1) = %.4f, log R(.05) >= %.4f (at y = .05), log 2 = %.4f\n", nm, R01, w, log(2)))
  ok(R01 <= bound + 1e-9 && bound < log(2) && w > log(2) && R005 >= w - 1e-9, sprintf("%s: .1 feasible, .05 infeasible -> disconnected identified set", nm))
}
ok(all(sapply(c(0.5, 1, 2), function(e) lsup(lg_l(ys + e) - lg_l(ys)) >= e - 1e-6)),
   "modulated logistic: log R(eps) >= |eps| (exponential-tail outer bound holds)")
ok(abs(lsup(lg_l(ys + 1) - lg_l(ys)) - 1) > 0.05, "modulated logistic: log R(1) != 1 (the closed form e^|eps| needs the Lipschitz condition)")
fs <- seq(0, 0.2, by = 0.0025); feas <- sapply(fs, function(e) log(0.5) + lsup(lg_c(ys + e) - lg_c(ys)) <= 0)
cat("   modulated Cauchy, feasible shifts in [0, .2] (p = .5):", paste(format(fs[feas]), collapse = " "), "\n")

cat("(7) Lipschitz neighbourhood\n")
L <- max(abs(diff(lg_c(ys)) / diff(ys)))
eps_ok <- log(2) / L
ok(all(sapply(c(-1, -0.5, 0.5, 1) * eps_ok, function(e) log(0.5) + lsup(lg_c(ys + e) - lg_c(ys)) <= 1e-9)),
   sprintf("modulated Cauchy: log f is %.2f-Lipschitz; every |eps| <= log(2)/L = %.4f is feasible at p = .5", L, eps_ok))

cat("(8) discrete maps (Corollary 4)\n")
feasible <- function(cmap, Fv, Qv, p) all(sapply(seq_along(Qv), function(a) p * Qv[a] <= sum(Fv[cmap == a]) + 1e-12))
ok(feasible(c(1, 2, 3), rep(1/3, 3), c(1, 0, 0), 1/3), "F2 uniform on {0,1,2}, p = 1/3, Q = point mass at 0: the identity is feasible, its image {0,1,2} strictly contains the reported {0}")
ok(feasible(c(1, 2), c(.5, .5), c(1, 0), .5), "F2 = (.5,.5), Q = (1,0), p = .5: identity feasible with an unreported category")
maps3 <- as.matrix(expand.grid(1:3, 1:3, 1:3)); maps3 <- maps3[apply(maps3, 1, function(v) all(diff(v) >= 0)), ]
feas3 <- maps3[apply(maps3, 1, function(v) feasible(v, rep(1/3, 3), c(2/3, 0, 1/3), 1)), , drop = FALSE]
ok(nrow(feas3) == 1 && all(feas3[1, ] == c(1, 1, 3)), "F2 uniform, p = 1, Q = (2/3, 0, 1/3): the unique feasible map is (0,0,2), a non-identity map")

cat("(9) repeated tenures with a general map (Corollary 5)\n")
for (mu in c(0, 1)) {
  lr9 <- log(0.1) + lsup(log(2) + dnorm(2 * ys, mu, 1, log = TRUE) - dnorm(ys, mu, 1, log = TRUE))
  ok(abs(lr9 - (log(0.1) + log(2) + mu^2 / 6)) < 1e-6 && lr9 < 0, sprintf("F = Q = N(%d,1), p = .1: c'(y) = 2y satisfies the constraint (p x max ratio = %.4f)", mu, exp(lr9)))
}
cat("   implied selection differentials of c'(y) = 2y: 0 and -1/2 (mean mu/2 against mu); the true identity has 0 and 0\n")

cat("(10) covariate predicting retention only (Proposition 1)\n")
hw_lap <- function(p) hw(dlap, p, seq(0, 3, by = 0.001))
ok(abs(hw_lap(0.5) - log(2)) < 0.0015 && abs(min(hw_lap(0.8), hw_lap(0.2)) - log(1.25)) < 0.0015,
   "pooled p = .5: half-width log 2; strata p = .8, .2: intersection half-width log 1.25")

cat("(11) bounded common loading (Theorem 4(a))\n")
Gs <- seq(0, 1, by = 0.001)
lap_ok <- function(t) abs(t) <= log(2) + 1e-12                      # episode 1: Laplace F = Q, p = .5, C1 = 0, N1 = 1
gau_ok <- function(t) abs(t - 10) <= 1e-12                          # episode 2: F = N(0,1), Q = N(10,1), p = .5 -> {10}
joint <- Gs[sapply(Gs, function(G) lap_ok(0 - G * 1) && gau_ok(10 - G * 1))]
ok(length(joint) == 1 && joint == 0, "common-loading feasible set {0}: tau_1 = C1 - Gamma N1 = 0, not the per-pair outer bound [-log 2, 0]")
cat("(12) marginal versus joint target-battery information (Theorem 4(a))\n")
b12 <- 1 / sqrt(2)                                                  # Laplace scale: variance 2 b^2 = 1, like X = +-1
dlapb <- function(y) -abs(y) / b12 - log(2 * b12)
px <- c(0.8, 0.2); wx <- c(0.5, 0.5); xv <- c(1, -1)               # retention by battery stratum, independent of Y*
p12 <- sum(wx * px)                                                 # pooled retention
N12 <- sum(wx * px * xv) / p12 - sum(wx * xv)                       # battery contrast: survivors' mean of X minus the cohort's
C12 <- 0                                                            # Y* independent of X, identity map: target contrast zero
ok(abs(p12 - 0.5) < 1e-15 && abs(N12 - 0.6) < 1e-15, "p = .5, C = 0, N = .6: the true loading is C / N = 0")
feas12 <- function(p, hs) log(p) + lsup(dlapb(ys + hs) - dlapb(ys)) <= 1e-9
Gs12 <- seq(0, 1, by = 0.0005)
marg12 <- Gs12[sapply(Gs12, function(G) feas12(p12, C12 - G * N12))]
cond12 <- Gs12[sapply(Gs12, function(G) all(sapply(px, function(p) feas12(p, C12 - G * N12))))]
gm <- b12 * log(2) / 0.6; gc <- b12 * log(1.25) / 0.6
ok(min(marg12) == 0 && max(marg12) <= gm && gm - max(marg12) < 5e-4,
   sprintf("marginal information: Gamma in [0, %.4f] on the grid; closed form [0, b log 2 / .6] = [0, %.6f]", max(marg12), gm))
ok(min(cond12) == 0 && max(cond12) <= gc && gc - max(cond12) < 5e-4,
   sprintf("joint information (both strata): Gamma in [0, %.4f] on the grid; closed form [0, b log 1.25 / .6] = [0, %.6f]", max(cond12), gc))
r_pool <- exp(log(0.5) + lsup(dlapb(ys - 0.3) - dlapb(ys))); r_hi <- exp(log(0.8) + lsup(dlapb(ys - 0.3) - dlapb(ys)))
ok(abs(r_pool - 0.5 * exp(0.3 / b12)) < 1e-9 && round(r_pool, 6) == 0.764233 && r_pool <= 1,
   sprintf("Gamma = .5 (shift -.3): pooled constraint p exp(.3/b) = %.6f <= 1", r_pool))
ok(round(r_hi, 6) == 1.222772 && r_hi > 1,
   sprintf("Gamma = .5: stratum X = 1 requires .8 exp(.3/b) = %.6f <= 1, which fails", r_hi))
att_hi <- min(exp(dlapb(ys)) - 0.8 * exp(dlapb(ys - 0.3)))           # implied attriter mass in stratum X = 1, times .2
ok(att_hi < 0, sprintf("so no completion keeps F(. | X = 1): the implied attriter density there has a negative part (min %.4f)", att_hi / 0.2))
cat("\nall checks passed\n")

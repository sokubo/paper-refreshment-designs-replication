# ============================================================
# Exact (enumerated, sampling-free) checks of Lemma 2, Theorem 2, Corollary 4 and Corollary 6 of the refreshment-designs paper
#   (A) Lemma 2 (pattern completion) and Theorem 2 (the longitudinal structure adds nothing):
#       a 4-wave cohort on a finite latent lattice, a binary entry-only negative control W,
#       interrupted response patterns allowed, selection depending on the WHOLE latent path
#       (including the never-reported coordinates).  For every candidate map c'_4 = c_4 + eps
#       that satisfies the wave-4 domination inequality, the completion of the proof is built
#       explicitly and shown to reproduce EVERY observed pattern sub-measure, the entry marginal
#       (including W), and the refreshment marginal F_4, with a valid kernel pi'(r | y*).
#       Candidates violating the inequality are shown to have a negative leftover L_4.
#       The perturbed structure moves tau(4) by eps and the survivors' latent mean by -eps while
#       leaving every observable (in particular the negative-control contrast) unchanged.
#   (A') the monotone-attrition special case, with the sequential kernel pi'_j = dN_{j+1}/dN_j.
#   (B) Corollary 4 (discrete outcomes): for m = 4 categories, enumerate all 35 weakly increasing
#       maps and compare the mass-domination criterion with an LP feasibility check of
#       "exists nu <= F_2 with c'_# nu = p Q".
#   (C) Corollary 6 (mean bounds for an unrestricted map): Gaussian trimmed-mean widths and an LP
#       check of the bathtub argument on a grid.
# ============================================================
options(width = 120)
tol <- 1e-10
s <- 4
pats <- cbind(r1 = 1L, as.matrix(expand.grid(r2 = 0:1, r3 = 0:1, r4 = 0:1)))   # response patterns, r_1 = 1
P <- nrow(pats); pat_lab <- apply(pats, 1, paste, collapse = "")

## ---- generic representation: a structure is (state table with columns w,y1..y4; mass vector F; kernel matrix pi [n x P]; shifts h)
## observables: for each pattern, the sub-measure of REPORTED values (w, y1, c_j(y_j) for reported j), keyed by strings
observables <- function(state, Fv, pi, h) {
  out <- vector("list", P)
  for (i in seq_len(P)) {
    r <- pats[i, ]; mass <- Fv * pi[, i]
    key_cols <- list(state$w, state$y1)
    for (j in 2:s) if (r[j] == 1) key_cols[[length(key_cols) + 1]] <- state[[paste0("y", j)]] + h[j]
    out[[i]] <- tapply(mass, do.call(paste, c(key_cols, sep = ",")), sum)
  }
  out
}
maxdiff_obs <- function(o1, o2) max(sapply(seq_len(P), function(i) {
  a <- o1[[i]]; b <- o2[[i]]; keys <- union(names(a), names(b)); x <- y <- setNames(numeric(length(keys)), keys)
  x[names(a)] <- a; y[names(b)] <- b; max(abs(x - y)) }))
marg <- function(state, Fv, cols) tapply(Fv, do.call(paste, c(state[cols], sep = ",")), sum)
maxdiff_named <- function(a, b) { keys <- union(names(a), names(b)); x <- y <- setNames(numeric(length(keys)), keys); x[names(a)] <- a; y[names(b)] <- b; max(abs(x - y)) }

## ---- the true structure ------------------------------------------------------------------------------
K <- 5; lat <- (-K):K; m <- length(lat)
st <- expand.grid(w = 0:1, y1 = lat, y2 = lat, y3 = lat, y4 = lat)
n_st <- nrow(st)
disc <- function(mu, sd = 1.2) { d <- dnorm(lat, mu, sd); d / sum(d) }
u_vals <- c(-1, 0, 1); pu <- c(.3, .4, .3)
Fv <- numeric(n_st)
for (k in seq_along(u_vals)) {
  py <- disc(0.8 * u_vals[k]); pw <- c(1 - plogis(0.9 * u_vals[k]), plogis(0.9 * u_vals[k]))
  Fv <- Fv + pu[k] * pw[st$w + 1] * py[match(st$y1, lat)] * py[match(st$y2, lat)] * py[match(st$y3, lat)] * py[match(st$y4, lat)]
}
stopifnot(abs(sum(Fv) - 1) < 1e-12)
h_true <- c(0, 1, 0, 1)                                   # true conditioning: integer shifts by tenure, c_1 = id
wgt <- sapply(seq_len(P), function(i) {                   # true kernel pi(r | y*), depends on the whole latent path
  r <- pats[i, ]
  a <- 0.7 * st$y2 * r[2] + 0.5 * st$y3 * r[3] + 0.9 * st$y4 * r[4] - 0.4 * st$w * (sum(r) - 1) + 0.3 * st$y1 * (sum(r) - 1) - 0.6 * (sum(r) - 1)
  e <- exp(a); e[r[4] == 1 & st$y4 == -K] <- 0; e         # nobody with the lowest latent y4 responds at wave 4
})
pi_true <- wgt / rowSums(wgt)
stopifnot(all(abs(rowSums(pi_true) - 1) < 1e-12))

lam <- observables(st, Fv, pi_true, h_true)
p_r <- setNames(sapply(lam, sum), pat_lab)
F1 <- marg(st, Fv, c("w", "y1")); F4 <- marg(st, Fv, "y4")
resp4 <- pats[, 4] == 1; p4 <- sum(p_r[resp4])
cat("(A) pattern probabilities:\n"); print(round(p_r, 4)); cat(sprintf("    response rate at wave 4: p_4 = %.4f\n", p4))

## ---- the completion of Lemma 2 for candidate shifts hp (hp[4] constrained by F_4; hp[2], hp[3] free) ----------
complete <- function(lam, hp, F1, F4, p4) {
  ## (1) pull every pattern measure back to the latent scale through the candidate maps
  lt <- vector("list", P); Q4 <- numeric(0)
  for (i in seq_len(P)) {
    r <- pats[i, ]; keys <- strsplit(names(lam[[i]]), ",")
    mat <- do.call(rbind, lapply(keys, as.numeric)); rep_j <- (2:s)[r[2:s] == 1]
    colnames(mat) <- c("w", "y1", if (length(rep_j)) paste0("y", rep_j) else character(0))
    for (j in rep_j) mat[, paste0("y", j)] <- mat[, paste0("y", j)] - hp[j]          # latent = report - h'_j
    lt[[i]] <- list(mat = mat, mass = as.numeric(lam[[i]]))
    if (r[4] == 1) { tb <- tapply(lt[[i]]$mass, mat[, "y4"], sum); Q4 <- c(Q4, tb) }
  }
  Q4 <- tapply(Q4, names(Q4), sum)                                                  # wave-4 respondents' pulled-back sub-measure
  ## (2) leftover at wave 4: L_4 = F_4 - Q4;  feasible iff L_4 >= 0 (support of Q4 inside that of F_4)
  keys4 <- union(names(F4), names(Q4)); L4 <- setNames(numeric(length(keys4)), keys4)
  L4[names(F4)] <- F4; L4[names(Q4)] <- L4[names(Q4)] - Q4
  if (min(L4) < -tol) return(list(feasible = FALSE, minL4 = min(L4)))
  L4 <- pmax(L4, 0); kappa4 <- L4 / sum(L4)                                          # = L_4 / (1 - p_4)
  stopifnot(abs(sum(L4) - (1 - p4)) < 1e-9)
  ## (3) unrestricted tenures: arbitrary kappa_j (here uniform on the values that occur)
  vals <- lapply(2:4, function(j) sort(unique(unlist(lapply(lt, function(z) if (paste0("y", j) %in% colnames(z$mat)) z$mat[, paste0("y", j)] else NULL)))))
  kappa <- list(NULL, setNames(rep(1 / length(vals[[1]]), length(vals[[1]])), vals[[1]]),
                setNames(rep(1 / length(vals[[2]]), length(vals[[2]])), vals[[2]]), kappa4)
  ## (4) nu_r = lambda~_r (x) kappa_j on unreported j;  F' = sum_r nu_r;  pi'(r | y*) = nu_r / F'
  rows <- list()
  for (i in seq_len(P)) {
    r <- pats[i, ]; mat <- lt[[i]]$mat; mass <- lt[[i]]$mass
    ex <- data.frame(idx = seq_len(nrow(mat)), wgt = mass)
    for (j in 2:s) {
      if (r[j] == 1) ex[[paste0("y", j)]] <- mat[ex$idx, paste0("y", j)]
      else {
        kj <- kappa[[j]]; ex <- ex[rep(seq_len(nrow(ex)), each = length(kj)), , drop = FALSE]
        ex[[paste0("y", j)]] <- rep(as.numeric(names(kj)), times = nrow(ex) / length(kj)); ex$wgt <- ex$wgt * rep(kj, times = nrow(ex) / length(kj))
      }
    }
    rows[[i]] <- data.frame(pat = i, w = mat[ex$idx, "w"], y1 = mat[ex$idx, "y1"], y2 = ex$y2, y3 = ex$y3, y4 = ex$y4, wgt = ex$wgt)
  }
  rows <- do.call(rbind, rows)
  key <- paste(rows$w, rows$y1, rows$y2, rows$y3, rows$y4, sep = ",")
  ukey <- unique(key); idx <- match(key, ukey)
  nu <- matrix(0, length(ukey), P)
  for (i in seq_len(P)) { sel <- rows$pat == i; nu[, i] <- tapply(c(rows$wgt[sel], rep(0, length(ukey))), c(idx[sel], seq_along(ukey)), sum) }
  state <- as.data.frame(do.call(rbind, lapply(strsplit(ukey, ","), as.numeric))); names(state) <- c("w", "y1", "y2", "y3", "y4")
  Fp <- rowSums(nu); keep <- Fp > 0                                                  # states of zero mass (kappa_4 = 0 there) are dropped
  nu <- nu[keep, , drop = FALSE]; state <- state[keep, , drop = FALSE]; Fp <- Fp[keep]; pip <- nu / Fp
  list(feasible = TRUE, minL4 = min(L4), state = state, Fp = Fp, pip = pip)
}

cat("\n(A) candidates h'_4 = h_4 + eps (with h'_2 = h_2 + 1 and h'_3 = h_3 - 1 as arbitrary intermediate maps):\n")
res <- data.frame()
m_true <- sum((Fv * rowSums(pi_true[, resp4, drop = FALSE])) * st$y4) / p4                     # survivors' latent mean at wave 4
nc_true <- sum((Fv * rowSums(pi_true[, resp4, drop = FALSE])) * st$w) / p4 - sum(Fv * st$w)  # negative-control contrast
for (eps in -3:3) {
  hp <- c(0, h_true[2] + 1, h_true[3] - 1, h_true[4] + eps)
  cp <- complete(lam, hp, F1, F4, p4)
  if (!cp$feasible) { res <- rbind(res, data.frame(eps = eps, feasible = FALSE, min_L4 = cp$minL4, max_obs_diff = NA, F1_diff = NA, F4_diff = NA, kernel_ok = NA, mean_shift = NA, nc_shift = NA)); next }
  ob <- observables(cp$state, cp$Fp, cp$pip, hp)
  kernel_ok <- all(cp$pip >= -tol & cp$pip <= 1 + tol) && all(abs(rowSums(cp$pip) - 1) < 1e-9)
  m_new  <- sum((cp$Fp * rowSums(cp$pip[, resp4, drop = FALSE])) * cp$state$y4) / p4
  nc_new <- sum((cp$Fp * rowSums(cp$pip[, resp4, drop = FALSE])) * cp$state$w) / p4 - sum(cp$Fp * cp$state$w)
  res <- rbind(res, data.frame(eps = eps, feasible = TRUE, min_L4 = cp$minL4, max_obs_diff = maxdiff_obs(lam, ob),
                               F1_diff = maxdiff_named(marg(cp$state, cp$Fp, c("w", "y1")), F1), F4_diff = maxdiff_named(marg(cp$state, cp$Fp, "y4"), F4),
                               kernel_ok = kernel_ok, mean_shift = m_new - m_true, nc_shift = nc_new - nc_true))
}
print(res, digits = 4)
cat("    feasible eps: every observed pattern measure, F_1 (with W) and F_4 are reproduced to machine precision, the kernel is valid,\n")
cat("    the survivors' latent mean moves by exactly -eps and the negative-control contrast by 0;\n")
cat("    infeasible eps: the leftover L_4 has a negative entry, i.e. the wave-4 domination inequality fails.\n")

## ---- (A') monotone attrition: sequential kernel pi'_j = dN_{j+1}/dN_j ----------------------------------------
cat("\n(A') monotone-attrition special case:\n")
mono <- apply(pats, 1, function(r) all(diff(r) <= 0))
pi_mono <- pi_true; pi_mono[, !mono] <- 0; pi_mono <- pi_mono / rowSums(pi_mono)
lam_m <- observables(st, Fv, pi_mono, h_true); p_rm <- sapply(lam_m, sum); p4m <- sum(p_rm[resp4])
hp <- c(0, 1, 0, h_true[4] + 1)
cp <- complete(lam_m, hp, F1, F4, p4m)
if (cp$feasible) {
  N <- sapply(1:s, function(j) rowSums(cp$pip[, pats[, j] == 1, drop = FALSE]) * cp$Fp)      # N_j = sub-measure still responding at j
  ok_nest <- all(N[, 1] - N[, 2] >= -tol & N[, 2] - N[, 3] >= -tol & N[, 3] - N[, 4] >= -tol)
  pij <- sapply(1:(s - 1), function(j) ifelse(N[, j] > 0, N[, j + 1] / N[, j], 0))
  ok_range <- all(pij >= -tol & pij <= 1 + tol)
  cum <- cp$Fp; ok_chain <- TRUE
  for (j in 1:(s - 1)) { cum <- cum * pij[, j]; ok_chain <- ok_chain && max(abs(cum - N[, j + 1])) < 1e-12 }
  ob <- observables(cp$state, cp$Fp, cp$pip, hp)
  cat(sprintf("    eps = +1 feasible (p_4 = %.4f); N_1 >= ... >= N_4: %s; pi'_j in [0,1]: %s; prod_{u<j} pi'_u F' = N_j: %s; max obs diff = %.2e\n",
              p4m, ok_nest, ok_range, ok_chain, maxdiff_obs(lam_m, ob)))
} else cat(sprintf("    eps = +1 infeasible in the monotone example (min L_4 = %.4f)\n", cp$minL4))

## ---- (B) Corollary 4: discrete outcomes ------------------------------------------------------------------------
cat("\n(B) Corollary 4, m = 4 categories: mass-domination criterion vs LP feasibility for all weakly increasing maps\n")
suppressPackageStartupMessages(library(lpSolve))
m4 <- 4
maps <- as.matrix(expand.grid(rep(list(1:m4), m4))); maps <- maps[apply(maps, 1, function(v) all(diff(v) >= 0)), ]
cat(sprintf("    number of weakly increasing maps: %d (= choose(2m-1, m) = %d)\n", nrow(maps), choose(2 * m4 - 1, m4)))
id_row <- which(apply(maps, 1, function(v) all(v == 1:m4)))
check_one <- function(F2, Q, p) {
  crit <- lp_ok <- logical(nrow(maps))
  for (i in seq_len(nrow(maps))) {
    cm <- maps[i, ]
    crit[i] <- all(sapply(1:m4, function(a) p * Q[a] <= sum(F2[cm == a]) + 1e-12))
    A <- t(sapply(1:m4, function(a) as.numeric(cm == a)))                            # (c'_# nu)(a) = sum_{y: c'(y) = a} nu(y)
    sol <- lp("min", rep(0, m4), rbind(A, diag(m4)), c(rep("=", m4), rep("<=", m4)), c(p * Q, F2))
    lp_ok[i] <- sol$status == 0
  }
  c(agree = all(crit == lp_ok), n_feasible = sum(crit), identity_feasible = crit[id_row])
}
set.seed(1)
for (trial in 1:5) {
  F2 <- as.numeric(rmultinom(1, 400, runif(m4))) / 400; Q <- as.numeric(rmultinom(1, 300, runif(m4))) / 300; p <- runif(1, .2, .7)
  r <- check_one(F2, Q, p)
  cat(sprintf("    random trial %d: p = %.2f; criterion == LP for all 35 maps: %s; feasible maps: %d; identity feasible: %s\n", trial, p, r["agree"] == 1, r["n_feasible"], r["identity_feasible"] == 1))
}
F2 <- c(.25, .25, .25, .25); Q <- c(.3, .2, .2, .3)
r <- check_one(F2, Q, .8); cat(sprintf("    all categories reported, p = .8 (p*Q <= F2 everywhere): feasible maps = %d (identity only)\n", r["n_feasible"]))
r <- check_one(F2, Q, .9); cat(sprintf("    same with p = .9 (p*Q(1) = .27 > .25): feasible maps = %d (the deterministic model is rejected)\n", r["n_feasible"]))
Q <- c(.3, .3, .4, 0)
for (p in c(.7, .6)) {
  feas <- apply(maps, 1, function(cm) all(sapply(1:m4, function(a) p * Q[a] <= sum(F2[cm == a]) + 1e-12)))
  cat(sprintf("    stayers never report category 4 (Q = .3,.3,.4,0), p = %.1f: feasible maps: %s\n", p,
              paste(apply(maps[feas, , drop = FALSE], 1, paste, collapse = ""), collapse = " ")))
}
cat("    (a map is written as the images of categories 1..4; 1233 merges category 4 into 3; the identity is 1234)\n")

## ---- (C) Corollary 6: trimmed-mean bounds --------------------------------------------------------------------
cat("\n(C) Corollary 6, Gaussian F_2 (sd 1): width of the unrestricted-map interval, 2*phi(Phi^{-1}(p))/p\n")
for (p in c(.9, .8, .6, .5)) {
  w_formula <- 2 * dnorm(qnorm(p)) / p
  mU <- integrate(function(y) y * dnorm(y), qnorm(1 - p), Inf)$value / p
  g <- seq(-6, 6, by = 0.01); f <- dnorm(g); f <- f / sum(f)
  sol <- lp("max", g, rbind(rep(1, length(g)), diag(length(g))), c("=", rep("<=", length(g))), c(1, f / p))   # max mean s.t. 0 <= nu <= F2/p, sum nu = 1
  cat(sprintf("    p = %.1f: formula width %.4f; numeric integration %.4f; LP maximum on the grid %.4f vs m_U = %.4f\n", p, w_formula, 2 * mU, sol$objval, mU))
}
cat("    -> at p = .8 the unrestricted-map interval has width .70 sd; in Simulation 1 the population location-shift identified sets at p = .8 have width 0 to .26, and the relaxed plug-in sets .25 to .32.\n")

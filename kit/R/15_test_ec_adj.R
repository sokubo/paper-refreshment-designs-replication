# ============================================================
# P1 JLPS kit — unit test of the entry-wave corrections in the bootstrap engine of 15_panelcond_designs.R
# (kit v1.2, 2026-09-30). Three exact population fixtures, each a set of equally likely types, are pushed through
# the engine body itself (extracted from the script by parsing, so that the tested code is the code that runs);
# no data file is read. Expected values are exact fractions.
#   1. changing-subgroup standardisation: no conditioning, equal selection terms within every covariate value,
#      raw correction exact; the covariate distributions of the survivors T and of the survivors with an entry
#      answer P differ. EC = 0 and EC-adj = 0 (kit v1.2); the rule of kit v1.1, which averaged the entry residual
#      over P and the fresh prediction over T, gives 1/12.
#   2. entry item completion depending on the answer: EC = EC-adj = -1/2 with no conditioning.
#   3. fresh item completion depending on the answer: EC = EC-adj = +1/6 with no conditioning.
# Usage: Rscript R/15_test_ec_adj.R [path/to/15_panelcond_designs.R]
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
script <- if (length(args)) args[1] else file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "15_panelcond_designs.R")
ast <- parse(script, keep.source = FALSE)
find_engine <- function(z) {
  if (missing(z)) return(NULL)                      # empty argument slots (e.g. x[i, ]) parse as missing symbols
  if (is.call(z) && identical(z[[1]], as.name("<-")) && identical(z[[2]], as.name("engine"))) return(z[[3]])
  if (is.call(z) || is.expression(z) || is.pairlist(z)) for (el in as.list(z)) { r <- find_engine(el); if (!is.null(r)) return(r) }
  NULL
}
engine_def <- find_engine(ast)
stopifnot(!is.null(engine_def))

run_engine <- function(y_old, y_new, entry, survival, x_old, x_new, adj = TRUE) {
  env <- new.env(parent = baseenv())
  env$Yo <- matrix(y_old, ncol = 1, dimnames = list(NULL, "item")); env$Yn <- matrix(y_new, ncol = 1, dimnames = list(NULL, "item"))
  env$Eo <- matrix(entry, ncol = 1, dimnames = list(NULL, "item"))
  env$Xo <- data.frame(x = x_old); env$Xn <- data.frame(x = x_new)
  env$S_next <- rep(1L, length(y_old)); env$Sm <- env$Sm1 <- rep(1L, length(y_new))
  env$JJ <- 1L; env$cols <- "item"; env$ec_ok_col <- TRUE; env$setNames <- stats::setNames
  env$cm <- function(Mat, mask) { Mm <- Mat; Mm[!mask] <- NA; s <- colSums(Mm, na.rm = TRUE); n <- colSums(!is.na(Mm)); r <- s / n; r[n == 0] <- NA; r }
  eng <- eval(engine_def, envir = env)
  eng(seq_along(y_old), seq_along(y_new), as.integer(survival), adj = adj)
}
ols_pred <- function(y, Xfit, Xpred) { ok <- !is.na(y); b <- stats::lm.fit(Xfit[ok, , drop = FALSE], y[ok])$coefficients; b[is.na(b)] <- 0; mean(Xpred %*% b) }

n_ok <- 0L; n_fail <- 0L
expect <- function(label, got, want, tol = 1e-12) {
  ok <- all(abs(got - want) < tol)
  cat(sprintf("%-72s %-22s %-14s %s\n", label, paste(format(round(got, 10)), collapse = " "), paste(format(want), collapse = " "), if (ok) "ok" else "FAIL"))
  if (ok) n_ok <<- n_ok + 1L else n_fail <<- n_fail + 1L
}
cat(sprintf("%-72s %-22s %-14s\n", "check", "got", "expected")); cat(strrep("-", 115), "\n")

## 1. sixteen types (X, U, V, W) Bernoulli(1/2): S = U 1{X = 1 or V = 1}; eligible at entry G_e = 1{X = 0 or W = 1};
##    Y_t = Y_t* = X U; Y_e = X U + X/2 where eligible; fresh cohort with the same law, everyone answers.
z <- expand.grid(x = 0:1, u = 0:1, v = 0:1, w = 0:1)
S <- z$u == 1 & (z$x == 1 | z$v == 1); Ge <- z$x == 0 | z$w == 1
yt <- z$x * z$u; ye <- z$x * z$u + z$x / 2; ye[!Ge] <- NA
a <- run_engine(yt, yt, ye, S, z$x, z$x)
expect("fixture 1: raw entry-wave correction (true effect 0)", a$ec, 0)
expect("fixture 1: EC-adj, one covariate distribution (true effect 0)", a$ec_adj, 0)
io <- S; ip <- io & !is.na(ye); Xo <- cbind(1, z$x)
old_rule <- mean(yt[io]) - (mean(ye[ip]) - ols_pred(ye, Xo, Xo[ip, , drop = FALSE])) - ols_pred(yt, Xo, Xo[io, , drop = FALSE])
expect("fixture 1: rule of kit v1.1 (entry residual over P, fresh prediction over T)", old_rule, 1 / 12)
expect("fixture 1: covariate shares P(X = 1 | T), P(X = 1 | P)", c(mean(z$x[io]), mean(z$x[ip])), c(2 / 3, 1 / 2))
expect("fixture 1: engine survival share equals 3/8", a$p_surv, 3 / 8)

## 2. entry item completion depends on the answer: Y_e = Y_t = Y, S = U independent of Y, entry answer observed only
##    when Y = U; fresh cohort answers completely.
z2 <- expand.grid(y = 0:1, u = 0:1, v = 0:1, w = 0:1); S2 <- z2$u == 1
ye2 <- z2$y; ye2[z2$y != z2$u] <- NA
b <- run_engine(z2$y, z2$y, ye2, S2, z2$w, z2$w)
expect("fixture 2: raw correction with answer-dependent entry completion", b$ec, -1 / 2)
expect("fixture 2: EC-adj with answer-dependent entry completion", b$ec_adj, -1 / 2)

## 3. fresh item completion depends on the answer: fresh entrants with Y = 1 answer with probability 1/2.
yn3 <- z2$y; yn3[z2$y == 1 & z2$v == 0] <- NA
c3 <- run_engine(z2$y, yn3, z2$y, S2, z2$w, z2$w)
expect("fixture 3: raw correction with answer-dependent fresh completion", c3$ec, 1 / 6)
expect("fixture 3: EC-adj with answer-dependent fresh completion", c3$ec_adj, 1 / 6)

## 4. with every survivor answering at entry (P = T) the two rules coincide: fixture 1 with G_e = 1 for all.
ye4 <- z$x * z$u + z$x / 2
d <- run_engine(yt, yt, ye4, S, z$x, z$x)
old4 <- mean(yt[io]) - (mean(ye4[io]) - ols_pred(ye4, Xo, Xo[io, , drop = FALSE])) - ols_pred(yt, Xo, Xo[io, , drop = FALSE])
expect("fixture 4: P = T, EC-adj equals the v1.1 rule", c(d$ec_adj, old4), c(0, 0))

cat("\n", n_ok + n_fail, " checks: ", n_ok, " ok, ", n_fail, " failed\n", sep = "")
if (n_fail > 0) quit(status = 1)

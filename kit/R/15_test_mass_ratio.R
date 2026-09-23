# ============================================================
# P1 JLPS kit — unit test of the mass-domination helper (15_mass_ratio.R). No data are read.
#   Rscript R/15_test_mass_ratio.R        (exits with status 1 on any failure)
# Cases: positive/zero, zero/zero (codebook category observed in neither group), a rare (sparse) category,
# an ordinary finite ratio, the population case Q = F2, and ineligible items.
# ============================================================
.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })
source(file.path(.here, "15_mass_ratio.R"))

ok <- 0L; bad <- 0L
expect <- function(what, got, want, tol = 1e-12) {
  hit <- length(got) == length(want) && all(is.na(got) == is.na(want)) &&
         all(abs(got[!is.na(got)] - want[!is.na(want)]) <= tol | got[!is.na(got)] == want[!is.na(want)])
  if (isTRUE(hit)) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("%-66s %-12s %-12s %s\n", what, paste(format(got), collapse = " "),
              paste(format(want), collapse = " "), if (isTRUE(hit)) "ok" else "FAIL"))
}
rp <- function(v, n) rep(v, n)

## rule shipped with v0.5 (for comparison only): positive/zero categories were set to NA and dropped
old_rule <- function(yo, yn, p) {
  vals <- sort(unique(c(yo, yn)))
  q <- table(factor(yo, levels = vals)) / length(yo); f2 <- table(factor(yn, levels = vals)) / length(yn)
  r <- as.numeric(p * q / f2); r[f2 == 0] <- NA; max(r, na.rm = TRUE)
}

## 1. positive/zero: 100 stayers (99 zeros, 1 one), 100 fresh respondents all zero, p = .5
yo <- c(rp(0, 99), 1); yn <- rp(0, 100); m <- mass_ratio_item(yo, yn, .5)
expect("positive/zero: v0.5 rule returned .495 (category 1 dropped)", old_rule(yo, yn, .5), .495)
expect("positive/zero: finite maximum over positive denominators", m$ratio_finite, .495)
expect("positive/zero: categories with zero fresh-cohort mass flagged", m$n_zero_denom, 1L)
expect("positive/zero: item is eligible (two categories)", c(m$eligible, m$ncat), c(TRUE, 2L))

## 2. zero/zero: five-point codebook, category 3 used by nobody in either group
yo <- c(rp(1, 20), rp(2, 30), rp(4, 30), rp(5, 20)); yn <- c(rp(1, 25), rp(2, 25), rp(4, 25), rp(5, 25))
m <- mass_ratio_item(yo, yn, .8, codebook = 1:5); m0 <- mass_ratio_item(yo, yn, .8)
expect("zero/zero: not a category (four observed categories)", m$ncat, 4L)
expect("zero/zero: recorded as an unobserved codebook category", m$n_zero_zero, 1L)
expect("zero/zero: no zero-denominator flag", m$n_zero_denom, 0L)
expect("zero/zero: ratio unchanged by the codebook (.8 * .30 / .25)", c(m$ratio_finite, m0$ratio_finite), c(.96, .96))

## 3. rare category: 2 of 100 fresh respondents in category 2 against 10 of 100 stayers, p = .8
yo <- c(rp(1, 90), rp(2, 10)); yn <- c(rp(1, 98), rp(2, 2)); m <- mass_ratio_item(yo, yn, .8)
expect("rare: finite maximum attained in the sparse category (.8 * .10 / .02)", m$ratio_finite, 4)
expect("rare: sparse categories with ratio above one", m$n_sparse_gt1, 1L)
expect("rare: maximum over supported categories (.8 * .90 / .98)", m$ratio_supported, .8 * .9 / .98)
expect("rare: no zero-denominator flag", m$n_zero_denom, 0L)

## 4. ordinary finite ratio: both categories supported, p = .9
yo <- c(rp(1, 40), rp(2, 60)); yn <- c(rp(1, 50), rp(2, 50)); m <- mass_ratio_item(yo, yn, .9)
expect("ordinary: finite = supported maximum (.9 * .60 / .50)", c(m$ratio_finite, m$ratio_supported), c(1.08, 1.08))
expect("ordinary: no sparse or zero-denominator flag", c(m$n_sparse_gt1, m$n_zero_denom), c(0L, 0L))

## 5. population case Q = F2 (identity map, attrition independent of the outcome): every ratio equals p
yo <- c(rp(1, 30), rp(2, 50), rp(3, 20)); yn <- yo; m <- mass_ratio_item(yo, yn, .7)
expect("Q = F2: maximum equals retention", c(m$ratio_finite, m$ratio_supported), c(.7, .7))

## 6. retention above one (an estimated subgroup share can undershoot): the ratio is still computed
m <- mass_ratio_item(c(rp(1, 50), rp(2, 50)), c(rp(1, 50), rp(2, 50)), 1.2)
expect("p = 1.2: ratio equals p (not rejected as ineligible)", c(m$eligible, m$ratio_finite), c(TRUE, 1.2))

## 7. ineligible items: one category, ten categories, an empty group
expect("one observed category: not eligible", mass_ratio_item(rp(1, 5), rp(1, 5), .5)$eligible, FALSE)
expect("ten observed categories: not eligible", mass_ratio_item(1:10, 1:10, .5)$eligible, FALSE)
expect("no fresh-cohort answer: not eligible", mass_ratio_item(1:3, c(NA, NA), .5)$eligible, FALSE)

cat(sprintf("\n%d checks: %d ok, %d failed\n", ok + bad, ok, bad))
if (bad > 0L) quit(status = 1L)

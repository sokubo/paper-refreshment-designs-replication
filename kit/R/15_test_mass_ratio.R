# ============================================================
# P1 JLPS kit — unit test of the mass-domination helper (15_mass_ratio.R). No data are read.
#   Rscript R/15_test_mass_ratio.R        (exits with status 1 on any failure)
# Cases: positive/zero, zero/zero (codebook category observed in neither group), a rare (sparse) category,
# an ordinary finite ratio, the population case Q = F2, ineligible items, and the missing-mass allocation
# (fresh item nonresponse) including a constructed population example.
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

## 8. fresh item nonresponse: a constructed population example. Both cohorts Y ~ Bernoulli(.5), no conditioning,
##    everyone eligible; panel survival .8 independent of Y and every survivor answers; among fresh entrants all
##    with Y = 0 answer and .2 of those with Y = 1 (fresh answerers 5/6 : 1/6). Complete-case ratio 2.4 at Y = 1,
##    population ratio .8 in both categories; missing mass .4 and needed mass .3: the identity map is compatible.
yo <- c(rp(0, 400), rp(1, 400)); yn <- c(rp(0, 500), rp(1, 100)); m <- mass_ratio_item(yo, yn, .8, n_new_reached = 1000)
expect("counterexample: complete-case ratio 2.4 (.8 * .5 / (1/6))", m$ratio_finite, 2.4, 1e-9)
expect("counterexample: missing mass .4 among fresh entrants who reached the item", m$missing_mass, .4)
expect("counterexample: mass needed to cover the stayers .3", m$needed_mass, .3, 1e-9)
expect("counterexample: identity map compatible once the missing mass is allocated", m$identity_feasible, TRUE)
## 9. no fresh nonresponse: an exceedance cannot be reconciled (needed > 0 = missing)
yo <- c(rp(1, 40), rp(2, 60)); yn <- c(rp(1, 50), rp(2, 50)); m <- mass_ratio_item(yo, yn, .9, n_new_reached = 100)
expect("no missing mass: needed mass .9 * .60 - .50 = .04, not reconcilable", c(m$missing_mass, m$needed_mass, m$identity_feasible), c(0, .04, FALSE), 1e-9)
## 10. positive/zero with enough missing mass: reconcilable; with too little: not
yo <- c(rp(0, 90), rp(1, 10)); yn <- rp(0, 80); m <- mass_ratio_item(yo, yn, .5, n_new_reached = 100)
expect("positive/zero, missing .2 >= needed .05: reconcilable", c(m$missing_mass, m$needed_mass, m$identity_feasible), c(.2, .05, TRUE), 1e-9)
m <- mass_ratio_item(yo, yn, .5, n_new_reached = 82)
expect("positive/zero, missing 2/82 < needed .05: not reconcilable", c(round(m$missing_mass, 4), m$needed_mass, m$identity_feasible), c(round(2/82, 4), .05, FALSE), 1e-9)
## 11. n_new_reached not supplied: the allocation fields are NA and everything else unchanged
m0 <- mass_ratio_item(yo, yn, .5); expect("no reach count: allocation fields NA", c(is.na(m0$missing_mass), is.na(m0$needed_mass), is.na(m0$identity_feasible)), c(TRUE, TRUE, TRUE))

## 12. shared reach of a nominal question (kit v1.1, 2026-09-29). A nominal question with categories A, B, C is
##     represented by three 0/1 indicators; 04 attaches its item-nonresponse (M) column to the FIRST indicator only.
##     The diagnostic must assess every indicator with the question's reach, and its scope and findings must not depend
##     on the order of the categories. Example: retention .8, stayer probabilities (.5, .4, .1), fresh
##     probabilities (.5, 0, .5), complete response: the B indicator has a positive/zero category (stayers report B,
##     no fresh entrant does) and the identity map is not compatible with it.
build_nominal <- function(order) {                     # order: permutation of the three categories in the columns
  n_old <- 1000L; n_new <- 200L
  S <- rep(c(1L, 0L), c(800L, 200L))                     # retention .8, independent of the answer
  ans_old <- c(rep(c(1L, 2L, 3L), c(400L, 320L, 80L)), rep(c(1L, 2L, 3L), c(100L, 80L, 20L)))   # survivors: (.5, .4, .1); attriters likewise
  ans_new <- rep(c(1L, 3L), c(100L, 100L))              # fresh: (.5, 0, .5)
  Yo <- sapply(order, function(k) as.numeric(ans_old == k)); Yn <- sapply(order, function(k) as.numeric(ans_new == k))
  Mo <- matrix(NA_real_, n_old, 3); Mn <- matrix(NA_real_, n_new, 3)
  Mo[, 1] <- 0; Mn[, 1] <- 0                            # everyone reached the question and answered: M = 0 on the first indicator only
  list(Yo = Yo, Yn = Yn, Mo = Mo, Mn = Mn, S = S, qgroup = rep("Q", 3), order = order)
}
d <- build_nominal(c(1L, 2L, 3L))
r <- mass_diagnostic_columns(d$Yo, d$Yn, d$Mo, d$Mn, d$S, d$qgroup)
expect("shared reach: every indicator of the question is assessed", r$assessed, c(TRUE, TRUE, TRUE))
expect("shared reach: the same reach and retention for every indicator", c(r$reach_new, r$funnel_p), c(1, 1, 1, .8, .8, .8), 1e-9)
expect("shared reach: category B (no fresh entrant) is flagged positive/zero and not reconcilable", c(r$funnel_zero_denom[2], r$funnel_identity_feasible[2]), c(1L, FALSE))
expect("shared reach: A has ratio .8 (no flag); C exceeds one in its complement category (.8 * .9 / .5 = 1.44)", c(r$funnel_ratio[1], r$funnel_ratio[3], r$funnel_zero_denom[1], r$funnel_zero_denom[3]), c(.8, 1.44, 0L, 0L), 1e-9)
## the rule before 2026-09-29 read each indicator's own M column: only the first indicator was assessed, so the
## finding depended on which category came first (A first: nothing flagged; B first: the flag appears)
own <- function(d) mass_diagnostic_columns(d$Yo, d$Yn, d$Mo, d$Mn, d$S, qgroup = rep("", 3))
expect("previous rule (own M column): only the first indicator assessed", own(d)$assessed, c(TRUE, FALSE, FALSE))
expect("previous rule: with A first no flag appears", own(d)$funnel_zero_denom[1], 0L)
expect("previous rule: with B first the flag appears (order-dependent)", own(build_nominal(c(2L, 1L, 3L)))$funnel_zero_denom[1], 1L)
## permutation invariance of the corrected rule: the results travel with the categories
r2 <- mass_diagnostic_columns(d$Yo, d$Yn, d$Mo, d$Mn, d$S, d$qgroup)   # baseline again
for (perm in list(c(2L, 1L, 3L), c(3L, 2L, 1L), c(3L, 1L, 2L))) {
  dp <- build_nominal(perm); rp_ <- mass_diagnostic_columns(dp$Yo, dp$Yn, dp$Mo, dp$Mn, dp$S, dp$qgroup)
  expect(sprintf("permutation %s: assessed set, ratios and flags identical up to the permutation", paste(perm, collapse = "")),
         c(rp_$assessed, rp_$funnel_ratio, rp_$funnel_zero_denom, as.numeric(rp_$funnel_identity_feasible)),
         c(r2$assessed[perm], r2$funnel_ratio[perm], r2$funnel_zero_denom[perm], as.numeric(r2$funnel_identity_feasible[perm])), 1e-9)
}
## a base item (no group) keeps its own reach; a question nobody reached is skipped
expect("reach source: base items map to themselves, indicators to the first of their group", mass_reach_source(c("", "Q", "Q", "", "R", "R", "R")), c(1L, 2L, 2L, 4L, 5L, 5L, 5L))
d0 <- d; d0$Mn[] <- NA_real_
expect("no fresh entrant reached the question: nothing assessed", mass_diagnostic_columns(d0$Yo, d0$Yn, d0$Mo, d0$Mn, d0$S, d0$qgroup)$assessed, c(FALSE, FALSE, FALSE))

cat(sprintf("\n%d checks: %d ok, %d failed\n", ok + bad, ok, bad))
if (bad > 0L) quit(status = 1L)

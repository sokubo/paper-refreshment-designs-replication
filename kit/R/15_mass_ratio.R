# ============================================================
# P1 JLPS kit — 15 helper: sample analogue of the mass-domination criterion (Corollary 4 of the
# refreshment-designs paper, eq. (mass)): for a discrete item, the population inequality is
#     p * Q({a}) <= F2({a})   for every category a,
# with Q the stayers' reported distribution, F2 the refreshment (fresh-cohort) distribution and p retention.
#
# Category set. The categories are the union of the categories observed among the stayers and among the
# fresh cohort. A codebook category observed in neither group (zero/zero) contributes 0 <= 0, which is no
# restriction; it is not a category here and does not count towards the two-to-nine-category eligibility.
# If `codebook` is supplied, such categories are counted in `n_zero_zero` for the record and otherwise ignored.
#
# Three kinds of category are kept apart (never silently discarded):
#   positive/zero  stayer mass > 0 and no fresh-cohort respondent: the sample ratio is infinite. Counted in
#                  `n_zero_denom`; excluded from both finite maxima below.
#   sparse         1 .. (min_fresh - 1) fresh-cohort respondents: the ratio is finite but its denominator rests
#                  on fewer than min_fresh persons. Counted in `n_sparse_gt1` when its ratio exceeds one.
#   supported      at least min_fresh fresh-cohort respondents.
# `ratio_finite` is the maximum over all categories with a positive denominator (sparse and supported);
# `ratio_supported` is the maximum over supported categories only.
# An infinite or large sample ratio is not evidence that a population mass is zero; the diagnostic is exploratory.
# p is the retention of the item's own population: for an item asked only of a subgroup G it is P(stayer and
# answer | G); the caller supplies it (15_panelcond_designs.R divides the cohort rate by the fresh share of G).
# ============================================================
mass_ratio_item <- function(yo, yn, p, max_cat = 9L, min_fresh = 10L, codebook = NULL) {
  yo <- as.numeric(yo[!is.na(yo)]); yn <- as.numeric(yn[!is.na(yn)])
  out <- list(eligible = FALSE, ncat = NA_integer_, ratio_finite = NA_real_, ratio_supported = NA_real_,
              n_zero_denom = NA_integer_, n_sparse_gt1 = NA_integer_, n_zero_zero = NA_integer_)
  ## p may exceed one in a sample (an estimated subgroup share can undershoot); the inequality is still evaluated
  if (!length(yo) || !length(yn) || !is.finite(p) || p <= 0) return(out)
  vals <- sort(unique(c(yo, yn)))                       # union of the observed supports
  if (length(vals) < 2L || length(vals) > max_cat) return(out)
  no <- tabulate(match(yo, vals), length(vals))         # stayer counts
  nn <- tabulate(match(yn, vals), length(vals))         # fresh-cohort counts
  q <- no / length(yo); f2 <- nn / length(yn)
  pos <- nn > 0L; sup <- nn >= min_fresh
  r <- rep(Inf, length(vals)); r[pos] <- p * q[pos] / f2[pos]
  r[!pos & no == 0L] <- NA_real_                        # cannot occur (union of supports); kept for safety
  list(eligible        = TRUE,
       ncat            = length(vals),
       ratio_finite    = if (any(pos)) max(r[pos]) else NA_real_,
       ratio_supported = if (any(sup)) max(r[sup]) else NA_real_,
       n_zero_denom    = sum(!pos & no > 0L),
       n_sparse_gt1    = sum(pos & !sup & r > 1),
       n_zero_zero     = if (is.null(codebook)) 0L else sum(!(unique(as.numeric(codebook)) %in% vals)))
}

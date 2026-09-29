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
#
# Item nonresponse in the fresh cohort (since 2026-09-24). F2 above is the distribution among fresh entrants who
# gave a substantive answer; the population inequality concerns the latent distribution of everyone in G. If
# answering is selective on the outcome, the complete-case ratio can exceed one with no conditioning at all. When
# `n_new_reached` (fresh entrants who reached the item, answering or not) is supplied, the helper also evaluates the
# identity map against the fresh distribution WITH its missing mass allocated freely: with
#     l_a = p * Q({a})                       (stayers who answered a, per member of G)
#     r_a = (1 - r_M) * F2({a})              (fresh entrants who answered a, per member of G)
#     r_M = 1 - n_answered / n_reached       (fresh entrants who reached the item but gave no substantive answer)
# a common latent distribution compatible with the identity map exists iff  sum_a max{l_a - r_a, 0} <= r_M
# (`needed_mass` <= `missing_mass`, `identity_feasible`). This is a population compatibility calculation, not a
# calibrated test; an item whose exceedance disappears once the missing mass is allocated is a diagnostic of item
# response or coding, not evidence against the population implication.
# ============================================================
mass_ratio_item <- function(yo, yn, p, max_cat = 9L, min_fresh = 10L, codebook = NULL, n_new_reached = NULL) {
  yo <- as.numeric(yo[!is.na(yo)]); yn <- as.numeric(yn[!is.na(yn)])
  out <- list(eligible = FALSE, ncat = NA_integer_, ratio_finite = NA_real_, ratio_supported = NA_real_,
              n_zero_denom = NA_integer_, n_sparse_gt1 = NA_integer_, n_zero_zero = NA_integer_,
              missing_mass = NA_real_, needed_mass = NA_real_, identity_feasible = NA)
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
  ## missing-mass allocation (identity map): needs the number of fresh entrants who reached the item
  mm <- nd <- NA_real_; feas <- NA
  if (!is.null(n_new_reached) && is.finite(n_new_reached) && n_new_reached >= length(yn)) {
    mm <- 1 - length(yn) / n_new_reached
    nd <- sum(pmax(p * q - (1 - mm) * f2, 0))
    feas <- nd <= mm + 1e-12
  }
  list(eligible        = TRUE,
       ncat            = length(vals),
       ratio_finite    = if (any(pos)) max(r[pos]) else NA_real_,
       ratio_supported = if (any(sup)) max(r[sup]) else NA_real_,
       n_zero_denom    = sum(!pos & no > 0L),
       n_sparse_gt1    = sum(pos & !sup & r > 1),
       n_zero_zero     = if (is.null(codebook)) 0L else sum(!(unique(as.numeric(codebook)) %in% vals)),
       missing_mass    = mm,
       needed_mass     = nd,
       identity_feasible = feas)
}

# ------------------------------------------------------------
# Reach of a question shared by its category indicators (since 2026-09-29; kit v1.1).
# 04_build_analysis.R attaches the item-nonresponse (M) column of a nominal question to its FIRST indicator only,
# so that the question's nonresponse is counted once in the B family. The mass diagnostic, however, needs the
# question's reach (entrants routed to it) for EVERY indicator: before this helper it read the reach from each
# indicator's own M column, so every indicator but the first had zero reach and was skipped, and which category
# was examined depended on the order of the codes. mass_reach_source() maps every column to the column that
# carries its question's M column (the first indicator for the members of a nominal group `qgroup`, the column
# itself otherwise); mass_diagnostic_columns() computes the diagnostic for every column with that shared reach,
# so the result does not depend on the order of the categories (unit tests in 15_test_mass_ratio.R).
# ------------------------------------------------------------
mass_reach_source <- function(qgroup) {
  qgroup <- as.character(qgroup); qgroup[is.na(qgroup)] <- ""
  src <- seq_along(qgroup)
  for (q in unique(qgroup[nzchar(qgroup)])) { i <- which(qgroup == q); src[i] <- i[1] }
  src
}

## Yo, Yn: item matrices (old cohort x items, fresh cohort x items; NA = no substantive answer)
## Mo, Mn: the item-nonresponse indicator matrices of the same items (NA = did not reach the question)
## S: survivor indicator of the old cohort (0/1); qgroup: "" for a base item, the question name for its indicators
mass_diagnostic_columns <- function(Yo, Yn, Mo, Mn, S, qgroup = rep("", ncol(Yo)), min_fresh = 10L, max_cat = 9L) {
  J <- ncol(Yo); S <- as.integer(S); src <- mass_reach_source(qgroup)
  IOf <- (S == 1L) & !is.na(Yo)
  reached_new <- colSums(!is.na(Mn[, src, drop = FALSE]))                 # fresh entrants who reached the question
  elig_new <- reached_new / nrow(Yn)
  reach_old_S <- colSums(!is.na(Mo[S == 1L, src, drop = FALSE])) / sum(S)   # survivors who reached the question
  ans_old_S <- colSums(IOf); ans_new <- colSums(!is.na(Yn)); p_surv <- colMeans(IOf)
  out <- data.frame(reach_new = NA_real_, reach_old_S = NA_real_, answer_new_given_reach = NA_real_, answer_old_S_given_reach = NA_real_,
                    funnel_p = NA_real_, funnel_ratio = NA_real_, funnel_ratio_supported = NA_real_, funnel_zero_denom = NA_integer_,
                    funnel_sparse_gt1 = NA_integer_, funnel_ncat = NA_integer_, funnel_missing_mass = NA_real_, funnel_needed_mass = NA_real_,
                    funnel_identity_feasible = NA, assessed = FALSE)[rep(1L, J), ]
  rownames(out) <- NULL
  for (j in seq_len(J)) {
    if (!(elig_new[j] > 0)) next
    out$funnel_p[j] <- p_surv[j] / elig_new[j]
    out$reach_new[j] <- elig_new[j]; out$reach_old_S[j] <- reach_old_S[j]
    out$answer_new_given_reach[j] <- ans_new[j] / reached_new[j]
    out$answer_old_S_given_reach[j] <- if (reach_old_S[j] > 0) ans_old_S[j] / (reach_old_S[j] * sum(S)) else NA_real_
    mr <- mass_ratio_item(Yo[IOf[, j], j], Yn[!is.na(Yn[, j]), j], out$funnel_p[j], max_cat = max_cat, min_fresh = min_fresh,
                          n_new_reached = reached_new[j])
    if (!mr$eligible) next
    out$assessed[j] <- TRUE
    out$funnel_ratio[j] <- mr$ratio_finite; out$funnel_ratio_supported[j] <- mr$ratio_supported
    out$funnel_zero_denom[j] <- mr$n_zero_denom; out$funnel_sparse_gt1[j] <- mr$n_sparse_gt1; out$funnel_ncat[j] <- mr$ncat
    out$funnel_missing_mass[j] <- mr$missing_mass; out$funnel_needed_mass[j] <- mr$needed_mass; out$funnel_identity_feasible[j] <- mr$identity_feasible
  }
  out
}

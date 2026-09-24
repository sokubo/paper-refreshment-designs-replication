# ============================================================
# P1 JLPS kit — response-style composites over a fixed list of rating-scale items
# Sourced by 04_build_analysis.R and 15_panelcond_designs.R (needs data.table; `.here` = the R/ folder).
#
# The person-level midpoint and extreme-category shares are computed over the items listed in
# 15_style_items.csv (set = "rating": agreement, satisfaction and evaluation scales with verbal anchors;
# set = "frequency": behaviour-frequency scales, used only in the sensitivity file). Before 2026-09-24 every
# item with four to seven consecutively numbered substantive codes entered the composites, whatever its
# meaning (marital status, occupational rank, education, ...); that rule is kept only as the "v0.7 rule"
# rows of 15_style_sensitivity.csv.
#
# For every listed item the codes found in the file are verified against the list: the substantive codes
# (the value labels minus the special codes: don't know, refusal, no answer, not applicable, and the
# NAP_OVERRIDE codes of 00_config.R) must be exactly 1..k, and an item declared to have a midpoint must
# have an odd k whose middle code is labelled as a neutral category (NEUTRAL_PAT). An item that fails the
# code check enters no composite; an item that fails only the midpoint check enters the extreme-category
# composite but not the midpoint composite. 15 writes the outcome for every listed item, with the labels
# read from the file, to 15_style_items_audit.csv.
# P1_STYLE_ITEMS names another list (the synthetic example uses synthetic/style_items_synthetic.csv).
# ============================================================

STYLE_ITEMS_FILE <- path.expand(Sys.getenv("P1_STYLE_ITEMS", file.path(.here, "15_style_items.csv")))
NEUTRAL_PAT <- "どちらともいえない|どちらでもない|ふつう|普通|変わらない"

read_style_items <- function(file = STYLE_ITEMS_FILE) {
  if (!file.exists(file)) stop("style item list not found: ", file, call. = FALSE)
  ## strict: a malformed line (e.g. an unquoted comma in a note) must stop the run, not truncate the list
  st <- withCallingHandlers(fread(file, encoding = "UTF-8", colClasses = list(character = c("var", "set", "battery", "note"))),
                            warning = function(w) stop("reading ", file, ": ", conditionMessage(w), call. = FALSE))
  n_lines <- sum(nzchar(trimws(readLines(file, encoding = "UTF-8", warn = FALSE)))) - 1L
  if (nrow(st) != n_lines) stop(sprintf("reading %s: %d rows read but %d data lines in the file", file, nrow(st), n_lines), call. = FALSE)
  st[, var := toupper(trimws(var))]; st[, k := as.integer(k)]; st[, midpoint := as.logical(midpoint)]
  stopifnot(!anyDuplicated(st$var), all(st$set %in% c("rating", "frequency")), all(st$k >= 3), all(st$k <= 9),
            all(!is.na(st$midpoint)), all(!st$midpoint | st$k %% 2L == 1L))
  st
}

## identity of the list actually used (recorded by 04 in the derived file and by 15 in its run record)
style_items_id <- function(file = STYLE_ITEMS_FILE) list(path = normalizePath(file), md5 = unname(tools::md5sum(file)))

## the observed substantive values of an item must lie in 1..k as well (labels alone do not guarantee it)
style_values_ok <- function(x, k) { v <- x[!is.na(x)]; !length(v) || all(v %in% seq_len(k)) }

## verify one item's value labels against its row of the list
style_verify <- function(vl, spec, k, midpoint, neutral_pat = NEUTRAL_PAT) {
  if (is.null(vl) || !length(vl)) return(list(ok = FALSE, mid_ok = FALSE, codes = "", labels = "", reason = "no value labels"))
  keep <- !(unname(vl) %in% spec)
  codes <- unname(vl)[keep]; labs <- names(vl)[keep]
  o <- order(codes); codes <- codes[o]; labs <- labs[o]
  codes_s <- paste(codes, collapse = " "); labs_s <- paste(labs, collapse = " | ")
  if (length(codes) != k || !isTRUE(all(codes == seq_len(k))))
    return(list(ok = FALSE, mid_ok = FALSE, codes = codes_s, labels = labs_s,
                reason = sprintf("substantive codes are not 1..%d", k)))
  mid_ok <- FALSE; reason <- ""
  if (isTRUE(midpoint)) {
    mid_ok <- grepl(neutral_pat, labs[(k + 1L) / 2L])
    if (!mid_ok) reason <- "middle category not labelled as neutral: excluded from the midpoint composite"
  }
  list(ok = TRUE, mid_ok = mid_ok, codes = codes_s, labels = labs_s, reason = reason)
}

## person-level shares over a matrix of substantive answers (columns named by item): the share of answered
## items on which the answer is an extreme code (lo or lo + k - 1) and, over the midpoint items, the middle code
style_shares <- function(Mat, ext_items, mid_items, k_of, lo_of = NULL) {
  n <- nrow(Mat); if (is.null(lo_of)) lo_of <- setNames(rep(1, length(k_of)), names(k_of))
  ext_cnt <- ext_n <- mid_cnt <- mid_n <- rep(0, n)
  for (v in ext_items) { x <- Mat[, v]; ok <- !is.na(x); lo <- lo_of[[v]]; k <- k_of[[v]]
    ext_cnt[ok] <- ext_cnt[ok] + (x[ok] == lo | x[ok] == lo + k - 1); ext_n[ok] <- ext_n[ok] + 1 }
  for (v in mid_items) { x <- Mat[, v]; ok <- !is.na(x); lo <- lo_of[[v]]; k <- k_of[[v]]
    mid_cnt[ok] <- mid_cnt[ok] + (x[ok] == lo + (k - 1) / 2); mid_n[ok] <- mid_n[ok] + 1 }
  list(mid = ifelse(mid_n > 0, mid_cnt / mid_n, NA_real_), ext = ifelse(ext_n > 0, ext_cnt / ext_n, NA_real_),
       n_mid = length(mid_items), n_ext = length(ext_items))
}

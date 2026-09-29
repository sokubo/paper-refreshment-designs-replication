# ============================================================
# P1 JLPS kit — the item specification table (15_item_scale.csv) and the helpers that apply it
# Sourced by 04_build_analysis.R and 15_panelcond_designs.R (needs data.table; `.here` = the R/ folder).
#
# Since v1.0 (2026-09-28) the scale of every wave-5 item, its special codes beyond the value-label rule,
# its recodes, the indicator construction of nominal single-choice questions, the derived items built from
# components (clock times, durations), the filters, and the entry-wave comparability class are declared
# item by item in R/15_item_scale.csv, which was checked against the 2011 and 2007 questionnaires
# (review_round06/SCALE_AUDIT_2011_ja.md, EC_ITEMS_REVIEW_ja.md, DECISIONS_2026-09-27_ja.md). Before v1.0
# the scale was inferred from the number of value labels alone and the special codes from label strings
# only, which let nominal codes and scale-external codes (e.g. "public sector", "other", "no parent")
# enter the means.
#
# Columns of 15_item_scale.csv
#   var             wave-5 variable (case-insensitive) or the name of a derived item
#   scale           binary | ordinal | continuous | nominal | exclude
#   codes_expected  substantive codes the file must show after the recode ("1;2;3;4"); empty = not checked
#   nap_add, dk_add, nr_add, ref_add   codes treated as not applicable / don't know / no answer / refusal in
#                   addition to the label rule (DK_PAT, NR_PAT, ...); applied to the entry-wave variable too
#   recode          "from=to;from=to" applied to the substantive codes at both waves ("." = missing)
#   construct       "indicators" (one 0/1 item per substantive code of a nominal question),
#                   "clock:X,Y,Z,kind[,leave_item]" (AM/PM, hour, minute components -> decimal hours),
#                   "months:Y,M" (years and months -> months)
#   filter          "VAR=code" keeps only respondents with VAR == code; "not:VAR=code" drops them
#   in_counts       TRUE/FALSE: whether the item (or its indicators) enters the multiplicity families and counts
#   ec_class        A/B/C/D or empty: comparability of the 2007 entry wave (EC_ITEMS_REVIEW_ja.md)
#   ec_subgroup     for class C: the 2007 subgroup within which the entry-wave term is computed
#                   (defined in 15_entry_subgroups.csv)
#   ec_main         TRUE/FALSE: whether the item is in the main entry-wave analysis (classes A, B, C)
#   note, source    free text
# P1_ITEM_SCALE names another table (the synthetic example uses one written by 15_make_synthetic_varmap.R).
# ============================================================

ITEM_SCALE_FILE <- path.expand(Sys.getenv("P1_ITEM_SCALE", file.path(.here, "15_item_scale.csv")))
ENTRY_SUBGROUPS_FILE <- file.path(.here, "15_entry_subgroups.csv")
CLOCK_SPECIAL <- c(88, 99, 888, 999)     # not applicable / no answer codes of the clock-time components

read_item_scale <- function(file = ITEM_SCALE_FILE) {
  if (!file.exists(file)) stop("item specification table not found: ", file, call. = FALSE)
  warns <- character(0)
  tb <- withCallingHandlers(fread(file, encoding = "UTF-8", colClasses = "character", na.strings = NULL),
                            warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (length(warns)) stop("reading ", file, ": ", paste(warns, collapse = "; "), call. = FALSE)
  need <- c("var", "scale", "codes_expected", "nap_add", "dk_add", "nr_add", "ref_add", "recode", "construct", "filter",
            "in_counts", "ec_class", "ec_subgroup", "ec_main")
  miss <- setdiff(need, names(tb)); if (length(miss)) stop("15_item_scale.csv lacks the column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  for (cc in names(tb)) { x <- tb[[cc]]; x[is.na(x)] <- ""; set(tb, j = cc, value = trimws(x)) }
  set(tb, j = "key", value = toupper(tb$var))
  fail <- function(what, bad) if (any(bad)) stop(sprintf("15_item_scale.csv: %s for %s", what, paste(unique(tb$var[bad]), collapse = ", ")), call. = FALSE)
  fail("duplicated item", duplicated(tb$key) | duplicated(tb$key, fromLast = TRUE))
  fail("unknown scale", !(tb$scale %in% c("binary", "ordinal", "continuous", "nominal", "exclude")))
  fail("in_counts is not TRUE or FALSE", !(tb$in_counts %in% c("TRUE", "FALSE")))
  fail("ec_main is not TRUE or FALSE", !(tb$ec_main %in% c("TRUE", "FALSE")))
  fail("ec_class is not A, B, C, D or empty", !(tb$ec_class %in% c("", "A", "B", "C", "D")))
  fail("construct must be indicators, clock:... or months:...", nzchar(tb$construct) & !grepl("^(indicators|clock:|months:)", tb$construct))
  fail("nominal items must construct indicators", tb$scale == "nominal" & tb$construct != "indicators")
  fail("indicators are only for nominal items", tb$construct == "indicators" & tb$scale != "nominal")
  fail("a derived item must be continuous", grepl("^(clock|months):", tb$construct) & tb$scale != "continuous")
  fail("filter must be VAR=code or not:VAR=code", nzchar(tb$filter) & !grepl("^(not:)?[A-Za-z0-9_]+=[-0-9.]+$", tb$filter))
  fail("recode must be from=to;from=to", nzchar(tb$recode) & !grepl("^[-0-9.]+=([-0-9.]+|\\.)(;[-0-9.]+=([-0-9.]+|\\.))*$", tb$recode))
  set(tb, j = "in_counts", value = tb$in_counts == "TRUE"); set(tb, j = "ec_main", value = tb$ec_main == "TRUE")
  set(tb, i = which(tb$scale == "exclude"), j = "in_counts", value = FALSE)
  tb
}

read_entry_subgroups <- function(file = ENTRY_SUBGROUPS_FILE) {
  if (!file.exists(file)) return(data.table(subgroup = character(0), w1_var = character(0), codes = character(0)))
  sg <- fread(file, encoding = "UTF-8", colClasses = "character", na.strings = NULL)
  set(sg, j = "subgroup", value = trimws(sg$subgroup)); set(sg, j = "w1_var", value = toupper(trimws(sg$w1_var))); set(sg, j = "codes", value = trimws(sg$codes))
  if (any(duplicated(sg$subgroup))) stop("15_entry_subgroups.csv: duplicated subgroup", call. = FALSE)
  sg
}

## identity of the table actually used (recorded by 04 in the derived file and by 15 in its run record)
item_scale_id <- function(file = ITEM_SCALE_FILE) list(path = file, md5 = unname(tools::md5sum(file)))

parse_codes <- function(s) { s <- trimws(s); if (!nzchar(s)) return(numeric(0)); as.numeric(strsplit(s, ";")[[1]]) }

## recode "from=to;..." on a numeric vector; "." on the right-hand side sets the value to missing
apply_recode <- function(x, recode) {
  recode <- trimws(recode); if (!nzchar(recode)) return(x)
  pairs <- strsplit(recode, ";")[[1]]
  from <- as.numeric(sub("=.*$", "", pairs)); to_s <- sub("^.*=", "", pairs)
  to <- suppressWarnings(as.numeric(to_s)); to[to_s == "."] <- NA
  hit <- match(x, from); out <- x; out[!is.na(hit)] <- to[hit[!is.na(hit)]]; out
}
## the same recode applied to a set of codes (label sets are compared after recoding)
recode_codes <- function(codes, recode) sort(unique(na.omit(apply_recode(codes, recode))))

## special codes of an item: label rule (patterns) plus the table's additions; returns a list by kind
special_codes_of <- function(vl, spec, pats) {
  codes_of <- function(pat) if (!is.null(vl)) unname(vl[grepl(pat, names(vl))]) else numeric(0)
  list(dk  = union(codes_of(pats$dk),  parse_codes(spec$dk_add)),
       ref = union(codes_of(pats$ref), parse_codes(spec$ref_add)),
       nr  = union(codes_of(pats$nr),  parse_codes(spec$nr_add)),
       nap = union(codes_of(pats$nap), parse_codes(spec$nap_add)))
}

## filter "VAR=code" / "not:VAR=code": returns a logical vector (TRUE = keep) given a column getter
filter_mask <- function(filter, n, getnum) {
  filter <- trimws(filter); if (!nzchar(filter)) return(rep(TRUE, n))
  neg <- startsWith(filter, "not:"); f <- sub("^not:", "", filter)
  v <- sub("=.*$", "", f); code <- as.numeric(sub("^.*=", "", f))
  x <- getnum(v); if (is.null(x)) stop("filter variable not in file: ", v, call. = FALSE)
  eq <- !is.na(x) & x == code
  if (neg) !eq else eq
}

## clock time in decimal hours from the AM/PM (1 = am, 2 = pm), hour (0-11) and minute components.
## kind: wake, leave — as recorded; return — +24 when earlier than the leave-home time of the same person
## (when that time is known); bed — hours 0-11 are after midnight (+24), 12 pm (Y = 0, X = 2) is read as
## midnight (24); bed_s12 — the same with 12 pm set to missing (sensitivity).
## (EC_ITEMS_REVIEW_ja.md §4.2; the rule is the same at both waves and for both cohorts.)
build_clock <- function(X, Y, Z, kind, leave = NULL) {
  bad <- function(v) is.na(v) | v %in% CLOCK_SPECIAL
  ok <- !(bad(X) | bad(Y) | bad(Z)) & X %in% c(1, 2) & Y >= 0 & Y <= 11 & Z >= 0 & Z <= 59
  h <- ifelse(ok, Y + 12 * (X == 2), NA_real_); t <- h + Z / 60
  if (kind == "return") {
    if (!is.null(leave)) { add <- !is.na(t) & !is.na(leave) & t < leave; t[add] <- t[add] + 24 }
  } else if (kind %in% c("bed", "bed_s12")) {
    early <- !is.na(h) & h <= 11; t[early] <- t[early] + 24
    noon <- !is.na(h) & h == 12
    if (kind == "bed") t[noon] <- 24 + Z[noon] / 60 else t[noon] <- NA_real_
  } else if (!(kind %in% c("wake", "leave"))) stop("unknown clock kind: ", kind, call. = FALSE)
  t
}

## months from a years and a months component (both must be present; special codes to missing)
build_months <- function(Y, M) {
  bad <- function(v) is.na(v) | v %in% CLOCK_SPECIAL | v < 0
  ifelse(bad(Y) | bad(M), NA_real_, Y * 12 + M)
}

## parse construct "clock:X,Y,Z,kind[,leave_item]" / "months:Y,M"
parse_construct <- function(s) {
  s <- trimws(s); if (!nzchar(s)) return(NULL)
  if (s == "indicators") return(list(type = "indicators"))
  typ <- sub(":.*$", "", s); parts <- trimws(strsplit(sub("^[a-z]+:", "", s), ",")[[1]])
  if (typ == "clock") { if (length(parts) < 4) stop("clock construct needs X,Y,Z,kind: ", s, call. = FALSE)
    return(list(type = "clock", X = parts[1], Y = parts[2], Z = parts[3], kind = parts[4], leave = if (length(parts) >= 5) parts[5] else NULL)) }
  if (typ == "months") { if (length(parts) != 2) stop("months construct needs Y,M: ", s, call. = FALSE); return(list(type = "months", Y = parts[1], M = parts[2])) }
  stop("unknown construct: ", s, call. = FALSE)
}

## 0/1 indicator matrix of a nominal item: one column per substantive code; missing where x is not a
## substantive code (special codes and true missing alike). Column names: <var>__<code>.
make_indicators <- function(x, codes, var) {
  out <- matrix(NA_real_, length(x), length(codes), dimnames = list(NULL, paste0(var, "__", codes)))
  ok <- !is.na(x) & x %in% codes
  for (j in seq_along(codes)) out[ok, j] <- as.numeric(x[ok] == codes[j])
  out
}

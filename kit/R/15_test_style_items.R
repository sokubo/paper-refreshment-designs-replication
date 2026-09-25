# ============================================================
# P1 JLPS kit — unit test of the response-style helpers (15_style_items.R). No data are read.
#   Rscript R/15_test_style_items.R       (exits with status 1 on any failure)
# Cases: the exclusion-reason update for zero, one and two items with entry-wave values outside 1..k (the
# statement 15 uses, applied to a constructed audit table); the reading of the item list, including malformed
# lists that must stop the run (k not a whole number, k out of range, duplicated item, unknown set, a midpoint
# declared on an even scale, an unreadable midpoint flag, an unquoted comma in a note, a missing file); the
# verification of codes and midpoint labels; the range check of observed values; and the person-level shares.
# ============================================================
suppressMessages(library(data.table))
.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("R") })
source(file.path(.here, "15_style_items.R"))

ok <- 0L; bad <- 0L
check <- function(what, cond) {
  hit <- isTRUE(tryCatch(cond, error = function(e) { cat("  error:", conditionMessage(e), "\n"); FALSE }))
  if (hit) ok <<- ok + 1L else bad <<- bad + 1L
  cat(sprintf("%-92s %s\n", what, if (hit) "ok" else "FAIL"))
}
## an error must occur, and its message must match the reason the case is about
stops <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  inherits(e, "error") && grepl(pattern, conditionMessage(e))
}
ADD <- "entry-wave values outside 1..k observed: outside the common battery"

## 1. the reason helper, vectorised
check("reason: no item -> empty result", identical(style_join_reason(character(0), ADD), character(0)))
check("reason: one item without an earlier reason", identical(style_join_reason("", ADD), ADD))
check("reason: one item with an earlier reason", identical(style_join_reason("x", ADD), paste0("x; ", ADD)))
check("reason: two items, one with an earlier reason", identical(style_join_reason(c("", "x"), ADD), c(ADD, paste0("x; ", ADD))))
check("reason: a missing earlier reason is treated as none", identical(style_join_reason(NA_character_, ADD), ADD))

## 2. the statement of 15 on a constructed audit table: zero, one and two items with invalid entry-wave values
audit_with <- function(bad_items) {
  audit <- data.table(var_file = c("DQ1", "DQ2", "DQ3"), reason = c("", "middle category not labelled as neutral", ""))
  st_ok <- data.table(var_file = c("DQ1", "DQ2", "DQ3"), w1_values_ok = !(c("DQ1", "DQ2", "DQ3") %in% bad_items))
  audit[var_file %in% st_ok[w1_values_ok == FALSE, var_file], reason := style_join_reason(reason, ADD)]
  audit$reason
}
check("audit: no invalid item leaves every reason unchanged",
      identical(audit_with(character(0)), c("", "middle category not labelled as neutral", "")))
check("audit: one invalid item", identical(audit_with("DQ1"), c(ADD, "middle category not labelled as neutral", "")))
check("audit: two invalid items (the case that stopped v0.8)",
      identical(audit_with(c("DQ1", "DQ2")), c(ADD, paste0("middle category not labelled as neutral; ", ADD), "")))
check("audit: three invalid items", identical(audit_with(c("DQ1", "DQ2", "DQ3"))[3], ADD))

## 3. reading the list
tmp <- tempfile(fileext = ".csv")
write_list <- function(lines) { writeLines(c("var,set,battery,k,midpoint,note", lines), tmp, useBytes = TRUE); tmp }
good <- c("DQ25A,rating,agree5,5,TRUE,agreement", "DQ04_1A,rating,agree4,4,FALSE,applies", "DQ07A,frequency,freq5,5,FALSE,frequency")
st <- read_style_items(write_list(good))
check("list: a well-formed list is read with integer k", nrow(st) == 3L && is.integer(st$k) && identical(st$k, c(5L, 4L, 5L)))
check("list: item names are trimmed and upper-cased", identical(read_style_items(write_list(c(" dq25a ,rating,agree5,5,TRUE,x")))$var, "DQ25A"))
check("list: k = 5.9 is rejected (not read as 5)", stops(read_style_items(write_list(c(good[1:2], "DQ07A,frequency,freq5,5.9,FALSE,f"))), "k is not a whole number for DQ07A"))
check("list: k = 5.0 is accepted as 5", read_style_items(write_list(c("DQ25A,rating,agree5,5.0,TRUE,x")))$k == 5L)
check("list: a non-numeric k is rejected", stops(read_style_items(write_list(c(good[1:2], "DQ07A,frequency,freq5,five,FALSE,f"))), "k is not a whole number for DQ07A"))
check("list: an empty k is rejected", stops(read_style_items(write_list(c(good[1:2], "DQ07A,frequency,freq5,,FALSE,f"))), "k is not a whole number for DQ07A"))
check("list: k = 2 is rejected", stops(read_style_items(write_list(c("DQ25A,rating,agree5,2,FALSE,x"))), "k is outside 3..9 for DQ25A"))
check("list: k = 10 is rejected", stops(read_style_items(write_list(c("DQ25A,rating,agree5,10,FALSE,x"))), "k is outside 3..9 for DQ25A"))
check("list: a duplicated item is rejected", stops(read_style_items(write_list(c(good, "DQ25A,rating,agree5,5,TRUE,again"))), "duplicated item for DQ25A"))
check("list: an unknown set is rejected", stops(read_style_items(write_list(c("DQ25A,ordinal,agree5,5,TRUE,x"))), "set is neither rating nor frequency for DQ25A"))
check("list: a midpoint on an even scale is rejected", stops(read_style_items(write_list(c("DQ04_1A,rating,agree4,4,TRUE,x"))), "a midpoint is declared on an even scale for DQ04_1A"))
check("list: an unreadable midpoint flag is rejected", stops(read_style_items(write_list(c("DQ25A,rating,agree5,5,maybe,x"))), "midpoint is not TRUE or FALSE for DQ25A"))
check("list: an unquoted comma in a note is rejected", stops(read_style_items(write_list(c(good[1], "DQ04_1A,rating,agree4,4,FALSE,applies, or not", good[3]))), "Stopped early|Expected|fields|rows read but"))
check("list: a quoted comma in a note is accepted", nrow(read_style_items(write_list(c(good[1], "DQ04_1A,rating,agree4,4,FALSE,\"applies, or not\"", good[3])))) == 3L)
check("list: a missing file is rejected", stops(read_style_items(file.path(tempdir(), "no_such_list.csv")), "style item list not found"))
shipped <- read_style_items(file.path(.here, "15_style_items.csv"))
check("list: the shipped list has 46 rating and 26 frequency items",
      nrow(shipped) == 72L && sum(shipped$set == "rating") == 46L && sum(shipped$set == "frequency") == 26L)
check("list: the shipped list declares 28 midpoints, all on odd scales",
      sum(shipped$midpoint) == 28L && all(shipped$k[shipped$midpoint] %% 2L == 1L))

## 4. verification of codes and midpoint labels (special codes 8, 9 removed first)
lab5 <- c("そう思う" = 1, "どちらかといえばそう思う" = 2, "どちらともいえない" = 3, "どちらかといえばそう思わない" = 4, "そう思わない" = 5,
          "わからない" = 8, "無回答" = 9)
v <- style_verify(lab5, spec = c(8, 9), k = 5L, midpoint = TRUE)
check("verify: codes 1..5 with a neutral middle label", v$ok && v$mid_ok)
v <- style_verify(lab5[c(1:4, 6, 7)], spec = c(8, 9), k = 5L, midpoint = TRUE)
check("verify: four codes where five are declared are rejected", !v$ok && !v$mid_ok)
lab_nn <- lab5; names(lab_nn)[3] <- "ややそう思う"
v <- style_verify(lab_nn, spec = c(8, 9), k = 5L, midpoint = TRUE)
check("verify: a middle label that is not neutral keeps the item for extremes only", v$ok && !v$mid_ok && nzchar(v$reason))
lab0 <- setNames(0:4, names(lab5)[1:5])
check("verify: codes 0..4 are rejected", !style_verify(lab0, spec = c(8, 9), k = 5L, midpoint = FALSE)$ok)
check("verify: an item without value labels is rejected", !style_verify(NULL, spec = c(8, 9), k = 5L, midpoint = FALSE)$ok)

## 5. observed values
check("values: 1..k with missing values pass", style_values_ok(c(1, 5, NA, 3), 5L))
check("values: a value outside 1..k fails", !style_values_ok(c(1, 7, NA), 5L))
check("values: an item with no observed value passes", style_values_ok(c(NA_real_, NA_real_), 5L))

## 6. person-level shares: denominators are the items a person answered
M <- cbind(A = c(1, 3, NA), B = c(5, 2, NA), C = c(4, NA, NA))
s <- style_shares(M, ext_items = c("A", "B", "C"), mid_items = c("A", "B"), k_of = list(A = 5L, B = 5L, C = 5L))
check("shares: extreme share over answered items", isTRUE(all.equal(s$ext, c(2/3, 0, NA))))
check("shares: midpoint share over answered midpoint items", isTRUE(all.equal(s$mid, c(0, 1/2, NA))))

cat(sprintf("\n%d checks: %d passed, %d failed\n", ok + bad, ok, bad))
if (bad > 0L) quit(status = 1L)

#!/usr/bin/env Rscript
# ============================================================
# 合成マスタ用の変数対応表(w1 ↔ w5)と項目指定表を作る。公開パイプラインを自己完結させるための補助。
# 実データでは、提供元の変数一覧から作った jlps_docs/varmap_w1-19.csv(この archive には同梱しない——第三者の
# 文書だから)と、著者が調査票と照合した R/15_item_scale.csv を使う。合成入力に対しては、同じ形の対応表と、
# 合成項目の行を加えた項目指定表をここで生成する。
# 実行: P1_DATA_DIR=<合成データの置き場> Rscript R/15_make_synthetic_varmap.R
#   → <P1_DATA_DIR>/varmap_synthetic.csv と <P1_DATA_DIR>/item_scale_synthetic.csv を書く。
#     以後 P1_VARMAP と P1_ITEM_SCALE でそれらを指す。
# ============================================================
suppressMessages({library(haven); library(data.table)})
.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })
DATA_DIR <- path.expand(Sys.getenv("P1_DATA_DIR", "/tmp/p1synth"))
f <- Sys.glob(file.path(DATA_DIR, "raw", "*.dta"))
if (!length(f)) stop("no synthetic .dta under ", file.path(DATA_DIR, "raw"))
d0 <- read_dta(f[1]); nm <- names(d0)
## --- variable map: 波接頭辞: Z = w1, A = w2, B = w3, C = w4, D = w5, ... (15_make_synthetic.R と同じ規約) ---
pre <- c(w1 = "Z", w2 = "A", w3 = "B", w4 = "C", w5 = "D", w6 = "E", w7 = "F", w8 = "G", w9 = "H", w10 = "I")
strip <- function(x, p) sub(paste0("^", p), "", x)
base <- unique(unlist(lapply(pre, function(p) strip(grep(paste0("^", p, "Q"), nm, value = TRUE), p))))
base <- base[nzchar(base)]
vm <- data.table(varname = paste0("Q", base))
for (w in names(pre)) set(vm, j = w, value = ifelse(paste0(pre[[w]], base) %in% nm, paste0(pre[[w]], base), ""))
for (w in paste0("w", 11:19)) set(vm, j = w, value = "")
set(vm, j = "exists_w1_19", value = "")
## items whose 2007 counterpart has a different name (as in the provider's map): the same pairs as in the real map
pairs <- list(c("DQ57B", "ZQ02B"), c("DQ57C", "ZQ02C"), c("DQ03_13", "ZQ05_5Y"), c("DQ54_1Y", "ZQ59_1Y"), c("DQ54_1M", "ZQ59_1M"),
              c("DQ03_1", "JC_1"), c("DQ03_2", "JC_2"))
for (p in pairs) if (all(p %in% nm)) { i <- which(vm$w5 == p[1]); if (length(i)) vm$w1[i] <- p[2] else vm <- rbind(vm, data.table(varname = p[1], w1 = p[2], w5 = p[1]), fill = TRUE) }
vm[is.na(vm)] <- ""
setcolorder(vm, c("varname", paste0("w", 1:19), "exists_w1_19"))
out <- file.path(DATA_DIR, "varmap_synthetic.csv")
fwrite(vm, out)
cat("wrote", out, ":", nrow(vm), "rows;", sum(vm$w1 != "" & vm$w5 != ""), "items present at both w1 and w5\n")

## --- item specification table for the synthetic input ------------------------------------------------------
## The rows of the real table whose items exist in the synthetic file (the items added to exercise the rules, the
## response markers and the derived items) are kept as they are; the generic synthetic items DQ101.. get a scale
## from their value labels (5 substantive labels: ordinal; 2: binary; none: continuous) and no entry-wave class
## (their entry wave is compared as in the real run whenever the codes match).
real <- fread(file.path(.here, "15_item_scale.csv"), encoding = "UTF-8", colClasses = "character", na.strings = NULL)
real[is.na(real)] <- ""
dq <- grep("^DQ[0-9]", nm, value = TRUE)
keep_real <- real[toupper(var) %in% toupper(dq) | grepl("^(clock|months):", construct)]
generic <- setdiff(dq, toupper(keep_real$var)); generic <- setdiff(generic, real$var)
spec_pat <- "わからない|分からない|ＤＫ|DK|答えたくない|回答したくない|無回答|不明|非該当"
gen_rows <- rbindlist(lapply(generic, function(v) {
  vl <- attr(d0[[v]], "labels", exact = TRUE)
  k <- if (is.null(vl)) 0L else sum(!grepl(spec_pat, names(vl)))
  sc <- if (k == 0L) "continuous" else if (k == 2L) "binary" else "ordinal"
  lab <- attr(d0[[v]], "label", exact = TRUE); if (is.null(lab)) lab <- ""
  data.table(var = v, label = as.character(lab), scale = sc, codes_expected = "", nap_add = "", dk_add = "", nr_add = "", ref_add = "", recode = "",
             construct = "", filter = "", in_counts = "TRUE", ec_class = "", ec_subgroup = "", ec_main = "TRUE",
             note = "synthetic item: scale from its value labels", source = "15_make_synthetic_varmap.R")
}))
## the synthetic items with an entry wave are treated as class A (comparable) so that the main entry-wave analysis is non-empty
gen_rows[toupper(var) %in% toupper(vm[w1 != "", w5]), ec_class := "A"]
## one synthetic item is put in class D (an entry wave that is not comparable): its entry-wave estimators are computed
## but stay outside the main entry-wave counts (they enter 15_ec_scope_sensitivity.csv)
gen_rows[var == "DQ103", `:=`(ec_class = "D", ec_main = "FALSE", note = "synthetic item: entry-wave class D (outside the main entry-wave analysis)")]
tab <- rbind(keep_real[, names(gen_rows), with = FALSE], gen_rows)
out2 <- file.path(DATA_DIR, "item_scale_synthetic.csv")
fwrite(tab, out2)
cat("wrote", out2, ":", nrow(tab), "rows (", nrow(keep_real), "from the real table,", nrow(gen_rows), "generic synthetic items )\n")

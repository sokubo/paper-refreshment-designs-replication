#!/usr/bin/env Rscript
# ============================================================
# 合成マスタ用の変数対応表(w1 ↔ w5)を作る。公開パイプラインを自己完結させるための補助。
# 実データでは、提供元の変数一覧から作った jlps_docs/varmap_w1-19.csv を使う(この archive には
# 同梱しない——第三者の文書だから)。合成入力に対しては、同じ形の対応表をここで生成する。
# 実行: P1_DATA_DIR=<合成データの置き場> Rscript R/15_make_synthetic_varmap.R
#   → <P1_DATA_DIR>/varmap_synthetic.csv を書く。以後 P1_VARMAP でそれを指す。
# ============================================================
suppressMessages({library(haven); library(data.table)})
DATA_DIR <- path.expand(Sys.getenv("P1_DATA_DIR", "/tmp/p1synth"))
f <- Sys.glob(file.path(DATA_DIR, "raw", "*.dta"))
if (!length(f)) stop("no synthetic .dta under ", file.path(DATA_DIR, "raw"))
nm <- names(read_dta(f[1], n_max = 0))
# 波接頭辞: Z = w1, A = w2, B = w3, C = w4, D = w5, ... (15_make_synthetic.R と同じ規約)
pre <- c(w1 = "Z", w2 = "A", w3 = "B", w4 = "C", w5 = "D", w6 = "E", w7 = "F", w8 = "G", w9 = "H", w10 = "I")
strip <- function(x, p) sub(paste0("^", p), "", x)
base <- unique(unlist(lapply(pre, function(p) strip(grep(paste0("^", p, "Q"), nm, value = TRUE), p))))
base <- base[nzchar(base)]
vm <- data.table(varname = paste0("Q", base))
for (w in names(pre)) vm[[w]] <- ifelse(paste0(pre[[w]], base) %in% nm, paste0(pre[[w]], base), "")
for (w in paste0("w", 11:19)) vm[[w]] <- ""
vm[, exists_w1_19 := ""]
setcolorder(vm, c("varname", paste0("w", 1:19), "exists_w1_19"))
out <- file.path(DATA_DIR, "varmap_synthetic.csv")
fwrite(vm, out)
cat("wrote", out, ":", nrow(vm), "rows;",
    sum(vm$w1 != "" & vm$w5 != ""), "items present at both w1 and w5\n")

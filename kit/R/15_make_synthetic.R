# ============================================================
# P1 JLPS kit — 15s: kit 15 の動作確認用 合成マスタ生成(個票なし・構造テスト専用)
# JLPS 統合 wide の構造(CN / sex / ybirth / 波接頭辞 Z=w1 … I=w10 / 値ラベル付き項目 /
# 応答マーカー / DQ74M 回答月 / 負対照 ZQ16-21,ZQ18_A-T ↔ DQ60-65,DQ62_A-T)を模した
# .dta を <P1_DATA_DIR>/raw/ に書く。数値は論文に使わない。
# 実行(クラウド/ローカルどちらでも):
#   P1_DATA_DIR=/tmp/p1synth P1_PROJECT_DIR=<P1ルート> Rscript analysis/R/15_make_synthetic.R
#   → 続けて 04 と 15 を同じ環境変数で実行すると E2E が回る。
# ============================================================
suppressMessages({library(haven); library(data.table)})
set.seed(20260915)
DATA_DIR <- path.expand(Sys.getenv("P1_DATA_DIR", "/tmp/p1synth"))
dir.create(file.path(DATA_DIR, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "derived_me"), showWarnings = FALSE)

n_cont <- 4800L; n_add <- 963L; n <- n_cont + n_add
cn <- c(rep(1L, n_cont), rep(2L, n_add))
sex <- sample(1:2, n, TRUE); yb <- sample(1966:1986, n, TRUE)
u <- rnorm(n)                                   # 固定形質(脱落と結果に共通)
educ <- sample(1:6, n, TRUE)

## 応答過程: 継続 w1=1、w2..w10 は形質依存(MNAR_trait 型)+わずかな状態依存
plogis_s <- function(a, b) plogis(a + b * u + 0.15 * rnorm(n))
R <- matrix(0L, n, 10)
R[cn == 1, 1] <- 1L
for (w in 2:10) { pr <- plogis_s(1.9, 0.45); R[, w] <- as.integer(runif(n) < pr) }
R[cn == 2, 1:4] <- 0L; R[cn == 2, 5] <- 1L     # 追加標本は w5 初回(全員回答)
## 継続の w5 回答者(=04 のリスク集合)は約 65%
cat("cont responded w1-5 all:", mean(rowSums(R[cn == 1, 1:5]) == 5), " w5:", mean(R[cn == 1, 5]), "\n")
cat("add responded w6-9 all:", mean(rowSums(R[cn == 2, 6:9]) == 4), "\n")

## 項目: 30 順序(5 件法)、5 二値、5 連続; うち 35 は w1 に同一コーディングの対応あり、5 は w5 のみ
## 条件付け効果(真値): 項目 1-3 に +0.30SD(継続 w5)、項目 4 に −0.30SD、他 0; DK 率は継続で低い
J_ord <- 30; J_bin <- 5; J_con <- 5
tau_true <- c(0.30, 0.30, 0.30, -0.30, rep(0, J_ord + J_bin + J_con - 4))
lab_ord <- c("そう思う" = 1, "どちらかといえばそう思う" = 2, "どちらともいえない" = 3,
             "どちらかといえばそう思わない" = 4, "そう思わない" = 5, "わからない" = 8, "無回答" = 9)
lab_bin <- c("はい" = 1, "いいえ" = 2, "無回答" = 9)
mk_ord <- function(base, tau, dk_rate) {
  z <- base + tau + rnorm(n) * 0.9
  v <- pmin(pmax(round(3 + z), 1), 5)
  dk <- runif(n) < dk_rate; nr <- runif(n) < 0.02
  v[dk] <- 8; v[nr] <- 9; v
}
d <- data.table(CN = cn, sex = sex, ybirth = yb)
items_w5 <- list(); items_w1 <- list()
for (j in seq_len(J_ord + J_bin + J_con)) {
  base <- 0.6 * u + 0.3 * rnorm(n)             # 項目固有の安定成分
  nm5 <- sprintf("DQ%03d", j + 100)             # DQ101.. (マーカー DQ02/DQ43、負対照 DQ60-65、除外パターン DQ59-74 と衝突しない)
  nm1 <- sprintf("ZQ%03d", j + 100)
  if (j <= J_ord) {
    v5 <- mk_ord(base, ifelse(cn == 1, tau_true[j], 0) * 1.1, dk_rate = ifelse(cn == 1, 0.03, 0.06))
    v1 <- mk_ord(base, 0, dk_rate = 0.05)
    items_w5[[nm5]] <- labelled(v5, lab_ord, label = paste0("w5問", j, "_意識項目", j))
    if (j <= 25) items_w1[[nm1]] <- labelled(v1, lab_ord, label = paste0("問", j, "_意識項目", j))
    if (j == 26) {  # コーディング不一致の例(w1 は 4 件法)
      items_w1[[nm1]] <- labelled(pmin(v1, 4), c("そう思う" = 1, "やや" = 2, "あまり" = 3, "思わない" = 4, "わからない" = 8, "無回答" = 9),
                                  label = paste0("問", j, "_意識項目", j))
    }
    # j 27-30: w1 対応なし
  } else if (j <= J_ord + J_bin) {
    p1 <- plogis(base + ifelse(cn == 1, tau_true[j], 0))
    v5 <- ifelse(runif(n) < p1, 1, 2); v5[runif(n) < 0.02] <- 9
    v1 <- ifelse(runif(n) < plogis(base), 1, 2); v1[runif(n) < 0.02] <- 9
    items_w5[[nm5]] <- labelled(v5, lab_bin, label = paste0("w5問", j, "_事実項目", j))
    items_w1[[nm1]] <- labelled(v1, lab_bin, label = paste0("問", j, "_事実項目", j))
  } else {
    v5 <- round(40 + 8 * base + ifelse(cn == 1, tau_true[j], 0) * 8 + rnorm(n) * 4); v5[runif(n) < 0.03] <- NA
    v1 <- round(40 + 8 * base + rnorm(n) * 4); v1[runif(n) < 0.03] <- NA
    items_w5[[nm5]] <- labelled(as.numeric(v5), NULL, label = paste0("w5問", j, "_労働時間", j))
    items_w1[[nm1]] <- labelled(as.numeric(v1), NULL, label = paste0("問", j, "_労働時間", j))
  }
}
for (nm in names(items_w5)) d[[nm]] <- items_w5[[nm]]
for (nm in names(items_w1)) d[[nm]] <- items_w1[[nm]]

## 負対照(時不変): 継続=Z 版、追加=D 版
nc_z <- c("ZQ16","ZQ17","ZQ19","ZQ21", paste0("ZQ18_", c("A","B","C","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T")))
nc_d <- c("DQ60","DQ61","DQ63","DQ65", paste0("DQ62_", c("A","B","C","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T")))
for (i in seq_along(nc_z)) {
  v <- sample(1:4, n, TRUE, prob = c(.4, .3, .2, .1))
  d[[nc_z[i]]] <- ifelse(cn == 1, v, NA_real_); d[[nc_d[i]]] <- ifelse(cn == 2, v, NA_real_)
}
## 応答マーカー(J3 strict と 04 の定義の両方を満たす): 各波の「就業有無」+ w5 の DQ74M/DQ43、w1 の ZQ23A/ w5 DQ69A(学歴)
mk_marker <- function(w, nm) { v <- ifelse(R[, w] == 1, sample(1:2, n, TRUE), NA_real_); d[[nm]] <<- labelled(v, c("いる" = 1, "いない" = 2), label = paste0("w", w, "_就業")) }
mk_marker(1, "ZQ03"); mk_marker(2, "AQ02"); mk_marker(3, "BQ02"); mk_marker(4, "CQ02"); mk_marker(5, "DQ02")
mk_marker(6, "EQ02"); mk_marker(7, "FQ02"); mk_marker(8, "GQ02"); mk_marker(9, "HQ02"); mk_marker(10, "IQ02")
d$DQ74M <- ifelse(R[, 5] == 1, sample(1:3, n, TRUE), NA_real_)
d$DQ43  <- labelled(ifelse(R[, 5] == 1, sample(1:2, n, TRUE), NA_real_), c("既婚" = 1, "未婚" = 2), label = "w5問43_婚姻")
d$ZQ50  <- labelled(ifelse(cn == 1, sample(1:2, n, TRUE), NA_real_), c("未婚" = 1, "既婚" = 2), label = "問50_婚姻")   # w1 は逆転コーディング
d$ZQ23A <- ifelse(cn == 1, educ, NA_real_); d$DQ69A <- ifelse(cn == 2, educ, NA_real_)
## 非回答者は w5 項目を欠測に(DQ02 系以外)
w5vars <- grep("^DQ", names(d), value = TRUE)
for (v in w5vars) { x <- d[[v]]; x[R[, 5] == 0] <- NA; d[[v]] <- x }
## 継続の w1 項目は全員観測(入所波)、追加は w1 なし
w1vars <- grep("^ZQ", names(d), value = TRUE)
for (v in w1vars) { x <- d[[v]]; x[cn == 2] <- NA; d[[v]] <- x }

f <- file.path(DATA_DIR, "raw", "SYNTH_ZQ100AQ100BQ100CQ100DQ100EQ100FQ100GQ100HQ100IQ100.dta")
write_dta(d, f)
## write_dta() stamps the Stata header with the current date and time (to the minute), so two runs a
## minute apart differ in four bytes and nothing else. Fix the stamp to that of the shipped file so that
## the output is byte-reproducible and its SHA-256 can be checked against synthetic/ (see KIT_README.md).
local({
  b <- readBin(f, "raw", file.info(f)$size); tag <- charToRaw("<timestamp>"); m <- length(tag)
  i <- which(vapply(seq_len(512L - m), function(k) all(b[k:(k + m - 1L)] == tag), logical(1)))[1]
  stopifnot(!is.na(i), as.integer(b[i + m]) == 17L)
  b[(i + m + 1L):(i + m + 17L)] <- charToRaw("21 Sep 2026 10:52")
  writeBin(b, f)
})
cat("wrote", f, ":", nrow(d), "x", ncol(d), "\n")
cat("true tau (items 1-4):", tau_true[1:4], " (others 0); DK rate cont .03 vs add .06\n")

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
## v1.0 (2026-09-28): 項目指定表(R/15_item_scale.csv)の規則を通す項目を、実データと同じ変数名で加える(数値は論文に使わない)
##   DQ15 (1・2 をまとめる再符号化) / DQ30 (名義 → 指標; 7=8 の再符号化; 10 = わからない) / DQ42_A, DQ42_B, DQ42_G (「わからない」の指標で欠損に)
##   DQ03_1・JC_1 (非該当追加 10;11 と 12=3) / DQ03_2・JC_2 (対応表の map:10=8; 区分 C の該当者 ZQ03 = 1)
##   DQ57A-D と構成要素 (時刻の派生; w1 は ZQ02A-D) / DQ54_1Y・M と ZQ59_1Y・M (月数の派生; 区分 C の該当者 ZQ58)
##   DQ03_5E・DQ03_5ES (フィルター) / DQ03_13 と ZQ05_5Y (区分 C の該当者 ZQ03 = 1) / DQ17A (尺度外コード 4 を非該当に)
lab_smoke <- c("喫煙したことがない" = 1, "禁煙した" = 2, "1～10本" = 3, "11～20本" = 4, "21本以上" = 5, "無回答" = 9)
mk_smoke <- function(tau) { z <- 0.5 * u + rnorm(n); v <- cut(z + tau, c(-Inf, -0.6, -0.1, 0.4, 1.0, Inf), labels = FALSE); v[runif(n) < 0.02] <- 9; v }
d$DQ15 <- labelled(as.numeric(mk_smoke(ifelse(cn == 1, 0.2, 0))), lab_smoke, label = "w5問15_喫煙")
d$ZQ15 <- labelled(as.numeric(mk_smoke(0)), lab_smoke, label = "問15_喫煙")
## party support: the real 9 codes; code 7 is a different party at the two waves (2011 みんなの党, 2007 新党日本), which
## the table folds into 8 = その他の政党 at both waves (recode 7=8)
lab_party <- c("自民党" = 1, "民主党" = 2, "公明党" = 3, "共産党" = 4, "社民党" = 5, "国民新党" = 6, "みんなの党" = 7,
               "その他の政党" = 8, "特に支持する政党はない" = 9, "わからない" = 10, "無回答" = 99)
lab_party_w1 <- lab_party; names(lab_party_w1)[7] <- "新党日本"
mk_party <- function(shift) { pr <- c(.25, .20, .05, .03, .02, .01, .04, .02, .38); v <- sample(1:9, n, TRUE, prob = pr); v[runif(n) < .08 + shift] <- 10; v[runif(n) < .02] <- 99; v }
d$DQ30 <- labelled(as.numeric(mk_party(-0.03)), lab_party, label = "w5問30_支持政党")
d$ZQ30 <- labelled(as.numeric(mk_party(0)), lab_party_w1, label = "問30_支持政党")
## work status (問3(1) / 2007 問4A(1)): the provider's 12 codes (10, 11 = not working, from the routing of 問2; 12 = students
## working non-regularly, split out of 3); the table sets 10, 11 to not applicable and 12 back to 3 (codes 1-9 as on the
## questionnaire). occupation (問3(2) / 2007 問4A(2)): 8 precodes; 2007 also has 9 = わからない and an after-code 10 = 農林,
## which 15_entry_overrides.csv folds into 8. The 2007 counterparts have the provider's names JC_1 and JC_2.
lab_work <- c("経営者、役員" = 1, "正社員・正職員" = 2, "パート・アルバイト・契約・臨時・嘱託" = 3, "派遣社員" = 4, "請負社員" = 5,
              "自営業主、自由業者" = 6, "家族従業者" = 7, "内職" = 8, "その他" = 9, "無職(学生は除く)" = 10, "学生(働いていない)" = 11,
              "学生(現在非正規で働いている)" = 12, "非該当" = 88, "無回答" = 99)
lab_work_w1 <- lab_work; names(lab_work_w1)[7] <- "家族従事者"
mk_work <- function(emp) {                       # emp: the wave's 就業 marker (1 = working, 2 = not, NA = not interviewed)
  v <- sample(c(1:9, 12), n, TRUE, prob = c(.03, .55, .18, .05, .02, .06, .02, .01, .02, .06))
  v[!is.na(emp) & emp == 2] <- sample(10:11, n, TRUE, prob = c(.75, .25))[!is.na(emp) & emp == 2]
  v[!is.na(emp) & runif(n) < .01] <- 99; v[is.na(emp)] <- NA; v
}
numv <- function(x) as.numeric(zap_labels(x))
d$DQ03_1 <- labelled(as.numeric(mk_work(numv(d$DQ02))), lab_work, label = "w5問3(1)_現職・働き方")
d$JC_1   <- labelled(as.numeric(mk_work(numv(d$ZQ03))), lab_work_w1, label = "問4A(1)_現職・働き方")
lab_occ <- c("専門職・技術職" = 1, "管理職" = 2, "事務職" = 3, "販売職" = 4, "サービス職" = 5, "生産現場職・技能職" = 6,
             "運輸・保安職" = 7, "その他" = 8, "非該当" = 88, "無回答" = 99)
lab_occ_w1 <- c(lab_occ[1:8], "わからない" = 9, "農林" = 10, "非該当" = 88, "無回答" = 99)
mk_occ <- function(emp, w1 = FALSE) {
  v <- sample(1:8, n, TRUE, prob = c(.20, .05, .25, .10, .12, .18, .05, .05))
  if (w1) { v[runif(n) < .02] <- 9; v[runif(n) < .01] <- 10 }
  v[!is.na(emp) & emp == 2] <- 88; v[!is.na(emp) & runif(n) < .01] <- 99; v[is.na(emp)] <- NA; v
}
d$DQ03_2 <- labelled(as.numeric(mk_occ(numv(d$DQ02))), lab_occ, label = "w5問3(2)_現職・職業―大分類(プリコード)")
d$JC_2   <- labelled(as.numeric(mk_occ(numv(d$ZQ03), w1 = TRUE)), lab_occ_w1, label = "問4A(2)_現職・職業―大分類(プリコード)")
lab_sel <- c("選択" = 1, "非選択" = 2, "無回答" = 9)
ins <- sample(1:3, n, TRUE, prob = c(.5, .4, .1))      # 1 = A, 2 = B, 3 = わからない(G)
d$DQ42_A <- labelled(as.numeric(ifelse(ins == 1, 1, 2)), lab_sel, label = "w5問42_健康保険―A")
d$DQ42_B <- labelled(as.numeric(ifelse(ins == 2, 1, 2)), lab_sel, label = "w5問42_健康保険―B")
d$DQ42_G <- labelled(as.numeric(ifelse(ins == 3, 1, 2)), lab_sel, label = "w5問42_健康保険―わからない")
## clock times: format code (1 = a time is given, 2 = not fixed, 3 = mainly at home for B and C) and the AM/PM, hour, minute components
lab_fmt2 <- c("だいたい午前／午後×時○分ころ" = 1, "特に決まっていない" = 2, "無回答" = 9)
lab_fmt3 <- c("だいたい午前／午後×時○分ころ" = 1, "特に決まっていない" = 2, "主に家にいる" = 3, "無回答" = 9)
lab_fmt3w1 <- c("だいたい×時○分ころ" = 1, "特に決まっていない" = 2, "主に家にいる" = 3, "無回答" = 9)
lab_ap <- c("午前" = 1, "午後" = 2, "非該当" = 8, "無回答" = 9); lab_hm <- c("非該当" = 88, "無回答" = 99)
mk_clock <- function(prefix, fmt_var, kind, w1 = FALSE, three = FALSE) {
  fmt <- sample(if (three) 1:3 else 1:2, n, TRUE, prob = if (three) c(.8, .12, .08) else c(.9, .1))
  if (w1 && !three) fmt[runif(n) < .01] <- 3                        # the w1 codebook artefact (a third code for A and D)
  h24 <- switch(kind, wake = round(rnorm(n, 6.5, 1)), leave = round(rnorm(n, 8, 1)), ret = round(rnorm(n, 19, 2.5)), bed = round(rnorm(n, 23.5, 1.2)))
  h24 <- ((h24 %% 24) + 24) %% 24
  X <- ifelse(h24 >= 12, 2, 1); Y <- h24 %% 12; Z <- sample(c(0, 15, 30, 45), n, TRUE)
  X[fmt != 1] <- 8; Y[fmt != 1] <- 88; Z[fmt != 1] <- 88
  nr <- runif(n) < .01; X[nr] <- 9; Y[nr] <- 99; Z[nr] <- 99
  d[[fmt_var]] <<- labelled(as.numeric(fmt), if (three) (if (w1) lab_fmt3w1 else lab_fmt3) else lab_fmt2, label = paste0(prefix, "_時刻の回答形式"))
  d[[paste0(fmt_var, "X")]] <<- labelled(as.numeric(X), lab_ap, label = paste0(prefix, "―午前午後"))
  d[[paste0(fmt_var, "Y")]] <<- labelled(as.numeric(Y), lab_hm, label = paste0(prefix, "―時"))
  d[[paste0(fmt_var, "Z")]] <<- labelled(as.numeric(Z), lab_hm, label = paste0(prefix, "―分"))
}
mk_clock("w5問57A_起床", "DQ57A", "wake"); mk_clock("w5問57B_家を出る", "DQ57B", "leave", three = TRUE)
mk_clock("w5問57C_帰宅", "DQ57C", "ret", three = TRUE); mk_clock("w5問57D_就寝", "DQ57D", "bed")
mk_clock("問2A_起床", "ZQ02A", "wake", w1 = TRUE); mk_clock("問2B_家を出る", "ZQ02B", "leave", w1 = TRUE, three = TRUE)
mk_clock("問2C_帰宅", "ZQ02C", "ret", w1 = TRUE, three = TRUE); mk_clock("問2D_就寝", "ZQ02D", "bed", w1 = TRUE)
## duration of the relationship (years and months) and the 2007 partner status (class C subgroup)
lab_partner <- c("婚約者がいる" = 1, "特定の交際相手がいる" = 2, "現在はいない" = 3, "無回答" = 9)
d$ZQ58 <- labelled(as.numeric(sample(1:3, n, TRUE, prob = c(.1, .3, .6))), lab_partner, label = "問58_現在交際している人はいるか")
mk_dur <- function() { yrs <- pmax(0, round(rnorm(n, 2, 1.5))); mos <- sample(0:11, n, TRUE); nap <- runif(n) < .6; yrs[nap] <- 88; mos[nap] <- 88; list(y = yrs, m = mos) }
du5 <- mk_dur(); du1 <- mk_dur()
d$DQ54_1Y <- labelled(as.numeric(du5$y), lab_hm, label = "w5問54(1)_交際期間―年"); d$DQ54_1M <- labelled(as.numeric(du5$m), lab_hm, label = "w5問54(1)_交際期間―月")
d$ZQ59_1Y <- labelled(as.numeric(du1$y), lab_hm, label = "問59(1)_交際期間―年"); d$ZQ59_1M <- labelled(as.numeric(du1$m), lab_hm, label = "問59(1)_交際期間―月")
## annual-salary amount recorded for everyone (a filter keeps those who chose that pay form)
d$DQ03_5E <- labelled(as.numeric(ifelse(runif(n) < .25, 1, 2)), lab_sel, label = "w5問3(5)-5_収入の形態―年俸")
d$DQ03_5ES <- labelled(as.numeric(ifelse(d$DQ03_5E == 1, round(300 + 80 * u + rnorm(n) * 60), 0)), NULL, label = "w5問3(5)-5_収入の金額―年俸")
## continuation of the current job (asked of workers only; 2007 also of the last job) and the scale-external code of DQ17A
lab_cont <- c("続けるつもり" = 1, "やめることを考えている" = 2, "すぐにやめる" = 3, "わからない" = 4, "非該当" = 8, "無回答" = 9)
mk_cont <- function() { v <- sample(1:3, n, TRUE, prob = c(.7, .2, .1)); v[runif(n) < .05] <- 4; v[runif(n) < .02] <- 9; v }
d$DQ03_13 <- labelled(as.numeric(ifelse(d$DQ02 == 1, mk_cont(), 8)), lab_cont, label = "w5問3(13)_現在の会社で当面仕事を続けるか")
d$ZQ05_5Y <- labelled(as.numeric(mk_cont()), lab_cont, label = "問5(5)_現在の会社での仕事や事業の継続")
lab_smk_par <- c("まったく吸ったことがない" = 1, "禁煙していた" = 2, "吸っていた" = 3, "その時父・母はいなかった" = 4, "無回答" = 9)
d$DQ17A <- labelled(as.numeric({ v <- sample(1:4, n, TRUE, prob = c(.3, .2, .45, .05)); v[runif(n) < .02] <- 9; v }), lab_smk_par, label = "w5問17A_中3時の父親の喫煙")

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

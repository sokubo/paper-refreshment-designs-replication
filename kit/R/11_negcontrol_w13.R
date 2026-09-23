# ============================================================
# P1 JLPS kit — 11: w13(2019年補充)エピソードの参入時負対照
#  TP2定理3の第2バッテリー実装: 2019年新規コホートの参入時負対照(LQ系)と、
#  w13残存の旧コホートの参入時測定値(2007=ZQ系 / 2011=DQ系)を比較し、
#  β̂_2(13)(2007コホート・テニュア13) と β̂_2(9)(2011コホート・テニュア9)
#  を項目別に推定する。対応表=results/03_negcontrol_sets_v1.csv(complete3=1)。
#  推定は09のw5版と同型: 標準化差d+コホートバンド×性別の層別逆分散合成+TOST。
#  w13残存の判定: MQ系(w13)ではなくLQ系…※w13本体の回答有無で判定
#   (L接頭辞はw13。継続コホートのw13回答者=LQ系に非欠測がある者)。
# 出力はすべて集計値(N<10抑制)。
# 実行: cd <P1ルート> && Rscript analysis/R/11_negcontrol_w13.R
# ============================================================

.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
suppressMessages({library(haven); library(data.table)})

SESOI_D <- 0.10
MIN_ANSWERS_W13 <- 3    # L接頭辞の非欠測がこれ以上 → w13参加とみなす

main <- function() {
  nc <- fread(file.path(RESULTS_DIR, "03_negcontrol_sets_v1.csv"))
  nc <- nc[complete3 == 1]
  cat("diag: complete3の負対照セット =", nrow(nc), "\n")

  # 2020ウェブ特別調査ファイル(JLPSYM_online_*)が同居してもマスタを確実に選ぶ
  f <- input_dta(); fp <- input_fingerprint(f)
  cat(sprintf("input: %s (%s bytes; md5 %s; sha256 %s)\n", fp$file, fp$bytes, fp$md5, fp$sha256))
  d <- as.data.table(read_dta(f))
  gc_ <- function(nm) { hit <- names(d)[toupper(trimws(names(d))) == toupper(nm)]
    if (!length(hit)) return(NULL); d[[hit[1]]] }
  num <- function(x) if (is.null(x)) NULL else suppressWarnings(as.numeric(zap_labels(x)))
  cn <- as.integer(zap_labels(gc_("CN")))

  # 特殊コードをNA化(値ラベルからDK/無回答/非該当系を検出——04と同じ規則)
  SPEC_PAT <- "わからない|分からない|ＤＫ|DK|答えたくない|回答したくない|無回答|不明|非該当"
  clean <- function(nm) {
    x <- gc_(nm); if (is.null(x)) return(NULL)
    vl <- attr(x, "labels", exact = TRUE)
    xv <- num(x)
    if (!is.null(vl)) { spec <- unname(vl[grepl(SPEC_PAT, names(vl))])
      xv[xv %in% spec] <- NA }
    xv
  }

  # w13参加(L接頭辞の非欠測数)
  lvars <- grep("^[Ll][Qq][0-9]", names(d), value = TRUE)
  lmat <- as.matrix(d[, lapply(.SD, function(x) suppressWarnings(as.numeric(zap_labels(x)))),
                      .SDcols = lvars])
  n_l <- rowSums(!is.na(lmat))
  in13 <- n_l >= MIN_ANSWERS_W13
  cat("diag: w13参加(L系非欠測>=", MIN_ANSWERS_W13, ") CN別:\n")
  print(table(cn, in13, useNA = "ifany"))
  cat("  (CN=1: 2007継続 / CN=2: 2011追加 / CN=3(想定): 2019新規——値はconsoleで要確認)\n")

  # コホート判定: CN=1(2007), CN=2(2011), それ以外のw13参加者=2019新規と仮定
  #  ※CNの実コードが異なる場合はこのdiag出力を見て修正
  grp <- rep(NA_character_, nrow(d))
  grp[cn == 1 & in13] <- "c2007_t13"
  grp[cn == 2 & in13] <- "c2011_t9"
  grp[!(cn %in% c(1, 2)) & in13] <- "c2019_t1"
  cat("diag: 群サイズ:", paste(names(table(grp)), table(grp), collapse = " / "), "\n")

  # 層別変数: 出生年×性別(04と同じ列名: sex / ybirth)。
  # v2修正(2026-08-25): 初版はZQ03Y/ZQ02を探して失敗し無層別で走った結果、
  #   出生年組成差(Δg——エアコン・パソコン等の年代依存の子ども時代資産)が
  #   βに混入して16/23棄却の見かけを生んだ。TP2の(A1)通り、
  #   **共通サポート(出生年の重なり)に制限+出生年バンド×性別で層別**する。
  yb <- num(gc_("ybirth")); sx <- num(gc_("sex"))
  if (is.null(yb) || is.null(sx))
    stop("sex/ybirth が見つかりません(04と同じ列名のはず)。names(d)を確認してください。")
  cat("diag: 出生年レンジ 群別:\n")
  for (g in c("c2007_t13", "c2011_t9", "c2019_t1"))
    cat(sprintf("  %s: [%s, %s]\n", g,
                min(yb[grp == g], na.rm = TRUE), max(yb[grp == g], na.rm = TRUE)))

  # 共通サポート: 比較する2群の出生年の重なりに制限する関数
  support_of <- function(gA, gB) {
    lo <- max(min(yb[grp == gA], na.rm = TRUE), min(yb[grp == gB], na.rm = TRUE))
    hi <- min(max(yb[grp == gA], na.rm = TRUE), max(yb[grp == gB], na.rm = TRUE))
    !is.na(yb) & yb >= lo & yb <= hi
  }

  # 層別逆分散合成の標準化差(09と同型; 層=出生年3年バンド×性別)
  ybband <- 3 * floor(yb / 3)
  strat_all <- interaction(ybband, sx, drop = TRUE)
  sdiff <- function(y, g1, g0) {
    strat <- strat_all
    est <- 0; wsum <- 0; nn1 <- 0; nn0 <- 0
    for (st in levels(factor(strat[g1 | g0]))) {
      a <- y[g1 & strat == st]; b <- y[g0 & strat == st]
      a <- a[!is.na(a)]; b <- b[!is.na(b)]
      if (length(a) < 10 || length(b) < 10) next
      sp <- sqrt((var(a) * (length(a) - 1) + var(b) * (length(b) - 1)) /
                 (length(a) + length(b) - 2))
      if (!is.finite(sp) || sp == 0) next
      dd <- (mean(a) - mean(b)) / sp
      vv <- 1 / length(a) + 1 / length(b)
      est <- est + dd / vv; wsum <- wsum + 1 / vv; nn1 <- nn1 + length(a); nn0 <- nn0 + length(b)
    }
    if (wsum == 0) return(NULL)
    d_hat <- est / wsum; se <- sqrt(1 / wsum)
    p_tost <- max(pnorm(-(d_hat + SESOI_D) / se), pnorm((d_hat - SESOI_D) / se))
    list(d = d_hat, se = se, p = 2 * pnorm(-abs(d_hat / se)), p_tost = p_tost,
         n1 = nn1, n0 = nn0)
  }

  # v3(2026-08-25): 初回実行で判明——**2019補充の出生年[1987,1998]は
  #  既存コホート[1966,1986]と共通サポートが空**(2019補充は年齢マッチ型でなく
  #  「新しい若年世代の投入」型)。よって
  #   (A) 主コントラスト = c2007_t13 vs c2011_t9(旧コホート同士: サポート完全重複)
  #       → NC版は Δβ^N(13,9) = テニュア13と9の選抜差 を点識別(TP2定理3の系)
  #   (B) 旧コホート vs 2019新規 = 世代合成(Δg+β)——βとして解釈不可。
  #       記述用に別ファイルへ出力(ラベル明示)。
  res <- list(); gen <- list()
  for (i in seq_len(nrow(nc))) {
    yZ <- clean(nc$w1_Z[i]); yD <- clean(nc$w5_D[i]); yL <- clean(nc$w13_L[i])
    # (A) Δβ^N(13,9): 2007コホート(Z測定) vs 2011コホート(D測定)、w13残存者同士
    if (!is.null(yZ) && !is.null(yD)) {
      y <- fifelse(grp == "c2007_t13", yZ, fifelse(grp == "c2011_t9", yD, NA_real_))
      sup <- support_of("c2007_t13", "c2011_t9")
      r <- sdiff(y, grp == "c2007_t13" & !is.na(y) & sup,
                 grp == "c2011_t9" & !is.na(y) & sup)
      if (!is.null(r)) res[[length(res) + 1]] <- data.table(
        stem = nc$stem[i], family = nc$family[i],
        contrast = "dbeta_t13_minus_t9_incumbents",
        n_1 = r$n1, n_0 = r$n0, d = round(r$d, 4), se = round(r$se, 4),
        p = r$p, p_tost = r$p_tost)
    }
    # (B) 世代合成(参考のみ): 旧コホート vs 2019新規(サポート外——Δg支配)
    if (!is.null(yL)) {
      for (gg in c("c2007_t13", "c2011_t9")) {
        yo <- if (gg == "c2007_t13") yZ else yD
        if (is.null(yo)) next
        y <- fifelse(grp == gg, yo, fifelse(grp == "c2019_t1", yL, NA_real_))
        r <- sdiff(y, grp == gg & !is.na(y), grp == "c2019_t1" & !is.na(y))
        if (!is.null(r)) gen[[length(gen) + 1]] <- data.table(
          stem = nc$stem[i], family = nc$family[i],
          contrast = paste0("GENERATIONAL_COMPOSITE_", gg, "_vs_c2019"),
          n_1 = r$n1, n_0 = r$n0, d = round(r$d, 4), se = round(r$se, 4))
      }
    }
  }
  if (length(res)) {
    out <- rbindlist(res)
    out[, q := p.adjust(p, "BH"), by = contrast]
    out[, equiv_tost := p_tost < 0.05]
    write_aggregate(out, "11_negcontrol_w13.csv")
    summ <- out[, .(n_items = .N, median_abs_d = round(median(abs(d)), 4),
                    max_abs_d = round(max(abs(d)), 4), n_reject_q10 = sum(q < 0.10),
                    n_equiv = sum(equiv_tost),
                    pooled_d_invvar = round(sum(d / se^2) / sum(1 / se^2), 4),
                    pooled_se = round(sqrt(1 / sum(1 / se^2)), 4)), by = contrast]
    write_aggregate(summ, "11_negcontrol_w13_summary.csv",
                    exempt = c("n_items", "n_reject_q10", "n_equiv"))
    cat("\ndiag: Δβ^N(13,9)の要約(旧コホート同士・共通サポート内; w5版の参照値: 中央値|d|=.031):\n")
    print(summ)
    cat("  解釈: ≈0なら「テニュア9→13で選抜組成は追加的に歪まない」——増分τ(13)−τ(9)の\n")
    cat("  負対照補正が小さいことを意味する(TP2定理3の系の実装)。\n")
  } else cat("\nWARN: 主コントラスト(A)が全滅——console上の群サイズ・サポートを確認。\n")
  if (length(gen)) {
    write_aggregate(rbindlist(gen), "11_generational_composites.csv")
    cat("\n注意: 11_generational_composites.csv は**世代合成(Δg+β)**であり選抜βでは\n")
    cat("ない(2019補充と旧コホートは出生年サポートが非重複)。エアコン・パソコン等が\n")
    cat("大きく出るのは世代効果そのもの。βとして引用しないこと。\n")
  }
  writeLines(c(paste("time:", format(Sys.time())), paste("R", getRversion()),
               paste("input:", fp$file), paste("input size:", fp$bytes),
               paste("input md5:", fp$md5), paste("input sha256:", fp$sha256)),
             file.path(RESULTS_DIR, "11_env.txt"))
  cat("\n== 11 v3 完了。共有してほしいもの ==\n")
  cat("コンソール出力全文 + 11_negcontrol_w13.csv / 11_negcontrol_w13_summary.csv / 11_env.txt\n")
}

run_guarded("11_negcontrol_w13", main)

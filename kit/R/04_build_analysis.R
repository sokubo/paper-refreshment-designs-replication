# ============================================================
# P1 JLPS kit — 04: w5コントラストの分析ファイル構築
# 内容:
#   (1) リスク集合: w5回答者(respM=DQ74M または婚姻DQ43の非欠測)、cn∈{1,2}
#   (2) 項目ユニバース: D系(w5)のうち両群で被覆のある実質項目
#       (回顧・パラデータ・文字列・識別子を除外)
#   (3) アウトカム構築: A族=実質値(標準化)、B族=項目無回答・DK、
#       グリッドstraightlining、人レベルDQ指標(欠測率・DK率・中間・極端・SL)
#   (4) バランシング: 参入時共変量(性・出生年・負対照23項目)で
#       IPW(継続→追加の分布に合わせるATC型)+バランス診断
#   (5) 派生データはローカル(DERIVED_DIR/analysis_w5.rds)に保存。
#       resultsへは集計診断のみ(ユニバース要約・バランス・重み診断)
# v4 2026-09-15: NAP_OVERRIDE(尺度付随の「いない」コード)を追加。
# 実行: cd <P1ルート> && Rscript analysis/R/04_build_analysis.R
# ============================================================

.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
suppressMessages({library(haven); library(data.table)})

RESP_W5 <- "DQ74M"; MAR_W5 <- "DQ43"
# 負対照23項目(03_negcontrol_sets_v1.csvで確定): 継続=Z版(w1)、追加=D版(w5)
NC_PAIRS <- list(
  c("ZQ16","DQ60"), c("ZQ17","DQ61"), c("ZQ19","DQ63"), c("ZQ21","DQ65"),
  c("ZQ18_A","DQ62_A"), c("ZQ18_B","DQ62_B"), c("ZQ18_C","DQ62_C"),
  c("ZQ18_E","DQ62_E"), c("ZQ18_F","DQ62_F"), c("ZQ18_G","DQ62_G"),
  c("ZQ18_H","DQ62_H"), c("ZQ18_I","DQ62_I"), c("ZQ18_J","DQ62_J"),
  c("ZQ18_K","DQ62_K"), c("ZQ18_L","DQ62_L"), c("ZQ18_M","DQ62_M"),
  c("ZQ18_N","DQ62_N"), c("ZQ18_O","DQ62_O"), c("ZQ18_P","DQ62_P"),
  c("ZQ18_Q","DQ62_Q"), c("ZQ18_R","DQ62_R"), c("ZQ18_S","DQ62_S"),
  c("ZQ18_T","DQ62_T"))
# 除外: 回顧ブロック(追加標本のみ回答)・回答月日・自由記述系
EXCL_PAT <- "^DQ(59|6[0-9]|7[0-2])|^DQ7[34]M$|^DQ7[34]D$|_u$|x$"
# 特殊コードの分類(値ラベル文字列で判定)——v2: DKだけでなく無回答・拒否・非該当を区別
DK_PAT   <- "わからない|分からない|ＤＫ|DK"
REF_PAT  <- "答えたくない|回答したくない"
NR_PAT   <- "無回答|不明"
NAP_PAT  <- "非該当"
# フィルター不一致の閾値: 回答ベース率の群間差がこれを超える項目は
# 設問ルーティング(調査票の版差・スキップ)の混入とみなし別分類
FILTER_GAP <- 0.15
# v4(2026-09-15): 尺度に付随する「該当者なし」コードを非該当として扱う項目別上書き
#   (満足度 DQ18A-E の 6=仕事/結婚/友人/親/子がいない、DQ04_1C の 5=部下はいない、DQ04_2 の 5=上司・同僚はいない)。
#   ラベル文字列では拾えない(「非該当」と書かれていない)ため明示する。kit 15 の EC 対応(w1 ZQ30A-E は 6=非該当)とも整合。
if (!exists("NAP_OVERRIDE")) NAP_OVERRIDE <- list(DQ18A = 6, DQ18B = 6, DQ18C = 6, DQ18D = 6, DQ18E = 6, DQ04_1C = 5, DQ04_2 = 5)   # 正本は 00_config.R

main <- function() {
  # 2020ウェブ特別調査ファイル(JLPSYM_online_*)が同居してもマスタを確実に選ぶ
  f <- input_dta(); fp <- input_fingerprint(f)
  cat(sprintf("input: %s (%s bytes; md5 %s; sha256 %s)\n", fp$file, fp$bytes, fp$md5, fp$sha256))
  d <- as.data.table(read_dta(f))
  getcol <- function(nm, required = TRUE) {
    hit <- names(d)[toupper(trimws(names(d))) == toupper(nm)]
    if (!length(hit)) { if (required) stop("column not found: ", nm, call. = FALSE); return(NULL) }
    d[[hit[1]]]
  }
  num <- function(v) { x <- getcol(v, required = FALSE); if (is.null(x)) return(NULL)
                       as.numeric(zap_labels(x)) }

  ## --- (1) リスク集合 ---------------------------------------------------------
  cn <- as.integer(zap_labels(getcol("CN")))
  r5 <- !is.na(num(RESP_W5)) | !is.na(num(MAR_W5))
  keep <- which(r5 & cn %in% c(1L, 2L))
  Tt <- as.integer(cn[keep] == 1L)               # 1=継続(経験) / 0=追加(新規)
  cat("diag: risk set n =", length(keep), " (T=1:", sum(Tt), ", T=0:", sum(1 - Tt), ")\n")

  ## --- (2) 項目ユニバース -----------------------------------------------------
  dvars <- grep("^[Dd][Qq][0-9]", names(d), value = TRUE)
  dvars <- dvars[!grepl(EXCL_PAT, dvars, ignore.case = TRUE)]
  is_chr <- vapply(d[, ..dvars], is.character, TRUE)
  dvars <- dvars[!is_chr]
  codes_of <- function(vl, pat) if (!is.null(vl)) unname(vl[grepl(pat, names(vl))]) else numeric(0)
  nap_codes_for <- function(v, vl) { x <- codes_of(vl, NAP_PAT); ov <- NAP_OVERRIDE[[toupper(v)]]; if (!is.null(ov)) x <- union(x, ov); x }
  meta <- rbindlist(lapply(dvars, function(v) {
    x <- d[[v]][keep]
    lab <- attr(d[[v]], "label", exact = TRUE); if (is.null(lab)) lab <- ""
    vl  <- attr(d[[v]], "labels", exact = TRUE)
    xv  <- as.numeric(zap_labels(x))
    dk_codes  <- codes_of(vl, DK_PAT);  ref_codes <- codes_of(vl, REF_PAT)
    nr_codes  <- codes_of(vl, NR_PAT);  nap_codes <- nap_codes_for(v, vl)
    spec <- c(dk_codes, ref_codes, nr_codes, nap_codes)
    subst <- xv; subst[subst %in% spec] <- NA
    n1 <- sum(!is.na(subst[Tt == 1])); n0 <- sum(!is.na(subst[Tt == 0]))
    # 回答ベース率(実質+特殊コード含む「その設問に到達した」割合)——フィルター診断
    base1 <- mean(!is.na(xv[Tt == 1]) & !(xv[Tt == 1] %in% nap_codes))
    base0 <- mean(!is.na(xv[Tt == 0]) & !(xv[Tt == 0] %in% nap_codes))
    k  <- if (!is.null(vl)) length(setdiff(unname(vl), spec)) else 0L
    un <- length(unique(na.omit(subst)))
    data.table(var = v, label = substr(as.character(lab), 1, 60),
               n_T1 = n1, n_T0 = n0, k_labels = k, n_unique = un,
               n_dk_codes = length(dk_codes), n_ref_codes = length(ref_codes),
               base_rate_T1 = round(base1, 3), base_rate_T0 = round(base0, 3))
  }))
  meta[, scale_type := fifelse(n_unique <= 1, "degenerate",
                        fifelse(k_labels == 2 | n_unique == 2, "binary",
                        fifelse(k_labels %in% 3:9, "ordinal", "continuous")))]
  meta[, filter_mismatch := abs(base_rate_T1 - base_rate_T0) > FILTER_GAP]
  meta[, in_universe := n_T1 >= 200 & n_T0 >= 50 & scale_type != "degenerate"]
  cat("diag: D系候補", nrow(meta), "項目 → ユニバース", sum(meta$in_universe),
      "項目(うちフィルター不一致フラグ", sum(meta$in_universe & meta$filter_mismatch), ")\n")
  print(meta[in_universe == TRUE, .N, by = scale_type])

  ## --- (3) アウトカム行列 -----------------------------------------------------
  uni <- meta[in_universe == TRUE, var]
  Ymat  <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  Mmat  <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  DKmat <- REFmat <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  mid_cnt <- ext_cnt <- mid_n <- ext_n <- rep(0L, length(keep))
  for (v in uni) {
    x  <- d[[v]][keep]; vl <- attr(d[[v]], "labels", exact = TRUE)
    xv <- as.numeric(zap_labels(x))
    dk_codes  <- codes_of(vl, DK_PAT);  ref_codes <- codes_of(vl, REF_PAT)
    nr_codes  <- codes_of(vl, NR_PAT);  nap_codes <- nap_codes_for(v, vl)
    spec  <- c(dk_codes, ref_codes, nr_codes, nap_codes)
    base  <- !is.na(xv) & !(xv %in% nap_codes)         # 設問に到達(非該当を除く)
    subst <- xv; subst[xv %in% spec] <- NA
    # B族は「設問到達者」の中で定義(フィルター構造と分離)。
    # 注: NAのままの真の項目無回答は base=FALSE 側に落ちる(NAと非該当を
    # データ上区別できないため)——無回答コード型のみをB_itemnonrespに使う。
    Mmat[base, v]   <- as.integer(xv[base] %in% nr_codes)
    DKmat[base, v]  <- as.integer(xv[base] %in% dk_codes)
    REFmat[base, v] <- as.integer(xv[base] %in% ref_codes)
    st <- meta[var == v, scale_type]
    if (st == "continuous") {                       # 1/99%ウィンザライズ
      qs <- quantile(subst, c(.01, .99), na.rm = TRUE)
      subst <- pmin(pmax(subst, qs[1]), qs[2])
    }
    Ymat[, v] <- subst
    # 中間・極端(実質値ラベルが4-7件の順序尺度; 値集合ベースで判定)
    k <- meta[var == v, k_labels]
    vals <- if (!is.null(vl)) sort(setdiff(unname(vl), spec)) else integer(0)
    if (k %in% 4:7 && length(vals) == k && all(diff(vals) == 1)) {
      ok <- !is.na(subst)
      if (k %% 2 == 1) { mid_cnt[ok] <- mid_cnt[ok] + as.integer(subst[ok] == vals[(k + 1) / 2])
                         mid_n[ok] <- mid_n[ok] + 1L }
      ext_cnt[ok] <- ext_cnt[ok] + as.integer(subst[ok] %in% c(vals[1], vals[k]))
      ext_n[ok] <- ext_n[ok] + 1L
    }
  }
  ## グリッド(同一問番号で4枝以上・同一尺度)のstraightlining
  meta[, qbase := sub("_[A-Za-z0-9]+$", "", var)]
  grids <- meta[in_universe == TRUE & grepl("_", var),
                .(nsub = .N, kk = uniqueN(k_labels)), by = qbase][nsub >= 4 & kk == 1, qbase]
  SL <- rep(0L, length(keep)); SLn <- rep(0L, length(keep))
  for (g in grids) {
    vs <- meta[in_universe == TRUE & qbase == g, var]
    sub <- Ymat[, vs, drop = FALSE]
    full <- rowSums(is.na(sub)) == 0
    sd0  <- apply(sub, 1, function(r) if (all(is.na(r))) NA else sd(r, na.rm = TRUE))
    SL[full]  <- SL[full] + as.integer(sd0[full] == 0)
    SLn[full] <- SLn[full] + 1L
  }
  cat("diag: グリッド数(SL計算対象) =", length(grids), "\n")
  ## v3: 無回答コード率の群間差が大きい項目(=更新設計・ルーティング差の疑い。
  ## 例: 結婚経緯ブロックは継続票では新婚のみ記入で他は「無回答」コード)を
  ## 人レベルmiss_shareから除外し、メタにフラグ(05でB族からも隔離)
  nr_gap <- vapply(uni, function(v) {
    m <- Mmat[, v]
    abs(mean(m[Tt == 1], na.rm = TRUE) - mean(m[Tt == 0], na.rm = TRUE))
  }, 0.0)
  meta[var %in% uni, nr_routing_flag := nr_gap[var] > 0.15]
  clean_m <- uni[!(nr_gap > 0.15)]
  cat("diag: 無回答ルーティング疑い項目 =", sum(nr_gap > 0.15),
      "(miss_share・B族主分類から隔離)\n")
  person_dq <- data.table(
    miss_share = rowMeans(Mmat[, clean_m, drop = FALSE], na.rm = TRUE),
    dk_share   = rowMeans(DKmat, na.rm = TRUE),
    ref_share  = rowMeans(REFmat, na.rm = TRUE),
    mid_share  = fifelse(mid_n > 0, mid_cnt / mid_n, NA_real_),
    ext_share  = fifelse(ext_n > 0, ext_cnt / ext_n, NA_real_),
    sl_share   = fifelse(SLn > 0, SL / SLn, NA_real_))
  cat("diag: 中間/極端の対象項目をもつ人 =", sum(mid_n > 0), "/", sum(ext_n > 0), "\n")

  ## --- (4) 共変量とIPW --------------------------------------------------------
  sex <- num("sex")[keep]; yb <- num("ybirth")[keep]
  X <- data.table(sex = sex, yb = yb)
  for (p in NC_PAIRS) {
    z <- num(p[1]); dd <- num(p[2])
    v <- rep(NA_real_, length(keep))
    if (!is.null(z))  v[Tt == 1] <- z[keep][Tt == 1]
    if (!is.null(dd)) v[Tt == 0] <- dd[keep][Tt == 0]
    X[[paste0("nc_", p[2])]] <- v
  }
  # 欠測はダミー法(中央値代入+欠測フラグ)
  for (cname in names(X)) {
    v <- X[[cname]]; mi <- is.na(v)
    if (any(mi)) { X[[paste0(cname, "_mi")]] <- as.integer(mi)
                   v[mi] <- median(v, na.rm = TRUE); X[[cname]] <- v }
  }
  X[, `:=`(yb2 = (yb - mean(yb))^2, yb3 = (yb - mean(yb))^3)]
  ps_fit <- glm(Tt ~ ., data = cbind(Tt = Tt, X), family = binomial())
  e <- pmin(pmax(fitted(ps_fit), 0.01), 0.99)
  w <- ifelse(Tt == 1, (1 - e) / e, 1)                   # ATC型: 継続→追加の分布へ
  w[Tt == 1] <- w[Tt == 1] / mean(w[Tt == 1])            # 安定化
  cap <- quantile(w[Tt == 1], .99); w[Tt == 1] <- pmin(w[Tt == 1], cap)
  ess <- sum(w[Tt == 1])^2 / sum(w[Tt == 1]^2)
  cat("diag: IPW ESS(T=1) =", round(ess), "/", sum(Tt), "; max w =", round(max(w), 2), "\n")

  ## バランス診断(標準化平均差、重み付け前後)
  bal <- rbindlist(lapply(setdiff(names(X), c("yb2", "yb3")), function(cname) {
    x <- X[[cname]]
    smd <- function(wt) { m1 <- weighted.mean(x[Tt == 1], wt[Tt == 1]); m0 <- mean(x[Tt == 0])
      s <- sqrt((var(x[Tt == 1]) + var(x[Tt == 0])) / 2); if (s == 0) return(NA)
      round((m1 - m0) / s, 4) }
    data.table(covariate = cname, smd_unw = smd(rep(1, length(x))), smd_ipw = smd(w))
  }))
  write_aggregate(bal, "04_balance.csv")
  write_aggregate(meta[, .(var, label, n_T1, n_T0, k_labels, scale_type,
                           base_rate_T1, base_rate_T0, filter_mismatch,
                           nr_routing_flag, in_universe)],
                  "04_item_meta.csv", exempt = c("k_labels"))
  write_aggregate(data.table(stat = c("n_risk", "n_T1", "n_T0", "n_universe", "n_grids",
                                      "ess_T1", "max_w"),
                             value = c(length(keep), sum(Tt), sum(1 - Tt),
                                       sum(meta$in_universe), length(grids),
                                       round(ess), round(max(w), 2))),
                  "04_summary.csv", exempt = "value")

  ## --- (5) 派生保存(ローカルのみ) --------------------------------------------
  idc <- names(d)[toupper(names(d)) %in% toupper(ID_COLS)]
  src <- list(input = fp, rows = as.integer(keep), cn = as.integer(cn[keep]),
              id = if (length(idc)) as.character(zap_labels(d[[idc[1]]]))[keep] else NULL)
  saveRDS(list(T = Tt, w = w, e = e, X = as.matrix(X), Y = Ymat, M = Mmat,
               DK = DKmat, REF = REFmat, person_dq = person_dq,
               meta = meta[in_universe == TRUE], src = src),
          file.path(DERIVED_DIR, "analysis_w5.rds"))
  cat("saved:", file.path(DERIVED_DIR, "analysis_w5.rds"), "\n")
  cat("\n== 04 完了。共有してほしいもの ==\n")
  cat("コンソール出力全文 + 04_summary.csv / 04_balance.csv / 04_item_meta.csv\n")
}

run_guarded("04_build_analysis", main)

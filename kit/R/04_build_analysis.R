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
# v5 2026-09-24: 中間・極端の人レベル指標は、R/15_style_items.csv に事前指定した評定尺度項目だけから作る
#   (R/15_style_items.R。以前は「実質コードが 4〜7 個で連番」の全項目が入り、婚姻状態・役職などの名義・分類項目も含まれていた)。
# v6 2026-09-28 (v1.0): 項目ごとの尺度・特殊コード・再符号化・指標化・派生・フィルターは R/15_item_scale.csv で明示する
#   (R/15_item_scale.R)。値ラベルの数だけで尺度を決める規則と、ラベル文字列だけで特殊コードを決める規則はやめた
#   (review_round06/SCALE_AUDIT_2011_ja.md)。表にない項目はユニバースに入れず 04_items_unlisted.csv に列挙する。
#   名義の単一選択設問は選択肢ごとの 0/1 指標(<var>__<code>)にし、時刻と交際期間は構成要素から派生項目を作る。
# v6.1 2026-09-29: 再符号化でまとめたコードの指標ラベルは、そのコード自身のラベル(あれば)にする(DQ30 の 7=8 で
#   「その他の政党」が「みんなの党」と出ていた)。表の変更(DQ03_1・DQ30)は R/15_item_scale.csv の note を見ること。
# 実行: cd <P1ルート> && Rscript analysis/R/04_build_analysis.R
# ============================================================

.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
suppressMessages({library(haven); library(data.table)})
source(file.path(.here, "15_style_items.R"))     # response-style composites: fixed item list + code verification
source(file.path(.here, "15_item_scale.R"))      # the item specification table (scales, special codes, recodes, indicators, derived items)

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

  ## --- (2) 項目ユニバース(v6: R/15_item_scale.csv による明示的な指定) ------------------------
  spec_tb <- read_item_scale(); sid_scale <- item_scale_id()
  cat("diag: item specification table:", basename(sid_scale$path), "md5", sid_scale$md5, ";", nrow(spec_tb), "rows\n")
  pats <- list(dk = DK_PAT, ref = REF_PAT, nr = NR_PAT, nap = NAP_PAT)
  dvars <- grep("^[Dd][Qq][0-9]", names(d), value = TRUE)
  dvars <- dvars[!grepl(EXCL_PAT, dvars, ignore.case = TRUE)]
  is_chr <- vapply(d[, ..dvars], is.character, TRUE)
  dvars <- dvars[!is_chr]
  spec_row <- function(v) { i <- match(toupper(v), spec_tb$key); if (is.na(i)) NULL else spec_tb[i] }
  unlisted <- dvars[vapply(dvars, function(v) is.null(spec_row(v)), TRUE)]
  if (length(unlisted)) cat("diag:", length(unlisted), "D-series variables are not listed in 15_item_scale.csv and are left out (04_items_unlisted.csv)\n")
  write_aggregate(data.table(var = unlisted, label = vapply(unlisted, function(v) { l <- attr(d[[v]], "label", exact = TRUE); if (is.null(l)) "" else substr(as.character(l), 1, 60) }, "")), "04_items_unlisted.csv")
  listed <- setdiff(dvars, unlisted)
  ## a candidate item: the numeric vector (special codes to missing, recoded, filtered), the reach indicator, its
  ## special-code kinds and substantive codes. Nominal items return one candidate per indicator.
  num_all <- function(nm) { x <- getcol(nm, required = FALSE); if (is.null(x)) NULL else as.numeric(zap_labels(x))[keep] }
  cand <- list(); vl_of <- list(); spec_list <- list(); parent_codes <- list()
  add_cand <- function(v, lab, subst, base, sc, scale, k, codes, labs, qgroup = "", parent = v, spec) {
    cand[[v]] <<- list(var = v, label = lab, subst = subst, base = base, sc = sc, scale = scale, k = k, codes = codes, labs = labs,
                       qgroup = qgroup, parent = parent, spec = spec)
  }
  for (v in listed) {
    sp <- spec_row(v)
    if (sp$scale == "exclude") { spec_list[[v]] <- sp; next }
    x <- d[[v]][keep]; vl <- attr(d[[v]], "labels", exact = TRUE)
    lab <- attr(d[[v]], "label", exact = TRUE); if (is.null(lab)) lab <- ""; lab <- substr(as.character(lab), 1, 60)
    xv <- as.numeric(zap_labels(x))
    sc <- special_codes_of(vl, sp, pats)
    ov <- NAP_OVERRIDE[[toupper(v)]]; if (!is.null(ov)) sc$nap <- union(sc$nap, ov)    # 00_config.R (kept for the entry-wave side of 15)
    spec <- unlist(sc)
    base <- !is.na(xv) & !(xv %in% sc$nap)
    subst <- xv; subst[subst %in% spec] <- NA
    subst <- apply_recode(subst, sp$recode)
    codes <- if (!is.null(vl)) recode_codes(setdiff(unname(vl), spec), sp$recode) else numeric(0)
    ## label of each (recoded) code: the label of the code itself when it survives the recode, else the first label merged into it
    labs  <- if (!is.null(vl)) vapply(codes, function(cd) { rc <- apply_recode(unname(vl), sp$recode); i <- which(!is.na(rc) & rc == cd & unname(vl) == cd)
                                                            if (!length(i)) i <- which(!is.na(rc) & rc == cd); if (!length(i)) "" else names(vl)[i[1]] }, "") else character(0)
    exp_codes <- parse_codes(sp$codes_expected)
    if (length(exp_codes) && !(length(codes) == length(exp_codes) && all(codes == exp_codes)))
      stop(sprintf("15_item_scale.csv: %s expects the substantive codes %s but the file shows %s", v, sp$codes_expected, paste(codes, collapse = ";")), call. = FALSE)
    if (sp$scale %in% c("binary", "ordinal", "nominal") && !length(codes))
      stop(sprintf("15_item_scale.csv: %s is declared %s but the file has no substantive value labels for it", v, sp$scale), call. = FALSE)
    if (nzchar(sp$filter)) { mk <- filter_mask(sp$filter, length(keep), num_all); subst[!mk] <- NA; base <- base & mk }
    vl_of[[v]] <- vl; spec_list[[v]] <- sp
    if (sp$scale == "nominal") {
      if (!length(codes)) stop("15_item_scale.csv: nominal item without value labels: ", v, call. = FALSE)
      ind <- make_indicators(subst, codes, v); parent_codes[[v]] <- codes
      for (j in seq_along(codes)) {
        nm <- colnames(ind)[j]
        ## the item-nonresponse and don't-know indicators of the question are attached to its first indicator only
        add_cand(nm, paste0(lab, ": ", labs[j]), ind[, j], if (j == 1L) base else rep(FALSE, length(base)), sc, "binary", 2L, c(0, 1), c("0", "1"), qgroup = v, parent = v, spec = sp)
      }
    } else {
      k <- if (sp$scale %in% c("binary", "ordinal")) length(codes) else 0L
      add_cand(v, lab, subst, base, sc, sp$scale, k, codes, labs, spec = sp)
    }
  }
  ## derived items (constructed from components read from the file)
  derived_tb <- spec_tb[grepl("^(clock|months):", construct)]
  for (i in seq_len(nrow(derived_tb))) {
    sp <- derived_tb[i]; cs <- parse_construct(sp$construct)
    comp <- if (cs$type == "clock") c(cs$X, cs$Y, cs$Z) else c(cs$Y, cs$M)
    vals <- lapply(comp, num_all)
    if (any(vapply(vals, is.null, TRUE))) { cat("diag: derived item", sp$var, "skipped: component(s) not in file:", paste(comp[vapply(vals, is.null, TRUE)], collapse = ", "), "\n"); next }
    if (cs$type == "clock") {
      leave <- if (!is.null(cs$leave)) { if (is.null(cand[[cs$leave]])) stop("clock construct of ", sp$var, " refers to ", cs$leave, " which is not built yet", call. = FALSE); cand[[cs$leave]]$subst } else NULL
      subst <- build_clock(vals[[1]], vals[[2]], vals[[3]], cs$kind, leave)
    } else subst <- build_months(vals[[1]], vals[[2]])
    if (nzchar(sp$filter)) { mk <- filter_mask(sp$filter, length(keep), num_all); subst[!mk] <- NA }
    base <- !is.na(subst)                              # reach = a usable time/duration (the format code is a separate item)
    add_cand(sp$var, sp$label, subst, base, list(dk = numeric(0), ref = numeric(0), nr = numeric(0), nap = numeric(0)), "continuous", 0L, numeric(0), character(0), spec = sp)
    spec_list[[sp$var]] <- sp
  }
  meta <- rbindlist(lapply(cand, function(cc) {
    n1 <- sum(!is.na(cc$subst[Tt == 1])); n0 <- sum(!is.na(cc$subst[Tt == 0]))
    base1 <- mean(cc$base[Tt == 1]); base0 <- mean(cc$base[Tt == 0])
    un <- length(unique(na.omit(cc$subst)))
    data.table(var = cc$var, label = cc$label, n_T1 = n1, n_T0 = n0, k_labels = as.integer(cc$k), n_unique = un,
               n_dk_codes = length(cc$sc$dk), n_ref_codes = length(cc$sc$ref),
               base_rate_T1 = round(base1, 3), base_rate_T0 = round(base0, 3),
               spec_scale = cc$spec$scale, qgroup = cc$qgroup, parent = cc$parent, in_counts = cc$spec$in_counts,
               ec_class = cc$spec$ec_class, ec_subgroup = cc$spec$ec_subgroup, ec_main = cc$spec$ec_main,
               recode = cc$spec$recode, construct = cc$spec$construct, filter = cc$spec$filter,
               scale_type = if (un <= 1) "degenerate" else cc$scale)
  }))
  meta[, filter_mismatch := abs(base_rate_T1 - base_rate_T0) > FILTER_GAP]
  meta[, in_universe := n_T1 >= 200 & n_T0 >= 50 & scale_type != "degenerate"]
  ## items of the table that are excluded or not built (for the record, outside the universe)
  n_excl <- sum(spec_tb$scale == "exclude" & spec_tb$key %in% toupper(listed))
  cat("diag: D系候補", length(listed), "項目(表にあるもの; 除外指定", n_excl, ") → 候補", nrow(meta), "列(名義の指標と派生項目を含む) → ユニバース", sum(meta$in_universe),
      "列(うちフィルター不一致フラグ", sum(meta$in_universe & meta$filter_mismatch), ")\n")
  print(meta[in_universe == TRUE, .N, by = scale_type])
  print(meta[in_universe == TRUE, .N, by = .(spec_scale, derived = nzchar(construct) & construct != "indicators")])

  ## --- (3) アウトカム行列 -----------------------------------------------------
  uni <- meta[in_universe == TRUE, var]
  Ymat  <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  Mmat  <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  DKmat <- REFmat <- matrix(NA_real_, length(keep), length(uni), dimnames = list(NULL, uni))
  spec_of_item <- list(); vals_of <- list(); labs_of <- list()   # special codes, substantive codes and their labels per universe item (used below and by 15)
  for (v in uni) {
    cc <- cand[[v]]; base <- cc$base; subst <- cc$subst
    ## B family: defined among those who reached the question (routing kept apart). A true item nonresponse recorded
    ## as NA cannot be told from routing; only the no-answer code enters B_itemnonresp.
    src_raw <- if (nzchar(cc$qgroup) || nzchar(cc$spec$construct)) NULL else as.numeric(zap_labels(d[[v]][keep]))
    if (!is.null(src_raw)) {
      Mmat[base, v]   <- as.integer(src_raw[base] %in% cc$sc$nr)
      DKmat[base, v]  <- as.integer(src_raw[base] %in% cc$sc$dk)
      REFmat[base, v] <- as.integer(src_raw[base] %in% cc$sc$ref)
    } else if (nzchar(cc$qgroup) && any(base)) {          # first indicator of a nominal question: the question's codes
      src_raw <- as.numeric(zap_labels(d[[cc$parent]][keep]))
      Mmat[base, v]   <- as.integer(src_raw[base] %in% cc$sc$nr)
      DKmat[base, v]  <- as.integer(src_raw[base] %in% cc$sc$dk)
      REFmat[base, v] <- as.integer(src_raw[base] %in% cc$sc$ref)
    }
    if (cc$scale == "continuous") {                 # 1/99%ウィンザライズ
      qs <- quantile(subst, c(.01, .99), na.rm = TRUE)
      subst <- pmin(pmax(subst, qs[1]), qs[2])
    }
    Ymat[, v] <- subst
    spec_of_item[[v]] <- unlist(cc$sc); vals_of[[v]] <- cc$codes; labs_of[[v]] <- cc$labs
  }
  ## value labels of the universe items (for the style-item verification below): the parent's labels for an indicator
  vl_uni <- lapply(uni, function(v) vl_of[[cand[[v]]$parent]]); names(vl_uni) <- uni
  vl_of <- vl_uni; spec_of <- spec_of_item
  ## 中間・極端の人レベル指標(v5): 事前指定リスト(評定尺度)の項目のうち、ユニバースにあり、実質コードが 1..k で
  ## 検証できたものだけ。中点指標は、中央のコードが中立ラベル(NEUTRAL_PAT)であることを確認できた奇数件法の項目だけ。
  st <- read_style_items()
  style <- rbindlist(lapply(seq_len(nrow(st)), function(i) {
    v <- names(vl_of)[toupper(names(vl_of)) == st$var[i]]
    if (!length(v)) return(data.table(st[i], var_file = NA_character_, in_universe = FALSE, ok = FALSE, mid_ok = FALSE, codes = "", labels = "", reason = "not in the item universe"))
    r <- style_verify(vl_of[[v]], spec_of[[v]], st$k[i], st$midpoint[i])
    if (r$ok && !style_values_ok(Ymat[, v], st$k[i])) { r$ok <- FALSE; r$mid_ok <- FALSE; r$reason <- sprintf("values outside 1..%d observed", st$k[i]) }
    data.table(st[i], var_file = v, in_universe = TRUE, ok = r$ok, mid_ok = r$mid_ok, codes = r$codes, labels = r$labels, reason = r$reason)
  }))
  k_of <- setNames(as.list(style$k), style$var_file)
  ext_items <- style[set == "rating" & ok == TRUE, var_file]; mid_items <- style[set == "rating" & mid_ok == TRUE, var_file]
  sh <- style_shares(Ymat, ext_items, mid_items, k_of)
  cat("diag: 様式指標(評定尺度リスト): 極端", length(ext_items), "項目 / 中点", length(mid_items), "項目; リストのうちユニバース外",
      sum(!style$in_universe), "、コード不一致", sum(style$in_universe & !style$ok), "\n")
  ## グリッド(同一問番号で4枝以上・同一尺度)のstraightlining
  meta[, qbase := sub("_[A-Za-z0-9]+$", "", var)]
  grids <- meta[in_universe == TRUE & grepl("_", var) & !nzchar(qgroup) & !nzchar(construct),
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
  ## v3: 無回答コード率の群間差が 15 ポイントを超える項目(継続票と追加票で経路や無回答の扱いが
  ## 違う疑い。実際の調査票の分岐は照合していないので、これは「比較可能性に疑いあり」の
  ## スクリーニング規則であって、経路差の記録ではない)を人レベル miss_share から除外し、
  ## メタにフラグ(05 で B 族からも隔離。15 は件数からも外し、閾値 .10/.20/なしでの感度を出す)
  nr_gap <- vapply(uni, function(v) {
    m <- Mmat[, v]
    abs(mean(m[Tt == 1], na.rm = TRUE) - mean(m[Tt == 0], na.rm = TRUE))
  }, 0.0)
  nr_gap[!is.finite(nr_gap)] <- 0        # derived items and the later indicators of a nominal question carry no item-nonresponse column
  meta[var %in% uni, nr_routing_flag := nr_gap[var] > 0.15]
  clean_m <- uni[!(nr_gap > 0.15)]
  cat("diag: 無回答ルーティング疑い項目 =", sum(nr_gap > 0.15),
      "(miss_share・B族主分類から隔離)\n")
  person_dq <- data.table(
    miss_share = rowMeans(Mmat[, clean_m, drop = FALSE], na.rm = TRUE),
    dk_share   = rowMeans(DKmat, na.rm = TRUE),
    ref_share  = rowMeans(REFmat, na.rm = TRUE),
    mid_share  = sh$mid,
    ext_share  = sh$ext,
    sl_share   = fifelse(SLn > 0, SL / SLn, NA_real_))
  cat("diag: 中間/極端の対象項目をもつ人 =", sum(!is.na(sh$mid)), "/", sum(!is.na(sh$ext)), "\n")

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
  write_aggregate(meta[, .(var, label, n_T1, n_T0, k_labels, scale_type, spec_scale, qgroup, in_counts, ec_class, ec_subgroup, ec_main,
                           recode, construct, filter,
                           base_rate_T1, base_rate_T0, filter_mismatch,
                           nr_routing_flag, in_universe)],
                  "04_item_meta.csv", exempt = c("k_labels"))
  ## the resolved specification of every universe item (an aggregate: codes and labels, no counts), for the author's
  ## review against the questionnaire and for freezing into codes_expected of 15_item_scale.csv
  write_aggregate(data.table(var = uni, label = meta$label[match(uni, meta$var)], scale = meta$scale_type[match(uni, meta$var)],
                             codes = vapply(uni, function(v) paste(vals_of[[v]], collapse = ";"), ""),
                             labels = vapply(uni, function(v) paste(labs_of[[v]], collapse = " | "), ""),
                             dk = vapply(uni, function(v) paste(cand[[v]]$sc$dk, collapse = ";"), ""),
                             nr = vapply(uni, function(v) paste(cand[[v]]$sc$nr, collapse = ";"), ""),
                             nap = vapply(uni, function(v) paste(cand[[v]]$sc$nap, collapse = ";"), ""),
                             ref = vapply(uni, function(v) paste(cand[[v]]$sc$ref, collapse = ";"), "")),
                  "04_item_scale_resolved.csv")
  write_aggregate(data.table(stat = c("n_risk", "n_T1", "n_T0", "n_universe", "n_grids",
                                      "ess_T1", "max_w"),
                             value = c(length(keep), sum(Tt), sum(1 - Tt),
                                       sum(meta$in_universe), length(grids),
                                       round(ess), round(max(w), 2))),
                  "04_summary.csv", exempt = "value")

  ## --- (5) 派生保存(ローカルのみ) --------------------------------------------
  idc <- names(d)[toupper(names(d)) %in% toupper(ID_COLS)]
  src <- list(input = fp, style_items = style_items_id(), item_scale = sid_scale, rows = as.integer(keep), cn = as.integer(cn[keep]),
              id = if (length(idc)) as.character(zap_labels(d[[idc[1]]]))[keep] else NULL)
  saveRDS(list(T = Tt, w = w, e = e, X = as.matrix(X), Y = Ymat, M = Mmat,
               DK = DKmat, REF = REFmat, person_dq = person_dq,
               meta = meta[in_universe == TRUE], src = src,
               style = style, codes = vals_of, code_labels = labs_of,   # v5: verified style items; substantive codes and labels of every universe item
               spec = spec_tb, parent_codes = parent_codes),            # v6: the item specification table as read; codes of the nominal questions
          file.path(DERIVED_DIR, "analysis_w5.rds"))
  cat("saved:", file.path(DERIVED_DIR, "analysis_w5.rds"), "\n")
  cat("\n== 04 完了。共有してほしいもの ==\n")
  cat("コンソール出力全文 + 04_summary.csv / 04_balance.csv / 04_item_meta.csv / 04_item_scale_resolved.csv / 04_items_unlisted.csv\n")
}

run_guarded("04_build_analysis", main)

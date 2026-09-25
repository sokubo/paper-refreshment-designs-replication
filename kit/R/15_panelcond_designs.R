# ============================================================
# P1 JLPS kit — 15: 補充標本設計の全項目適用(panelcond: naive / SM / SSM / EC / EC-adj + 診断)
# 目的(構成 I の論文 E、ウォーターフォール D4 SSM / D5 EC の素材):
#   04 の派生(analysis_w5.rds: w5 リスク集合 × 全項目ユニバース)と統合マスタ(.dta)から、
#   継続コホート(CN=1、全 4,800 人)と 2011 追加コホート(CN=2)を panelcond の入力形に組み、
#   項目ごとに naive / SM / SSM / EC / EC-adj、入所波選抜差分 δ、診断 T1(SM−EC: 非定常性)・
#   T2(SSM−SM: 状態依存)、漏斗比 sup_y p·q(y)/f2(y)(TP2 §6 の項目別診断)を推定する。
#   A 族(実質値)、B 族(項目無回答・DK)、人レベル様式指標(欠測率・DK 率・中点・極端)に同じ枠を適用。
#   中点・極端の指標(2026-09-24 以降): R/15_style_items.csv の評定尺度項目のうち入所波(w1)にも同じコードで存在する
#   項目(共通バッテリー)で、w5・w1 とも同じ項目集合から計算する。項目構成の感度(全評定項目、頻度尺度を加えた集合、
#   v0.7 までの規則「実質コードが 4〜7 個で連番の全項目」)は 15_style_sensitivity.csv に出す(R/15_style_items.R)。
#   推論: (i) panelcond::pc_estimate の解析 SE、(ii) 項目横断の**同時人ブートストラップ**(B 回、
#   コホート内で人を再抽出、全推定量と T1/T2 差を反復内で再計算)。(ii) はベクトル化した自前エンジンで、
#   フルデータ上の点推定が pc_estimate と一致することを 15_validation.csv で検証する。
# 用量の定義: DOSE_DEF = "exact"(継続は w1–w5 全波回答 = 用量 4 ちょうど; J3 の規約)
#              または "any"(w5 回答者全員 = 04 のリスク集合; 用量 ≤ 4 の混合)。既定は両方。
# 入所波(w1)の対応変数: jlps_docs/varmap_w1-19.csv の w1/w5 列 + 15_entry_overrides.csv(手動)。
#   値ラベル集合(特殊コードを除く)が w1 と w5 で一致する項目だけ EC を計算(不一致は 15_ec_unavailable.csv に列挙)。
# 出力(すべて集計値、N<10 抑制): 15_designs_items.csv / 15_tests_items.csv / 15_detection_counts.csv /
#   15_diagnostics_summary.csv / 15_mass_diagnostic_summary.csv / 15_mass_diagnostic_flags.csv / 15_routing_sensitivity.csv /
#   15_style_sensitivity.csv / 15_style_items_audit.csv / 15_core32.csv /
#   15_arms.csv / 15_validation.csv / 15_ec_unavailable.csv /
#   15_env.txt(R・パッケージ・入力ファイルの版)
# 実行: cd <P1ルート> && Rscript analysis/R/15_panelcond_designs.R [--B 500] [--seed 20260915] [--dose exact,any]
#   必要: haven, data.table, panelcond(papers/J3_dose_response/panelcond から R CMD INSTALL、または
#         remotes::install_github("sokubo/panelcond"))。所要: 実データで約 13 分(B=500、2 用量; 2026-09-15 実測)。
# ============================================================

.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
source(file.path(.here, "15_mass_ratio.R"))          # mass-domination helper (tested by 15_test_mass_ratio.R)
suppressMessages({library(haven); library(data.table); library(panelcond)})
source(file.path(.here, "15_style_items.R"))         # response-style composites over the fixed item list

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1] }
B_BOOT   <- as.integer(getarg("--B", "500"))
SEED15   <- as.integer(getarg("--seed", "20260915"))
DOSE_DEFS <- strsplit(getarg("--dose", "exact,any"), ",")[[1]]
K_DOSE   <- 4L                       # w5 の継続コホートの用量(w1–w4 の 4 回)
T_WAVE   <- 5L
RESP_W5 <- "DQ74M"; MAR_W5 <- "DQ43"  # 04 と同じリスク集合の定義
VARMAP  <- path.expand(Sys.getenv("P1_VARMAP", file.path(PROJECT_DIR, "..", "..", "jlps_docs", "varmap_w1-19.csv")))
OVERRIDES <- file.path(.here, "15_entry_overrides.csv")
EXCLUDE   <- file.path(.here, "15_item_exclude.csv")     # 名義尺度など平均対比が無意味な項目(出力には残しフラグを付ける)
CORE32  <- path.expand(Sys.getenv("P1_CORE32", file.path(PROJECT_DIR, "..", "J3_dose_response", "analysis", "j3_items.csv")))
SESOI_D <- 0.10; SESOI_PP <- 0.02; Q_CUT <- 0.10
DK_PAT <- "わからない|分からない|ＤＫ|DK"; REF_PAT <- "答えたくない|回答したくない"; NR_PAT <- "無回答|不明"; NAP_PAT <- "非該当"

## J3 の strict 応答マーカー(各波で全員に聞く項目; 小文字。w1–w10 まで)
MARKERS <- list(
  w1 = c("zq03","zq13_1","zq30d","zq50"), w2 = c("aq02","aq11_1","aq32d","aq52"), w3 = c("bq02","bq13_1","bq23f","bq42"),
  w4 = c("cq02","cq17_1","cq20f","cq45"), w5 = c("dq02","dq10_1","dq18f","dq43"), w6 = c("eq02","eq12_1","eq18f","eq41"),
  w7 = c("fq02","fq10_1","fq21f","fq44"), w8 = c("gq02","gq13_1","gq19f","gq43"), w9 = c("hq02","hq13_1","hq22f","hq44"),
  w10 = c("iq02","iq13_1","iq21f","iq45"))
WAVE_PREFIX <- c("z","a","b","c","d","e","f","g","h","i")
NA_MARKER_CODES <- c(88, 99, 888, 999, 8888, 9999, 88888, 99999)

main <- function() {
  t0 <- Sys.time()
  ## ---- 0. 入力 -------------------------------------------------------------------
  a <- readRDS(file.path(DERIVED_DIR, "analysis_w5.rds"))
  Tt <- a$T; Y <- a$Y; M <- a$M; DK <- a$DK; meta <- as.data.table(a$meta); pdq <- as.data.table(a$person_dq)
  stopifnot(ncol(Y) == nrow(meta), length(Tt) == nrow(Y))
  f <- input_dta(); fp <- input_fingerprint(f)
  cat(sprintf("input: %s (%s bytes; md5 %s; sha256 %s)\n", fp$file, fp$bytes, fp$md5, fp$sha256))
  if (is.null(a$src)) stop("analysis_w5.rds has no input record (built by an older 04); rerun 04_build_analysis.R first", call. = FALSE)
  if (is.null(a$style) || is.null(a$codes)) stop("analysis_w5.rds has no style-item record (built by a 04 older than 2026-09-24); rerun 04_build_analysis.R first", call. = FALSE)
  style04 <- as.data.table(a$style); codes_of_var <- a$codes
  sid <- style_items_id()
  if (is.null(a$src$style_items) || !identical(a$src$style_items$md5, sid$md5))
    stop(sprintf("04 used the style-item list %s (md5 %s) but 15 finds %s (md5 %s); rerun 04 with the current list",
                 if (is.null(a$src$style_items)) "<none recorded>" else a$src$style_items$path, if (is.null(a$src$style_items)) "" else a$src$style_items$md5, sid$path, sid$md5), call. = FALSE)
  if (!identical(a$src$input$md5, fp$md5))
    stop(sprintf("04 read %s (md5 %s) but 15 is reading %s (md5 %s); rerun 04 on the same input", a$src$input$file, a$src$input$md5, fp$file, fp$md5), call. = FALSE)
  d <- as.data.table(read_dta(f))
  nm_low <- tolower(trimws(names(d)))
  getcol <- function(nm, required = TRUE) { hit <- which(nm_low == tolower(nm)); if (!length(hit)) { if (required) stop("column not found: ", nm, call. = FALSE); return(NULL) }; d[[hit[1]]] }
  num <- function(nm, required = FALSE) { x <- getcol(nm, required); if (is.null(x)) return(NULL); suppressWarnings(as.numeric(zap_labels(x))) }
  cn <- as.integer(zap_labels(getcol("CN")))
  r5 <- !is.na(num(RESP_W5)) | !is.na(num(MAR_W5))
  keep <- which(r5 & cn %in% c(1L, 2L))
  stopifnot(length(keep) == length(Tt))
  ## the risk set must be the same persons, in the same order, as in 04 (row positions, cohort, and the identifier if present)
  if (!identical(as.integer(keep), a$src$rows) || !identical(as.integer(cn[keep]), a$src$cn))
    stop("the w5 risk set differs from the one saved by 04 (rows or cohorts); rerun 04 on the same input", call. = FALSE)
  idc <- names(d)[toupper(names(d)) %in% toupper(ID_COLS)]
  if (!is.null(a$src$id) && (!length(idc) || !identical(as.character(zap_labels(d[[idc[1]]]))[keep], a$src$id)))
    stop("respondent identifiers of the risk set differ from those saved by 04", call. = FALSE)
  cat(sprintf("diag: risk set verified against 04 (%d persons; row positions and cohorts%s)\n", length(keep),
              if (!is.null(a$src$id)) " and identifiers" else ""))
  cat("diag: master", nrow(d), "rows;", "risk set", length(keep), "; cont total", sum(cn == 1L), "; add total", sum(cn == 2L), "\n")

  ## ---- 1. 波別応答フラグ(strict) と生存指標 ------------------------------------------
  resp <- sapply(1:10, function(w) {
    cols <- intersect(c(MARKERS[[paste0("w", w)]], paste0(WAVE_PREFIX[w], "q02")), nm_low)
    if (!length(cols)) return(rep(NA_integer_, nrow(d)))
    ok <- Reduce(`|`, lapply(cols, function(v) { x <- num(v); !is.na(x) & !(x %in% NA_MARKER_CODES) }))
    as.integer(ok)
  })
  colnames(resp) <- paste0("r_w", 1:10)
  cat("diag: マーカー被覆(波ごとの検出列数):", sapply(1:10, function(w) length(intersect(c(MARKERS[[paste0("w", w)]], paste0(WAVE_PREFIX[w], "q02")), nm_low))), "\n")
  cat("diag: 継続の応答率 w1..w10:", round(colMeans(resp[cn == 1L, ], na.rm = TRUE), 3), "\n")
  cat("diag: 追加の応答率 w5..w10:", round(colMeans(resp[cn == 2L, 5:10], na.rm = TRUE), 3), "\n")
  surv_all <- function(rows, from, to) { m <- resp[rows, from:to, drop = FALSE]; m[is.na(m)] <- 0L; as.integer(rowSums(m) == (to - from + 1)) }
  io_rows <- which(cn == 1L); in_rows <- which(cn == 2L)
  n_old <- length(io_rows); n_new <- length(in_rows)
  ## 継続: リスク集合(keep)内での位置
  pos_old <- match(io_rows, keep); pos_new <- match(in_rows, keep)   # NA = w5 非回答(継続)/ リスク集合外
  cat("diag: 継続で w5 リスク集合に入る人:", sum(!is.na(pos_old)), "/", n_old, "; 追加:", sum(!is.na(pos_new)), "/", n_new, "\n")
  S_exact <- surv_all(io_rows, 1, T_WAVE); S_any <- as.integer(!is.na(pos_old))
  S_next  <- resp[io_rows, T_WAVE + 1]                   # w6 回答(NA なら SSM 不能)
  Sm  <- surv_all(in_rows, T_WAVE + 1, T_WAVE + K_DOSE)  # 追加: w6–w9
  Sm1 <- surv_all(in_rows, T_WAVE + 1, T_WAVE + K_DOSE + 1)  # w6–w10
  cat("diag: S_exact(w1-5 完走)=", sum(S_exact), " S_any(w5 回答)=", sum(S_any), " 追加 s_m(w6-9)=", sum(Sm), " s_m1(w6-10)=", sum(Sm1), "\n")
  if (sum(S_exact & is.na(pos_old)) > 0) cat("diag: 注意: strict マーカーで w1-5 完走だが 04 リスク集合に入らない継続者", sum(S_exact & is.na(pos_old)), "人(y は NA 扱い)\n")

  ## ---- 2. 共変量(EC-adj 用: 男性・出生年(中心化)・入所時学歴+欠測フラグ) -------------------
  sex <- num("sex"); yb <- num("ybirth")
  educ_old <- num("zq23a"); educ_new <- num("dq69a")
  educ <- rep(NA_real_, nrow(d)); if (!is.null(educ_old)) educ[cn == 1L] <- educ_old[cn == 1L]; if (!is.null(educ_new)) educ[cn == 2L] <- educ_new[cn == 2L]
  educ[educ %in% c(88, 99)] <- NA
  has_educ <- any(!is.na(educ))
  if (!has_educ) cat("diag: 学歴変数(zq23a/dq69a)が見つからない → EC-adj は性・出生年のみ\n")
  mkX <- function(rows) {
    X <- data.frame(male = as.numeric(sex[rows] == 1), yb_c = yb[rows] - mean(yb, na.rm = TRUE))
    X$yb_c[is.na(X$yb_c)] <- 0; X$male[is.na(X$male)] <- mean(X$male, na.rm = TRUE)
    if (has_educ) { e <- educ[rows]; X$educ_mi <- as.numeric(is.na(e)); e[is.na(e)] <- median(educ, na.rm = TRUE); X$educ_c <- e - median(educ, na.rm = TRUE) }
    X
  }
  Xo <- mkX(io_rows); Xn <- mkX(in_rows); ADJ <- names(Xo)

  ## ---- 3. 入所波(w1)の対応変数と結果行列 ----------------------------------------------
  vm <- if (file.exists(VARMAP)) fread(VARMAP, encoding = "UTF-8") else { cat("diag: varmap not found:", VARMAP, "\n"); NULL }
  w1_of <- character(0)
  if (!is.null(vm) && all(c("w1", "w5") %in% names(vm))) {
    vv <- vm[nzchar(w1) & nzchar(w5), .(w5 = toupper(trimws(w5)), w1 = toupper(trimws(w1)))]
    vv <- vv[!duplicated(w5)]; w1_of <- setNames(vv$w1, vv$w5)
  }
  if (file.exists(OVERRIDES)) {
    ov <- fread(OVERRIDES, encoding = "UTF-8"); ov <- ov[nzchar(w5_var)]
    for (i in seq_len(nrow(ov))) w1_of[toupper(ov$w5_var[i])] <- toupper(ov$w1_var[i])
    ov_recode <- setNames(ov$recode_w1, toupper(ov$w5_var))
  } else ov_recode <- character(0)
  cat("diag: w5→w1 対応表:", length(w1_of), "件(varmap+overrides)\n")
  codes_of <- function(vl, pat) if (!is.null(vl)) unname(vl[grepl(pat, names(vl))]) else numeric(0)
  spec_codes <- function(vl) c(codes_of(vl, DK_PAT), codes_of(vl, REF_PAT), codes_of(vl, NR_PAT), codes_of(vl, NAP_PAT))
  vars <- meta$var; J <- length(vars)
  E_A <- matrix(NA_real_, n_old, J, dimnames = list(NULL, vars))   # 入所波の実質値(継続、全員)
  E_M <- E_DK <- matrix(NA_real_, n_old, J, dimnames = list(NULL, vars))
  ec_status <- data.table(var = vars, w1_var = NA_character_, ec_ok = FALSE, reason = "no w1 counterpart")
  dk_cnt <- dk_n <- miss_cnt <- miss_n <- rep(0, n_old)
  for (j in seq_len(J)) {
    v <- vars[j]; v1 <- w1_of[toupper(v)]
    if (is.na(v1) || !nzchar(v1)) next
    x1 <- getcol(v1, required = FALSE)
    if (is.null(x1)) { ec_status[j, `:=`(w1_var = v1, reason = "w1 variable not in file")]; next }
    x5 <- getcol(v, required = FALSE)
    vl1 <- attr(x1, "labels", exact = TRUE); vl5 <- attr(x5, "labels", exact = TRUE)
    xv1 <- suppressWarnings(as.numeric(zap_labels(x1)))[io_rows]
    rc <- ov_recode[toupper(v)]
    if (!is.na(rc) && nzchar(rc)) {
      if (rc == "reverse12") xv1 <- ifelse(xv1 == 1, 2, ifelse(xv1 == 2, 1, xv1))
      else if (grepl("^map:", rc)) { pairs <- strsplit(sub("^map:", "", rc), ";")[[1]]; from <- as.numeric(sub("=.*", "", pairs)); to <- as.numeric(sub(".*=", "", pairs)); xv1 <- ifelse(xv1 %in% from, to[match(xv1, from)], xv1) }
    }
    sp1 <- spec_codes(vl1); sp5 <- spec_codes(vl5)
    if (exists("NAP_OVERRIDE") && !is.null(NAP_OVERRIDE[[toupper(v)]])) sp5 <- union(sp5, NAP_OVERRIDE[[toupper(v)]])   # 04 と同じ非該当上書き
    vals1 <- if (!is.null(vl1)) sort(setdiff(unname(vl1), sp1)) else numeric(0); vals5 <- if (!is.null(vl5)) sort(setdiff(unname(vl5), sp5)) else numeric(0)
    if (!is.na(rc) && nzchar(rc) && grepl("^map:", rc)) {                       # 対応表の再符号化をラベル集合にも適用してから比較
      pairs <- strsplit(sub("^map:", "", rc), ";")[[1]]; from <- as.numeric(sub("=.*", "", pairs)); to <- as.numeric(sub(".*=", "", pairs))
      vals1 <- sort(unique(ifelse(vals1 %in% from, to[match(vals1, from)], vals1)))
    }
    st <- meta$scale_type[j]
    same <- (st == "continuous" && length(vals1) == 0) || (length(vals1) > 0 && length(vals1) == length(vals5) && all(vals1 == vals5))
    if (!same && !(st == "continuous")) { ec_status[j, `:=`(w1_var = v1, reason = sprintf("coding mismatch (w1 %d vs w5 %d substantive labels)", length(vals1), length(vals5)))]; next }
    dk1 <- codes_of(vl1, DK_PAT); nr1 <- codes_of(vl1, NR_PAT); nap1 <- codes_of(vl1, NAP_PAT); ref1 <- codes_of(vl1, REF_PAT)
    base1 <- !is.na(xv1) & !(xv1 %in% nap1)
    subst <- xv1; subst[xv1 %in% c(dk1, nr1, nap1, ref1)] <- NA
    if (st == "continuous") { qs <- quantile(subst, c(.01, .99), na.rm = TRUE); subst <- pmin(pmax(subst, qs[1]), qs[2]) }
    E_A[, j] <- subst
    E_M[base1, j] <- as.integer(xv1[base1] %in% nr1); E_DK[base1, j] <- as.integer(xv1[base1] %in% dk1)
    dk_cnt[base1] <- dk_cnt[base1] + (xv1[base1] %in% dk1); dk_n[base1] <- dk_n[base1] + 1
    miss_cnt[base1] <- miss_cnt[base1] + (xv1[base1] %in% nr1); miss_n[base1] <- miss_n[base1] + 1
    ec_status[j, `:=`(w1_var = v1, ec_ok = TRUE, reason = "")]
  }
  cat("diag: EC 可能項目:", sum(ec_status$ec_ok), "/", J, "\n")
  print(ec_status[ec_ok == FALSE, .N, by = reason])
  write_aggregate(ec_status[ec_ok == FALSE, .(var, label = meta$label[match(var, meta$var)], w1_var, reason)], "15_ec_unavailable.csv")

  ## ---- 4. 解析行列(継続 = 全入所者、追加 = 全員; 列 = 族×項目) -----------------------------
  fill_old <- function(Mat) { out <- matrix(NA_real_, n_old, ncol(Mat), dimnames = list(NULL, colnames(Mat))); ok <- !is.na(pos_old); out[ok, ] <- Mat[pos_old[ok], , drop = FALSE]; out }
  fill_new <- function(Mat) { out <- matrix(NA_real_, n_new, ncol(Mat), dimnames = list(NULL, colnames(Mat))); ok <- !is.na(pos_new); out[ok, ] <- Mat[pos_new[ok], , drop = FALSE]; out }
  ## ---- 3b. 人レベル様式指標(欠測率・DK 率は 04 の person_dq / 上の w1 計算; 中点・極端は事前指定リストから) -----------
  ## Main composites (since 2026-09-24): midpoint and extreme-category shares over the COMMON battery = the rating-scale
  ## items of R/15_style_items.csv whose codes were verified at w5 (04) and which have an entry-wave counterpart with
  ## the same codes (ec_ok), so that the entry-wave and comparison-wave composites, and all five designs, use one item
  ## set. The audit table lists every listed item with its codes and labels and the reason for any exclusion.
  st_ok <- style04[in_universe == TRUE & ok == TRUE]
  st_ok[, ec_ok := ec_status$ec_ok[match(var_file, ec_status$var)]]
  st_ok[, w1_var := ec_status$w1_var[match(var_file, ec_status$var)]]
  ## the entry-wave values of a common-battery item must lie in 1..k as well (the label sets are identical by ec_ok)
  st_ok[, w1_values_ok := vapply(seq_len(.N), function(i) !ec_ok[i] || style_values_ok(E_A[, var_file[i]], k[i]), logical(1))]
  st_ok[ec_ok & !w1_values_ok, ec_ok := FALSE]
  k_of <- setNames(as.list(st_ok$k), st_ok$var_file)
  sets <- list(
    rating_all      = list(ext = st_ok[set == "rating", var_file],                 mid = st_ok[set == "rating" & mid_ok, var_file], common = FALSE),
    rating_common   = list(ext = st_ok[set == "rating" & ec_ok, var_file],         mid = st_ok[set == "rating" & mid_ok & ec_ok, var_file], common = TRUE),
    bipolar_common  = list(ext = st_ok[battery %in% c("agree5", "satis5") & ec_ok, var_file], mid = st_ok[battery %in% c("agree5", "satis5") & mid_ok & ec_ok, var_file], common = TRUE),
    ## the agreement items alone: asked of every respondent (no routing by employment or marital status), so every
    ## person's composite is taken over the same items apart from item nonresponse
    agree_common    = list(ext = st_ok[battery == "agree5" & ec_ok, var_file],               mid = st_ok[battery == "agree5" & mid_ok & ec_ok, var_file], common = TRUE),
    ratingfreq_all  = list(ext = st_ok[, var_file],                                mid = st_ok[mid_ok == TRUE, var_file], common = FALSE),
    ratingfreq_common = list(ext = st_ok[ec_ok == TRUE, var_file],                 mid = st_ok[mid_ok & ec_ok, var_file], common = TRUE))
  ## the rule used up to v0.7 (2026-09-24): every universe item with four to seven consecutively numbered substantive
  ## codes, whatever its meaning; kept for the sensitivity file only. "v07_asrun" pairs the w5 composite over all such
  ## items with the w1 composite over those with an entry-wave counterpart (the pairing of the earlier runs);
  ## "v07_common" uses the common battery at both waves.
  leg <- vars[vapply(vars, function(v) { cv <- codes_of_var[[v]]; length(cv) %in% 4:7 && all(diff(cv) == 1) }, logical(1))]
  leg_k <- setNames(as.list(vapply(leg, function(v) length(codes_of_var[[v]]), 1L)), leg)
  leg_lo <- setNames(as.list(vapply(leg, function(v) codes_of_var[[v]][1], 1)), leg)
  leg_mid <- leg[vapply(leg, function(v) leg_k[[v]] %% 2L == 1L, logical(1))]
  leg_ec <- leg[ec_status$ec_ok[match(leg, ec_status$var)]]; leg_mid_ec <- intersect(leg_mid, leg_ec)
  comp5 <- function(ext, mid, kk = k_of, lo = NULL) style_shares(Y, ext, mid, kk, lo)
  comp1 <- function(ext, mid, kk = k_of, lo = NULL) style_shares(E_A, ext, mid, kk, lo)
  main5 <- comp5(sets$rating_common$ext, sets$rating_common$mid); main1 <- comp1(sets$rating_common$ext, sets$rating_common$mid)
  ## consistency with 04: its person_dq composites are the "rating_all" composites recomputed here
  chk5 <- comp5(sets$rating_all$ext, sets$rating_all$mid)
  stopifnot(isTRUE(all.equal(chk5$mid, pdq$mid_share)), isTRUE(all.equal(chk5$ext, pdq$ext_share)))
  P5 <- as.matrix(pdq[, .(pdq_miss_share = miss_share, pdq_dk_share = dk_share)]); P5 <- cbind(P5, pdq_mid_share = main5$mid, pdq_ext_share = main5$ext)
  P1 <- cbind(pdq_miss_share = ifelse(miss_n > 0, miss_cnt / miss_n, NA), pdq_dk_share = ifelse(dk_n > 0, dk_cnt / dk_n, NA),
              pdq_mid_share = main1$mid, pdq_ext_share = main1$ext)
  ## sensitivity composites (family P_style_sens; estimated and bootstrapped with everything else, reported apart)
  S5 <- NULL; S1 <- NULL; sens_meta <- list()
  add_sens <- function(name, c5, c1, n5, n1) {
    S5 <<- cbind(S5, c5$mid, c5$ext); S1 <<- cbind(S1, if (is.null(c1)) rep(NA_real_, n_old) else c1$mid, if (is.null(c1)) rep(NA_real_, n_old) else c1$ext)
    sens_meta[[length(sens_meta) + 1]] <<- data.table(col = paste0("sty_", name, c("_mid", "_ext")), battery = name, indicator = c("mid", "ext"),
                                                      n_items_w5 = c(n5$n_mid, n5$n_ext), n_items_w1 = if (is.null(n1)) NA_integer_ else c(n1$n_mid, n1$n_ext),
                                                      has_entry = !is.null(c1))
  }
  for (nm in names(sets)) { ss <- sets[[nm]]; c5 <- comp5(ss$ext, ss$mid); c1 <- if (ss$common) comp1(ss$ext, ss$mid) else NULL; add_sens(nm, c5, c1, c5, c1) }
  v07_5 <- comp5(leg, leg_mid, leg_k, leg_lo); v07_1as <- comp1(leg_ec, leg_mid_ec, leg_k, leg_lo)
  add_sens("v07_asrun", v07_5, v07_1as, v07_5, v07_1as)
  v07c5 <- comp5(leg_ec, leg_mid_ec, leg_k, leg_lo); add_sens("v07_common", v07c5, v07_1as, v07c5, v07_1as)
  sens_meta <- rbindlist(sens_meta); colnames(S5) <- colnames(S1) <- sens_meta$col
  sens_has_entry <- sens_meta$has_entry
  cat(sprintf("diag: 様式指標: 評定尺度 %d 項目(中点 %d)、共通バッテリー %d 項目(中点 %d); v0.7 規則 %d 項目(w1 対応 %d)\n",
              length(sets$rating_all$ext), length(sets$rating_all$mid), length(sets$rating_common$ext), length(sets$rating_common$mid), length(leg), length(leg_ec)))
  ## audit table: every listed item, its codes and labels as read from the file, and where it entered
  audit <- copy(style04)
  audit[, `:=`(w1_var = ec_status$w1_var[match(var_file, ec_status$var)], entry_wave_ok = var_file %in% st_ok[ec_ok == TRUE, var_file])]
  audit[var_file %in% st_ok[w1_values_ok == FALSE, var_file], reason := style_join_reason(reason, "entry-wave values outside 1..k observed: outside the common battery")]
  audit[, `:=`(in_main = var_file %in% sets$rating_common$ext, in_main_midpoint = var_file %in% sets$rating_common$mid,
               in_rating_all = var_file %in% sets$rating_all$ext, in_v07_rule = var_file %in% leg)]
  audit[, label := meta$label[match(var_file, meta$var)]]
  audit <- audit[, .(var, label, set, battery, k, midpoint_declared = midpoint, in_universe, codes_verified = ok, midpoint_verified = mid_ok,
                     codes, labels, w1_var, entry_wave_ok, in_main, in_main_midpoint, in_rating_all, in_v07_rule, reason, note)]
  ## items the v0.7 rule admitted that are not on the list (nominal codes, classifications, quantity bands, ...)
  extra <- setdiff(leg, style04$var_file)
  if (length(extra)) {
    code_labels <- a$code_labels
    audit <- rbind(audit, data.table(var = toupper(extra), label = meta$label[match(extra, meta$var)], set = "", battery = "", k = vapply(extra, function(v) leg_k[[v]], 1L),
                                     midpoint_declared = NA, in_universe = TRUE, codes_verified = NA, midpoint_verified = NA,
                                     codes = vapply(extra, function(v) paste(codes_of_var[[v]], collapse = " "), ""),
                                     labels = vapply(extra, function(v) paste(code_labels[[v]], collapse = " | "), ""),
                                     w1_var = ec_status$w1_var[match(extra, ec_status$var)], entry_wave_ok = ec_status$ec_ok[match(extra, ec_status$var)] %in% TRUE,
                                     in_main = FALSE, in_main_midpoint = FALSE, in_rating_all = FALSE, in_v07_rule = TRUE,
                                     reason = "not on the list: entered the composites under the v0.7 rule only", note = ""), fill = TRUE)
  }
  fam <- c(rep("A_substantive", J), rep("B_itemnonresp", J), rep("B_dk", J), rep("P_style", 4), rep("P_style_sens", ncol(S5)))
  cols <- c(vars, paste0(vars, "__M"), paste0(vars, "__DK"), colnames(P5), colnames(S5))
  Yo <- cbind(fill_old(Y), fill_old(M), fill_old(DK), fill_old(P5), fill_old(S5)); colnames(Yo) <- cols
  Yn <- cbind(fill_new(Y), fill_new(M), fill_new(DK), fill_new(P5), fill_new(S5)); colnames(Yn) <- cols
  Eo <- cbind(E_A, E_M, E_DK, P1, S1); colnames(Eo) <- cols
  ec_ok_col <- c(ec_status$ec_ok, ec_status$ec_ok, ec_status$ec_ok, rep(TRUE, 4), sens_has_entry)
  Eo[, !ec_ok_col] <- NA
  binary_col <- c(meta$scale_type == "binary", rep(TRUE, 2 * J), rep(FALSE, 4 + ncol(S5)))
  JJ <- ncol(Yo)

  ## ---- 5. ベクトル化エンジン(pc_point と同じ定義) ----------------------------------------
  cm <- function(Mat, mask) { Mm <- Mat; Mm[!mask] <- NA; s <- colSums(Mm, na.rm = TRUE); n <- colSums(!is.na(Mm)); r <- s / n; r[n == 0] <- NA; r }
  engine <- function(ic, jf, S, adj = TRUE) {
    yo <- Yo[ic, , drop = FALSE]; yn <- Yn[jf, , drop = FALSE]; eo <- Eo[ic, , drop = FALSE]
    Sv <- S[ic]; S1v <- S_next[ic]; Smv <- Sm[jf]; Sm1v <- Sm1[jf]
    IO <- (Sv == 1L) & !is.na(yo)                                   # 項目別: 生存 & y 観測
    m_oS <- cm(yo, IO)
    m_n  <- colMeans(yn, na.rm = TRUE)
    m_nSm <- cm(yn, matrix(Smv == 1L, nrow(yn), JJ))
    naive <- m_oS - m_n; sm <- m_oS - m_nSm
    ## SSM
    IO1 <- IO & !is.na(S1v) & S1v == 1L
    m_oS1 <- cm(yo, IO1); m_nSm1 <- cm(yn, matrix(Smv == 1L & !is.na(Sm1v) & Sm1v == 1L, nrow(yn), JJ))
    ssm <- m_oS1 - m_nSm1
    ## EC
    ES <- eo; ES[!IO] <- NA
    delta <- colMeans(ES, na.rm = TRUE) - colMeans(eo, na.rm = TRUE)
    delta[colSums(!is.na(ES)) == 0] <- NA
    ec <- (m_oS - delta) - m_n
    ## EC-adj(回帰標準化; pc_point の ols_pred と同じ)
    ec_adj <- setNames(rep(NA_real_, JJ), cols)
    if (adj) {
      Xo_b <- cbind(1, as.matrix(Xo[ic, , drop = FALSE])); Xn_b <- cbind(1, as.matrix(Xn[jf, , drop = FALSE]))
      for (j in which(ec_ok_col)) {
        e_j <- eo[, j]; ok_e <- !is.na(e_j); ip <- IO[, j] & ok_e
        if (sum(ok_e) <= ncol(Xo_b) || !any(ip)) next
        b1 <- stats::lm.fit(Xo_b[ok_e, , drop = FALSE], e_j[ok_e])$coefficients; b1[is.na(b1)] <- 0
        m_c <- mean(Xo_b[ip, , drop = FALSE] %*% b1)
        y_j <- yn[, j]; ok_y <- !is.na(y_j)
        if (sum(ok_y) <= ncol(Xn_b) || !any(IO[, j])) next
        b0 <- stats::lm.fit(Xn_b[ok_y, , drop = FALSE], y_j[ok_y])$coefficients; b0[is.na(b0)] <- 0
        m_0 <- mean(Xo_b[IO[, j], , drop = FALSE] %*% b0)
        ec_adj[j] <- m_oS[j] - (mean(e_j[ip]) - m_c) - m_0
      }
    }
    list(naive = naive, sm = sm, ssm = ssm, ec = ec, ec_adj = ec_adj, delta = delta, d1 = sm - ec, d2 = ssm - sm,
         n_oS = colSums(IO), n_oS1 = colSums(IO1), p_surv = colMeans(IO))
  }

  ## ---- 6. 用量定義ごとに実行 -------------------------------------------------------------
  res_all <- list(); tests_all <- list(); val_all <- list(); arms_all <- list(); sens_all <- list(); style_sens_all <- list()
  for (dose_def in DOSE_DEFS) {
    S <- if (dose_def == "exact") S_exact else S_any
    cat(sprintf("\n== dose_def = %s: old survivors (S=1) = %d / %d; new = %d (s_m=1: %d, s_m1=1: %d)\n", dose_def, sum(S), n_old, n_new, sum(Sm), sum(Sm1)))
    ## (a) panelcond 点推定 + 解析 SE(項目ごと)
    pc <- lapply(seq_len(JJ), function(j) {
      old <- data.frame(y = Yo[, j], y_entry = Eo[, j], s_prior = S, s_next = S_next)
      new <- data.frame(y = Yn[, j], s_m = Sm, s_m1 = Sm1)
      if (sum(S == 1L & !is.na(old$y)) < 30 || sum(!is.na(new$y)) < 30) return(NULL)
      use_adj <- ec_ok_col[j] && sum(!is.na(Eo[, j])) > 50
      fit <- tryCatch(suppressMessages(suppressWarnings(pc_estimate(cbind(old, Xo), cbind(new, Xn), k = K_DOSE, m = K_DOSE,
                                                                     adjust = if (use_adj) ADJ else NULL, nboot = 0))), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      e <- fit$estimates; tt <- fit$tests
      data.table(col = cols[j], family = fam[j], var = sub("__(M|DK)$", "", cols[j]), estimator = e$estimator, estimate = e$estimate, se_analytic = e$se,
                 n_old = e$n_old, n_new = e$n_new, delta_entry = fit$delta_entry, n_pairs = fit$info$n_pairs, p_survive = fit$info$p_survive,
                 T1_stat = tt$statistic[1], T1_p = tt$p[1], T2_stat = tt$statistic[2], T2_p = tt$p[2])
    })
    pc <- rbindlist(pc); pc <- pc[estimator != "ipw" & !is.na(estimate)]     # EC 不能項目の ec/ec_adj 行は落とす
    cat("diag: panelcond 推定完了:", uniqueN(pc$col), "列 ×", uniqueN(pc$estimator), "推定量 (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min )\n")
    ## (b) フルデータでのエンジン照合
    full <- engine(seq_len(n_old), seq_len(n_new), S)
    val <- rbindlist(lapply(c("naive", "sm", "ssm", "ec", "ec_adj"), function(es) {
      a_ <- pc[estimator == es, .(col, estimate)]; b_ <- full[[es]][a_$col]
      rel <- abs(a_$estimate - b_) / pmax(1, abs(a_$estimate))      # 円単位の項目(年収 3e7 等)は絶対差で見ると浮動小数点誤差が出るので相対差
      data.table(dose_def = dose_def, estimator = es, n_compared = sum(!is.na(a_$estimate) & !is.na(b_)),
                 max_abs_diff = max(abs(a_$estimate - b_), na.rm = TRUE), max_rel_diff = max(rel, na.rm = TRUE), n_na_mismatch = sum(is.na(a_$estimate) != is.na(b_)))
    }))
    print(val)
    if (any(val$max_rel_diff > 1e-8, na.rm = TRUE)) warning("engine と pc_estimate の点推定が一致しません(15_validation.csv 参照)")
    val_all[[dose_def]] <- val
    ## (c) 同時人ブートストラップ
    set.seed(SEED15)
    keys <- c("naive", "sm", "ssm", "ec", "ec_adj", "delta", "d1", "d2")
    boot <- array(NA_real_, c(B_BOOT, JJ, length(keys)), dimnames = list(NULL, cols, keys))
    for (b in seq_len(B_BOOT)) {
      ic <- sample.int(n_old, n_old, TRUE); jf <- sample.int(n_new, n_new, TRUE)
      r <- engine(ic, jf, S)
      for (kk in keys) boot[b, , kk] <- r[[kk]]
      if (b %% 50 == 0) cat(sprintf("   replicate %d / %d (%.1f min)\n", b, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    sdb <- function(x) if (sum(!is.na(x)) >= 2) sd(x, na.rm = TRUE) else NA_real_
    se_boot <- apply(boot, c(2, 3), sdb)
    ## (d) 標準化 SD(フルデータ: 生存者 vs 新規のプール SD)、漏斗比
    IOf <- (S == 1L) & !is.na(Yo)
    sd_pool <- sapply(seq_len(JJ), function(j) { yo <- Yo[IOf[, j], j]; yn <- Yn[!is.na(Yn[, j]), j]; sqrt((var(yo) + var(yn)) / 2) })
    ## Mass-domination diagnostic (sample analogue of Corollary 4): max over categories a of p * Q(a) / F2(a),
    ## computed by mass_ratio_item() (15_mass_ratio.R), which keeps positive/zero categories (flagged, never
    ## dropped), sparse categories (fewer than MIN_CELL fresh-cohort respondents) and supported ones apart.
    ## (Before 2026-09-23 positive/zero categories were silently discarded.)
    ## Retention for the item's own population: an item asked only of a subgroup G (e.g. those with a partner) compares
    ## the stayers' sub-measure with the refreshment distribution WITHIN G, so p must be P(stayer and answer | G) =
    ## [answering stayers / cohort size] / P(G), with P(G) estimated by the share of fresh entrants who reached the
    ## question (not routed past it: column M of 04 is defined exactly for them). For items asked of everyone P(G) = 1.
    ## (Before 2026-09-24 the cohort-wide rate was used for every item, which diluted the ratio of subgroup items.)
    ## Fresh item nonresponse (since 2026-09-24): the fresh distribution above is that of entrants who gave a
    ## substantive answer, while p refers to everyone in G. The reach and answer shares of both arms are recorded
    ## (disclosure-safe shares), and the identity map is also evaluated with the fresh missing mass allocated
    ## freely (mass_ratio_item(): needed_mass <= missing_mass). The observed reach rate stands in for eligibility;
    ## a true item nonresponse recorded as NA rather than with a no-answer code cannot be told from routing (04).
    reached_new <- colSums(!is.na(Yn[, J + seq_len(J), drop = FALSE]))
    elig_new <- reached_new / n_new
    reach_old_S  <- colSums(!is.na(Yo[S == 1L, J + seq_len(J), drop = FALSE])) / sum(S)          # survivors who reached the item
    ans_old_S    <- colSums(IOf[, seq_len(J), drop = FALSE])                                        # survivors who answered
    ans_new      <- colSums(!is.na(Yn[, seq_len(J), drop = FALSE]))                                 # fresh entrants who answered
    funnel <- funnel_sup <- p_item <- rep(NA_real_, JJ); funnel_zd <- funnel_sp <- funnel_ncat <- rep(NA_integer_, JJ)
    funnel_miss <- funnel_need <- rep(NA_real_, JJ); funnel_feas <- rep(NA, JJ)
    reach_new_v <- reach_oldS_v <- ans_new_v <- ans_oldS_v <- rep(NA_real_, JJ)
    for (j in seq_len(J)) {                                          # A 族の離散項目のみ
      if (!(elig_new[j] > 0)) next
      p_item[j] <- full$p_surv[j] / elig_new[j]
      reach_new_v[j] <- elig_new[j]; reach_oldS_v[j] <- reach_old_S[j]
      ans_new_v[j] <- ans_new[j] / reached_new[j]
      ans_oldS_v[j] <- if (reach_old_S[j] > 0) ans_old_S[j] / (reach_old_S[j] * sum(S)) else NA_real_
      mr <- mass_ratio_item(Yo[IOf[, j], j], Yn[!is.na(Yn[, j]), j], p_item[j], max_cat = 9L, min_fresh = MIN_CELL,
                            n_new_reached = reached_new[j])
      if (!mr$eligible) next
      funnel[j] <- mr$ratio_finite; funnel_sup[j] <- mr$ratio_supported
      funnel_zd[j] <- mr$n_zero_denom; funnel_sp[j] <- mr$n_sparse_gt1; funnel_ncat[j] <- mr$ncat
      funnel_miss[j] <- mr$missing_mass; funnel_need[j] <- mr$needed_mass; funnel_feas[j] <- mr$identity_feasible
    }
    ## (e) 表の組み立て
    pc[, `:=`(dose_def = dose_def)]
    pc[, se_boot := se_boot[cbind(match(col, cols), match(estimator, keys))]]
    pc[, sd_pool := sd_pool[match(col, cols)]]
    pc[, d_std := estimate / sd_pool]
    pc[, funnel_ratio := funnel[match(col, cols)]]
    pc[, funnel_ratio_supported := funnel_sup[match(col, cols)]]
    pc[, funnel_p := p_item[match(col, cols)]]
    pc[, funnel_zero_denom := funnel_zd[match(col, cols)]]
    pc[, funnel_sparse_gt1 := funnel_sp[match(col, cols)]]
    pc[, funnel_ncat := funnel_ncat[match(col, cols)]]
    pc[, funnel_missing_mass := funnel_miss[match(col, cols)]]
    pc[, funnel_needed_mass := funnel_need[match(col, cols)]]
    pc[, funnel_identity_feasible := funnel_feas[match(col, cols)]]
    pc[, reach_new := reach_new_v[match(col, cols)]]
    pc[, reach_old_S := reach_oldS_v[match(col, cols)]]
    pc[, answer_new_given_reach := ans_new_v[match(col, cols)]]
    pc[, answer_old_S_given_reach := ans_oldS_v[match(col, cols)]]
    pc[, se := fifelse(is.na(se_boot), se_analytic, se_boot)]
    pc[, z := estimate / se]; pc[, p := 2 * pnorm(-abs(z))]
    pc[, binary := binary_col[match(col, cols)]]
    pc[, sesoi := fifelse(binary | family == "P_style", SESOI_PP, SESOI_D * sd_pool)]
    pc[, p_tost := pmax(pnorm(-(estimate + sesoi) / se), pnorm((estimate - sesoi) / se))]
    ## items listed in 15_item_exclude.csv (nominal or paradata codes: a mean contrast is not meaningful) stay in the
    ## item files with exclude_flag, but are outside every multiplicity family and every count (since 2026-09-23;
    ## before, they were excluded from the counts but still entered the BH adjustment of the other items)
    if (file.exists(EXCLUDE)) { ex <- fread(EXCLUDE, encoding = "UTF-8"); pc[, exclude_flag := toupper(var) %in% toupper(ex$var)] } else pc[, exclude_flag := FALSE]
    ## 04 のメタにある衛生フラグ(調査票の版差・ルーティング疑い)を項目に付ける(B 族の解釈に必須)
    if ("filter_mismatch" %in% names(meta)) pc[, filter_mismatch := meta$filter_mismatch[match(var, meta$var)]]
    if ("nr_routing_flag" %in% names(meta)) pc[, nr_routing_flag := meta$nr_routing_flag[match(var, meta$var)]] else pc[, nr_routing_flag := NA]
    ## items outside the counts (since 2026-09-24): those listed in 15_item_exclude.csv (nominal codes, date and
    ## clock-time components, duplicate recodes) and those flagged as potentially incomparable between the arms
    ## because of differential item nonresponse or routing (04's nr_routing_flag: the item-nonresponse rate differs
    ## between continuing respondents and entrants by more than the threshold; the questionnaire filters behind the
    ## difference have not been verified, so the flag is a screening rule, not a documented routing difference).
    ## A follow-up component of a flagged question (the same variable name plus a component suffix, e.g. DQ46Y,
    ## years of premarital cohabitation, after DQ46) inherits the flag: its non-routed respondents are coded missing
    ## or not applicable rather than with the no-answer code, so 04's rate rule does not see it (since 2026-09-24).
    ## classify() applies the exclusion, the BH families and the three-way classification for a given routing
    ## flag; the main run uses 04's flag, and the threshold sensitivity below re-applies it with other thresholds.
    classify <- function(pc, routing) {
      pc <- copy(pc); pc[, nr_routing_flag := routing]
      rf <- unique(toupper(pc[nr_routing_flag %in% TRUE, var]))
      fu <- function(v) { if (!length(rf)) return(FALSE); v <- toupper(v); any(startsWith(v, rf) & grepl("^[A-Z_][A-Z0-9_]*$", substring(v, nchar(rf) + 1L)) & nchar(v) > nchar(rf)) }
      pc[, routing_followup := !(nr_routing_flag %in% TRUE) & vapply(var, fu, logical(1))]
      pc[, count_exclude := exclude_flag | (nr_routing_flag %in% TRUE) | routing_followup]
      pc[, `:=`(q = NA_real_, q_tost = NA_real_)]
      pc[count_exclude == FALSE, q := p.adjust(p, "BH"), by = .(family, estimator)]
      pc[count_exclude == FALSE, q_tost := p.adjust(p_tost, "BH"), by = .(family, estimator)]
      pc[, class3 := fifelse(is.na(se) | se <= 0, "n/a (no variation)", fifelse(q < Q_CUT, "affected", fifelse(q_tost < Q_CUT, "equivalent", "undetermined")))]
      pc[exclude_flag == TRUE & !(is.na(se) | se <= 0), class3 := "excluded (listed code)"]
      pc[exclude_flag == FALSE & count_exclude == TRUE & !(is.na(se) | se <= 0), class3 := "excluded (routing flag)"]
      pc
    }
    pc <- classify(pc, pc$nr_routing_flag)
    ## the sensitivity composites leave the item tables and the counts here (their own file below)
    ps <- pc[family == "P_style_sens"]; pc <- pc[family != "P_style_sens"]
    ps <- merge(ps[, .(dose_def, col, estimator, estimate, se_boot, sd_pool, d_std, p, n_old, n_new)], sens_meta[, .(col, battery, indicator, n_items_w5, n_items_w1)], by = "col")
    main_rows <- pc[family == "P_style" & var %in% c("pdq_mid_share", "pdq_ext_share"), .(dose_def, col, estimator, estimate, se_boot, sd_pool, d_std, p, n_old, n_new)]
    main_rows[, `:=`(battery = "rating_common (main)", indicator = fifelse(col == "pdq_mid_share", "mid", "ext"))]
    main_rows[, n_items_w5 := fifelse(indicator == "mid", length(sets$rating_common$mid), length(sets$rating_common$ext))]; main_rows[, n_items_w1 := n_items_w5]
    style_sens_all[[dose_def]] <- rbind(main_rows, ps)[order(battery != "rating_common (main)", battery, indicator, match(estimator, c("naive", "sm", "ssm", "ec", "ec_adj")))]
    if (any(pc$routing_followup)) cat("diag: follow-ups of routing-flagged questions (outside the counts):", paste(unique(pc[routing_followup == TRUE, var]), collapse = ", "), "\n")
    pc[, label := meta$label[match(var, meta$var)]]
    res_all[[dose_def]] <- pc
    ## Routing-threshold sensitivity (since 2026-09-24): the same counts with the flag recomputed at other
    ## thresholds of the item-nonresponse gap (04 uses .15; the gap is recomputed here over the same risk set
    ## and compared with 04's flag at .15), and with no routing exclusion at all.
    Mo <- Yo[, J + seq_len(J), drop = FALSE]; Mn <- Yn[, J + seq_len(J), drop = FALSE]
    nr_gap <- abs(colMeans(Mo, na.rm = TRUE) - colMeans(Mn, na.rm = TRUE)); nr_gap[!is.finite(nr_gap)] <- 0
    names(nr_gap) <- vars
    gap_of <- nr_gap[match(pc$var, vars)]; gap_of[is.na(gap_of)] <- 0
    sens <- rbindlist(lapply(c(0.10, 0.15, 0.20, Inf), function(thr) {
      pcs <- classify(pc, if (is.finite(thr)) (gap_of > thr) else rep(FALSE, nrow(pc)))
      A <- pcs[family == "A_substantive" & class3 != "n/a (no variation)" & !count_exclude]
      cnt <- A[, .(n_items = .N, n_affected = sum(class3 == "affected")), by = estimator]
      common <- A[estimator == "ec", unique(var)]
      cc <- sapply(c("naive", "sm", "ssm", "ec", "ec_adj"), function(e) A[estimator == e & var %in% common & class3 == "affected", uniqueN(var)])
      five <- Reduce(intersect, lapply(c("naive", "sm", "ssm", "ec", "ec_adj"), function(e) A[estimator == e & class3 == "affected", unique(var)]))
      u <- unique(pcs[family == "A_substantive", .(var, nr_routing_flag, routing_followup, exclude_flag)])
      data.table(dose_def = dose_def, threshold = thr,
                 n_routing = sum(u$nr_routing_flag %in% TRUE & !u$exclude_flag), n_followup = sum(u$routing_followup & !u$exclude_flag),
                 n_disagree_with_04 = if (is.finite(thr) && abs(thr - 0.15) < 1e-9) sum((gap_of > thr) != (pc$nr_routing_flag %in% TRUE)) else NA_integer_,
                 n_items_naive = cnt[estimator == "naive", n_items], n_items_sm = cnt[estimator == "sm", n_items], n_items_ec = cnt[estimator == "ec", n_items],
                 aff_naive = cnt[estimator == "naive", n_affected], aff_sm = cnt[estimator == "sm", n_affected], aff_ssm = cnt[estimator == "ssm", n_affected],
                 aff_ec = cnt[estimator == "ec", n_affected], aff_ec_adj = cnt[estimator == "ec_adj", n_affected],
                 common_naive = cc[["naive"]], common_sm = cc[["sm"]], common_ssm = cc[["ssm"]], common_ec = cc[["ec"]], common_ec_adj = cc[["ec_adj"]],
                 n_all5 = length(five), all5_items = paste(sort(five), collapse = "; "))
    }))
    sens_all[[dose_def]] <- sens
    ## 診断表(項目ごと; ブートストラップ SE 版の T1/T2)
    te <- unique(pc[, .(dose_def, col, family, var, label, T1_stat, T1_p, T2_stat, T2_p, delta_entry, n_pairs, p_survive, funnel_ratio, funnel_ratio_supported, funnel_p, funnel_zero_denom, funnel_sparse_gt1, funnel_ncat, funnel_missing_mass, funnel_needed_mass, funnel_identity_feasible, reach_new, reach_old_S, answer_new_given_reach, answer_old_S_given_reach, exclude_flag, count_exclude, nr_routing_flag)])
    te[, `:=`(d1 = full$d1[match(col, cols)], d2 = full$d2[match(col, cols)],
              d1_se_boot = se_boot[cbind(match(col, cols), match("d1", keys))], d2_se_boot = se_boot[cbind(match(col, cols), match("d2", keys))],
              delta_se_boot = se_boot[cbind(match(col, cols), match("delta", keys))])]
    te[, `:=`(T1_p_boot = 2 * pnorm(-abs(d1 / d1_se_boot)), T2_p_boot = 2 * pnorm(-abs(d2 / d2_se_boot)))]
    tests_all[[dose_def]] <- te
    arms_all[[dose_def]] <- data.table(dose_def = dose_def, n_old_total = n_old, n_old_S = sum(S), n_old_S_next = sum(S == 1L & S_next %in% 1L),
                                       n_new_total = n_new, n_new_sm = sum(Sm), n_new_sm1 = sum(Sm1), n_ec_items = sum(ec_status$ec_ok), B = B_BOOT, seed = SEED15)
  }
  res <- rbindlist(res_all); tests <- rbindlist(tests_all); arms <- rbindlist(arms_all)

  ## ---- 7. 出力 -----------------------------------------------------------------------------
  out_items <- res[, .(dose_def, family, var, label, estimator, estimate = round(estimate, 5), se_analytic = round(se_analytic, 5), se_boot = round(se_boot, 5),
                       d_std = round(d_std, 4), p = signif(p, 4), q = signif(q, 4), p_tost = signif(p_tost, 4), q_tost = signif(q_tost, 4), class3,
                       n_old, n_new, n_pairs, funnel_ratio = round(funnel_ratio, 3), funnel_ratio_supported = round(funnel_ratio_supported, 3),
                       funnel_zero_denom, funnel_sparse_gt1, funnel_ncat,
                       exclude_flag, count_exclude,
                       filter_mismatch = if ("filter_mismatch" %in% names(res)) filter_mismatch else NA,
                       nr_routing_flag = if ("nr_routing_flag" %in% names(res)) nr_routing_flag else NA,
                       routing_followup)]
  write_aggregate(out_items, "15_designs_items.csv")
  write_aggregate(tests[, .(dose_def, family, var, label, delta_entry = round(delta_entry, 5), delta_se_boot = round(delta_se_boot, 5), n_pairs, p_survive = round(p_survive, 4),
                            T1_stat = round(T1_stat, 3), T1_p = signif(T1_p, 4), T1_p_boot = signif(T1_p_boot, 6), T2_stat = round(T2_stat, 3), T2_p = signif(T2_p, 4), T2_p_boot = signif(T2_p_boot, 6),
                            funnel_ratio = round(funnel_ratio, 3), funnel_ratio_supported = round(funnel_ratio_supported, 3), funnel_p = round(funnel_p, 4),
                            funnel_zero_denom, funnel_sparse_gt1, funnel_ncat, funnel_missing_mass = round(funnel_missing_mass, 4), funnel_needed_mass = round(funnel_needed_mass, 4),
                            funnel_identity_feasible, reach_new = round(reach_new, 4), reach_old_S = round(reach_old_S, 4),
                            answer_new_given_reach = round(answer_new_given_reach, 4), answer_old_S_given_reach = round(answer_old_S_given_reach, 4),
                            exclude_flag, count_exclude, nr_routing_flag)], "15_tests_items.csv")
  write_aggregate(rbindlist(sens_all), "15_routing_sensitivity.csv", exempt = grep("^(n_|aff_|common_)", names(rbindlist(sens_all)), value = TRUE))
  ## response-style composites under alternative item sets (item counts are numbers of items, not persons)
  sty <- rbindlist(style_sens_all)
  sty <- sty[, .(dose_def, battery, indicator, estimator, n_items_w5, n_items_w1, estimate = round(estimate, 5), se_boot = round(se_boot, 5),
                 sd_pool = round(sd_pool, 5), d_std = round(d_std, 4), p = signif(p, 4), n_old, n_new)]
  write_aggregate(sty, "15_style_sensitivity.csv", exempt = c("n_items_w5", "n_items_w1"))
  write_aggregate(audit, "15_style_items_audit.csv")
  det <- res[class3 != "n/a (no variation)" & !count_exclude, .(n_items = .N, n_affected = sum(class3 == "affected", na.rm = TRUE), n_equivalent = sum(class3 == "equivalent", na.rm = TRUE),
                 n_undetermined = sum(class3 == "undetermined", na.rm = TRUE),
                 median_abs_d = round(median(abs(d_std), na.rm = TRUE), 4), share_positive = round(mean(estimate > 0, na.rm = TRUE), 3)), by = .(dose_def, family, estimator)]
  write_aggregate(det, "15_detection_counts.csv", exempt = c("n_items", "n_affected", "n_equivalent", "n_undetermined"))
  ## the diagnostics are contrasts of means: the same items as in the detection counts; n_T1/n_T2 are the numbers
  ## of items on which each statistic exists (T1 = SM - EC needs the entry-wave arm)
  diag <- tests[count_exclude == FALSE, .(n_items = .N, n_T1_tests = sum(!is.na(T1_p_boot)), n_T2_tests = sum(!is.na(T2_p_boot)),
                    n_T1_reject = sum(T1_p_boot < .05, na.rm = TRUE), n_T2_reject = sum(T2_p_boot < .05, na.rm = TRUE), share_T1_p05 = round(mean(T1_p_boot < .05, na.rm = TRUE), 3), share_T2_p05 = round(mean(T2_p_boot < .05, na.rm = TRUE), 3),
                    share_T1_positive = round(mean(d1 > 0, na.rm = TRUE), 3), share_T2_positive = round(mean(d2 > 0, na.rm = TRUE), 3),
                    median_delta_entry = round(median(delta_entry, na.rm = TRUE), 4), share_funnel_gt1 = round(mean(funnel_ratio > 1, na.rm = TRUE), 3),
                    median_funnel = round(median(funnel_ratio, na.rm = TRUE), 3)), by = .(dose_def, family)]
  write_aggregate(diag, "15_diagnostics_summary.csv", exempt = c("n_items", "n_T1_tests", "n_T2_tests", "n_T1_reject", "n_T2_reject"))
  ## mass-domination diagnostic: eligible items (two to nine observed categories), finite exceedances, exceedances
  ## in supported categories, positive/zero flags, and the finite maxima. Item counts, not person counts.
  md <- tests[family == "A_substantive" & !is.na(funnel_ncat)]
  mx <- function(x) if (any(!is.na(x))) round(max(x, na.rm = TRUE), 3) else NA_real_
  mass_sum <- md[, .(n_items_eligible = .N,
                     n_finite_gt1 = sum(funnel_ratio > 1, na.rm = TRUE),
                     max_finite_ratio = mx(funnel_ratio),
                     n_supported_gt1 = sum(funnel_ratio_supported > 1, na.rm = TRUE),
                     max_supported_ratio = mx(funnel_ratio_supported),
                     n_sparse_only_gt1 = sum(funnel_ratio > 1 & (is.na(funnel_ratio_supported) | funnel_ratio_supported <= 1), na.rm = TRUE),
                     n_zero_denom_items = sum(funnel_zero_denom > 0),
                     n_zero_denom_categories = sum(funnel_zero_denom),
                     n_flagged_any = sum(funnel_ratio > 1 | funnel_zero_denom > 0, na.rm = TRUE),
                     ## flags that survive the allocation of the fresh missing mass (identity map infeasible even then)
                     n_flagged_not_reconcilable = sum((funnel_ratio > 1 | funnel_zero_denom > 0) & funnel_identity_feasible %in% FALSE, na.rm = TRUE),
                     n_finite_gt1_not_reconcilable = sum(funnel_ratio > 1 & funnel_identity_feasible %in% FALSE, na.rm = TRUE),
                     n_zero_denom_not_reconcilable = sum(funnel_zero_denom > 0 & funnel_identity_feasible %in% FALSE, na.rm = TRUE),
                     n_supported_gt1_not_reconcilable = sum(funnel_ratio_supported > 1 & funnel_identity_feasible %in% FALSE, na.rm = TRUE),
                     median_fresh_missing_mass = round(median(funnel_missing_mass, na.rm = TRUE), 4),
                     min_fresh = MIN_CELL), by = dose_def]
  write_aggregate(mass_sum, "15_mass_diagnostic_summary.csv",
                  exempt = grep("^n_", names(mass_sum), value = TRUE))
  write_aggregate(md[funnel_ratio > 1 | funnel_zero_denom > 0,
                     .(dose_def, var, label, funnel_ncat, funnel_p = round(funnel_p, 4), funnel_ratio = round(funnel_ratio, 3),
                       funnel_ratio_supported = round(funnel_ratio_supported, 3), funnel_zero_denom, funnel_sparse_gt1,
                       funnel_missing_mass = round(funnel_missing_mass, 4), funnel_needed_mass = round(funnel_needed_mass, 4), funnel_identity_feasible,
                       reach_new = round(reach_new, 4), reach_old_S = round(reach_old_S, 4),
                       answer_new_given_reach = round(answer_new_given_reach, 4), answer_old_S_given_reach = round(answer_old_S_given_reach, 4),
                       nr_routing_flag)][
                       order(dose_def, -funnel_zero_denom, -funnel_ratio)],
                  "15_mass_diagnostic_flags.csv")
  write_aggregate(arms, "15_arms.csv", exempt = c("B", "seed", "n_ec_items"))
  write_aggregate(rbindlist(val_all), "15_validation.csv", exempt = c("n_compared", "n_na_mismatch"))
  ## 中核 32 項目(J3 の j3_items.csv の w5_cont を var に対応)
  if (file.exists(CORE32)) {
    it <- fread(CORE32, encoding = "UTF-8")
    core <- it[nzchar(w5_cont), .(key, label_en, class, w5_var = toupper(trimws(w5_cont)))]
    core <- core[!grepl("\\+", w5_var)]
    hit <- res[toupper(var) %in% core$w5_var & family == "A_substantive"]
    hit[, key := core$key[match(toupper(var), core$w5_var)]]
    cat("diag: 中核項目の一致:", uniqueN(hit$var), "/", nrow(core), "\n")
    write_aggregate(hit[, .(dose_def, key, var, label, estimator, estimate = round(estimate, 5), se_boot = round(se_boot, 5), d_std = round(d_std, 4), q = signif(q, 4), class3, n_old, n_new)], "15_core32.csv")
  } else cat("diag: j3_items.csv not found (", CORE32, ") — 15_core32.csv は省略\n")
  ## 環境
  si <- capture.output(sessionInfo())
  writeLines(c(paste("time:", format(Sys.time())), paste("input:", basename(f)), paste("input size:", file.size(f)), paste("input mtime:", format(file.mtime(f))),
               paste("input md5:", fp$md5), paste("input sha256:", fp$sha256),
               paste("derived:", file.path(DERIVED_DIR, "analysis_w5.rds"), format(file.mtime(file.path(DERIVED_DIR, "analysis_w5.rds")))),
               paste("varmap:", VARMAP), paste("style items:", basename(sid$path), "md5:", sid$md5), paste("B:", B_BOOT, "seed:", SEED15, "dose:", paste(DOSE_DEFS, collapse = ",")),
               paste("panelcond:", as.character(packageVersion("panelcond"))), "", si), file.path(RESULTS_DIR, "15_env.txt"))
  cat("\n== 15 完了 (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ). 共有してほしいもの ==\n")
  cat("コンソール出力全文 + results/15_*.csv + 15_env.txt(すべて集計値)\n")
  print(det[estimator %in% c("naive", "sm", "ssm", "ec", "ec_adj") & family == "A_substantive"])
}

run_guarded("15_panelcond_designs", main)

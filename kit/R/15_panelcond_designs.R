# ============================================================
# P1 JLPS kit — 15: 補充標本設計の全項目適用(panelcond: naive / SM / SSM / EC / EC-adj + 診断)
# 目的(構成 I の論文 E、ウォーターフォール D4 SSM / D5 EC の素材):
#   04 の派生(analysis_w5.rds: w5 リスク集合 × 全項目ユニバース)と統合マスタ(.dta)から、
#   継続コホート(CN=1、全 4,800 人)と 2011 追加コホート(CN=2)を panelcond の入力形に組み、
#   項目ごとに naive / SM / SSM / EC / EC-adj、入所波選抜差分 δ、診断 T1(SM−EC: 非定常性)・
#   T2(SSM−SM: 状態依存)、漏斗比 sup_y p·q(y)/f2(y)(TP2 §6 の項目別診断)を推定する。
#   A 族(実質値)、B 族(項目無回答・DK)、人レベル様式指標(欠測率・DK 率・中点・極端)に同じ枠を適用。
#   推論: (i) panelcond::pc_estimate の解析 SE、(ii) 項目横断の**同時人ブートストラップ**(B 回、
#   コホート内で人を再抽出、全推定量と T1/T2 差を反復内で再計算)。(ii) はベクトル化した自前エンジンで、
#   フルデータ上の点推定が pc_estimate と一致することを 15_validation.csv で検証する。
# 用量の定義: DOSE_DEF = "exact"(継続は w1–w5 全波回答 = 用量 4 ちょうど; J3 の規約)
#              または "any"(w5 回答者全員 = 04 のリスク集合; 用量 ≤ 4 の混合)。既定は両方。
# 入所波(w1)の対応変数: jlps_docs/varmap_w1-19.csv の w1/w5 列 + 15_entry_overrides.csv(手動)。
#   値ラベル集合(特殊コードを除く)が w1 と w5 で一致する項目だけ EC を計算(不一致は 15_ec_unavailable.csv に列挙)。
# 出力(すべて集計値、N<10 抑制): 15_designs_items.csv / 15_tests_items.csv / 15_detection_counts.csv /
#   15_diagnostics_summary.csv / 15_mass_diagnostic_summary.csv / 15_mass_diagnostic_flags.csv / 15_core32.csv /
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
  mid_cnt <- ext_cnt <- mid_n <- ext_n <- rep(0, n_old); dk_cnt <- dk_n <- miss_cnt <- miss_n <- rep(0, n_old)
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
    k <- length(vals1)
    if (k %in% 4:7 && all(diff(vals1) == 1)) {
      ok <- !is.na(subst)
      if (k %% 2 == 1) { mid_cnt[ok] <- mid_cnt[ok] + (subst[ok] == vals1[(k + 1) / 2]); mid_n[ok] <- mid_n[ok] + 1 }
      ext_cnt[ok] <- ext_cnt[ok] + (subst[ok] %in% c(vals1[1], vals1[k])); ext_n[ok] <- ext_n[ok] + 1
    }
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
  ## 人レベル様式指標(w5 は 04 の person_dq; w1 は上で計算)
  P5 <- as.matrix(pdq[, .(pdq_miss_share = miss_share, pdq_dk_share = dk_share, pdq_mid_share = mid_share, pdq_ext_share = ext_share)])
  P1 <- cbind(pdq_miss_share = ifelse(miss_n > 0, miss_cnt / miss_n, NA), pdq_dk_share = ifelse(dk_n > 0, dk_cnt / dk_n, NA),
              pdq_mid_share = ifelse(mid_n > 0, mid_cnt / mid_n, NA), pdq_ext_share = ifelse(ext_n > 0, ext_cnt / ext_n, NA))
  fam <- c(rep("A_substantive", J), rep("B_itemnonresp", J), rep("B_dk", J), rep("P_style", 4))
  cols <- c(vars, paste0(vars, "__M"), paste0(vars, "__DK"), colnames(P5))
  Yo <- cbind(fill_old(Y), fill_old(M), fill_old(DK), fill_old(P5)); colnames(Yo) <- cols
  Yn <- cbind(fill_new(Y), fill_new(M), fill_new(DK), fill_new(P5)); colnames(Yn) <- cols
  Eo <- cbind(E_A, E_M, E_DK, P1); colnames(Eo) <- cols
  ec_ok_col <- c(ec_status$ec_ok, ec_status$ec_ok, ec_status$ec_ok, rep(TRUE, 4))
  Eo[, !ec_ok_col] <- NA
  binary_col <- c(meta$scale_type == "binary", rep(TRUE, 2 * J), rep(FALSE, 4))
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
  res_all <- list(); tests_all <- list(); val_all <- list(); arms_all <- list()
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
    elig_new <- colMeans(!is.na(Yn[, J + seq_len(J), drop = FALSE]))
    funnel <- funnel_sup <- p_item <- rep(NA_real_, JJ); funnel_zd <- funnel_sp <- funnel_ncat <- rep(NA_integer_, JJ)
    for (j in seq_len(J)) {                                          # A 族の離散項目のみ
      if (!(elig_new[j] > 0)) next
      p_item[j] <- full$p_surv[j] / elig_new[j]
      mr <- mass_ratio_item(Yo[IOf[, j], j], Yn[!is.na(Yn[, j]), j], p_item[j], max_cat = 9L, min_fresh = MIN_CELL)
      if (!mr$eligible) next
      funnel[j] <- mr$ratio_finite; funnel_sup[j] <- mr$ratio_supported
      funnel_zd[j] <- mr$n_zero_denom; funnel_sp[j] <- mr$n_sparse_gt1; funnel_ncat[j] <- mr$ncat
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
    ## clock-time components, duplicate recodes) and those whose routing differs between the cohorts' questionnaires
    ## (04's nr_routing_flag: the continuing cohort answers only if newly married, so the arms are different subgroups)
    ## a follow-up component of a routing-flagged question (the same variable name plus a component suffix, e.g. DQ46Y,
    ## years of premarital cohabitation, after DQ46) inherits the flag: its non-routed respondents are coded missing
    ## or not applicable rather than with the no-answer code, so 04's rate rule does not see it (since 2026-09-24)
    rf <- unique(toupper(pc[nr_routing_flag %in% TRUE, var]))
    fu <- function(v) { if (!length(rf)) return(FALSE); v <- toupper(v); any(startsWith(v, rf) & grepl("^[A-Z_][A-Z0-9_]*$", substring(v, nchar(rf) + 1L)) & nchar(v) > nchar(rf)) }
    pc[, routing_followup := !(nr_routing_flag %in% TRUE) & vapply(var, fu, logical(1))]
    if (any(pc$routing_followup)) cat("diag: follow-ups of routing-flagged questions (outside the counts):", paste(unique(pc[routing_followup == TRUE, var]), collapse = ", "), "\n")
    pc[, count_exclude := exclude_flag | (nr_routing_flag %in% TRUE) | routing_followup]
    pc[, `:=`(q = NA_real_, q_tost = NA_real_)]
    pc[count_exclude == FALSE, q := p.adjust(p, "BH"), by = .(family, estimator)]
    pc[count_exclude == FALSE, q_tost := p.adjust(p_tost, "BH"), by = .(family, estimator)]
    pc[, class3 := fifelse(is.na(se) | se <= 0, "n/a (no variation)", fifelse(q < Q_CUT, "affected", fifelse(q_tost < Q_CUT, "equivalent", "undetermined")))]
    pc[exclude_flag == TRUE & !(is.na(se) | se <= 0), class3 := "excluded (listed code)"]
    pc[exclude_flag == FALSE & count_exclude == TRUE & !(is.na(se) | se <= 0), class3 := "excluded (routing differs)"]
    pc[, label := meta$label[match(var, meta$var)]]
    res_all[[dose_def]] <- pc
    ## 診断表(項目ごと; ブートストラップ SE 版の T1/T2)
    te <- unique(pc[, .(dose_def, col, family, var, label, T1_stat, T1_p, T2_stat, T2_p, delta_entry, n_pairs, p_survive, funnel_ratio, funnel_ratio_supported, funnel_p, funnel_zero_denom, funnel_sparse_gt1, funnel_ncat, exclude_flag, count_exclude)])
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
                            funnel_zero_denom, funnel_sparse_gt1, funnel_ncat, exclude_flag, count_exclude)], "15_tests_items.csv")
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
                     min_fresh = MIN_CELL), by = dose_def]
  write_aggregate(mass_sum, "15_mass_diagnostic_summary.csv",
                  exempt = grep("^n_", names(mass_sum), value = TRUE))
  write_aggregate(md[funnel_ratio > 1 | funnel_zero_denom > 0,
                     .(dose_def, var, label, funnel_ncat, funnel_p = round(funnel_p, 4), funnel_ratio = round(funnel_ratio, 3),
                       funnel_ratio_supported = round(funnel_ratio_supported, 3), funnel_zero_denom, funnel_sparse_gt1)][
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
               paste("varmap:", VARMAP), paste("B:", B_BOOT, "seed:", SEED15, "dose:", paste(DOSE_DEFS, collapse = ",")),
               paste("panelcond:", as.character(packageVersion("panelcond"))), "", si), file.path(RESULTS_DIR, "15_env.txt"))
  cat("\n== 15 完了 (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ). 共有してほしいもの ==\n")
  cat("コンソール出力全文 + results/15_*.csv + 15_env.txt(すべて集計値)\n")
  print(det[estimator %in% c("naive", "sm", "ssm", "ec", "ec_adj") & family == "A_substantive"])
}

run_guarded("15_panelcond_designs", main)

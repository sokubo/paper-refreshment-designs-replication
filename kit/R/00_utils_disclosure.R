# ============================================================
# P1 JLPS kit — 開示ユーティリティ(N<10セル抑制)
# すべての results/ への書き出しは write_aggregate() を通すこと。
# ============================================================

# 度数系の列(n, n_で始まる, count系)の小セルを抑制
# v2修正(2026-08-24): 旧正規表現 "^(n|n_|N_|count|freq)" は "n" で始まる全列
#   (nonmiss_rate 等の率列も!)に誤マッチし、02_covariate_coverage.csv の
#   sex/ybirth 行を誤抑制していた。度数列を「n / N 単独、または n_/N_/count/freq
#   で始まる列」に限定。率・割合列(*_rate, *_share, prop_*)は抑制対象外
#   (率の開示は分母セルの抑制で守る設計)。
# exempt: 人数でないメタデータ度数列(例: n_waves=波数, n_vars=変数数,
#   n_value_labels=値ラベル数)を明示的に抑制対象から外す。既定は空(安全側)。
suppress_counts <- function(df, min_cell = MIN_CELL, exempt = character(0)) {
  cnt_cols <- grep("^(n_|N_|count|freq)", names(df), value = TRUE)
  cnt_cols <- union(cnt_cols, names(df)[names(df) %in% c("n", "N", "count", "freq")])
  cnt_cols <- setdiff(cnt_cols, exempt)
  for (cc in cnt_cols) {
    if (is.numeric(df[[cc]])) {
      small <- !is.na(df[[cc]]) & df[[cc]] > 0 & df[[cc]] < min_cell
      if (any(small)) {
        df[[cc]][small] <- NA
        df$suppressed <- if (is.null(df$suppressed)) small else (df$suppressed | small)
      }
    }
  }
  df
}

# 集計値の唯一の書き出し口
write_aggregate <- function(df, filename, min_cell = MIN_CELL, exempt = character(0)) {
  stopifnot(is.data.frame(df))
  df <- suppress_counts(df, min_cell, exempt = exempt)
  # 個票らしき列が紛れていないかの防波堤
  banned <- c("PanelID", "panelid", "pid", "zip", "zip7", "zipcode")
  hit <- intersect(tolower(names(df)), tolower(banned))
  if (length(hit) > 0) stop(sprintf("write_aggregate: 個票識別子らしき列があります: %s", paste(hit, collapse = ",")), call. = FALSE)
  path <- file.path(RESULTS_DIR, filename)
  write.csv(df, path, row.names = FALSE)
  cat("wrote:", path, "(", nrow(df), "rows )\n")
  invisible(path)
}

# 度数表(小セル抑制つき)
safe_freq <- function(x, varname = "value", min_cell = MIN_CELL) {
  tb <- as.data.frame(table(x, useNA = "ifany"), stringsAsFactors = FALSE)
  names(tb) <- c(varname, "n")
  suppress_counts(tb, min_cell)
}

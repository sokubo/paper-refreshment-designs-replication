# ============================================================
# P1 JLPS kit — 15b: 項目の衛生フラグ(04 のメタ)を集計として書き出す(数秒)
#   15 の結果(特に B 族=項目無回答)を解釈するには、調査票の版差・ルーティング疑いの
#   フラグ(04 の filter_mismatch / nr_routing_flag)が要る。15 v1 の出力には無かったので別出し。
# 実行: cd <P1ルート> && Rscript analysis/R/15b_item_flags.R  → results/15_item_flags.csv
# ============================================================
.here <- local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) dirname(normalizePath(a[1])) else file.path("analysis", "R") })  # this script's folder (kit: R/)
source(file.path(.here, "00_config.R"))
source(file.path(.here, "00_utils_disclosure.R"))
suppressMessages(library(data.table))
main <- function() {
  a <- readRDS(file.path(DERIVED_DIR, "analysis_w5.rds")); meta <- as.data.table(a$meta)
  out <- meta[, .(var, label, scale_type, k_labels, n_dk_codes, n_ref_codes, filter_mismatch, nr_routing_flag,
                  base_rate_T1 = round(base_rate_T1, 3), base_rate_T0 = round(base_rate_T0, 3))]
  cat("diag: universe", nrow(out), "; filter_mismatch", sum(out$filter_mismatch, na.rm = TRUE), "; nr_routing_flag", sum(out$nr_routing_flag, na.rm = TRUE), "\n")
  write_aggregate(out, "15_item_flags.csv", exempt = c("k_labels", "n_dk_codes", "n_ref_codes"))
}
run_guarded("15b_item_flags", main)

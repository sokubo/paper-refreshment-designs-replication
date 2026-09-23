# ============================================================
# Common settings of the JLPS pipeline (public kit version).
# The licensed microdata are read from a directory OUTSIDE the project tree, named by the environment
# variable P1_DATA_DIR; derived person-level files are written to <P1_DATA_DIR>/derived_me/, and only
# aggregates (cells below MIN_CELL suppressed) are written to P1_RESULTS_DIR (default <P1_PROJECT_DIR>/results).
#   P1_DATA_DIR=<data root> P1_PROJECT_DIR=<this kit directory> Rscript R/15_panelcond_designs.R
# Both variables are required: this version has no default location (see KIT_README.md).
# ============================================================
.need_env <- function(v) {
  x <- Sys.getenv(v, "")
  if (!nzchar(x)) stop(sprintf("set the environment variable %s (see KIT_README.md)", v), call. = FALSE)
  path.expand(x)
}
DATA_DIR    <- .need_env("P1_DATA_DIR")
RAW_DIR     <- file.path(DATA_DIR, "raw")
ZIP_DIR     <- file.path(DATA_DIR, "zipcode")
DERIVED_DIR <- file.path(DATA_DIR, "derived_me")   # person-level derived files: outside the project tree

# project side: the only exit for outputs (aggregates only)
PROJECT_DIR <- .need_env("P1_PROJECT_DIR")
RESULTS_DIR <- path.expand(Sys.getenv("P1_RESULTS_DIR", file.path(PROJECT_DIR, "results")))

## --- the licensed input: one explicitly named file for every script ------------------------------------
## P1_INPUT_DTA names the integrated master file (a path, or a file name inside RAW_DIR). If it is unset, the
## single master .dta in RAW_DIR is used (files whose names contain "online" are the 2020 web supplement and are
## ignored); the run stops if there is more than one candidate instead of choosing by a naming heuristic.
input_dta <- function() {
  f <- Sys.getenv("P1_INPUT_DTA", "")
  if (nzchar(f)) {
    if (!file.exists(f)) f <- file.path(RAW_DIR, f)
    if (!file.exists(f)) stop("P1_INPUT_DTA does not name an existing file: ", Sys.getenv("P1_INPUT_DTA"), call. = FALSE)
    return(normalizePath(f))
  }
  f_all <- list.files(RAW_DIR, pattern = "\\.dta$", full.names = TRUE)
  f_all <- f_all[!grepl("online", basename(f_all), ignore.case = TRUE)]
  if (length(f_all) != 1)
    stop(sprintf("RAW_DIR holds %d candidate master .dta files (%s); set P1_INPUT_DTA to the one to analyse",
                 length(f_all), paste(basename(f_all), collapse = ", ")), call. = FALSE)
  normalizePath(f_all)
}
## fingerprint of the input file: size, MD5 (base R) and SHA-256 when a system tool is available
input_fingerprint <- function(f) {
  sha <- tryCatch({
    out <- suppressWarnings(system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE, stderr = FALSE))
    if (!length(out)) out <- suppressWarnings(system2("sha256sum", shQuote(f), stdout = TRUE, stderr = FALSE))
    if (length(out)) sub("\\s.*$", "", out[1]) else NA_character_
  }, error = function(e) NA_character_)
  list(file = basename(f), bytes = file.size(f), md5 = unname(tools::md5sum(f)), sha256 = sha)
}
## respondent identifier column of the master file, if present (used only to verify alignment across scripts)
ID_COLS <- c("PanelID", "PANELID", "ID", "CASEID")

MIN_CELL  <- 10          # 秘匿閾値(N<10セルは抑制)
# 尺度に付随する「該当者なし」コードを非該当として扱う項目別上書き(04 と 15 で共有; 2026-09-16)
#   満足度 DQ18A-E の 6(仕事/結婚/友人/親/子がいない)、DQ04_1C の 5(部下はいない)、DQ04_2 の 5(上司・同僚はいない)。
NAP_OVERRIDE <- list(DQ18A = 6, DQ18B = 6, DQ18C = 6, DQ18D = 6, DQ18E = 6, DQ04_1C = 5, DQ04_2 = 5)
SEED      <- 20260824
REFRESH_WAVE_HYP <- 5L   # リフレッシュ標本投入波の作業仮説(2011=w5)。01-02で実証確認する
MAX_WAVE  <- 19L

## --- ガード -----------------------------------------------------------------
stopifnot_dir <- function(p, what) {
  if (!dir.exists(p)) stop(sprintf("%s が見つかりません: %s", what, p), call. = FALSE)
}
# DATA_DIR がプロジェクト内を指していたら停止(個票混入ガード)
if (grepl(normalizePath(PROJECT_DIR, mustWork = FALSE), normalizePath(DATA_DIR, mustWork = FALSE), fixed = TRUE))
  stop("P1_DATA_DIR resolves inside the project tree; restricted data must stay outside it.", call. = FALSE)
stopifnot_dir(RAW_DIR, "RAW_DIR(個票フォルダ)")
dir.create(DERIVED_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## --- 共通ランナー: main()+tryCatch+エラーファイル ---------------------------
run_guarded <- function(script_name, main_fn) {
  err_file <- file.path(RESULTS_DIR, paste0(script_name, "_ERROR.txt"))
  if (file.exists(err_file)) file.remove(err_file)
  tryCatch(main_fn(), error = function(e) {
    writeLines(c(paste("time:", format(Sys.time())),
                 paste("script:", script_name),
                 paste("error:", conditionMessage(e)),
                 capture.output(print(sys.calls()))), err_file)
    stop(e)
  })
}

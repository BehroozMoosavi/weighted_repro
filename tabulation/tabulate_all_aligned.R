## ===========================================================================
##  tabulate.R
##
##  Produces every table and figure of the paper from the simulation output
##  under codes/.  Base R only; no packages.
##
##  Usage
##  -----
##      setwd("<project root>")          # the folder that will hold results/
##      source("codes/tabulation/tabulate.R")
##
##  or, to load the functions without writing anything,
##
##      options(tabulate.autorun = FALSE)
##      source("codes/tabulation/tabulate.R")
##      tab_init()                       # resolve the repository root
##      tab_exp1_main()                  # then call whichever you want
##
##  Output
##  ------
##      results/tables/*.tex     booktabs fragments and figure wrappers,
##                               each a complete float with caption and label
##      results/figures/*.pdf    vector figures
##
##  One float per distinct set of numbers: four tables and seven figures, so
##  no table restates a matrix that a figure already displays.  Section 8
##  lists the writers that are defined but deliberately not registered.
##
##  The figure wrappers emit \includegraphics{results/figures/...}, so the
##  paths resolve on Overleaf when the main document sits beside results/.
##
##  Where the numbers come from
##  ---------------------------
##  Every value for the proposed procedure is read from a results file under
##  codes/; each reader in Section 4 names its inputs and checks the columns
##  it needs, stopping with the names actually found rather than indexing
##  positionally.  Section 3 holds the comparison values of Awan and Wang
##  (2025) and of Wang, Chang and Awan, which exist only in those papers and
##  are therefore transcribed, with the source recorded beside each object.
##  A missing input is reported as [skip] and the corresponding table or
##  figure is omitted, so a partial checkout still produces what it can.
##
##  Figure conventions
##  ------------------
##  Grey-scale heat maps with the value printed in each cell and a colour
##  break at the nominal level, so the threshold is legible without the key.
##  Panel titles are the method name alone; attribution is made in the
##  caption.  The shared key is a horizontal strip of fixed height beneath
##  the panels, and axis titles are drawn only on the outer panels of a grid.
##
##  Contents
##  --------
##      1  repository root and result directories
##      2  reference levels, formatting, LaTeX writers
##      3  comparison values from the literature
##      4  readers
##      5  figure geometry and line plots
##      6  table writers
##      7  figure writers
##      8  driver
## ===========================================================================


## ===========================================================================
## 1.  REPOSITORY ROOT AND RESULT DIRECTORIES
## ===========================================================================
## ---------------------------------------------------------------------------
## 1.  Repository root
##
##  Walks up from the working directory until a folder containing "codes" is
##  found, so the scripts run from any depth.  Override with
##      options(repro.root = "/path/to/asymptotic_optimal_repro")
##  or the environment variable REPRO_ROOT.
## ---------------------------------------------------------------------------

repo_root <- function() {
  
  opt <- getOption("repro.root", default = NULL)
  if (!is.null(opt)) {
    if (!dir.exists(file.path(opt, "codes")))
      stop("options(repro.root=) is set but ", opt, "/codes does not exist")
    return(normalizePath(opt))
  }
  
  env <- Sys.getenv("REPRO_ROOT", unset = "")
  if (nzchar(env)) {
    if (!dir.exists(file.path(env, "codes")))
      stop("REPRO_ROOT is set but ", env, "/codes does not exist")
    return(normalizePath(env))
  }
  
  here <- normalizePath(getwd())
  for (i in 1:8) {
    if (dir.exists(file.path(here, "codes"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) break
    here <- parent
  }
  stop("could not locate the repository root: no ancestor of '", getwd(),
       "' contains a 'codes' directory. Set options(repro.root=).")
}

ROOT  <- repo_root()
CODES <- file.path(ROOT, "codes")

## result directories, as laid out in the repository map
dir_list <- function(CODES) list(
  exp1       = file.path(CODES, "exp1_location_scale_normal"),
  exp1_cov   = file.path(CODES, "exp1_location_scale_normal", "results_coverage"),
  exp1_apply = file.path(CODES, "exp1_location_scale_normal", "results_apply"),
  
  exp2       = file.path(CODES, "exp2_linear_regression"),
  exp2_eff   = file.path(CODES, "exp2_linear_regression", "results_eff"),
  exp2_sigma = file.path(CODES, "exp2_linear_regression",
                         "linear_sigma_penalized_efficient"),
  
  exp3       = file.path(CODES, "exp3_objective_perturbation"),
  exp3_res   = file.path(CODES, "exp3_objective_perturbation", "results"),
  
  sens       = file.path(CODES, "sensitivity_analysis"),
  
  ## present in some checkouts only; every use is guarded
  pbadi      = file.path(CODES, "PB_ADI_Failure")
)

## ---------------------------------------------------------------------------
## Initialisation
##
##  Sourcing this file only defines things.  tab_init() resolves the
##  repository root, builds the directory list and creates the output
##  folders; every writer calls it if it has not run yet.  Pass a path, or
##  set options(repro.root=) or REPRO_ROOT, to override the search.
## ---------------------------------------------------------------------------

DIR_TAB <- "results/tables"
DIR_FIG <- "results/figures"

tab_init <- function(root = NULL) {
  
  if (!is.null(root))
    options(repro.root = normalizePath(root, mustWork = TRUE))
  
  ROOT  <<- repo_root()
  CODES <<- file.path(ROOT, "codes")
  DIR   <<- dir_list(CODES)
  
  dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)
  dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
  
  message("repository root : ", ROOT)
  message("tables to       : ", normalizePath(DIR_TAB, mustWork = FALSE))
  message("figures to      : ", normalizePath(DIR_FIG, mustWork = FALSE))
  invisible(ROOT)
}

tab_need_init <- function()
  if (!exists("ROOT", envir = globalenv(), inherits = FALSE)) tab_init()


## ===========================================================================
## 2.  REFERENCE LEVELS, FORMATTING, LATEX WRITERS
## ===========================================================================

## The exact rank cut-off is a = floor(alpha (R+1)) + 1; the size bound is
## (a-1)/(R+1) and the coverage bound its complement.  The bound is an upper
## bound on size and a lower bound on coverage, attained when ties in the
## depth have probability zero.  These are computed, never written in, so a
## caption cannot quote the wrong level.

rank_cutoff <- function(R, alpha) floor(alpha * (R + 1)) + 1
size_bound  <- function(R, alpha) (rank_cutoff(R, alpha) - 1) / (R + 1)
cov_bound   <- function(R, alpha) 1 - size_bound(R, alpha)
mc_se       <- function(p, n) sqrt(p * (1 - p) / n)

ALPHA_12 <- 0.05     # Experiments 1 and 2, and their sweeps
ALPHA_3  <- 0.10     # Experiment 3 runs at 90 per cent confidence
R_SYN    <- 200

BND_12 <- cov_bound(R_SYN, ALPHA_12)     # 0.95025
SZB_12 <- size_bound(R_SYN, ALPHA_12)    # 0.04975
BND_3  <- cov_bound(R_SYN, ALPHA_3)      # 0.90050

E3_BOX  <- c(-10, 10)                    # search box for beta1
E3_BOXW <- diff(E3_BOX)
E3_SAT  <- 0.90 * E3_BOXW                # width above which an interval is
# box-limited and uninformative
## ---------------------------------------------------------------------------
## 4.  Input helpers
## ---------------------------------------------------------------------------

## TRUE with a note if absent, so a caller can skip a whole section
have <- function(path, what = NULL) {
  ok <- file.exists(path)
  if (!ok) {
    msg <- if (is.null(what)) "" else paste0(" (", what, ")")
    message("  [skip] not found", msg, ": ", sub(ROOT, "<root>", path, fixed = TRUE))
  }
  ok
}

## read.csv with the conventions these results files use: R's write.csv
## quotes headers and character cells, and several headers repeat, so names
## are kept verbatim and duplicates are left for positional access.
read_results <- function(path, row.names = NULL, check.names = FALSE) {
  if (!file.exists(path)) stop("missing results file: ", path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE,
                       check.names = check.names, row.names = row.names)
  if (nrow(x) == 0L) warning("results file has no rows: ", path)
  x
}

## first match of a glob inside a directory, or NA
find_one <- function(dir, pattern) {
  if (!dir.exists(dir)) return(NA_character_)
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) return(NA_character_)
  if (length(f) > 1L)
    message("  [note] ", length(f), " files match '", pattern, "' in ",
            basename(dir), "; using ", basename(f[1]))
  f[1]
}

## pull n / epsilon / nSIM out of a file name such as
## coverage_study_n100_ep1_nSIM1000.csv
parse_stem <- function(path) {
  b <- basename(path)
  g <- function(re) {
    m <- regmatches(b, regexpr(re, b))
    if (length(m) == 0L) return(NA_real_)
    as.numeric(sub(re, "\\1", m))
  }
  list(n    = g("n([0-9]+)"),
       ep   = g("ep([0-9.]+)"),
       nSIM = g("nSIM([0-9]+)"))
}

## strip the quotes write.csv leaves inside character cells, if any survive
unquote <- function(x) gsub('^"|"$', "", trimws(as.character(x)))
## ---------------------------------------------------------------------------
## Number formatting
## ---------------------------------------------------------------------------

fnum <- function(x, d = 3, dash = "--")
  ifelse(is.na(x) | !is.finite(x), dash, formatC(x, format = "f", digits = d))

## "0.991 (0.003)"
fse <- function(x, se, d = 3)
  ifelse(is.na(x), "--", paste0(fnum(x, d), "\\,(", fnum(se, d), ")"))

## exact 0 and 1 without decimals, as in the published grids
fcell <- function(x, d = 3)
  ifelse(is.na(x), "--", ifelse(x == 1, "1", ifelse(x == 0, "0", fnum(x, d))))

fpct <- function(x, d = 1)
  ifelse(is.na(x), "--", paste0(fnum(x, d), "\\,\\%"))
## ---------------------------------------------------------------------------
## Grey ramps.  Both are two-segment: light below the nominal level, dark
## above it, so the threshold is visible in the plot itself.  The return value
## is percent black.
## ---------------------------------------------------------------------------

grey_reject <- function(v, alpha = ALPHA_12)
  ifelse(is.na(v), NA_real_,
         ifelse(v <= alpha, 18 * v / alpha,
                18 + 62 * (v - alpha) / (1 - alpha)))

grey_cov <- function(v, nominal = 1 - ALPHA_12)
  ifelse(is.na(v), NA_real_,
         ifelse(v <= nominal, 20 * v / nominal,
                20 + 60 * (v - nominal) / (1 - nominal)))

## percent black -> a grey colour
pct_to_col <- function(g) ifelse(is.na(g), "white", grey(1 - g / 100))
## ---------------------------------------------------------------------------
## LaTeX writers
## ---------------------------------------------------------------------------

## body: character matrix of formatted cells; header: character vector;
## align: e.g. "lcccc"; span/cmid: optional raw lines placed above and below
## the header; note: optional footnote paragraph.
latex_tabular <- function(body, header, align, span = NULL, cmid = NULL,
                          note = NULL, size = "\\small", rowsep = NULL) {
  body <- as.matrix(body)
  if (ncol(body) != length(header))
    stop("body has ", ncol(body), " columns but header has ", length(header))
  if (nchar(align) != length(header))
    stop("align has ", nchar(align), " entries but header has ", length(header))
  
  out <- c(size, paste0("\\begin{tabular}{", align, "}"), "\\toprule")
  if (!is.null(span)) out <- c(out, span)
  out <- c(out, paste0(paste(header, collapse = " & "), " \\\\"))
  if (!is.null(cmid)) out <- c(out, cmid)
  out <- c(out, "\\midrule")
  for (i in seq_len(nrow(body))) {
    ## a row separator has to follow the row terminator, never sit in a cell
    term <- if (!is.null(rowsep) && i %in% rowsep) " \\\\[2pt]" else " \\\\"
    out <- c(out, paste0(paste(body[i, ], collapse = " & "), term))
  }
  out <- c(out, "\\bottomrule", "\\end{tabular}")
  if (!is.null(note))
    out <- c(out, "", "\\vspace{2pt}",
             "\\begin{minipage}{.94\\linewidth}\\footnotesize",
             note, "\\end{minipage}")
  out
}

write_table <- function(inner, caption, label, file, placement = "!htbp") {
  lines <- c(paste0("\\begin{table}[", placement, "]"), "\\centering",
             paste0("\\caption{", caption, "}"),
             paste0("\\label{", label, "}"),
             inner, "\\end{table}")
  path <- file.path(DIR_TAB, file)
  writeLines(lines, path)
  message("  [table] ", path)
  invisible(path)
}

## a figure float that includes a PDF this script produced, so the .tex and
## the .pdf stay in step
write_figure_wrapper <- function(pdfname, caption, label, width = 0.9,
                                 file = NULL, placement = "!htbp") {
  if (is.null(file)) file <- sub("\\.pdf$", ".tex", pdfname)
  lines <- c(paste0("\\begin{figure}[", placement, "]"), "\\centering",
             sprintf("\\includegraphics[width=%.2f\\linewidth]{%s/%s}",
                     width, DIR_FIG, pdfname),
             paste0("\\caption{", caption, "}"),
             paste0("\\label{", label, "}"),
             "\\end{figure}")
  path <- file.path(DIR_TAB, file)
  writeLines(lines, path)
  message("  [figure wrapper] ", path)
  invisible(path)
}


## ===========================================================================
## 3.  COMPARISON VALUES FROM THE LITERATURE
##
##  These do not come from our simulation output and are therefore stored
##  explicitly, with the source of each array recorded below.
##
##  Awan and Wang (2025), Figure 5:
##    Repro linear-regression rejection probabilities.
##  Awan and Wang (2025), Figure 6:
##    Repro and parametric-bootstrap clamping sensitivity.
##  Wang, Chang and Awan (2026), Figure 5:
##    PB-naive and PB-ADI linear-regression rejection probabilities.
##    For PB-ADI we use:
##      PB (Indirect estimator + approximate-pivot)
##  Wang, Chang and Awan (2026), Figure 6:
##    Repro, PB-ADI and PB-naive logistic-regression
##    90% CI coverage and average width.
##
##  A note on the Experiment 3 level.  Wang, Chang and Awan describe the
##  quantities in their Figure 6 as 90% confidence intervals, in the figure
##  and in its caption, notwithstanding an apparent alpha = 0.05 in the
##  surrounding setup text.  We follow the figure, so ALPHA_3 = 0.10 and the
##  applicable coverage bound is 0.90050.
##
##  Grids are aligned to our design by row and column value rather than by
##  position, so a change of design cannot silently misalign a comparison; a
##  cell absent from a source grid is left as NA and rendered as a dash.
## ===========================================================================

LIT_E2_N    <- c(100, 200, 300, 400, 500, 1000, 2000, 5000)
LIT_E2_BETA <- c(0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

lit_grid <- function(v, rows = LIT_E2_BETA, cols = LIT_E2_N, src = NA_character_) {
  m <- matrix(v, nrow = length(rows), ncol = length(cols), byrow = TRUE)
  dimnames(m) <- list(as.character(rows), as.character(cols))
  attr(m, "source") <- src
  m
}

## [AW] Figure 5, panel "Repro Sample"
LIT_E2_REPRO <- lit_grid(c(
  0.000, 0.000, 0.001, 0.000, 0.000, 0.001, 0.000, 0.002,
  0.000, 0.012, 0.023, 0.051, 0.130, 0.718, 1.000, 1.000,
  0.007, 0.086, 0.278, 0.562, 0.836, 1.000, 1.000, 1.000,
  0.074, 0.666, 0.984, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.281, 0.978, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.566, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.762, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
  src = "Awan and Wang (2025), Fig. 5, Repro")

## [WCA] Figure 5, panel "PB (Indirect estimator + approximate-pivot)"
LIT_E2_ADI <- lit_grid(c(
  0.060, 0.053, 0.052, 0.033, 0.051, 0.043, 0.042, 0.042,
  0.073, 0.185, 0.282, 0.438, 0.592, 0.981, 1.000, 1.000,
  0.151, 0.508, 0.793, 0.942, 0.990, 1.000, 1.000, 1.000,
  0.492, 0.967, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.779, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.924, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.966, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
  src = "Wang, Chang and Awan (2026), Fig. 5, PB (indirect + approximate pivot)")

## [WCA] Figure 5, panel "PB (Naive estimator + F-statistic)"
LIT_E2_NAIVE <- lit_grid(c(
  0.002, 0.020, 0.035, 0.035, 0.048, 0.040, 0.077, 0.117,
  0.003, 0.106, 0.194, 0.329, 0.463, 0.925, 1.000, 1.000,
  0.018, 0.390, 0.734, 0.904, 0.988, 1.000, 1.000, 1.000,
  0.164, 0.913, 0.997, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.468, 0.897, 0.984, 0.997, 1.000, 1.000, 1.000, 1.000,
  0.594, 0.843, 0.951, 0.988, 0.993, 1.000, 1.000, 1.000,
  0.630, 0.780, 0.898, 0.961, 0.982, 1.000, 1.000, 1.000),
  src = "Wang, Chang and Awan (2026), Fig. 5, PB-naive")

## ---------------------------------------------------------------------------
## Experiment 2, clamping sweep: [AW] Figure 6.
## Rows = clamping range Delta, columns = sample size.
## ---------------------------------------------------------------------------

LIT_DELTA <- c(0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0)

lit_clamp <- function(v, src) lit_grid(v, rows = LIT_DELTA, cols = LIT_E2_N, src = src)

## [AW] Fig. 6, "Repro Sample (beta1 = 0)"  -- type I error
LIT_CLAMP_REPRO_H0 <- lit_clamp(c(
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.001,
  0.000, 0.000, 0.000, 0.000, 0.001, 0.001, 0.000, 0.002,
  0.000, 0.001, 0.000, 0.000, 0.001, 0.001, 0.002, 0.002,
  0.000, 0.000, 0.002, 0.000, 0.001, 0.001, 0.001, 0.002,
  0.000, 0.000, 0.001, 0.000, 0.000, 0.001, 0.000, 0.002,
  0.001, 0.000, 0.000, 0.000, 0.001, 0.001, 0.001, 0.000,
  0.003, 0.001, 0.000, 0.000, 0.003, 0.000, 0.001, 0.000),
  "Awan and Wang (2025), Fig. 6, Repro, beta1 = 0")

## [AW] Fig. 6, "Repro Sample (beta1 = 1)"  -- power
LIT_CLAMP_REPRO_H1 <- lit_clamp(c(
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.998, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.967, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.762, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.005, 0.079, 0.238, 0.519, 0.824, 1.000, 1.000, 1.000,
  0.001, 0.001, 0.001, 0.006, 0.018, 0.128, 0.833, 1.000),
  "Awan and Wang (2025), Fig. 6, Repro, beta1 = 1")

## [AW] Fig. 6, "Parametric Bootstrap (beta1 = 0)"  -- type I error
LIT_CLAMP_PB_H0 <- lit_clamp(c(
  0.007, 0.018, 0.019, 0.020, 0.029, 0.015, 0.009, 0.005,
  0.014, 0.037, 0.053, 0.075, 0.089, 0.114, 0.199, 0.357,
  0.017, 0.045, 0.068, 0.107, 0.118, 0.186, 0.361, 0.674,
  0.011, 0.034, 0.047, 0.071, 0.083, 0.130, 0.236, 0.475,
  0.002, 0.019, 0.035, 0.036, 0.047, 0.037, 0.078, 0.107,
  0.000, 0.000, 0.000, 0.004, 0.007, 0.035, 0.048, 0.043,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.004, 0.042),
  "Awan and Wang (2025), Fig. 6, parametric bootstrap, beta1 = 0")

## [AW] Fig. 6, "Parametric Bootstrap (beta1 = 1)"  -- power
LIT_CLAMP_PB_H1 <- lit_clamp(c(
  0.987, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.883, 0.989, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000,
  0.823, 0.957, 0.996, 0.999, 1.000, 1.000, 1.000, 1.000,
  0.718, 0.842, 0.957, 0.987, 0.993, 1.000, 1.000, 1.000,
  0.629, 0.780, 0.898, 0.961, 0.982, 1.000, 1.000, 1.000,
  0.000, 0.030, 0.322, 0.550, 0.672, 0.829, 0.969, 1.000,
  0.000, 0.000, 0.000, 0.000, 0.000, 0.124, 0.699, 0.902),
  "Awan and Wang (2025), Fig. 6, parametric bootstrap, beta1 = 1")

## ---------------------------------------------------------------------------
## Experiment 3: logistic regression, 90% confidence.  [WCA] Figure 6.
## Rows = epsilon, columns = sample size.
## ---------------------------------------------------------------------------

LIT_E3_N   <- c(100, 200, 500, 1000, 2000)
LIT_E3_EPS <- c(0.1, 0.3, 1.0, 3.0, 10.0)

lit_e3 <- function(v, src) lit_grid(v, rows = LIT_E3_EPS, cols = LIT_E3_N, src = src)

LIT_E3_COV_REPRO <- lit_e3(c(
  0.94, 0.95, 0.97, 0.98, 0.99,
  0.96, 0.96, 0.97, 0.98, 0.97,
  0.96, 0.97, 0.98, 0.98, 0.98,
  0.97, 0.97, 0.98, 0.99, 0.98,
  0.97, 0.97, 0.98, 0.99, 0.98), "Wang, Chang and Awan (2026), Fig. 6, Repro, coverage")

LIT_E3_WID_REPRO <- lit_e3(c(
  18.59, 18.75, 18.05, 15.48, 9.75,
  18.42, 15.86, 10.21,  3.54, 1.28,
  11.40,  6.89,  1.76,  0.98, 0.62,
  4.45,  2.01,  0.96,  0.76, 0.54,
  2.66,  1.70,  1.04,  0.74, 0.53), "Wang, Chang and Awan (2026), Fig. 6, Repro, width")

LIT_E3_COV_ADI <- lit_e3(c(
  0.99, 0.96, 0.98, 0.97, 0.90,
  0.95, 0.95, 0.94, 0.91, 0.85,
  0.91, 0.94, 0.89, 0.89, 0.87,
  0.95, 0.89, 0.90, 0.87, 0.89,
  0.89, 0.89, 0.88, 0.89, 0.88), "Wang, Chang and Awan (2026), Fig. 6, PB-ADI, coverage")

LIT_E3_WID_ADI <- lit_e3(c(
  9.51, 10.14, 11.09, 8.76, 3.69,
  6.50,  6.06,  3.89, 1.59, 0.70,
  3.15,  2.37,  1.06, 0.60, 0.38,
  1.74,  1.29,  0.73, 0.51, 0.34,
  1.59,  1.15,  0.68, 0.49, 0.34), "Wang, Chang and Awan (2026), Fig. 6, PB-ADI, width")

LIT_E3_COV_NAIVE <- lit_e3(c(
  0.38, 0.30, 0.20, 0.15, 0.13,
  0.26, 0.21, 0.16, 0.16, 0.18,
  0.25, 0.28, 0.33, 0.45, 0.54,
  0.58, 0.62, 0.69, 0.69, 0.69,
  0.74, 0.73, 0.73, 0.68, 0.71), "Wang, Chang and Awan (2026), Fig. 6, PB-naive, coverage")

LIT_E3_WID_NAIVE <- lit_e3(c(
  3.51, 2.85, 1.67, 0.97, 0.58,
  2.38, 1.49, 0.77, 0.47, 0.29,
  1.21, 0.79, 0.46, 0.32, 0.23,
  0.95, 0.69, 0.44, 0.32, 0.23,
  1.12, 0.75, 0.46, 0.33, 0.23), "Wang, Chang and Awan (2026), Fig. 6, PB-naive, width")

## ---------------------------------------------------------------------------
## Align a published grid to the keys of one of ours.  Missing keys give NA,
## which the heat maps render as a dashed empty cell rather than silently
## dropping a row or column.
## ---------------------------------------------------------------------------

lit_subset <- function(m, rows, cols) {
  r <- as.character(rows); cc <- as.character(cols)
  out <- matrix(NA_real_, nrow = length(r), ncol = length(cc),
                dimnames = list(r, cc))
  ri <- intersect(r, rownames(m)); ci <- intersect(cc, colnames(m))
  if (length(ri) && length(ci)) out[ri, ci] <- m[ri, ci]
  miss_r <- setdiff(r, rownames(m)); miss_c <- setdiff(cc, colnames(m))
  if (length(miss_r))
    message("  [note] published grid has no rows ", paste(miss_r, collapse = ", "))
  if (length(miss_c))
    message("  [note] published grid has no columns ", paste(miss_c, collapse = ", "))
  attr(out, "source") <- attr(m, "source")
  out
}

## ===========================================================================
## 4.  READERS
##
##  One function per results object.  Each locates its file under codes/,
##  verifies the columns it needs, and returns the object in the shape the
##  table and figure writers expect: a data frame for a sweep, a matrix with
##  dimnames for a grid.  A missing file returns NULL, and the caller skips
##  that table or figure rather than failing.
##
##  Column names are checked rather than assumed.  A file whose layout has
##  changed stops with the names actually found, which is a diagnosable
##  failure; indexing positionally would return a plausible wrong number.
## ===========================================================================

## helper: stop with a useful message when a column is absent
need_cols <- function(x, cols, what) {
  miss <- setdiff(cols, names(x))
  if (length(miss))
    stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
         "\n  columns found: ", paste(names(x), collapse = ", "))
  invisible(TRUE)
}

## helper: long format (row key, col key, value) -> matrix with dimnames
long_to_matrix <- function(d, rowvar, colvar, valvar) {
  rk <- sort(unique(d[[rowvar]]))
  ck <- sort(unique(d[[colvar]]))
  m <- matrix(NA_real_, length(rk), length(ck),
              dimnames = list(as.character(rk), as.character(ck)))
  for (i in seq_len(nrow(d)))
    m[as.character(d[[rowvar]][i]), as.character(d[[colvar]][i])] <-
    d[[valvar]][i]
  m
}

## ---------------------------------------------------------------------------
## 4.1  Experiment 1
## ---------------------------------------------------------------------------

read_exp1 <- function() {
  
  f <- find_one(DIR$exp1_cov,
                "^coverage_study_n[0-9]+_ep[0-9.]+_nSIM[0-9]+\\.csv$")
  if (is.na(f)) f <- file.path(DIR$exp1_cov, "coverage_all_results.csv")
  if (!have(f, "Experiment 1 summary")) return(NULL)
  
  S <- read_results(f)
  need_cols(S, c("method",
                 "coverage_mu", "coverage_mu_se",
                 "coverage_sigma", "coverage_sigma_se",
                 "coverage_joint", "coverage_joint_se",
                 "width_mu", "width_mu_se",
                 "width_sigma", "width_sigma_se",
                 "area", "area_se", "failure_rate"),
            "Experiment 1 summary")
  
  meta <- parse_stem(f)
  reps <- if (!is.na(meta$nSIM)) meta$nSIM else
    if ("N" %in% names(S)) S$N[1] else NA
  n_obs <- if (!is.na(meta$n)) meta$n else NA
  
  list(S = S, method = relabel_method(S$method), n = n_obs, reps = reps,
       source = f)
}

read_exp1_targeted <- function() {
  
  f <- file.path(DIR$exp1_cov, "coverage_targeted_efficient.csv")
  if (!have(f, "Experiment 1 targeted comparison")) return(NULL)
  T1 <- read_results(f)
  
  pick <- function(...) {
    for (nm in c(...)) if (nm %in% names(T1)) return(nm)
    stop("targeted file has none of: ", paste(c(...), collapse = ", "),
         "\n  columns found: ", paste(names(T1), collapse = ", "))
  }
  data.frame(
    target = unquote(T1[[pick("target", "parameter")]]),
    mah    = T1[[pick("mahalanobis_width", "mah_width")]],
    pw     = T1[[pick("efficient_width", "eff_width",
                      "penalized_wald_width")]],
    reduce = T1[[pick("width_reduction_percent", "reduction",
                      "width_reduction")]],
    cov    = T1[[pick("coverage")]],
    cov_se = T1[[pick("coverage_se")]],
    stringsAsFactors = FALSE)
}

## Diagnostics the summary does not record: how often the full-box
## convention fired, and whether the summary follows from the replications.
read_exp1_raw <- function(mu_box = c(-2, 5), sigma_box = c(0.1, 5),
                          blocks = c("Mahalanobis", "Penalized Wald (mu)",
                                     "Penalized Wald (sigma)", "PB-ADI"),
                          per_block = 7L) {
  
  f <- find_one(DIR$exp1_cov,
                "^coverage_study_raw_n[0-9]+_ep[0-9.]+_nSIM[0-9]+\\.csv$")
  if (is.na(f)) { message("  [skip] no Experiment 1 replication file"); return(NULL) }
  
  RW <- read_results(f)
  if (ncol(RW) != per_block * length(blocks)) {
    message(sprintf(paste0("  [skip] replication file has %d columns, ",
                           "expected %d; the block layout has changed and ",
                           "positional reading would be unsafe"),
                    ncol(RW), per_block * length(blocks)))
    return(NULL)
  }
  
  num <- function(v) suppressWarnings(as.numeric(v))
  out <- NULL
  for (b in seq_along(blocks)) {
    o <- (b - 1L) * per_block
    w_mu <- num(RW[[o + 4L]]); w_sig <- num(RW[[o + 5L]]); fail <- num(RW[[o + 7L]])
    out <- rbind(out, data.frame(
      method   = blocks[b],
      full_mu  = sum(abs(w_mu  - diff(mu_box))    < 1e-6, na.rm = TRUE),
      full_sig = sum(abs(w_sig - diff(sigma_box)) < 1e-6, na.rm = TRUE),
      max_mu   = max(w_mu,  na.rm = TRUE),
      max_sig  = max(w_sig, na.rm = TRUE),
      failure  = mean(fail, na.rm = TRUE),
      stringsAsFactors = FALSE))
  }
  attr(out, "reps") <- nrow(RW)
  attr(out, "box")  <- c(diff(mu_box), diff(sigma_box))
  out
}

## ---------------------------------------------------------------------------
## 4.2  Experiment 2
## ---------------------------------------------------------------------------

## Efficient_INT.csv is written wide, with the true slope as row names and
## the sample sizes as column names, so the header cell is empty.
read_exp2_grid <- function(file = "Efficient_INT.csv", what = "power grid") {
  
  f <- file.path(DIR$exp2_eff, file)
  if (!have(f, paste("Experiment 2", what))) return(NULL)
  G <- as.matrix(read_results(f, row.names = 1, check.names = FALSE))
  storage.mode(G) <- "double"
  G
}

## ---------------------------------------------------------------------------
## 4.3  Experiment 3
## ---------------------------------------------------------------------------

read_exp3 <- function() {
  
  f <- file.path(DIR$exp3_res, "heatmap_data.csv")
  if (!have(f, "Experiment 3 grid")) return(NULL)
  H <- read_results(f)
  need_cols(H, c("n", "epsilon", "coverage", "width_all"),
            "heatmap_data.csv")
  list(cov = long_to_matrix(H, "epsilon", "n", "coverage"),
       wid = long_to_matrix(H, "epsilon", "n", "width_all"),
       empty = if ("n_no_accepted_point" %in% names(H))
         long_to_matrix(H, "epsilon", "n", "n_no_accepted_point") else NULL,
       reps = if ("reps" %in% names(H)) H$reps[1] else NA)
}

## ---------------------------------------------------------------------------
## 4.4  Sensitivity sweeps
##
##  Discovered by scanning rather than listed, because the sweeps present
##  differ between checkouts.  The design column is the one that varies.
## ---------------------------------------------------------------------------

SWEEP_KEYS <- c("lambda_n", "R_aux", "R_synthetic", "Delta", "epsilon")

read_sweeps <- function() {
  
  if (!dir.exists(DIR$sens)) {
    message("  [skip] no sensitivity_analysis/ directory")
    return(list())
  }
  files <- list.files(DIR$sens, pattern = "_results\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("_raw\\.csv$", files)]
  if (length(files) == 0L) {
    message("  [skip] no *_results.csv under sensitivity_analysis/")
    return(list())
  }
  
  out <- list()
  for (f in files) {
    S <- read_results(f)
    cand <- intersect(SWEEP_KEYS, names(S))
    cand <- cand[vapply(cand, function(k) length(unique(S[[k]])) > 1L,
                        logical(1))]
    if (length(cand) == 0L) {
      message("  [skip] no varying design column in ", basename(f))
      next
    }
    key <- cand[1]
    if (length(cand) > 1L)
      message("  [note] ", basename(f), " varies ",
              paste(cand, collapse = " and "), "; using ", key)
    S <- S[order(S[[key]]), , drop = FALSE]
    attr(S, "key")  <- key
    attr(S, "file") <- f
    out[[key]] <- S
    message("  [read]  ", key, " sweep: ", nrow(S), " grid points")
  }
  out
}

## the clamping sweep is long-format and needs reshaping into two grids
read_delta <- function() {
  
  if (!dir.exists(DIR$sens)) return(NULL)
  f <- list.files(DIR$sens, pattern = "delta_sensitivity_long\\.csv$",
                  recursive = TRUE, full.names = TRUE)
  if (length(f) == 0L) {
    message("  [skip] no delta_sensitivity_long.csv"); return(NULL)
  }
  D <- read_results(f[1])
  need_cols(D, c("beta1_true", "Delta", "n", "rejection_probability"),
            basename(f[1]))
  split_by <- function(b) {
    d <- D[D$beta1_true == b, , drop = FALSE]
    long_to_matrix(d, "Delta", "n", "rejection_probability")
  }
  list(h0 = split_by(0), h1 = split_by(1), raw = D)
}

## ---------------------------------------------------------------------------
## 4.5  Sigma-target study
## ---------------------------------------------------------------------------

read_sigma <- function() {
  
  g <- function(nm) {
    f <- file.path(DIR$exp2_sigma, nm)
    if (file.exists(f)) read_results(f) else NULL
  }
  out <- list(summary  = g("Table_Exp2_sigma_summary.csv"),
              audit    = g("exp2_sigma_search_audit_summary.csv"),
              levelset = g("exp2_sigma_empirical_levelset_summary.csv"),
              settings = g("exp2_sigma_settings.csv"))
  if (all(vapply(out, is.null, logical(1)))) {
    message("  [skip] no sigma-target output in ",
            sub(ROOT, "<root>", DIR$exp2_sigma, fixed = TRUE))
    return(NULL)
  }
  out
}

## ---------------------------------------------------------------------------
## 4.6  Clamping-boundary studies
## ---------------------------------------------------------------------------

## master_grid_summary.csv is long: one row per (mu_star, epsilon, method)
read_boundary_grid <- function() {
  
  f <- file.path(DIR$pbadi, "results_grid_study_v1", "master_grid_summary.csv")
  if (!have(f, "boundary grid")) return(NULL)
  B <- read_results(f)
  need_cols(B, c("mu_star", "epsilon", "method"), basename(f))
  B$method <- unquote(B$method)
  
  key <- function(m) {
    if (grepl("Mahal", m, ignore.case = TRUE))   return("mah")
    if (grepl("Effic|Wald", m, ignore.case = TRUE)) return("eff")
    if (grepl("PB|ADI",  m, ignore.case = TRUE)) return("pb")
    NA_character_
  }
  B$key <- vapply(B$method, key, character(1))
  if (any(is.na(B$key)))
    message("  [note] unrecognised method label(s): ",
            paste(unique(B$method[is.na(B$key)]), collapse = ", "))
  
  grab <- function(col) {
    if (!(col %in% names(B))) return(NULL)
    setNames(lapply(c("mah", "eff", "pb"), function(k)
      long_to_matrix(B[B$key == k, , drop = FALSE], "mu_star", "epsilon", col)),
      c("mah", "eff", "pb"))
  }
  list(joint  = grab("cov_joint"),
       cov_mu = grab("cov_mu"),
       wid_mu = grab("width_mu"),
       wid_sig = grab("width_sigma"),
       raw = B)
}

read_boundary_soft <- function() {
  
  f <- file.path(DIR$pbadi, "results_softclamp_study_v2",
                 "master_softclamp_summary.csv")
  if (!have(f, "soft-clamp study")) return(NULL)
  S <- read_results(f)
  need_cols(S, c("clamp_type", "epsilon", "method"), basename(f))
  S$clamp_type <- unquote(S$clamp_type)
  S$method     <- unquote(S$method)
  
  key <- function(m) {
    if (grepl("Mahal", m, ignore.case = TRUE))      return("mah")
    if (grepl("Effic|Wald", m, ignore.case = TRUE)) return("eff")
    if (grepl("PB|ADI",  m, ignore.case = TRUE))    return("pb")
    NA_character_
  }
  S$key <- vapply(S$method, key, character(1))
  eps <- sort(unique(S$epsilon))
  
  out <- list()
  for (k in c("mah", "eff", "pb"))
    for (ct in c("hard", "soft")) {
      d <- S[S$key == k & S$clamp_type == ct, , drop = FALSE]
      d <- d[order(d$epsilon), , drop = FALSE]
      pull <- function(col) if (col %in% names(d)) d[[col]] else rep(NA_real_, nrow(d))
      out[[paste0(k, ".", ct)]] <- list(
        joint    = pull("cov_joint"),
        cov_mu   = pull("cov_mu"),
        width_mu = pull("width_mu"),
        area     = pull("area"))
    }
  attr(out, "eps") <- eps
  out
}


## ===========================================================================
## 5.  FIGURE GEOMETRY
## ===========================================================================


## ===========================================================================
## 4.  FIGURE WRITERS
##
##  Four problems with an earlier version of these figures are corrected here,
##  and each correction is marked in the code below.
##
##  (F1) Panel titles carried attributions such as "(this paper)" and
##       "(Wang, Chang and Awan)".  A panel title is a label, not a citation:
##       titles are now the method name alone, identical in every figure, and
##       attribution is made once in the caption.
##
##  (F2) The shared colour key was a tall narrow column beside the panels.
##       That leaves a large empty region, places the key title far from the
##       bar, and crowds the tick labels.  The key is now a horizontal strip
##       beneath the panels at a fixed physical height, so it cannot be
##       squeezed by the device size.
##
##  (F3) Row labels and the axis title were clipped on the left, because a
##       2 x 2 grid repeated both axis titles on all four panels and the left
##       margin was sized for neither.  Axis titles are now drawn only on the
##       outer panels -- the y title on the left column, the x title on the
##       bottom row -- which removes the clutter and leaves room for them.
##
##  (F4) The width panels carried a legend inside the plotting region, over
##       the curves.  It is now a single horizontal key beneath the row.
##
##  Heat maps are drawn cell by cell with rect() rather than with image(),
##  because image() interpolates the colour scale linearly and would lose the
##  break at the nominal level.  Symbols are plotmath throughout; note that
##  c("text", expression(x)) coerces the expression to its source text, so
##  titles are carried in a list and indexed with [[ ]], and legend entries
##  that mix words and symbols are built in one call to parse().
## ===========================================================================

## ---------------------------------------------------------------------------
## 4.1  Labels, defined once
## ---------------------------------------------------------------------------

LAB_N      <- expression("Sample size " * italic(n))
LAB_N_S    <- expression(italic(n))
LAB_BETA1  <- expression("True " * beta[1]^"*")
LAB_EPS    <- expression("Privacy budget " * epsilon)
LAB_DELTA  <- expression("Clamping bound " * Delta)
LAB_MUSTAR <- expression("True mean " * mu^"*")
LAB_LAMBDA <- expression(lambda[n])
LAB_WIDTH  <- "Average interval width"
LAB_CIW    <- "CI width"

## (F1) method labels: the name only, no attribution
M_PW  <- "Penalized Wald"
M_RP  <- "Repro"
M_ADI <- "PB-ADI"
M_NV  <- "PB-naive"
M_MAH <- "Mahalanobis"

## ---------------------------------------------------------------------------
## 4.2  One heat-map panel
##
##  Margins are in lines.  The left margin is sized for the row labels plus
##  an axis title at 2.6 lines; when no y title is drawn the margin shrinks,
##  which is what lets the panels fill the width.
## ---------------------------------------------------------------------------

heat_panel <- function(m, greyfn, main = "", xlab = NULL, ylab = NULL,
                       dec = 3, cex_cell = 0.70, cex_ax = 0.78,
                       mar = c(3.1, 3.9, 1.8, 0.6)) {
  
  m  <- as.matrix(m)
  nr <- nrow(m); nc <- ncol(m)
  
  ## (F3) margins are uniform across panels, so every panel is the same
  ## size and the cells are directly comparable; the axis titles are drawn
  ## only on the outer panels, which is what removes the clutter
  par(mar = mar, mgp = c(1.9, 0.40, 0), tcl = -0.18)
  
  plot(NA, xlim = c(0, nc), ylim = c(0, nr), axes = FALSE,
       xlab = "", ylab = "", main = "", xaxs = "i", yaxs = "i")
  
  for (i in seq_len(nr)) {
    y <- nr - i
    for (j in seq_len(nc)) {
      v <- m[i, j]
      g <- greyfn(v)
      rect(j - 1, y, j, y + 1, col = pct_to_col(g),
           border = "grey85", lwd = 0.4)
      text(j - 0.5, y + 0.5,
           if (is.na(v)) "--" else fcell(v, dec), cex = cex_cell,
           col = if (is.na(g)) "grey55" else if (g >= 46) "white" else "black")
    }
  }
  
  axis(1, at = seq_len(nc) - 0.5, labels = colnames(m),
       tick = FALSE, line = -0.65, cex.axis = cex_ax)
  axis(2, at = rev(seq_len(nr)) - 0.5, labels = rownames(m),
       tick = FALSE, line = -0.55, las = 1, cex.axis = cex_ax)
  box(col = "grey70", lwd = 0.5)
  
  ## (M1) the x title sits directly under the tick labels, not 25pt below;
  ## (M4) the y title clears the row labels without touching the edge
  if (!identical(main, "")) mtext(main, side = 3, line = 0.35, cex = 0.84)
  if (!is.null(xlab))       mtext(xlab, side = 1, line = 1.75, cex = 0.82)
  if (!is.null(ylab))       mtext(ylab, side = 2, line = 2.55, cex = 0.82)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## 4.3  (F2) Horizontal colour key
##
##  The bar is inset at each end and the value scale mapped onto the inset
##  range, so a tick label at 0 or 1 stays inside the device while the key
##  keeps the same margins as the panels above it.
## ---------------------------------------------------------------------------

colour_key <- function(greyfn, ticks, label, nseg = 160, pad = 0.03,
                       mar = c(1.9, 3.9, 0.2, 0.6)) {
  
  at_x <- function(v) pad + (1 - 2 * pad) * v
  
  par(mar = mar, mgp = c(1.4, 0.25, 0), tcl = -0.13)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (k in seq_len(nseg)) {
    v0 <- (k - 1) / nseg; v1 <- k / nseg
    rect(at_x(v0), 0.45, at_x(v1), 0.95,
         col = pct_to_col(greyfn((v0 + v1) / 2)), border = NA)
  }
  rect(at_x(0), 0.45, at_x(1), 0.95, border = "grey60", lwd = 0.5)
  axis(1, at = at_x(ticks),
       labels = formatC(ticks, format = "f", digits = 2),
       pos = 0.45, cex.axis = 0.72, lwd = 0, lwd.ticks = 0.6)
  mtext(label, side = 1, line = 1.05, cex = 0.80)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## 4.4  Grid of heat-map panels over a shared key
##
##  nrow x ncol panels, filled row-wise.  Axis titles are drawn on the outer
##  panels only.  The key occupies a fixed 2.2 cm row, so shrinking the
##  device takes space from the panels and never collapses the bar.
## ---------------------------------------------------------------------------

fig_heat_grid <- function(panels, titles, greyfn, ticks, key_label,
                          file, xlab, ylab, nrow, ncol,
                          dec = NULL, width, height,
                          cex_cell = 0.70, cex_ax = 0.78,
                          row_labels = NULL, device = TRUE, key = TRUE) {
  
  k <- length(panels)
  stopifnot(is.list(titles), length(titles) == k, k <= nrow * ncol)
  if (is.null(dec)) dec <- rep(3, k)
  
  path <- file.path(DIR_FIG, file)
  if (device) {
    pdf(path, width = width, height = height)
    on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
            add = TRUE)
  }
  
  ## 0 is layout()'s code for a blank cell, and the figure indices must run
  ## contiguously from 1, so the key is always index k + 1
  idx <- rep(0L, nrow * ncol)
  idx[seq_len(k)] <- seq_len(k)
  lay <- matrix(idx, nrow = nrow, ncol = ncol, byrow = TRUE)
  if (key) {
    lay <- rbind(lay, rep(k + 1L, ncol))    # key row spans the full width
    layout(lay, heights = c(rep(1, nrow), lcm(1.9)))
  } else layout(lay, heights = rep(1, nrow))
  
  for (i in seq_len(k)) {
    r <- ((i - 1) %/% ncol) + 1L
    c_ <- ((i - 1) %% ncol) + 1L
    ## (F3) outer panels only
    is_bottom <- (r == nrow) || (i + ncol > k)
    is_left   <- (c_ == 1L)
    ## a row descriptor sits further out in the left margin than the y title
    extra <- if (is_left && !is.null(row_labels)) 1.5 else 0
    heat_panel(panels[[i]], greyfn, main = titles[[i]],
               xlab = if (is_bottom) xlab else NULL,
               ylab = if (is_left)   ylab else NULL,
               dec = dec[i], cex_cell = cex_cell, cex_ax = cex_ax,
               mar = c(3.1, 3.9 + extra, 1.8, 0.6))
    if (is_left && !is.null(row_labels) && r <= length(row_labels))
      mtext(row_labels[[r]], side = 2, line = 4.0, cex = 0.86, font = 1)
  }
  if (key) colour_key(greyfn, ticks, key_label)
  invisible(path)
}

## ---------------------------------------------------------------------------
## 4.5  (F4) Horizontal key for the line plots
## ---------------------------------------------------------------------------

line_key <- function(labels, cols, ltys, pchs, mar = c(0.2, 3.9, 0.2, 0.6)) {
  n <- length(labels)
  par(mar = mar)
  plot(NA, xlim = c(0, n), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) {
    x0 <- i - 1 + 0.06; x1 <- x0 + 0.26
    segments(x0, 0.5, x1, 0.5, col = cols[i], lty = ltys[i], lwd = 1.7)
    if (!is.na(pchs[i]))
      points((x0 + x1) / 2, 0.5, col = cols[i], pch = pchs[i], cex = 0.9)
    text(x1 + 0.07, 0.5, labels[[i]], adj = c(0, 0.5), cex = 0.78)
  }
  invisible(NULL)
}

## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 5.5  Line plots, parameterised by the data the readers supply
## ---------------------------------------------------------------------------

fig_exp2_size_plot <- function(v_pw, v_rp, v_adi, v_nv, n_all) {
  
  path <- file.path(DIR_FIG, "fig_exp2_size.pdf")
  pdf(path, width = 6.3, height = 3.6)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  ## (F4) plot on top, key beneath
  layout(matrix(1:2, nrow = 2), heights = c(1, lcm(1.0)))
  
  ymax <- max(c(v_pw, v_rp, v_adi, v_nv), na.rm = TRUE) * 1.08
  par(mar = c(3.3, 4.1, 0.7, 0.8), mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
  plot(NA, xlim = range(n_all), ylim = c(0, ymax), log = "x",
       xlab = "", ylab = "", axes = FALSE)
  axis(1, at = n_all, labels = n_all, cex.axis = 0.80)
  axis(2, cex.axis = 0.80)
  box(col = "grey60", lwd = 0.6)
  mtext(LAB_N, side = 1, line = 2.0, cex = 0.90)
  mtext(expression("Empirical size at " * beta[1]^"*" == 0),
        side = 2, line = 2.8, cex = 0.90, las = 0)
  abline(h = pretty(c(0, ymax)), col = "grey93", lwd = 0.6)
  abline(h = ALPHA_12, lty = 2, col = "grey25", lwd = 1.1)
  abline(h = SZB_12,   lty = 3, col = "grey50", lwd = 1.1)
  
  cols <- c("black", "grey40", "grey55", "grey72")
  ltys <- c(1, 1, 2, 4)
  pchs <- c(19, 15, 17, 18)
  vals <- list(v_pw, v_rp, v_adi, v_nv)
  for (i in 1:4) {
    keep <- !is.na(vals[[i]])
    lines(n_all[keep], vals[[i]][keep], col = cols[i], lty = ltys[i], lwd = 1.7)
    points(n_all[keep], vals[[i]][keep], col = cols[i], pch = pchs[i], cex = 0.85)
  }
  
  line_key(list(M_PW, M_RP, M_ADI, M_NV), cols, ltys, pchs)
  
}

fig_exp3_width_plot <- function(W_pw, W_rp, W_adi, eps, ns) {
  
  path <- file.path(DIR_FIG, "fig_exp3_width.pdf")
  pdf(path, width = 8.0, height = 3.4)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  ## (F4) three panels on top, one key beneath, no legend over the curves
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
         heights = c(1, lcm(1.0)))
  
  panels <- list(W_pw, W_rp, W_adi)
  titles <- list(M_PW, M_RP, M_ADI)
  cols <- c("black", "grey40", "grey55", "grey72")
  ltys <- c(1, 2, 3, 4)
  pchs <- c(19, 15, 17, 18)
  yt   <- c(0.5, 1, 2, 5, 10, 20)
  
  for (p in seq_along(panels)) {
    W <- panels[[p]]
    ## (F3) y title on the first panel only
    par(mar = c(3.3, if (p == 1L) 4.1 else 2.4, 1.8, 0.7),
        mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
    plot(NA, xlim = range(ns), ylim = c(0.35, 30), log = "xy",
         xlab = "", ylab = "", axes = FALSE)
    axis(1, at = ns, labels = ns, cex.axis = 0.78)
    axis(2, at = yt, labels = yt, cex.axis = 0.78)
    box(col = "grey60", lwd = 0.6)
    mtext(titles[[p]], side = 3, line = 0.45, cex = 0.82)
    mtext(LAB_N_S, side = 1, line = 1.9, cex = 0.82)
    if (p == 1L) mtext(LAB_CIW, side = 2, line = 2.7, cex = 0.82, las = 0)
    abline(h = yt, col = "grey93", lwd = 0.5)
    abline(h = E3_BOXW, lty = 2, col = "grey55", lwd = 1.1)
    for (i in seq_along(eps)) {
      lines(ns, W[i, ], col = cols[i], lty = ltys[i], lwd = 1.7)
      points(ns, W[i, ], col = cols[i], pch = pchs[i], cex = 0.85)
    }
  }
  
  line_key(as.list(parse(text = c(sprintf("epsilon == %s", eps),
                                  '"search box"'))),
           c(cols[seq_along(eps)], "grey55"),
           c(ltys[seq_along(eps)], 2),
           c(pchs[seq_along(eps)], NA),
           mar = c(0.3, 4.4, 0.3, 0.9))
  
}

fig_lambda_plot <- function(LAMBDA) {
  
  path <- file.path(DIR_FIG, "fig_sens_lambda.pdf")
  pdf(path, width = 6.1, height = 3.5)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  layout(matrix(1:2, nrow = 2), heights = c(1, lcm(1.0)))
  
  ## the "selected" column is optional; fall back to the minimum objective
  sel <- if ("selected" %in% names(LAMBDA))
    as.logical(unquote(LAMBDA$selected)) else
      seq_len(nrow(LAMBDA)) == which.min(
        LAMBDA$width_mu / min(LAMBDA$width_mu) +
          LAMBDA$width_sigma / min(LAMBDA$width_sigma))
  x  <- seq_along(LAMBDA$lambda_n) - 1        # geometric grid, ordinal spacing
  yl <- range(c(LAMBDA$width_mu, LAMBDA$width_sigma)) * c(0.90, 1.08)
  
  par(mar = c(3.3, 4.1, 0.7, 0.8), mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
  plot(NA, xlim = range(x), ylim = yl, axes = FALSE, xlab = "", ylab = "")
  axis(1, at = x, labels = formatC(LAMBDA$lambda_n, format = "g"), cex.axis = 0.80)
  axis(2, cex.axis = 0.80)
  box(col = "grey60", lwd = 0.6)
  mtext(LAB_LAMBDA, side = 1, line = 2.0, cex = 0.90)
  mtext(LAB_WIDTH,  side = 2, line = 2.8, cex = 0.90, las = 0)
  abline(h = pretty(yl), col = "grey93", lwd = 0.6)
  abline(v = x[sel], lty = 3, col = "grey45", lwd = 1.2)
  
  lines(x, LAMBDA$width_mu,   col = "black",  lwd = 1.7)
  points(x, LAMBDA$width_mu,  col = "black",  pch = 19, cex = 0.85)
  lines(x, LAMBDA$width_sigma,  col = "grey45", lwd = 1.7, lty = 2)
  points(x, LAMBDA$width_sigma, col = "grey45", pch = 15, cex = 0.85)
  
  line_key(as.list(parse(text = c('mu * " target"', 'sigma * " target"',
                                  '"selected"'))),
           c("black", "grey45", "grey45"), c(1, 2, 3), c(19, 15, NA))
  
}

fig_boundary_soft_plot <- function(S2) {
  
  path <- file.path(DIR_FIG, "fig_pbadi_soft.pdf")
  pdf(path, width = 6.4, height = 4.0)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  layout(matrix(1:2, nrow = 2), heights = c(1, lcm(1.6)))
  
  eps <- as.numeric(attr(S2, "eps"))
  par(mar = c(3.3, 4.3, 0.7, 0.8), mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
  plot(NA, xlim = range(eps), ylim = c(0.10, 1.04), log = "x",
       axes = FALSE, xlab = "", ylab = "")
  axis(1, at = eps, labels = format(eps), cex.axis = 0.80)
  axis(2, at = seq(0.2, 1.0, by = 0.2), cex.axis = 0.80)
  box(col = "grey60", lwd = 0.6)
  mtext(LAB_EPS, side = 1, line = 2.0, cex = 0.90)
  mtext(expression("Joint coverage at " * mu^"*" == 3),
        side = 2, line = 3.0, cex = 0.90, las = 0)
  abline(h = seq(0.2, 1.0, by = 0.2), col = "grey93", lwd = 0.6)
  abline(h = BND_12, lty = 2, col = "grey40", lwd = 1.1)
  
  mm   <- c("mah", "eff", "pb")
  cols <- c(mah = "black", eff = "grey42", pb = "grey62")
  pchs <- c(mah = 19, eff = 15, pb = 17)
  for (m in mm) {
    H <- S2[[paste0(m, ".hard")]]$joint
    S <- S2[[paste0(m, ".soft")]]$joint
    lines(eps, H, col = cols[[m]], lwd = 1.7, lty = 1)
    points(eps, H, col = cols[[m]], pch = pchs[[m]], cex = 0.9)
    lines(eps, S, col = cols[[m]], lwd = 1.7, lty = 2)
    points(eps, S, col = cols[[m]], pch = pchs[[m]], cex = 0.9)
  }
  
  ## two rows of three entries: method by clamp type
  par(mar = c(0.3, 4.8, 0.3, 1.0))
  plot(NA, xlim = c(0, 3), ylim = c(0, 2), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  lab <- c(M_MAH, M_PW, M_ADI)
  for (j in seq_along(mm)) for (r in 1:2) {
    y  <- if (r == 1L) 1.45 else 0.45
    x0 <- j - 1 + 0.04; x1 <- x0 + 0.24
    segments(x0, y, x1, y, col = cols[[mm[j]]],
             lty = if (r == 1L) 1 else 2, lwd = 1.7)
    points((x0 + x1) / 2, y, col = cols[[mm[j]]], pch = pchs[[mm[j]]], cex = 0.9)
    text(x1 + 0.06, y,
         paste0(lab[j], if (r == 1L) ", hard" else ", soft"),
         adj = c(0, 0.5), cex = 0.76)
  }
  
}


## ---------------------------------------------------------------------------
## 5.6  Figures the simulation scripts produce directly
##
##  Three Experiment-1 figures and three from the sigma-target study are
##  drawn by the simulation code, not here.  They are copied into
##  results/figures unchanged and given a wrapper, so the whole figure set
##  lives in one place and the paths resolve the same way.
## ---------------------------------------------------------------------------

copy_figure <- function(src_dir, pattern, out_name, caption, label,
                        width = 0.80) {
  
  f <- find_one(src_dir, pattern)
  if (is.na(f)) {
    message("  [skip] no match for ", pattern, " in ",
            sub(ROOT, "<root>", src_dir, fixed = TRUE))
    return(invisible(NULL))
  }
  dest <- file.path(DIR_FIG, out_name)
  ok <- file.copy(f, dest, overwrite = TRUE)
  if (!ok) { message("  [ERROR] could not copy ", f); return(invisible(NULL)) }
  message("  [copy]  ", dest)
  write_figure_wrapper(out_name, caption, label, width = width)
}

fig_exp1_extra <- function() {
  
  tab_need_init()
  
  copy_figure(DIR$exp1_apply, "^comparison_.*_region\\.pdf$",
              "exp1_comparison_region.pdf",
              paste0("Experiment 1, single draw. Boundaries of the joint acceptance ",
                     "region for Mahalanobis depth and for the two Penalized Wald targets, ",
                     "together with the \\textsc{pb-adi} bootstrap ellipse, all computed ",
                     "from the same released summary, inference cloud and auxiliary cloud, ",
                     "so the regions differ only through the statistic."),
              "fig:exp1-region", width = 0.72)
  
  copy_figure(DIR$exp1_apply, "^marginal_CI_mu_.*\\.pdf$",
              "exp1_marginal_CI_mu.pdf",
              paste0("Experiment 1, single draw. The four $95\\%$ marginal intervals ",
                     "for $\\mu$ from the same released summary, with the true value ",
                     "marked."),
              "fig:exp1-mci-mu", width = 0.72)
  
  copy_figure(DIR$exp1_apply, "^marginal_CI_sigma_.*\\.pdf$",
              "exp1_marginal_CI_sigma.pdf",
              paste0("Experiment 1, single draw. The four $95\\%$ marginal intervals ",
                     "for $\\sigma$ from the same released summary."),
              "fig:exp1-mci-sigma", width = 0.72)
}

fig_sigma_extra <- function() {
  
  tab_need_init()
  
  copy_figure(DIR$exp2_sigma, "^Figure_Exp2_sigma_performance\\.pdf$",
              "Figure_Exp2_sigma_performance.pdf",
              paste0("Experiment 2, $\\sigma$ target. Empirical coverage with an ",
                     "exact Clopper--Pearson Monte Carlo interval against the nominal level ",
                     "and the finite-$R$ bound, and the distribution of interval widths ",
                     "over replications yielding a finite interval."),
              "fig:exp2-sigma-performance", width = 0.92)
  
  copy_figure(DIR$exp2_sigma, "^Figure_Exp2_sigma_empirical_geometry\\.pdf$",
              "Figure_Exp2_sigma_empirical_geometry.pdf",
              paste0("Experiment 2, $\\sigma$ target. Monte Carlo mean level sets in ",
                     "the $(\\beta_1,\\beta_0)$ nuisance section, with the remaining ",
                     "coordinates fixed at their true values: the nuisance-orthogonal Wald ",
                     "component, the Mahalanobis penalty, and the combined statistic."),
              "fig:exp2-sigma-geometry", width = 0.92)
  
  copy_figure(DIR$exp2_sigma,
              "^Figure_Exp2_sigma_empirical_Repro_region\\.pdf$",
              "Figure_Exp2_sigma_empirical_Repro_region.pdf",
              paste0("Experiment 2, $\\sigma$ target. Empirical inclusion probability ",
                     "on the $(\\beta_1,\\beta_0)$ grid. At the true parameter this ",
                     "estimates the coverage of the rule; away from it the surface is a ",
                     "repeated-sampling inclusion probability and not a confidence region, ",
                     "so its $0.95$ contour is a level set rather than the boundary of any ",
                     "realised set."),
              "fig:exp2-sigma-region", width = 0.68)
}

## ---------------------------------------------------------------------------
## 5.7  Sweep figures for the two cloud sizes
##
##  Same layout as the penalty-weight figure.  The R sweep additionally
##  carries the moving coverage bound, which is drawn as a step reference.
## ---------------------------------------------------------------------------

fig_sweep_plot <- function(S, key, xlab, file, ordinal = TRUE,
                           bound_col = NULL) {
  
  path <- file.path(DIR_FIG, file)
  pdf(path, width = 6.1, height = 3.5)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  layout(matrix(1:2, nrow = 2), heights = c(1, lcm(1.0)))
  
  v  <- S[[key]]
  x  <- if (ordinal) seq_along(v) - 1 else v
  yl <- range(c(S$width_mu, S$width_sigma)) * c(0.92, 1.08)
  
  par(mar = c(3.3, 4.1, 0.7, 0.8), mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
  plot(NA, xlim = range(x), ylim = yl, axes = FALSE, xlab = "", ylab = "")
  axis(1, at = x, labels = format(v, trim = TRUE), cex.axis = 0.80)
  axis(2, cex.axis = 0.80)
  box(col = "grey60", lwd = 0.6)
  mtext(xlab,      side = 1, line = 2.0, cex = 0.90)
  mtext(LAB_WIDTH, side = 2, line = 2.8, cex = 0.90, las = 0)
  abline(h = pretty(yl), col = "grey93", lwd = 0.6)
  
  lines(x, S$width_mu,     col = "black",  lwd = 1.7)
  points(x, S$width_mu,    col = "black",  pch = 19, cex = 0.85)
  lines(x, S$width_sigma,  col = "grey45", lwd = 1.7, lty = 2)
  points(x, S$width_sigma, col = "grey45", pch = 15, cex = 0.85)
  
  line_key(as.list(parse(text = c('mu * " target"', 'sigma * " target"'))),
           c("black", "grey45"), c(1, 2), c(19, 15))
  invisible(path)
}

fig_sweeps_clouds <- function(SW = read_sweeps()) {
  
  A <- SW[["R_aux"]]
  if (!is.null(A)) {
    fig_sweep_plot(A, "R_aux", expression(R[aux]), "fig_sens_raux.pdf")
    write_figure_wrapper("fig_sens_raux.pdf",
                         sprintf(paste0("Sensitivity to the auxiliary cloud size ",
                                        "$R_{\\mathrm{aux}}$. Over a %g-fold range the $\\mu$-target ",
                                        "width varies by $%.4f$ and the $\\sigma$-target width by $%.4f$, ",
                                        "in both cases well inside one Monte Carlo standard error of the ",
                                        "width itself. $R_{\\mathrm{aux}}$ never enters the acceptance ",
                                        "rule, so the coverage bound is $%.5f$ throughout."),
                                 max(A$R_aux) / min(A$R_aux), diff(range(A$width_mu)),
                                 diff(range(A$width_sigma)), BND_12),
                         "fig:sens-raux", width = 0.76)
  }
  
  B <- SW[["R_synthetic"]]
  if (!is.null(B)) {
    fig_sweep_plot(B, "R_synthetic", expression(R), "fig_sens_R.pdf")
    bnd_txt <- if ("finite_R_reference" %in% names(B))
      sprintf(paste0(" The coverage bound moves with $R$, from $%.5f$ at ",
                     "$R=%g$ to $%.5f$ at $R=%g$, and is reported per row in ",
                     "Table~\\ref{tab:sens-R}."),
              B$finite_R_reference[1], B$R_synthetic[1],
              B$finite_R_reference[nrow(B)], B$R_synthetic[nrow(B)]) else ""
    write_figure_wrapper("fig_sens_R.pdf",
                         paste0("Sensitivity to the inference cloud size $R$: average interval ",
                                "width for each target. Grid points are drawn at equal ordinal ",
                                "spacing because the design grid is geometric.", bnd_txt),
                         "fig:sens-R", width = 0.76)
  }
  invisible(NULL)
}


## ---------------------------------------------------------------------------
## 5.8  Combined figures
##
##  Experiment 3 and the clamping sweep each produced three floats from one
##  matrix in an earlier version.  Each is now a single figure: rows carry
##  the two quantities and columns the methods, which is also the layout the
##  source figures use.
## ---------------------------------------------------------------------------

## Experiment 3: coverage on the first row, width on the second, four methods
fig_exp2_sigma_geometry <- function() {
  
  tab_need_init()
  copy_figure(DIR$exp2_sigma, "^Figure_Exp2_sigma_empirical_geometry\\.pdf$",
              "fig_exp2_sigma_geometry.pdf",
              paste0("Experiment 2, $\\sigma$ target. Monte Carlo mean level sets in ",
                     "the $(\\beta_1,\\beta_0)$ nuisance section, with the remaining ",
                     "coordinates fixed at their true values: the nuisance-orthogonal Wald ",
                     "component, the Mahalanobis penalty, and the combined statistic. The ",
                     "figure carries information the summary table does not, namely the ",
                     "shape of the accepted set in the nuisance directions and the relative ",
                     "contribution of the two components."),
              "fig:exp2-sigma-geometry", width = 0.92)
}


fig_exp3_comparison <- function(E3 = read_exp3()) {
  
  tab_need_init()
  if (is.null(E3)) return(invisible(NULL))
  
  eps <- rownames(E3$cov); ns <- colnames(E3$cov)
  nsn <- as.numeric(ns)
  gc3 <- function(v) grey_cov(v, nominal = 1 - ALPHA_3)
  
  cov_p <- list(E3$cov,
                lit_subset(LIT_E3_COV_REPRO, eps, ns),
                lit_subset(LIT_E3_COV_ADI,   eps, ns),
                lit_subset(LIT_E3_COV_NAIVE, eps, ns))
  wid_p <- list(E3$wid,
                lit_subset(LIT_E3_WID_REPRO, eps, ns),
                lit_subset(LIT_E3_WID_ADI,   eps, ns),
                lit_subset(LIT_E3_WID_NAIVE, eps, ns))
  ttl   <- list(M_PW, M_RP, M_ADI, M_NV)
  dec   <- c(3, 2, 2, 2)
  
  path <- file.path(DIR_FIG, "fig_exp3_comparison.pdf")
  pdf(path, width = 8.6, height = 6.4)
  on.exit({ while (dev.cur() > 1L) dev.off(); message("  [pdf]   ", path) },
          add = TRUE)
  
  ## coverage panels, their key, width panels, their key
  layout(matrix(c(1, 2, 3, 4,
                  5, 5, 5, 5,
                  6, 7, 8, 9,
                  10, 10, 10, 10), nrow = 4, byrow = TRUE),
         heights = c(1, lcm(1.75), 1, lcm(1.0)))
  
  for (i in 1:4)
    heat_panel(cov_p[[i]], gc3, main = ttl[[i]],
               xlab = LAB_N, ylab = if (i == 1L) LAB_EPS else NULL,
               dec = dec[i], mar = c(3.1, if (i == 1L) 3.9 else 2.3, 1.8, 0.6))
  colour_key(gc3, c(0, 0.25, 0.50, 0.75, 1 - ALPHA_3, 1), "Coverage",
             mar = c(1.9, 3.9, 0.2, 0.6))
  
  cols <- c("black", "grey40", "grey55", "grey72")
  ltys <- c(1, 2, 3, 4); pchs <- c(19, 15, 17, 18)
  yt   <- c(0.5, 1, 2, 5, 10, 20)
  for (i in 1:4) {
    W <- wid_p[[i]]
    par(mar = c(3.3, if (i == 1L) 4.1 else 2.4, 1.6, 0.7),
        mgp = c(2.0, 0.5, 0), tcl = -0.22, las = 1)
    plot(NA, xlim = range(nsn), ylim = c(0.35, 30), log = "xy",
         xlab = "", ylab = "", axes = FALSE)
    axis(1, at = nsn, labels = ns, cex.axis = 0.76)
    axis(2, at = yt, labels = yt, cex.axis = 0.76)
    box(col = "grey60", lwd = 0.6)
    mtext(LAB_N_S, side = 1, line = 1.9, cex = 0.80)
    if (i == 1L) mtext(LAB_CIW, side = 2, line = 2.7, cex = 0.80, las = 0)
    abline(h = yt, col = "grey93", lwd = 0.5)
    abline(h = E3_BOXW, lty = 2, col = "grey55", lwd = 1.1)
    for (j in seq_along(eps)) {
      keep <- !is.na(W[j, ])
      if (!any(keep)) next
      lines(nsn[keep], W[j, keep], col = cols[j], lty = ltys[j], lwd = 1.6)
      points(nsn[keep], W[j, keep], col = cols[j], pch = pchs[j], cex = 0.8)
    }
  }
  line_key(as.list(parse(text = c(sprintf("epsilon == %s", eps),
                                  '"search-box width"'))),
           c(cols[seq_along(eps)], "grey55"),
           c(ltys[seq_along(eps)], 2),
           c(pchs[seq_along(eps)], NA))
  
  sat <- sum(E3$wid >= E3_SAT, na.rm = TRUE)
  write_figure_wrapper("fig_exp3_comparison.pdf",
                       sprintf(paste0("Experiment 3. Coverage of the $90\\%%$ interval for ",
                                      "$\\beta_1^\\ast$ (first row) and its average width (second row), for ",
                                      "the four procedures. Rows of each heat map are the privacy budget and ",
                                      "columns the sample size; the grey scale breaks at the nominal level ",
                                      "$%.2f$, so cells below it are light. The finite-$R$ coverage bound is ",
                                      "$%.5f$. In the width panels the dashed line marks the total width ",
                                      "$%g$ of the search region $[%g,%g]$: an interval approaching it is ",
                                      "box-limited and therefore uninformative about efficiency, and we draw ",
                                      "no comparison from such a cell. This is a property of the interval ",
                                      "and not of the search: %d of %d of our cells are box-limited, and the ",
                                      "implementation does not substitute the search region when a search ",
                                      "fails to resolve. Comparison values are those reported by Wang, Chang ",
                                      "and Awan (2026, Fig.~6)."),
                               1 - ALPHA_3, BND_3, E3_BOXW, E3_BOX[1], E3_BOX[2],
                               sat, length(E3$wid)),
                       "fig:exp3-comparison", width = 0.98)
  invisible(path)
}

## Clamping sweep: size on the first row, power on the second, three methods
fig_delta_sensitivity <- function(D = read_delta()) {
  
  tab_need_init()
  if (is.null(D)) return(invisible(NULL))
  
  dl <- rownames(D$h0); ns <- colnames(D$h0)
  panels <- list(D$h0,
                 lit_subset(LIT_CLAMP_REPRO_H0, dl, ns),
                 lit_subset(LIT_CLAMP_PB_H0,    dl, ns),
                 D$h1,
                 lit_subset(LIT_CLAMP_REPRO_H1, dl, ns),
                 lit_subset(LIT_CLAMP_PB_H1,    dl, ns))
  ttl <- list(M_PW, M_RP, "Parametric bootstrap", "", "", "")
  
  fig_heat_grid(
    panels = panels, titles = ttl, greyfn = grey_reject,
    ticks = c(0, ALPHA_12, 0.25, 0.50, 0.75, 1),
    key_label = "Rejection probability",
    file = "fig_delta_sensitivity.pdf",
    xlab = LAB_N, ylab = LAB_DELTA,
    nrow = 2, ncol = 3, width = 7.8, height = 5.6,
    row_labels = list(expression("Size at " * beta[1]^"*" == 0),
                      expression("Power at " * beta[1]^"*" == 1)))
  
  write_figure_wrapper("fig_delta_sensitivity.pdf",
                       sprintf(paste0("Sensitivity to the clamping range. The first row is ",
                                      "empirical size at $\\beta_1^\\ast=0$ and the second is power at ",
                                      "$\\beta_1^\\ast=1$; rows of each heat map are the clamping bound and ",
                                      "columns the sample size. The grey scale breaks at $%.2f$, so in the ",
                                      "first row any dark cell is a violation of the nominal level, and the ",
                                      "exact finite-$R$ size bound is $%.5f$. A small clamping bound ",
                                      "discards information while a large one inflates the privacy noise, ",
                                      "whose scale is of order $\\Delta^2/n$ in the second-moment ",
                                      "coordinates, so power is attainable only in between. Comparison ",
                                      "values are those reported by Awan and Wang (2025, Fig.~6), ",
                                      "restricted to the cells of our design."), ALPHA_12, SZB_12),
                       "fig:delta-sensitivity", width = 0.95)
}


## ===========================================================================
## 6.  TABLE WRITERS
##
##  Each takes the object its reader returned, or returns silently if that
##  object is NULL.  The reference levels quoted in a caption are recomputed
##  from R and alpha, never written in, so a caption cannot quote 0.95 where
##  the applicable bound is 0.95025 or 0.90050.
## ===========================================================================

tab_exp1_main <- function(E = read_exp1(), T1 = read_exp1_targeted()) {
  
  if (is.null(E)) return(invisible(NULL))
  S <- E$S; lab <- E$method
  bnd <- BND_12; se <- mc_se(bnd, if (is.na(E$reps)) 1000 else E$reps)
  
  ## the reduction against Mahalanobis, matched to the row that targets it
  red <- rep("", nrow(S))
  if (!is.null(T1)) {
    tgt <- ifelse(T1$target == "mu", "\\mu", "\\sigma")
    for (j in seq_len(nrow(T1))) {
      hit <- grep(paste0("\\$", tgt[j], "\\$"), lab, fixed = FALSE)
      if (length(hit) == 1L) red[hit] <- fpct(T1$reduce[j], 1)
    }
  }
  
  body <- NULL
  for (i in seq_len(nrow(S))) {
    body <- rbind(body,
                  c(sprintf("\\multirow{2}{*}{%s}", lab[i]), "Coverage",
                    fse(S$coverage_mu[i],    S$coverage_mu_se[i]),
                    fse(S$coverage_sigma[i], S$coverage_sigma_se[i]),
                    fse(S$coverage_joint[i], S$coverage_joint_se[i]), "", ""),
                  c("", "Width",
                    fse(S$width_mu[i],    S$width_mu_se[i]),
                    fse(S$width_sigma[i], S$width_sigma_se[i]), "",
                    fse(S$area[i],        S$area_se[i]), red[i]))
  }
  
  fail <- if (all(S$failure_rate == 0))
    "Numerical failures were zero for every method." else
      paste0("Failure rates: ",
             paste(sprintf("%s %s", lab, fnum(S$failure_rate, 3)),
                   collapse = "; "), ".")
  
  inner <- latex_tabular(
    body   = body,
    header = c("Method", "", "$\\mu$", "$\\sigma$", "Joint", "Area",
               "Reduction"),
    align  = "llccccc",
    span   = "& & \\multicolumn{2}{c}{Marginal} & & & \\\\",
    cmid   = "\\cmidrule(lr){3-4}",
    rowsep = seq(2, nrow(body) - 2, by = 2),
    note   = paste("Notes: the two Penalized Wald rows are separate runs,",
                   "each orthogonalised against the other coordinate, so the off-target",
                   "column records the price of that orthogonalisation and the pair does",
                   "not carry simultaneous coverage. Area is reported for completeness",
                   "but is not the relevant criterion for a targeted statistic. The last",
                   "column gives the width of the targeted interval relative to the",
                   "Mahalanobis interval for the same parameter.", fail))
  
  write_table(inner,
              sprintf(paste0("Experiment 1. Empirical coverage and average width of the ",
                             "nominal $95\\%%$ intervals, and average area of the joint region, over ",
                             "$%s$ replications at $n=%s$. Monte Carlo standard errors in ",
                             "parentheses. Joint coverage is obtained by evaluating the acceptance ",
                             "rule at the true parameter, so it involves no nuisance optimisation ",
                             "and is a direct check on the implemented statistic; its reference ",
                             "value is the finite-$R$ bound $%.5f$, with a Monte Carlo standard ",
                             "error of $%.4f$."),
                      format(E$reps), format(E$n), bnd, se),
              "tab:exp1-main", "tab_exp1_main.tex")
}

tab_exp1_targeted <- function(T1 = read_exp1_targeted()) {
  
  if (is.null(T1)) return(invisible(NULL))
  lab <- ifelse(T1$target == "mu", "$\\mu$",
                ifelse(T1$target == "sigma", "$\\sigma$", tex_escape(T1$target)))
  body <- cbind(lab, fnum(T1$mah, 3), fnum(T1$pw, 3),
                fpct(T1$reduce, 1), fse(T1$cov, T1$cov_se))
  
  write_table(latex_tabular(body,
                            header = c("Target", "Mahalanobis width", "Penalized Wald width",
                                       "Reduction", "Coverage"),
                            align  = "ccccc"),
              sprintf(paste0("Experiment 1. Width of the targeted Penalized Wald ",
                             "interval relative to the Mahalanobis interval for the same parameter, ",
                             "with the coverage of the targeted interval. The finite-$R$ coverage ",
                             "bound is $%.5f$."), BND_12),
              "tab:exp1-targeted", "tab_exp1_targeted.tex")
}

tab_exp1_reproduction <- function(D = read_exp1_raw()) {
  
  if (is.null(D)) return(invisible(NULL))
  reps <- attr(D, "reps"); box <- attr(D, "box")
  total <- sum(D$full_mu) + sum(D$full_sig)
  
  note <- if (total == 0)
    sprintf(paste0("No interval was returned at a search-box boundary in any ",
                   "replication: the largest widths recorded are $%.3f$ for $\\mu$ and ",
                   "$%.3f$ for $\\sigma$, against box widths of $%.1f$ and $%.1f$. The ",
                   "conservative full-box convention was therefore never invoked, and no ",
                   "reported width is affected by it."),
            max(D$max_mu), max(D$max_sig), box[1], box[2]) else
              sprintf(paste0("The full-box convention was invoked %d time(s). Such a ",
                             "replication is not flagged as a failure, because a full-box interval ",
                             "has finite endpoints, so its width enters the average and its ",
                             "coverage indicator is one; the counts are given per method so the ",
                             "affected averages can be discounted."), total)
  
  body <- cbind(relabel_method(D$method),
                format(D$full_mu), format(D$full_sig),
                fnum(D$max_mu, 3), fnum(D$max_sig, 3), fnum(D$failure, 3))
  
  write_table(latex_tabular(body,
                            header = c("Method", "Full box $\\mu$", "Full box $\\sigma$",
                                       "Max width $\\mu$", "Max width $\\sigma$", "Failure rate"),
                            align  = "lccccc", note = note),
              sprintf(paste0("Experiment 1, recomputed from the replication-level file ",
                             "over $%d$ replications. The two full-box columns count replications ",
                             "whose interval equalled the search box, which is the convention the ",
                             "simulation uses when the nuisance search does not resolve."), reps),
              "tab:exp1-reproduction", "tab_exp1_reproduction.tex")
}

tab_exp2_power <- function(G = read_exp2_grid()) {
  
  if (is.null(G)) return(invisible(NULL))
  body <- cbind(rownames(G), matrix(fcell(G, 3), nrow = nrow(G)))
  write_table(latex_tabular(body,
                            header = c("$\\beta_1^\\ast$", paste0("$n=", colnames(G), "$")),
                            align  = paste(rep("c", 1 + ncol(G)), collapse = "")),
              sprintf(paste0("Experiment 2. Rejection probability of $H_0:\\beta_1=0$ ",
                             "at the $%.2f$ level. The first row is empirical size and should be ",
                             "read against the exact finite-$R$ bound $%.5f$; the remaining rows are ",
                             "power. Unresolved replications count as non-rejections, so this is the ",
                             "conservative power."), ALPHA_12, SZB_12),
              "tab:exp2-power", "tab_exp2_power.tex")
}

tab_exp3_grid <- function(E3 = read_exp3()) {
  
  if (is.null(E3)) return(invisible(NULL))
  eps <- rownames(E3$cov); ns <- colnames(E3$cov)
  rp  <- lit_subset(LIT_E3_COV_REPRO, eps, ns)
  adi <- lit_subset(LIT_E3_COV_ADI,   eps, ns)
  nv  <- lit_subset(LIT_E3_COV_NAIVE, eps, ns)
  wrp <- lit_subset(LIT_E3_WID_REPRO, eps, ns)
  wad <- lit_subset(LIT_E3_WID_ADI,   eps, ns)
  sat <- E3$wid >= E3_SAT
  
  body <- NULL
  for (i in seq_along(eps))
    for (j in seq_along(ns))
      body <- rbind(body, c(
        if (j == 1L) eps[i] else "", ns[j],
        fnum(E3$cov[i, j], 3), fnum(rp[i, j], 2),
        fnum(adi[i, j], 2), fnum(nv[i, j], 2),
        paste0(fnum(E3$wid[i, j], 3),
               if (isTRUE(sat[i, j])) "$^{\\dagger}$" else ""),
        fnum(wrp[i, j], 2), fnum(wad[i, j], 2)))
  
  write_table(latex_tabular(body,
                            header = c("$\\varepsilon$", "$n$",
                                       "\\textsc{pw}", "\\textsc{repro}", "\\textsc{pb-adi}",
                                       "\\textsc{pb-naive}",
                                       "\\textsc{pw}", "\\textsc{repro}", "\\textsc{pb-adi}"),
                            align  = "ccccccccc",
                            span   = paste0("& & \\multicolumn{4}{c}{Coverage} & ",
                                            "\\multicolumn{3}{c}{Width} \\\\"),
                            cmid   = "\\cmidrule(lr){3-6}\\cmidrule(lr){7-9}",
                            note   = sprintf(paste0("$^{\\dagger}$ width at or above $%.1f$, that ",
                                                    "is at least $90\\%%$ of the search box $[%g,%g]$: the interval is ",
                                                    "box-limited and its coverage follows from truncation rather than ",
                                                    "from the procedure. Comparison values are those reported by ",
                                                    "\\citet{wang2025optimal}."), E3_SAT, E3_BOX[1], E3_BOX[2])),
              sprintf(paste0("Experiment 3. Coverage and average width of the $90\\%%$ ",
                             "interval for $\\beta_1^\\ast$. The finite-$R$ coverage bound is ",
                             "$%.5f$, not $0.95$: this experiment runs at $\\alpha=%.2f$, and the ",
                             "Monte Carlo standard error at that level is $%.4f$."),
                      BND_3, ALPHA_3, mc_se(BND_3, 1000)),
              "tab:exp3-grid", "tab_exp3_grid.tex")
}

## One table per cloud-size sweep.  The penalty-weight sweep is presented as
## a figure only, so no table is written for it.
tab_sweeps <- function(SW = read_sweeps(), keys = c("R_aux", "R_synthetic")) {
  
  if (length(SW) == 0L) return(invisible(NULL))
  SYM <- c(lambda_n = "$\\lambda_n$", R_aux = "$R_{\\mathrm{aux}}$",
           R_synthetic = "$R$", Delta = "$\\Delta$",
           epsilon = "$\\varepsilon$")
  NAME <- c(lambda_n = "penalty weight", R_aux = "auxiliary cloud size",
            R_synthetic = "inference cloud size",
            Delta = "clamping range", epsilon = "privacy budget")
  
  for (key in intersect(keys, names(SW))) {
    S <- SW[[key]]
    two <- all(c("coverage_mu", "width_mu",
                 "coverage_sigma", "width_sigma") %in% names(S))
    if (!two) next
    
    body <- cbind(format(S[[key]], trim = TRUE),
                  fse(S$coverage_mu, S$coverage_mu_se),
                  fse(S$width_mu, S$width_mu_se),
                  fse(S$coverage_sigma, S$coverage_sigma_se),
                  fse(S$width_sigma, S$width_sigma_se))
    header <- c(SYM[[key]], "Coverage", "Width", "Coverage", "Width")
    align  <- "ccccc"
    cmid   <- "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}"
    span   <- paste0("& \\multicolumn{2}{c}{$\\mu$ target} & ",
                     "\\multicolumn{2}{c}{$\\sigma$ target} \\\\")
    
    if (all(c("failure_rate_mu", "failure_rate_sigma") %in% names(S))) {
      body   <- cbind(body, paste0(fnum(S$failure_rate_mu, 3), " / ",
                                   fnum(S$failure_rate_sigma, 3)))
      header <- c(header, "Failure $\\mu$ / $\\sigma$")
      align  <- paste0(align, "c")
      span   <- sub("\\\\\\\\$", "& \\\\\\\\", span)
    }
    ## the R sweep records the moving cut-off and bound
    moving <- all(c("rank_cutoff", "finite_R_reference") %in% names(S))
    if (moving) {
      body   <- cbind(body[, 1, drop = FALSE],
                      format(S$rank_cutoff, trim = TRUE),
                      fnum(S$finite_R_reference, 5),
                      body[, -1, drop = FALSE])
      header <- c(header[1], "$a_{R,\\alpha}$", "Bound", header[-1])
      align  <- paste0("ccc", substring(align, 2))
      cmid   <- "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}"
      span   <- paste0("& & & \\multicolumn{2}{c}{$\\mu$ target} & ",
                       "\\multicolumn{2}{c}{$\\sigma$ target}",
                       if (ncol(body) > 7) " &" else "", " \\\\")
    }
    
    bnd_txt <- if (moving)
      paste0("The acceptance rule depends on this constant, so the coverage ",
             "bound moves with it and is given in its own column.") else
               sprintf(paste0("The acceptance rule does not depend on this constant, ",
                              "so the coverage bound is $%.5f$ in every row and the ",
                              "sweep measures width alone."), BND_12)
    
    write_table(latex_tabular(body, header, align, span = span, cmid = cmid,
                              note = "Monte Carlo standard errors in parentheses."),
                paste0("Sensitivity to the ", NAME[[key]], " ", SYM[[key]],
                       ", $1000$ replications per grid point. ", bnd_txt),
                paste0("tab:sens-", key), paste0("tab_sens_", key, ".tex"))
  }
}

tab_delta <- function(D = read_delta()) {
  
  if (is.null(D)) return(invisible(NULL))
  dl <- rownames(D$h0); ns <- colnames(D$h0)
  body <- cbind(dl, matrix(fnum(D$h0, 3), nrow = nrow(D$h0)),
                matrix(fnum(D$h1, 3), nrow = nrow(D$h1)))
  k <- length(ns)
  
  write_table(latex_tabular(body,
                            header = c("$\\Delta$", paste0("$n{=}", ns, "$"),
                                       paste0("$n{=}", ns, "$")),
                            align  = paste(rep("c", 1 + 2 * k), collapse = ""),
                            span   = sprintf(paste0("& \\multicolumn{%d}{c}{Size at ",
                                                    "$\\beta_1^\\ast=0$} & \\multicolumn{%d}{c}{Power at ",
                                                    "$\\beta_1^\\ast=1$} \\\\"), k, k),
                            cmid   = sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}",
                                             1 + k, 2 + k, 1 + 2 * k),
                            note   = sprintf(paste0("Size is at or below the exact bound $%.5f$ in ",
                                                    "every cell."), SZB_12)),
              paste0("Sensitivity to the clamping range $\\Delta$ in Experiment 2. ",
                     "$\\Delta$ changes the model rather than the rule: a small ",
                     "$\\Delta$ discards information while a large one inflates the ",
                     "privacy noise, whose scale is of order $\\Delta^2/n$ in the ",
                     "second-moment coordinates."),
              "tab:sens-delta", "tab_sens_delta.tex")
}

tab_exp2_sigma <- function(G = read_sigma()) {
  
  if (is.null(G)) return(invisible(NULL))
  
  if (!is.null(G$summary)) {
    S <- G$summary
    g <- function(nm, d = 3) if (nm %in% names(S)) fnum(S[[nm]][1], d) else "--"
    gv <- function(nm) if (nm %in% names(S)) S[[nm]][1] else NA
    bnd <- if (!is.na(gv("R"))) cov_bound(gv("R"), ALPHA_12) else BND_12
    rows <- rbind(
      c("Sample size $n$",                      g("n", 0)),
      c("Privacy ($\\mu$-GDP)",                 g("GDP", 2)),
      c("Clamping range $\\Delta$",             g("Delta", 2)),
      c("Inference cloud $R$",                  g("R", 0)),
      c("Auxiliary cloud $R_{\\mathrm{aux}}$",  g("R_aux", 0)),
      c("Penalty weight $\\lambda_n$",          g("Lambda", 6)),
      c("Replications",                         g("N", 0)),
      c("Empirical coverage",
        paste0(g("Coverage", 4), "\\,(", g("MC_SE_Coverage", 4), ")")),
      c("Clopper--Pearson $95\\%$ interval",
        paste0("[", g("CP_Lower", 4), ", ", g("CP_Upper", 4), "]")),
      c("Finite-$R$ coverage bound",            fnum(bnd, 6)),
      c("Mean width",
        paste0(g("Mean_Width", 4), "\\,(", g("MC_SE_Width", 4), ")")),
      c("Median width",                         g("Median_Width", 4)),
      c("Width quartiles",
        paste0("[", g("Width_Q25", 4), ", ", g("Width_Q75", 4), "]")),
      c("Failure rate",                         g("Failure_Rate", 4)),
      c("Disconnected-set rate",                g("Disconnected_Rate", 4)))
    
    ## two figures from the companion studies, folded in so that the
    ## sigma-target evidence sits in a single table
    if (!is.null(G$audit) && "Overturn_Rate" %in% names(G$audit))
      rows <- rbind(rows, c("Search-audit overturn rate",
                            fnum(G$audit$Overturn_Rate[1], 4)))
    if (!is.null(G$levelset) && "Inclusion_At_Truth" %in% names(G$levelset))
      rows <- rbind(rows, c("Inclusion probability at the truth",
                            fnum(G$levelset$Inclusion_At_Truth[1], 4)))
    
    write_table(latex_tabular(rows, c("Quantity", "Value"), "ll",
                              note = paste("Replications whose accepted set on the $\\sigma$ grid",
                                           "is empty or has more than one connected run are scored as",
                                           "non-coverage and excluded from the width distribution; both rates",
                                           "are reported so the frequency of that convention is visible.",
                                           "Since the hull of a disconnected set would be a valid and",
                                           "conservative alternative, this convention costs coverage rather",
                                           "than inflating it. The overturn rate is the fraction of rejected",
                                           "$\\sigma$ values that a stronger nuisance search reverses, which",
                                           "is the only channel by which the numerical profiling can push",
                                           "coverage below its bound; the inclusion probability at the truth",
                                           "involves no search and so checks the statistic directly.")),
                paste0("Experiment 2 with $\\sigma$ as the target and ",
                       "$(\\beta_1,\\beta_0,\\mathbb{E}X,\\sqrt{\\operatorname{Var}X})$ ",
                       "profiled out: settings and finite-sample performance."),
                "tab:exp2-sigma-summary", "tab_exp2_sigma_summary.tex")
  }
  
}


tab_boundary <- function(B = read_boundary_grid()) {
  
  if (is.null(B) || is.null(B$joint)) return(invisible(NULL))
  mu <- rownames(B$joint$mah); eps <- colnames(B$joint$mah)
  body <- NULL
  for (i in seq_along(mu))
    for (j in seq_along(eps))
      body <- rbind(body, c(
        if (j == 1L) mu[i] else "", eps[j],
        fnum(B$joint$mah[i, j], 2), fnum(B$joint$eff[i, j], 2),
        fnum(B$joint$pb[i, j], 2),
        fnum(B$cov_mu$mah[i, j], 2), fnum(B$cov_mu$eff[i, j], 2),
        fnum(B$cov_mu$pb[i, j], 2),
        fnum(B$wid_mu$mah[i, j], 3), fnum(B$wid_mu$eff[i, j], 3),
        fnum(B$wid_mu$pb[i, j], 3)))
  
  write_table(latex_tabular(body,
                            header = c("$\\mu^\\ast$", "$\\varepsilon$",
                                       "\\textsc{mah}", "\\textsc{pw}", "\\textsc{adi}",
                                       "\\textsc{mah}", "\\textsc{pw}", "\\textsc{adi}",
                                       "\\textsc{mah}", "\\textsc{pw}", "\\textsc{adi}"),
                            align  = "ccccccccccc",
                            span   = paste0("& & \\multicolumn{3}{c}{Joint coverage} & ",
                                            "\\multicolumn{3}{c}{Coverage of $\\mu$} & ",
                                            "\\multicolumn{3}{c}{Width for $\\mu$} \\\\"),
                            cmid   = paste0("\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
                                            "\\cmidrule(lr){9-11}"),
                            note   = sprintf(paste0("Each cell uses $100$ replications, so the ",
                                                    "binomial standard error at a coverage of $0.95$ is $%.3f$ and ",
                                                    "differences below about $0.04$ should not be read as real. At ",
                                                    "$\\varepsilon=0.1$ the $\\sigma$-widths of both repro-family ",
                                                    "procedures are box-saturated, so that column carries no information ",
                                                    "about efficiency in the scale parameter."), mc_se(0.95, 100))),
              paste0("Clamping-boundary study: the true mean is moved towards the ",
                     "upper clamp at four privacy budgets, with $\\sigma^\\ast=1$ ",
                     "fixed. The coverage bound is $", fnum(BND_12, 5), "$."),
              "tab:pbadi-grid", "tab_pbadi_grid.tex")
}

tab_boundary_soft <- function(S = read_boundary_soft()) {
  
  if (is.null(S)) return(invisible(NULL))
  eps <- attr(S, "eps")
  lab <- c(mah = "Mahalanobis (\\textsc{repro})",
           eff = "Penalized Wald ($\\mu$)", pb = "\\textsc{pb-adi}")
  body <- NULL
  for (m in c("mah", "eff", "pb")) {
    H <- S[[paste0(m, ".hard")]]; So <- S[[paste0(m, ".soft")]]
    for (j in seq_along(eps))
      body <- rbind(body, c(
        if (j == 1L) lab[[m]] else "", format(eps[j]),
        fnum(H$joint[j], 2),    fnum(So$joint[j], 2),
        fnum(H$cov_mu[j], 2),   fnum(So$cov_mu[j], 2),
        fnum(H$width_mu[j], 3), fnum(So$width_mu[j], 3),
        fnum(H$area[j], 2),     fnum(So$area[j], 2)))
  }
  
  write_table(latex_tabular(body,
                            header = c("Method", "$\\varepsilon$", "hard", "soft", "hard", "soft",
                                       "hard", "soft", "hard", "soft"),
                            align  = "lccccccccc",
                            span   = paste0("& & \\multicolumn{2}{c}{Joint coverage} & ",
                                            "\\multicolumn{2}{c}{Coverage of $\\mu$} & ",
                                            "\\multicolumn{2}{c}{Width for $\\mu$} & ",
                                            "\\multicolumn{2}{c}{Area} \\\\"),
                            cmid   = paste0("\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
                                            "\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}"),
                            rowsep = c(length(eps), 2 * length(eps)),
                            note   = paste("The soft clamp is the logistic map into $(L,U)$ with",
                                           "$k=1.5$. Because the transformed data still lie in $[L,U]$ the",
                                           "sensitivities are unchanged, so the two regimes are compared at",
                                           "exactly the same privacy cost, and the clamp type is applied to the",
                                           "release, the region grid, the indirect estimator and the resampler",
                                           "alike.")),
              paste0("Hard against soft clamping at $\\mu^\\ast=3$, $100$ replications ",
                     "per cell."),
              "tab:pbadi-soft", "tab_pbadi_soft.tex")
}


## ===========================================================================
## 7.  FIGURE WRITERS THAT READ
##
##  The geometry is in Section 5; these supply the data and the captions.
## ===========================================================================

fig_exp2_comparison <- function(G = read_exp2_grid()) {
  
  if (is.null(G)) return(invisible(NULL))
  b <- rownames(G); ns <- colnames(G)
  
  fig_heat_grid(
    panels = list(G,
                  lit_subset(LIT_E2_REPRO, b, ns),
                  lit_subset(LIT_E2_ADI,   b, ns),
                  lit_subset(LIT_E2_NAIVE, b, ns)),
    titles = list(M_PW, M_RP, M_ADI, M_NV),
    greyfn = grey_reject,
    ticks  = c(0, ALPHA_12, 0.25, 0.50, 0.75, 1),
    key_label = "Rejection probability",
    file = "fig_exp2_grid.pdf",
    xlab = LAB_N, ylab = LAB_BETA1,
    nrow = 2, ncol = 2, width = 6.9, height = 5.7)
  
  write_figure_wrapper("fig_exp2_grid.pdf",
                       sprintf(paste0("Experiment 2. Rejection probability of $H_0:\\beta_1=0$ ",
                                      "at the $%.2f$ level. Rows are the true slope, columns the sample size. ",
                                      "The grey scale runs from $0$ to $1$ and breaks at $%.2f$, so cells at ",
                                      "or below the nominal level are light and cells above it dark: the ",
                                      "bottom row of each panel is empirical size and the rest is power. ",
                                      "\\textsc{repro} is the procedure of \\citet{awan2025simulation}; ",
                                      "\\textsc{pb-adi} and \\textsc{pb-naive} are the debiased and naive ",
                                      "parametric bootstraps as reported by \\citet{wang2025optimal}."),
                               ALPHA_12, ALPHA_12),
                       "fig:exp2-grid", width = 0.95)
  invisible(NULL)
}


fig_lambda_sensitivity <- function(SW = read_sweeps()) {
  
  L <- SW[["lambda_n"]]
  if (is.null(L)) return(invisible(NULL))
  fig_lambda_plot(L)
  
  sel <- if ("selected" %in% names(L)) as.logical(unquote(L$selected)) else
    seq_len(nrow(L)) == which.min(L$width_mu / min(L$width_mu) +
                                    L$width_sigma / min(L$width_sigma))
  write_figure_wrapper("fig_sens_lambda.pdf",
                       sprintf(paste0("Sensitivity to the penalty weight $\\lambda_n$: average ",
                                      "interval width for each target. Grid points are drawn at equal ",
                                      "ordinal spacing because the design grid is geometric; the dotted line ",
                                      "marks the selected value $\\lambda_n=%g$. The $\\mu$-target width ",
                                      "falls from $%.3f$ at $\\lambda_n=%g$ to $%.3f$ at the selected value ",
                                      "and then rises slowly. Coverage clears the Monte-Carlo-adjusted floor ",
                                      "at every grid point, so the figure concerns width alone."),
                               L$lambda_n[sel][1], L$width_mu[1], L$lambda_n[1], L$width_mu[sel][1]),
                       "fig:sens-lambda", width = 0.76)
  invisible(NULL)
}

fig_boundary_joint <- function(B = read_boundary_grid()) {
  
  if (!is.null(B) && !is.null(B$joint)) {
    fig_heat_grid(
      panels = list(B$joint$mah, B$joint$eff, B$joint$pb),
      titles = list(M_MAH, M_PW, M_ADI),
      greyfn = grey_cov,
      ticks = c(0, 0.25, 0.50, 0.75, 0.95),
      key_label = "Joint coverage",
      file = "fig_pbadi_joint.pdf",
      xlab = LAB_EPS, ylab = LAB_MUSTAR,
      nrow = 1, ncol = 3, dec = rep(2, 3),
      width = 7.8, height = 3.2)
    
    write_figure_wrapper("fig_pbadi_joint.pdf",
                         sprintf(paste0("Joint coverage of the nominal $95\\%%$ region for ",
                                        "$(\\mu,\\sigma)$ as the true mean approaches the clamping bound. ",
                                        "Rows are $\\mu^\\ast$, columns the privacy budget. The grey scale ",
                                        "breaks at $0.95$, so light cells are under-coverage. Each cell uses ",
                                        "$100$ replications, so the Monte Carlo standard error near $0.95$ is ",
                                        "$%.3f$."), mc_se(0.95, 100)),
                         "fig:pbadi-joint", width = 0.95)
  }
  invisible(NULL)
}


fig_boundary_soft <- function(S = read_boundary_soft()) {
  
  tab_need_init()
  if (is.null(S)) return(invisible(NULL))
  fig_boundary_soft_plot(S)
  write_figure_wrapper("fig_pbadi_soft.pdf",
                       sprintf(paste0("Joint coverage at $\\mu^\\ast=3$, the cell in which the ",
                                      "clamp binds on half the sample, under hard (solid) and soft (dashed) ",
                                      "clamping. Because the smooth map sends the data into $(L,U)$ the ",
                                      "replace-one sensitivities are unchanged, so the two regimes carry the ",
                                      "same privacy cost and the only difference between them is the ",
                                      "differentiability of the clamp. The dashed horizontal line is the ",
                                      "finite-$R$ bound $%.5f$. Each cell uses $100$ replications, so the ",
                                      "Monte Carlo standard error near $0.95$ is $%.3f$."),
                               BND_12, mc_se(0.95, 100)),
                       "fig:boundary-soft", width = 0.78)
  invisible(NULL)
}


## ===========================================================================
## 8.  DRIVER
## ===========================================================================

## One float per distinct set of numbers.  Writers that are defined but not
## registered here produce material that duplicates a registered float or
## belongs in an appendix; enable one by adding its name to a list below.
##
##   defined, not registered:
##     tab_exp1_targeted        folded into tab_exp1_main
##     tab_exp1_reproduction    appendix diagnostic on the replication file
##     tab_exp2_power           the same matrix as fig_exp2_comparison
##     tab_exp3_grid            the same matrices as fig_exp3_comparison
##     tab_delta                the same matrices as fig_delta_sensitivity
##     tab_boundary             the same matrix as fig_boundary_joint
##     tab_boundary_soft        the same matrix as fig_boundary_soft
##     fig_exp1_extra           single-draw geometry from the simulation
##     fig_sigma_extra          sigma-target performance and inclusion surface
##     fig_sweeps_clouds        the cloud sweeps, which are tabulated instead

TABLE_WRITERS <- c(
  "tab_exp1_main",
  "tab_exp2_sigma",
  "tab_sweeps"          # writes tab_sens_R_aux and tab_sens_R_synthetic
)

FIGURE_WRITERS <- c(
  "fig_exp2_comparison",      # Penalized Wald, Repro, PB-ADI, PB-naive
  "fig_exp2_sigma_geometry",  # geometry of the accepted set
  "fig_exp3_comparison",      # coverage and width in one figure
  "fig_lambda_sensitivity",
  "fig_delta_sensitivity",    # size and power in one figure
  "fig_boundary_joint",
  "fig_boundary_soft"
)

## Produces everything.  A writer that stops with an error is reported and
## the run continues, so one unreadable input cannot cost the whole set; any
## open graphics device is closed on the way out.
make_all <- function(what = c("both", "tables", "figures"), root = NULL) {
  
  what <- match.arg(what)
  tab_init(root)
  
  todo <- switch(what, both = c(TABLE_WRITERS, FIGURE_WRITERS),
                 tables = TABLE_WRITERS, figures = FIGURE_WRITERS)
  failed <- character(0)
  
  for (f in todo) {
    message("\n---- ", f, "()")
    ok <- tryCatch({ do.call(f, list()); TRUE },
                   error = function(e) {
                     message("  [ERROR] ", conditionMessage(e)); FALSE })
    if (!ok) {
      failed <- c(failed, f)
      while (dev.cur() > 1L) dev.off()
    }
  }
  
  message("\n===========================================================")
  message("  tables  ", normalizePath(DIR_TAB, mustWork = FALSE), ": ",
          length(list.files(DIR_TAB, pattern = "\\.tex$")))
  message("  figures ", normalizePath(DIR_FIG, mustWork = FALSE), ": ",
          length(list.files(DIR_FIG, pattern = "\\.pdf$")))
  if (length(failed))
    message("  did not complete: ", paste(failed, collapse = ", "))
  message("===========================================================")
  invisible(failed)
}

## ---------------------------------------------------------------------------
## Autorun.  Sourcing the file produces everything; set
##     options(tabulate.autorun = FALSE)
## beforehand to load the functions without writing anything.
## ---------------------------------------------------------------------------

if (isTRUE(getOption("tabulate.autorun", default = TRUE))) make_all()
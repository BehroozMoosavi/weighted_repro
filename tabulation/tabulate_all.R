## ===========================================================================
##  results_tables_and_plots.R
##
##  Produces every table and figure of the numerical study.
##
##  Self-contained: all values are written into section 2 of this file, so the
##  script runs on its own and does not read the simulation output.  Values
##  for the Penalized Wald procedure are our own runs; the Repro and
##  parametric bootstrap grids are the published ones and are attributed in
##  section 2 and in every caption where they appear.
##
##  Base R only.  No packages.
##
##  Usage
##  -----
##      source("results_tables_and_plots.R")
##
##  or
##
##      Rscript results_tables_and_plots.R
##
##  Writes, relative to the working directory:
##
##      tables/   *.tex   booktabs fragments, each a complete float with a
##                        caption and a label, ready to \input
##      figures/  *.pdf   vector figures
##
##  Figure conventions follow Awan and Wang (2025): grey-scale cells with the
##  value printed inside, and a two-segment colour ramp whose break sits at
##  the nominal level, so the significance or coverage threshold is legible
##  without consulting the key.
##
##  Contents
##  --------
##      1  configuration, reference levels, helpers
##      2  data
##      3  table writers
##      4  figure writers
##      5  driver
## ===========================================================================


## ===========================================================================
## 1.  CONFIGURATION, REFERENCE LEVELS, HELPERS
## ===========================================================================

DIR_TAB <- "tables"
DIR_FIG <- "figures"
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## Reference levels.  The exact rank cut-off is a = floor(alpha (R+1)) + 1;
## the size bound is (a-1)/(R+1) and the coverage bound its complement.  The
## bound is an upper bound on size and a lower bound on coverage, attained
## when ties in the depth have probability zero.
## ---------------------------------------------------------------------------

rank_cutoff <- function(R, alpha) floor(alpha * (R + 1)) + 1
size_bound  <- function(R, alpha) (rank_cutoff(R, alpha) - 1) / (R + 1)
cov_bound   <- function(R, alpha) 1 - size_bound(R, alpha)
mc_se       <- function(p, n) sqrt(p * (1 - p) / n)

ALPHA_12 <- 0.05   # Experiments 1 and 2, and their sweeps
ALPHA_3  <- 0.10   # Experiment 3 runs at 90 per cent confidence
R_SYN    <- 200

BND_12   <- cov_bound(R_SYN, ALPHA_12)    # 0.95025
SZB_12   <- size_bound(R_SYN, ALPHA_12)   # 0.04975
BND_3    <- cov_bound(R_SYN, ALPHA_3)     # 0.90050

E3_BOX   <- c(-10, 10)                    # search box for beta1
E3_BOXW  <- diff(E3_BOX)
E3_SAT   <- 0.90 * E3_BOXW                # width above which an interval is
# box-limited and uninformative

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

## ---------------------------------------------------------------------------
## Matrix constructor used by the data section
## ---------------------------------------------------------------------------

mk <- function(v, rows, cols)
  matrix(v, nrow = length(rows), ncol = length(cols), byrow = TRUE,
         dimnames = list(rows, cols))


## ===========================================================================
## 2.  DATA
##
##  Our own runs
##  ------------
##  Experiment 1 coverage, width, area and targeted reduction
##      exp1_location_scale_normal/exp1_coverage.R
##  Experiment 1 sweeps in lambda, R_aux and R
##      sensitivity_analysis/exp1_normal_sensitivity/
##  Experiment 2 size and power
##      exp2_linear_regression/linearreg_eff.R
##  Experiment 2 clamping sweep
##      sensitivity_analysis/exp2_linear_regression_sensitivity/
##  Experiment 3 coverage and width
##      exp3_objective_perturbation/eff.R
##  Clamping-boundary studies
##      PB_adi_failure/failure_v1.R and failure_v2.R
##
##  Published comparison values
##  ---------------------------
##  E2_REPRO            Awan and Wang (2025), Fig. 5, "Repro Sample"
##  E2_ADI              Wang, Chang and Awan, Fig. 5, indirect + pivot
##  E2_NAIVE            Wang, Chang and Awan, Fig. 5, naive + F
##  DEL_REPRO_*         Awan and Wang (2025), Fig. 6, Repro panels
##  DEL_PB_*            Awan and Wang (2025), Fig. 6, bootstrap panels
##  E3_*_REPRO          Wang, Chang and Awan, Fig. 6, Repro Samples
##  E3_*_ADI            Wang, Chang and Awan, Fig. 6, Debiased Estimator
##  E3_*_NAIVE          Wang, Chang and Awan, Fig. 6, Naive Estimator
##
##  Published grids are shown only on the cells of our own design, except
##  DEL_*_FULL, which are reproduced in full for context.  Rows are the first
##  factor and columns the second in every matrix; dimnames carry the keys, so
##  a subset by name cannot silently misalign.
## ===========================================================================

## --- Experiment 1: our run ------------------------------------------------
EXP1 <- data.frame(
  method = c("Mahalanobis", "Penalized Wald (mu)",
             "Penalized Wald (sigma)", "PB-ADI"),
  cov_mu     = c(0.991, 0.966, 1, 0.959),
  cov_mu_se  = c(0.002986, 0.005731, 0, 0.00627),
  cov_sig    = c(0.985, 1, 0.964, 0.96),
  cov_sig_se = c(0.003844, 0, 0.005891, 0.006197),
  cov_jt     = c(0.963, 0.956, 0.956, 0.967),
  cov_jt_se  = c(0.005969, 0.006486, 0.006486, 0.005649),
  w_mu       = c(0.599854, 0.504201, 1.07633, 0.462868),
  w_mu_se    = c(0.003245, 0.002664, 0.004628, 0.002585),
  w_sig      = c(0.757885, 1.73005, 0.613595, 0.576474),
  w_sig_se   = c(0.004228, 0.009236, 0.003437, 0.003291),
  area       = c(0.341496, 0.65695, 0.502823, 0.337577),
  area_se    = c(0.00366, 0.007249, 0.004925, 0.003881),
  stringsAsFactors = FALSE)

EXP1_TGT <- data.frame(
  target  = c("mu", "sigma"),
  mah     = c(0.599854, 0.757885),
  pw      = c(0.504201, 0.613595),
  reduce  = c(15.946, 19.039),
  cov     = c(0.966, 0.964),
  cov_se  = c(0.005731, 0.005891),
  stringsAsFactors = FALSE)

## --- Experiment 1: sweeps -------------------------------------------------
LAMBDA <- data.frame(
  lambda  = c(0, 0.1, 0.5, 1, 2, 4, 8),
  cov_mu     = c(0.994, 0.98, 0.979, 0.981, 0.986, 0.983, 0.977),
  cov_mu_se  = c(0.002442, 0.004427, 0.004534, 0.004317, 0.003715, 0.004088, 0.00474),
  w_mu       = c(2.21253, 0.592582, 0.58722, 0.604612, 0.623857, 0.641953, 0.656998),
  w_mu_se    = c(0.039306, 0.002798, 0.002921, 0.002982, 0.003047, 0.003121, 0.00316),
  cov_sig    = c(0.989, 0.976, 0.975, 0.977, 0.978, 0.979, 0.982),
  cov_sig_se = c(0.003298, 0.00484, 0.004937, 0.00474, 0.004639, 0.004534, 0.004204),
  w_sig      = c(0.851937, 0.663229, 0.688129, 0.712397, 0.745553, 0.773708, 0.791345),
  w_sig_se   = c(0.006358, 0.003486, 0.003591, 0.003622, 0.003848, 0.003943, 0.00402),
  fail_mu    = c(0.001, 0.01, 0.01, 0.009, 0.007, 0.011, 0.018),
  fail_sig   = c(0, 0.008, 0.01, 0.011, 0.011, 0.011, 0.015),
  obj        = c(2.526, 1.005, 1.019, 1.052, 1.093, 1.13, 1.156),
  selected = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)
LAMBDA_EVAL <- list(lambda = 0.1, cov_mu = 0.971, cov_mu_se = 0.005307, w_mu = 0.596281,
                    w_mu_se = 0.003022, cov_sig = 0.979, cov_sig_se = 0.004534,
                    w_sig = 0.665381, w_sig_se = 0.003811, fail_mu = 0.011, fail_sig = 0.004)

RAUX <- data.frame(
  R_aux   = c(50, 100, 200, 400, 800),
  cov_mu     = c(0.979, 0.979, 0.978, 0.981, 0.983),
  cov_mu_se  = c(0.004534, 0.004534, 0.004639, 0.004317, 0.004088),
  w_mu       = c(0.581795, 0.58221, 0.582468, 0.581742, 0.581491),
  w_mu_se    = c(0.002826, 0.002851, 0.002851, 0.002847, 0.002859),
  cov_sig    = c(0.976, 0.976, 0.977, 0.976, 0.974),
  cov_sig_se = c(0.00484, 0.00484, 0.00474, 0.00484, 0.005032),
  w_sig      = c(0.668259, 0.669025, 0.669722, 0.66904, 0.669032),
  w_sig_se   = c(0.00349, 0.003537, 0.003543, 0.003532, 0.00354),
  fail_mu    = c(0.005, 0.006, 0.01, 0.005, 0.005),
  fail_sig   = c(0.008, 0.005, 0.006, 0.007, 0.005),
  stringsAsFactors = FALSE)

RPC <- data.frame(
  R       = c(50, 100, 200, 400, 800),
  cutoff  = c(3, 6, 11, 21, 41),
  bound   = c(0.960784, 0.950495, 0.950249, 0.950125, 0.950062),
  cov_mu     = c(0.984, 0.979, 0.981, 0.975, 0.975),
  cov_mu_se  = c(0.003968, 0.004534, 0.004317, 0.004937, 0.004937),
  w_mu       = c(0.641526, 0.587342, 0.581742, 0.57764, 0.574313),
  w_mu_se    = c(0.004673, 0.003047, 0.002847, 0.002713, 0.002571),
  cov_sig    = c(0.986, 0.98, 0.976, 0.974, 0.976),
  cov_sig_se = c(0.003715, 0.004427, 0.00484, 0.005032, 0.00484),
  w_sig      = c(0.738688, 0.679219, 0.66904, 0.662394, 0.662361),
  w_sig_se   = c(0.00493, 0.003903, 0.003532, 0.003255, 0.003242),
  fail_mu    = c(0.003, 0.006, 0.005, 0.012, 0.014),
  fail_sig   = c(0.002, 0.001, 0.007, 0.011, 0.008),
  stringsAsFactors = FALSE)

## --- Experiment 2 ---------------------------------------------------------
E2_N <- c(100, 200, 300, 400, 500, 1000)
E2_B <- c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")
E2_OURS <- mk(c(0.028, 0.03, 0.032, 0.026, 0.041, 0.028, 0.108, 0.439, 0.751, 0.923, 0.989, 1, 0.399, 0.958, 0.999, 1, 1, 1, 0.731, 0.999, 1, 1, 1, 1, 0.902, 1, 1, 1, 1, 1, 0.97, 1, 1, 1, 1, 1),
              rows = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cols = c("100", "200", "300", "400", "500", "1000"))
E2_REPRO <- mk(c(0, 0, 0.001, 0, 0, 0.001, 0.007, 0.086, 0.278, 0.562, 0.836, 1, 0.074, 0.666, 0.984, 1, 1, 1, 0.281, 0.978, 1, 1, 1, 1, 0.566, 0.999, 1, 1, 1, 1, 0.762, 1, 1, 1, 1, 1),
               rows = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cols = c("100", "200", "300", "400", "500", "1000"))
E2_ADI <- mk(c(0.06, 0.053, 0.052, 0.033, 0.051, 0.043, 0.151, 0.508, 0.793, 0.942, 0.99, 1, 0.492, 0.967, 1, 1, 1, 1, 0.779, 0.999, 1, 1, 1, 1, 0.924, 1, 1, 1, 1, 1, 0.966, 1, 1, 1, 1, 1),
             rows = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cols = c("100", "200", "300", "400", "500", "1000"))
E2_NAIVE <- mk(c(0.002, 0.02, 0.035, 0.035, 0.048, 0.04, 0.018, 0.39, 0.734, 0.904, 0.988, 1, 0.164, 0.913, 0.997, 1, 1, 1, 0.468, 0.897, 0.984, 0.997, 1, 1, 0.594, 0.843, 0.951, 0.988, 0.993, 1, 0.63, 0.78, 0.898, 0.961, 0.982, 1),
               rows = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"), cols = c("100", "200", "300", "400", "500", "1000"))
E2_N8 <- c(100, 200, 300, 400, 500, 1000, 2000, 5000)
T1_OURS <- c(0.028, 0.03, 0.032, 0.026, 0.041, 0.028, NA, NA)
T1_REPRO <- c(0, 0, 0.001, 0, 0, 0.001, 0, 0.002)
T1_ADI <- c(0.06, 0.053, 0.052, 0.033, 0.051, 0.043, 0.042, 0.042)
T1_NAIVE <- c(0.002, 0.02, 0.035, 0.035, 0.048, 0.04, 0.077, 0.117)

## --- Experiment 2: clamping sweep ----------------------------------------
D_DELTA <- c("0.5", "1", "2", "5", "10")
D_N <- c(100, 200, 500, 1000)
DEL_OURS_H0 <- mk(c(0.01, 0.019, 0.036, 0.02, 0.015, 0.021, 0.033, 0.037, 0.026, 0.027, 0.036, 0.041, 0.028, 0.022, 0.034, 0.031, 0.021, 0.017, 0.029, 0.029),
                  rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DEL_OURS_H1 <- mk(c(1, 1, 1, 1, 1, 1, 1, 1, 0.968, 1, 1, 1, 0.091, 0.328, 0.984, 1, 0.023, 0.035, 0.152, 0.529),
                  rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DEL_REPRO_H0 <- mk(c(0, 0, 0, 0, 0, 0.001, 0.001, 0.001, 0, 0, 0, 0.001, 0.001, 0, 0.001, 0.001, 0.003, 0.001, 0.003, 0),
                   rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DEL_REPRO_H1 <- mk(c(1, 1, 1, 1, 0.998, 1, 1, 1, 0.762, 1, 1, 1, 0.005, 0.079, 0.824, 1, 0.001, 0.001, 0.018, 0.128),
                   rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DEL_PB_H0 <- mk(c(0.007, 0.018, 0.029, 0.015, 0.017, 0.045, 0.118, 0.186, 0.002, 0.019, 0.047, 0.037, 0, 0, 0.007, 0.035, 0, 0, 0, 0),
                rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DEL_PB_H1 <- mk(c(0.987, 1, 1, 1, 0.823, 0.957, 1, 1, 0.629, 0.78, 0.982, 1, 0, 0.03, 0.672, 0.829, 0, 0, 0, 0.124),
                rows = c("0.5", "1", "2", "5", "10"), cols = c("100", "200", "500", "1000"))
DF_DELTA <- c("0.5", "0.8", "1.0", "1.5", "2.0", "5.0", "10.0")
DF_N <- c(100, 200, 300, 400, 500, 1000, 2000, 5000)
DEL_PB_H0_FULL <- mk(c(0.007, 0.018, 0.019, 0.02, 0.029, 0.015, 0.009, 0.005, 0.014, 0.037, 0.053, 0.075, 0.089, 0.114, 0.199, 0.357, 0.017, 0.045, 0.068, 0.107, 0.118, 0.186, 0.361, 0.674, 0.011, 0.034, 0.047, 0.071, 0.083, 0.13, 0.236, 0.475, 0.002, 0.019, 0.035, 0.036, 0.047, 0.037, 0.078, 0.107, 0, 0, 0, 0.004, 0.007, 0.035, 0.048, 0.043, 0, 0, 0, 0, 0, 0, 0.004, 0.042),
                     rows = c("0.5", "0.8", "1.0", "1.5", "2.0", "5.0", "10.0"), cols = c("100", "200", "300", "400", "500", "1000", "2000", "5000"))
DEL_REPRO_H0_FULL <- mk(c(0, 0, 0, 0, 0, 0, 0, 0.001, 0, 0, 0, 0, 0.001, 0.001, 0, 0.002, 0, 0.001, 0, 0, 0.001, 0.001, 0.002, 0.002, 0, 0, 0.002, 0, 0.001, 0.001, 0.001, 0.002, 0, 0, 0.001, 0, 0, 0.001, 0, 0.002, 0.001, 0, 0, 0, 0.001, 0.001, 0.001, 0, 0.003, 0.001, 0, 0, 0.003, 0, 0.001, 0),
                        rows = c("0.5", "0.8", "1.0", "1.5", "2.0", "5.0", "10.0"), cols = c("100", "200", "300", "400", "500", "1000", "2000", "5000"))

## --- Experiment 3 ---------------------------------------------------------
E3_N <- c(100, 200, 500, 1000)
E3_EPS <- c("0.1", "0.3", "1", "3")
E3_COV_OURS <- mk(c(0.999, 0.999, 1, 0.999, 1, 0.999, 0.998, 0.98, 1, 0.997, 0.983, 0.946, 0.993, 0.978, 0.945, 0.937),
                  rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_WID_OURS <- mk(c(19.967, 19.963, 19.86, 18.243, 19.912, 19.585, 13.751, 3.796, 17.366, 9.658, 1.648, 0.739, 5.151, 1.872, 0.885, 0.575),
                  rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_COV_REPRO <- mk(c(0.94, 0.95, 0.97, 0.98, 0.96, 0.96, 0.97, 0.98, 0.96, 0.97, 0.98, 0.98, 0.97, 0.97, 0.98, 0.99),
                   rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_WID_REPRO <- mk(c(18.59, 18.75, 18.05, 15.48, 18.42, 15.86, 10.21, 3.54, 11.4, 6.89, 1.76, 0.98, 4.45, 2.01, 0.96, 0.76),
                   rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_COV_ADI <- mk(c(0.99, 0.96, 0.98, 0.97, 0.95, 0.95, 0.94, 0.91, 0.91, 0.94, 0.89, 0.89, 0.95, 0.89, 0.9, 0.87),
                 rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_WID_ADI <- mk(c(9.51, 10.14, 11.09, 8.76, 6.5, 6.06, 3.89, 1.59, 3.15, 2.37, 1.06, 0.6, 1.74, 1.29, 0.73, 0.51),
                 rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_COV_NAIVE <- mk(c(0.38, 0.3, 0.2, 0.15, 0.26, 0.21, 0.16, 0.16, 0.25, 0.28, 0.33, 0.45, 0.58, 0.62, 0.69, 0.69),
                   rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))
E3_WID_NAIVE <- mk(c(3.51, 2.85, 1.67, 0.97, 2.38, 1.49, 0.77, 0.47, 1.21, 0.79, 0.46, 0.32, 0.95, 0.69, 0.44, 0.32),
                   rows = c("0.1", "0.3", "1", "3"), cols = c("100", "200", "500", "1000"))

## --- Clamping-boundary studies -------------------------------------------
PB_MU  <- c("1", "1.5", "2", "2.5", "3")
PB_EPS <- c("0.1", "0.2", "0.5", "1")
G1_JOINT <- list(
  mah = mk(c(0.94, 0.94, 0.95, 0.98, 0.96, 0.95, 0.97, 0.98, 0.96, 0.96, 0.98, 0.98, 0.96, 0.96, 0.97, 0.97, 0.97, 0.96, 0.97, 0.97), rows = PB_MU, cols = PB_EPS),
  eff = mk(c(0.94, 0.95, 0.97, 0.99, 0.96, 0.97, 0.99, 0.98, 0.99, 0.97, 0.97, 0.96, 0.97, 0.97, 0.97, 0.92, 0.95, 0.94, 0.9, 0.87), rows = PB_MU, cols = PB_EPS),
  pb = mk(c(0.54, 0.81, 0.93, 0.94, 0.6, 0.9, 0.96, 0.94, 0.57, 0.83, 0.94, 0.93, 0.36, 0.61, 0.87, 0.9, 0.19, 0.4, 0.73, 0.82), rows = PB_MU, cols = PB_EPS))
G1_COVMU <- list(
  mah = mk(c(0.98, 0.99, 0.99, 1, 0.98, 0.99, 0.99, 1, 1, 0.99, 0.99, 0.99, 1, 1, 1, 1, 1, 1, 1, 1), rows = PB_MU, cols = PB_EPS),
  eff = mk(c(0.97, 0.95, 0.99, 0.99, 0.96, 0.99, 0.99, 0.99, 0.99, 0.98, 0.99, 0.99, 0.98, 0.98, 0.98, 0.96, 0.96, 0.95, 0.93, 0.88), rows = PB_MU, cols = PB_EPS),
  pb = mk(c(0.82, 0.93, 0.95, 0.96, 0.92, 0.99, 0.97, 0.97, 0.91, 0.97, 0.95, 0.97, 0.74, 0.79, 0.88, 0.93, 0.55, 0.49, 0.73, 0.82), rows = PB_MU, cols = PB_EPS))
G1_WIDMU <- list(
  mah = mk(c(4.313, 2.788, 0.928, 0.673, 4.892, 2.385, 0.756, 0.61, 4.301, 2.817, 0.992, 0.691, 3.467, 2.757, 1.831, 0.968, 3.09, 2.456, 2.176, 1.77), rows = PB_MU, cols = PB_EPS),
  eff = mk(c(3.972, 1.96, 0.758, 0.568, 4.417, 1.649, 0.633, 0.522, 3.932, 2.039, 0.806, 0.59, 3.097, 2.503, 1.365, 0.795, 2.574, 2.295, 1.729, 1.443), rows = PB_MU, cols = PB_EPS),
  pb = mk(c(3.988, 2.157, 0.657, 0.467, 5.317, 1.75, 0.531, 0.437, 6.825, 3.299, 0.677, 0.478, 7.712, 5.769, 1.677, 0.678, 7.823, 7.06, 4.471, 1.771), rows = PB_MU, cols = PB_EPS))
G1_WIDSIG <- list(
  mah = mk(c(4.736, 4.187, 1.56, 0.827, 4.728, 3.947, 1.353, 0.75, 4.738, 4.2, 1.623, 0.853, 4.751, 4.475, 2.612, 1.184, 4.646, 4.246, 3.172, 2.081), rows = PB_MU, cols = PB_EPS),
  eff = mk(c(4.879, 4.827, 3.606, 1.783, 4.9, 4.89, 3.777, 1.712, 4.9, 4.88, 3.554, 1.806, 4.896, 4.778, 3.291, 2.069, 4.881, 4.533, 3.586, 2.493), rows = PB_MU, cols = PB_EPS),
  pb = mk(c(3.548, 2.381, 1.132, 0.584, 3.491, 2.179, 0.984, 0.529, 3.766, 2.328, 1.166, 0.597, 4.31, 2.979, 1.662, 0.837, 3.937, 3.572, 2.277, 1.539), rows = PB_MU, cols = PB_EPS))

S2 <- list(
  "mah.hard" = list(joint = c(0.97, 0.96, 0.97, 0.97), cov_mu = c(1, 1, 1, 1),
                    width_mu = c(3.09, 2.456, 2.176, 1.77), area = c(10.753, 5.62, 2.446, 1.332)),
  "mah.soft" = list(joint = c(0.97, 0.97, 0.96, 0.96), cov_mu = c(1, 1, 1, 1),
                    width_mu = c(3.312, 2.8, 2.205, 1.412), area = c(11.962, 6.707, 2.52, 0.97)),
  "eff.hard" = list(joint = c(0.95, 0.94, 0.9, 0.87), cov_mu = c(0.96, 0.95, 0.93, 0.88),
                    width_mu = c(2.574, 2.295, 1.729, 1.443), area = c(11.708, 9.324, 4.259, 2.162)),
  "eff.soft" = list(joint = c(0.95, 0.95, 0.95, 0.96), cov_mu = c(0.97, 0.97, 0.97, 0.98),
                    width_mu = c(3.23, 2.757, 1.958, 1.154), area = c(14.677, 11.16, 4.983, 1.776)),
  "pb.hard" = list(joint = c(0.19, 0.4, 0.73, 0.82), cov_mu = c(0.55, 0.49, 0.73, 0.82),
                   width_mu = c(7.823, 7.06, 4.471, 1.771), area = c(31.371, 26.606, 12.971, 3.373)),
  "pb.soft" = list(joint = c(0.26, 0.53, 0.76, 0.87), cov_mu = c(0.55, 0.61, 0.78, 0.88),
                   width_mu = c(7.842, 6.809, 2.943, 1.01), area = c(32.496, 23.528, 7.598, 1.349)))


## ===========================================================================
## 3.  TABLE WRITERS
## ===========================================================================

MLAB <- c("Mahalanobis"            = "Mahalanobis",
          "Penalized Wald (mu)"    = "Penalized Wald ($\\mu$)",
          "Penalized Wald (sigma)" = "Penalized Wald ($\\sigma$)",
          "PB-ADI"                 = "\\textsc{pb-adi}")

## ---------------------------------------------------------------------------
## 3.1  Experiment 1: coverage, width, area
## ---------------------------------------------------------------------------

tab_exp1_coverage <- function() {
  
  S <- EXP1
  lab <- MLAB[S$method]
  body <- NULL
  for (i in seq_len(nrow(S))) {
    body <- rbind(body,
                  c(sprintf("\\multirow{2}{*}{%s}", lab[i]), "Coverage",
                    fse(S$cov_mu[i],  S$cov_mu_se[i]),
                    fse(S$cov_sig[i], S$cov_sig_se[i]),
                    fse(S$cov_jt[i],  S$cov_jt_se[i]), ""),
                  c("", "Width",
                    fse(S$w_mu[i],  S$w_mu_se[i]),
                    fse(S$w_sig[i], S$w_sig_se[i]), "",
                    fse(S$area[i],  S$area_se[i])))
  }
  
  note <- paste(
    "Notes: the two Penalized Wald rows are separate runs, each orthogonalised",
    "against the other coordinate, so the off-target column records the price",
    "of that orthogonalisation and the pair does not carry simultaneous",
    "coverage. Area is reported for completeness but is not the relevant",
    "criterion for a targeted statistic. Numerical failures were zero for",
    "every method, and no interval was returned at the search-box boundary in",
    "any replication.")
  
  inner <- latex_tabular(
    body   = body,
    header = c("Method", "", "$\\mu$", "$\\sigma$", "Joint", "Area"),
    align  = "llcccc",
    span   = "& & \\multicolumn{2}{c}{Marginal} & & \\\\",
    cmid   = "\\cmidrule(lr){3-4}",
    rowsep = seq(2, nrow(body) - 2, by = 2),
    note   = note)
  
  write_table(inner,
              sprintf(paste0("Experiment 1. Empirical coverage and average width of the ",
                             "nominal $95\\%%$ intervals, and average area of the joint ",
                             "region, over $1000$ replications at $n=100$. Monte Carlo ",
                             "standard errors in parentheses. Joint coverage is obtained ",
                             "by evaluating the acceptance rule at the true parameter, so ",
                             "it involves no nuisance optimisation and is a direct check ",
                             "on the implemented statistic; its reference value is the ",
                             "finite-$R$ bound $%.5f$, with a Monte Carlo standard error ",
                             "of $%.4f$."), BND_12, mc_se(BND_12, 1000)),
              "tab:exp1-coverage", "tab_exp1_coverage.tex")
}

## ---------------------------------------------------------------------------
## 3.2  Experiment 1: targeted width reduction
## ---------------------------------------------------------------------------

tab_exp1_targeted <- function() {
  
  T1 <- EXP1_TGT
  body <- cbind(paste0("$\\", T1$target, "$"),
                fnum(T1$mah, 3), fnum(T1$pw, 3),
                fpct(T1$reduce, 1), fse(T1$cov, T1$cov_se))
  
  inner <- latex_tabular(body,
                         header = c("Target", "Mahalanobis width", "Penalized Wald width",
                                    "Reduction", "Coverage"),
                         align  = "ccccc")
  
  write_table(inner,
              sprintf(paste0("Experiment 1. Width of the targeted Penalized Wald ",
                             "interval relative to the Mahalanobis interval for the same ",
                             "parameter, with the coverage of the targeted interval. The ",
                             "finite-$R$ coverage bound is $%.5f$."), BND_12),
              "tab:exp1-targeted", "tab_exp1_targeted.tex")
}

## ---------------------------------------------------------------------------
## 3.3  Experiment 1: the three sweeps
## ---------------------------------------------------------------------------

tab_lambda <- function() {
  
  S <- LAMBDA
  body <- cbind(
    paste0(formatC(S$lambda, format = "g"),
           ifelse(S$selected, "\\;$^{\\star}$", "")),
    fse(S$cov_mu, S$cov_mu_se), fse(S$w_mu, S$w_mu_se),
    fse(S$cov_sig, S$cov_sig_se), fse(S$w_sig, S$w_sig_se),
    paste0(fnum(S$fail_mu, 3), " / ", fnum(S$fail_sig, 3)),
    fnum(S$obj, 3))
  
  E <- LAMBDA_EVAL
  body <- rbind(body, c(
    paste0(formatC(E$lambda, format = "g"), "$^{\\dagger}$"),
    fse(E$cov_mu, E$cov_mu_se), fse(E$w_mu, E$w_mu_se),
    fse(E$cov_sig, E$cov_sig_se), fse(E$w_sig, E$w_sig_se),
    paste0(fnum(E$fail_mu, 3), " / ", fnum(E$fail_sig, 3)), ""))
  
  floor_c <- 0.95 - 1.96 * sqrt(0.95 * 0.05 / 1000)
  note <- sprintf(paste0(
    "$^{\\star}$ selected. $^{\\dagger}$ independent re-evaluation of the ",
    "selected value on a fresh seed stream, $1000$ further replications. ",
    "\\textsc{obj} is the width objective, the average of the two widths each ",
    "normalised by its grid minimum; the selection keeps grid values whose ",
    "coverage clears the Monte-Carlo-adjusted floor $%.4f$ for both targets, ",
    "which all seven do, so the choice is driven by width alone. Two caveats: ",
    "the sweep locates endpoints to a bisection tolerance of $0.1$ against ",
    "$10^{-3}$ in Table~\\ref{tab:exp1-coverage}, so its absolute widths are ",
    "inflated although the ordering in $\\lambda_n$ is not affected; and the ",
    "selection uses coverage of the known truth, so the chosen value is an ",
    "oracle choice."), floor_c)
  
  inner <- latex_tabular(body,
                         header = c("$\\lambda_n$", "Coverage", "Width", "Coverage", "Width",
                                    "Failure $\\mu$ / $\\sigma$", "\\textsc{obj}"),
                         align  = "ccccccc",
                         span   = paste0("& \\multicolumn{2}{c}{$\\mu$ target} & ",
                                         "\\multicolumn{2}{c}{$\\sigma$ target} & & \\\\"),
                         cmid   = "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
                         rowsep = nrow(body) - 1,
                         note   = note)
  
  write_table(inner,
              paste0("Sensitivity to the penalty weight $\\lambda_n$, $1000$ ",
                     "replications per grid point. The acceptance rule does not depend ",
                     "on $\\lambda_n$, so the coverage bound is $",
                     fnum(BND_12, 5), "$ in every row and the sweep measures width."),
              "tab:sens-lambda", "tab_sens_lambda.tex")
}

tab_raux <- function() {
  
  S <- RAUX
  body <- cbind(format(S$R_aux),
                fse(S$cov_mu, S$cov_mu_se), fse(S$w_mu, S$w_mu_se),
                fse(S$cov_sig, S$cov_sig_se), fse(S$w_sig, S$w_sig_se),
                paste0(fnum(S$fail_mu, 3), " / ", fnum(S$fail_sig, 3)))
  
  rng_mu  <- diff(range(S$w_mu))
  rng_sig <- diff(range(S$w_sig))
  note <- sprintf(paste0(
    "Across a sixteen-fold range the $\\mu$-target width varies by $%.4f$ and ",
    "the $\\sigma$-target width by $%.4f$, in both cases well inside one Monte ",
    "Carlo standard error of the width itself. The direction is therefore ",
    "already stable at $R_{\\mathrm{aux}}=50$ in this model, and this is the ",
    "one design constant that can be reduced to save computation without a ",
    "measurable cost."), rng_mu, rng_sig)
  
  inner <- latex_tabular(body,
                         header = c("$R_{\\mathrm{aux}}$", "Coverage", "Width", "Coverage", "Width",
                                    "Failure $\\mu$ / $\\sigma$"),
                         align  = "cccccc",
                         span   = paste0("& \\multicolumn{2}{c}{$\\mu$ target} & ",
                                         "\\multicolumn{2}{c}{$\\sigma$ target} & \\\\"),
                         cmid   = "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
                         note   = note)
  
  write_table(inner,
              paste0("Sensitivity to the auxiliary cloud size $R_{\\mathrm{aux}}$, ",
                     "$1000$ replications per grid point. $R_{\\mathrm{aux}}$ enters only ",
                     "the direction and never the acceptance rule, so the coverage bound ",
                     "is $", fnum(BND_12, 5), "$ throughout."),
              "tab:sens-raux", "tab_sens_raux.tex")
}

tab_rpc <- function() {
  
  S <- RPC
  body <- cbind(format(S$R), format(S$cutoff), fnum(S$bound, 5),
                fse(S$cov_mu, S$cov_mu_se), fnum(S$w_mu, 3),
                fse(S$cov_sig, S$cov_sig_se), fnum(S$w_sig, 3),
                paste0(fnum(S$fail_mu, 3), " / ", fnum(S$fail_sig, 3)))
  
  note <- sprintf(paste0(
    "Width falls monotonically in $R$, from $%.3f$ to $%.3f$ for $\\mu$, a ",
    "gain of about $%.0f\\%%$ over a sixteen-fold increase in cost, while the ",
    "failure rate rises from $%.3f$ to $%.3f$ because a finer rank grid places ",
    "more candidates near the cut-off. $R=200$ is the point past which the ",
    "bound is flat to four decimals and the width gain is under $1.5\\%%$ per ",
    "doubling."),
    S$w_mu[1], S$w_mu[nrow(S)],
    100 * (S$w_mu[1] - S$w_mu[nrow(S)]) / S$w_mu[1],
    S$fail_mu[1], S$fail_mu[nrow(S)])
  
  inner <- latex_tabular(body,
                         header = c("$R$", "$a_{R,\\alpha}$", "Bound", "Coverage", "Width",
                                    "Coverage", "Width", "Failure $\\mu$ / $\\sigma$"),
                         align  = "cccccccc",
                         span   = paste0("& & & \\multicolumn{2}{c}{$\\mu$ target} & ",
                                         "\\multicolumn{2}{c}{$\\sigma$ target} & \\\\"),
                         cmid   = "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
                         note   = note)
  
  write_table(inner,
              paste0("Sensitivity to the inference cloud size $R$, $1000$ replications ",
                     "per grid point. The acceptance rule depends on $R$ through the ",
                     "cut-off $a_{R,\\alpha}=\\lfloor\\alpha(R+1)\\rfloor+1$, so the ",
                     "coverage bound moves with it and is given in its own column; ",
                     "coverage must be read against the bound in the same row."),
              "tab:sens-R", "tab_sens_R.tex")
}

## ---------------------------------------------------------------------------
## 3.4  Experiment 2: size and power, and the clamping sweep
## ---------------------------------------------------------------------------

tab_exp2_power <- function() {
  
  body <- cbind(E2_B, matrix(fcell(E2_OURS, 3), nrow = nrow(E2_OURS)))
  inner <- latex_tabular(body,
                         header = c("$\\beta_1^\\ast$", paste0("$n=", E2_N, "$")),
                         align  = paste(rep("c", 1 + length(E2_N)), collapse = ""))
  
  write_table(inner,
              sprintf(paste0("Experiment 2. Rejection probability of $H_0:\\beta_1=0$ at ",
                             "the $%.2f$ level, $1000$ replications per cell. The first ",
                             "row is empirical size and should be read against the exact ",
                             "finite-$R$ bound $%.5f$; the remaining rows are power. ",
                             "Unresolved replications count as non-rejections, so this ",
                             "is the conservative power."), ALPHA_12, SZB_12),
              "tab:exp2-power", "tab_exp2_power.tex")
}

tab_delta <- function() {
  
  body <- cbind(D_DELTA,
                matrix(fnum(DEL_OURS_H0, 3), nrow = nrow(DEL_OURS_H0)),
                matrix(fnum(DEL_OURS_H1, 3), nrow = nrow(DEL_OURS_H1)))
  
  inner <- latex_tabular(body,
                         header = c("$\\Delta$", paste0("$n{=}", D_N, "$"), paste0("$n{=}", D_N, "$")),
                         align  = paste(rep("c", 1 + 2 * length(D_N)), collapse = ""),
                         span   = sprintf("& \\multicolumn{%d}{c}{Size at $\\beta_1^\\ast=0$} & \\multicolumn{%d}{c}{Power at $\\beta_1^\\ast=1$} \\\\",
                                          length(D_N), length(D_N)),
                         cmid   = sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}",
                                          1 + length(D_N), 2 + length(D_N), 1 + 2 * length(D_N)),
                         note   = paste("All forty cells resolved; the failure rate was zero",
                                        "throughout. Size is at or below the exact bound",
                                        paste0("$", fnum(SZB_12, 5), "$"),
                                        "in every cell."))
  
  write_table(inner,
              paste0("Sensitivity to the clamping range $\\Delta$ in Experiment 2, ",
                     "$1000$ replications per cell. $\\Delta$ changes the model rather ",
                     "than the rule: a small $\\Delta$ destroys information while a ",
                     "large one inflates the privacy noise, which scales as $\\Delta^2$ ",
                     "for the second-moment coordinates."),
              "tab:sens-delta", "tab_sens_delta.tex")
}

## ---------------------------------------------------------------------------
## 3.5  Experiment 3
## ---------------------------------------------------------------------------

tab_exp3 <- function() {
  
  sat <- E3_WID_OURS >= E3_SAT
  body <- NULL
  for (i in seq_along(E3_EPS))
    for (j in seq_along(E3_N))
      body <- rbind(body, c(
        if (j == 1L) E3_EPS[i] else "",
        format(E3_N[j]),
        fnum(E3_COV_OURS[i, j], 3), fnum(E3_COV_REPRO[i, j], 2),
        fnum(E3_COV_ADI[i, j], 2),  fnum(E3_COV_NAIVE[i, j], 2),
        paste0(fnum(E3_WID_OURS[i, j], 3),
               if (sat[i, j]) "$^{\\dagger}$" else ""),
        fnum(E3_WID_REPRO[i, j], 2), fnum(E3_WID_ADI[i, j], 2)))
  
  note <- sprintf(paste0(
    "$^{\\dagger}$ width at or above $%.1f$, that is at least $90\\%%$ of the ",
    "search box $[%g,%g]$: the interval is box-limited and its coverage ",
    "follows from truncation rather than from the procedure, so those cells ",
    "carry no information about relative efficiency. Published values are from ",
    "Wang, Chang and Awan (Fig.~6)."), E3_SAT, E3_BOX[1], E3_BOX[2])
  
  inner <- latex_tabular(body,
                         header = c("$\\varepsilon$", "$n$",
                                    "\\textsc{pw}", "\\textsc{repro}", "\\textsc{pb-adi}",
                                    "\\textsc{pb-naive}",
                                    "\\textsc{pw}", "\\textsc{repro}", "\\textsc{pb-adi}"),
                         align  = "ccccccccc",
                         span   = "& & \\multicolumn{4}{c}{Coverage} & \\multicolumn{3}{c}{Width} \\\\",
                         cmid   = "\\cmidrule(lr){3-6}\\cmidrule(lr){7-9}",
                         note   = note)
  
  write_table(inner,
              sprintf(paste0("Experiment 3. Coverage and average width of the $90\\%%$ ",
                             "interval for $\\beta_1^\\ast$, $1000$ replications per cell. ",
                             "The finite-$R$ coverage bound is $%.5f$, not $0.95$: this ",
                             "experiment runs at $\\alpha=%.2f$, and the Monte Carlo ",
                             "standard error at that level is $%.4f$."),
                      BND_3, ALPHA_3, mc_se(BND_3, 1000)),
              "tab:exp3-grid", "tab_exp3_grid.tex")
}

## ---------------------------------------------------------------------------
## 3.6  Clamping-boundary studies
## ---------------------------------------------------------------------------

tab_pbadi_grid <- function() {
  
  body <- NULL
  for (i in seq_along(PB_MU))
    for (j in seq_along(PB_EPS))
      body <- rbind(body, c(
        if (j == 1L) PB_MU[i] else "", PB_EPS[j],
        fnum(G1_JOINT$mah[i, j], 2), fnum(G1_JOINT$eff[i, j], 2),
        fnum(G1_JOINT$pb[i, j], 2),
        fnum(G1_COVMU$mah[i, j], 2), fnum(G1_COVMU$eff[i, j], 2),
        fnum(G1_COVMU$pb[i, j], 2),
        fnum(G1_WIDMU$mah[i, j], 3), fnum(G1_WIDMU$eff[i, j], 3),
        fnum(G1_WIDMU$pb[i, j], 3)))
  
  note <- sprintf(paste0(
    "These studies use $100$ replications per cell, not $1000$, so the ",
    "binomial standard error at a coverage of $0.95$ is $%.3f$ and differences ",
    "smaller than about $0.04$ should not be read as real. The search box for ",
    "$\\sigma$ is $[0.1,5]$, so a $\\sigma$-width of $4.9$ is the whole box; at ",
    "$\\varepsilon=0.1$ the $\\sigma$-widths of both repro-family methods are ",
    "box-saturated, and that column carries no information about efficiency in ",
    "the scale parameter."), mc_se(0.95, 100))
  
  inner <- latex_tabular(body,
                         header = c("$\\mu^\\ast$", "$\\varepsilon$",
                                    "\\textsc{rp}", "\\textsc{pw}", "\\textsc{adi}",
                                    "\\textsc{rp}", "\\textsc{pw}", "\\textsc{adi}",
                                    "\\textsc{rp}", "\\textsc{pw}", "\\textsc{adi}"),
                         align  = "ccccccccccc",
                         span   = paste0("& & \\multicolumn{3}{c}{Joint coverage} & ",
                                         "\\multicolumn{3}{c}{Coverage of $\\mu$} & ",
                                         "\\multicolumn{3}{c}{Width for $\\mu$} \\\\"),
                         cmid   = "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
                         note   = note)
  
  write_table(inner,
              paste0("Study 1 of the clamping-boundary experiment: the true mean is ",
                     "moved towards the upper clamp $U=3$ at four privacy budgets, with ",
                     "$\\sigma^\\ast=1$ fixed. \\textsc{rp} is Mahalanobis depth, ",
                     "\\textsc{pw} the Penalized Wald statistic targeting $\\mu$, and ",
                     "\\textsc{adi} the debiased parametric bootstrap. The coverage ",
                     "bound is $", fnum(BND_12, 5), "$."),
              "tab:pbadi-grid", "tab_pbadi_grid.tex")
}

tab_pbadi_soft <- function() {
  
  key <- function(m, ct) S2[[paste0(m, ".", ct)]]
  body <- NULL
  for (m in c("mah", "eff", "pb")) {
    lab <- c(mah = "Mahalanobis (\\textsc{repro})",
             eff = "Penalized Wald ($\\mu$)",
             pb  = "\\textsc{pb-adi}")[m]
    H <- key(m, "hard"); S <- key(m, "soft")
    for (j in seq_along(PB_EPS))
      body <- rbind(body, c(
        if (j == 1L) lab else "", PB_EPS[j],
        fnum(H$joint[j], 2),    fnum(S$joint[j], 2),
        fnum(H$cov_mu[j], 2),   fnum(S$cov_mu[j], 2),
        fnum(H$width_mu[j], 3), fnum(S$width_mu[j], 3),
        fnum(H$area[j], 2),     fnum(S$area[j], 2)))
  }
  
  note <- paste(
    "The soft clamp is the logistic map into $(L,U)$ with $k=1.5$. Because the",
    "transformed data still lie in $[L,U]$ the sensitivities are unchanged, so",
    "the two regimes are compared at exactly the same privacy cost, and the",
    "clamp type is threaded through the release, the region grid, the indirect",
    "estimator and the bootstrap resampler, so no method is solved under a",
    "mismatched model.")
  
  inner <- latex_tabular(body,
                         header = c("Method", "$\\varepsilon$", "hard", "soft", "hard", "soft",
                                    "hard", "soft", "hard", "soft"),
                         align  = "lccccccccc",
                         span   = paste0("& & \\multicolumn{2}{c}{Joint coverage} & ",
                                         "\\multicolumn{2}{c}{Coverage of $\\mu$} & ",
                                         "\\multicolumn{2}{c}{Width for $\\mu$} & ",
                                         "\\multicolumn{2}{c}{Area} \\\\"),
                         cmid   = paste0("\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
                                         "\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}"),
                         rowsep = c(length(PB_EPS), 2 * length(PB_EPS)),
                         note   = note)
  
  write_table(inner,
              paste0("Study 2 of the clamping-boundary experiment: hard against soft ",
                     "clamping at $\\mu^\\ast=3$, $100$ replications per cell. Softening ",
                     "the clamp restores the Penalized Wald statistic to the bound and ",
                     "improves the bootstrap without repairing it, while leaving ",
                     "Mahalanobis depth unchanged, which is what isolates the cause to ",
                     "the non-differentiability of the hard clamp."),
              "tab:pbadi-soft", "tab_pbadi_soft.tex")
}


## ===========================================================================
## 4.  FIGURE WRITERS
##
##  Base R graphics.  Each heat map is drawn cell by cell with rect() so that
##  the two-segment grey ramp and the printed value are under direct control;
##  image() would interpolate the colour scale linearly and lose the break at
##  the nominal level, which is the whole point of the convention.
## ===========================================================================

## ---------------------------------------------------------------------------
## 4.1  One heat-map panel
## ---------------------------------------------------------------------------

heat_panel <- function(m, greyfn, main = "", xlab = "", ylab = "",
                       dec = 3, cex_cell = 0.62, cex_ax = 0.72) {
  
  m  <- as.matrix(m)
  nr <- nrow(m); nc <- ncol(m)
  
  par(mar = c(3.1, 3.3, 1.9, 0.6), mgp = c(2.0, 0.55, 0), tcl = -0.2)
  plot(NA, xlim = c(0, nc), ylim = c(0, nr), axes = FALSE,
       xlab = xlab, ylab = ylab, main = main, cex.main = 0.92,
       cex.lab = 0.82, font.main = 1, xaxs = "i", yaxs = "i")
  
  for (i in seq_len(nr)) {
    y <- nr - i
    for (j in seq_len(nc)) {
      v <- m[i, j]
      g <- greyfn(v)
      rect(j - 1, y, j, y + 1, col = pct_to_col(g),
           border = "grey85", lwd = 0.4)
      if (!is.na(v))
        text(j - 0.5, y + 0.5, fcell(v, dec), cex = cex_cell,
             col = if (!is.na(g) && g >= 46) "white" else "black")
      else
        text(j - 0.5, y + 0.5, "--", cex = cex_cell, col = "grey55")
    }
  }
  axis(1, at = seq_len(nc) - 0.5, labels = colnames(m),
       tick = FALSE, line = -0.7, cex.axis = cex_ax)
  axis(2, at = rev(seq_len(nr)) - 0.5, labels = rownames(m),
       tick = FALSE, line = -0.7, las = 1, cex.axis = cex_ax)
  box(col = "grey70", lwd = 0.5)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## 4.2  Vertical colour key, drawn as its own panel
## ---------------------------------------------------------------------------

colour_key <- function(greyfn, ticks, label, nseg = 64) {
  
  par(mar = c(3.1, 0.4, 1.9, 3.4), mgp = c(2.0, 0.5, 0), tcl = -0.18)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (k in seq_len(nseg)) {
    v <- (k - 0.5) / nseg
    rect(0.15, (k - 1) / nseg, 0.85, k / nseg,
         col = pct_to_col(greyfn(v)), border = NA)
  }
  rect(0.15, 0, 0.85, 1, border = "grey60", lwd = 0.5)
  axis(4, at = ticks, labels = formatC(ticks, format = "f", digits = 2),
       las = 1, cex.axis = 0.66, pos = 0.85)
  mtext(label, side = 4, line = 2.2, cex = 0.68)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## 4.3  A 2x2 block of heat maps with a shared key
## ---------------------------------------------------------------------------

fig_heat_2x2 <- function(panels, titles, greyfn, ticks, key_label,
                         file, xlab, ylab, dec = rep(3, 4),
                         width = 7.4, height = 5.2) {
  
  stopifnot(length(panels) == 4, length(titles) == 4)
  path <- file.path(DIR_FIG, file)
  pdf(path, width = width, height = height)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  layout(matrix(c(1, 2, 5,
                  3, 4, 5), nrow = 2, byrow = TRUE),
         widths = c(1, 1, 0.30))
  for (i in 1:4)
    heat_panel(panels[[i]], greyfn, main = titles[i],
               xlab = xlab, ylab = ylab, dec = dec[i])
  colour_key(greyfn, ticks, key_label)
  invisible(path)
}

## ---------------------------------------------------------------------------
## 4.4  A 1x3 row of heat maps with a shared key
## ---------------------------------------------------------------------------

fig_heat_1x3 <- function(panels, titles, greyfn, ticks, key_label,
                         file, xlab, ylab, dec = rep(3, 3),
                         width = 8.6, height = 2.9) {
  
  stopifnot(length(panels) == 3, length(titles) == 3)
  path <- file.path(DIR_FIG, file)
  pdf(path, width = width, height = height)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.34))
  for (i in 1:3)
    heat_panel(panels[[i]], greyfn, main = titles[i],
               xlab = xlab, ylab = ylab, dec = dec[i])
  colour_key(greyfn, ticks, key_label)
  invisible(path)
}

## ---------------------------------------------------------------------------
## 4.5  Experiment 2: the four-panel rejection grid, and the size plot
## ---------------------------------------------------------------------------

fig_exp2_grid <- function() {
  fig_heat_2x2(
    panels = list(E2_OURS, E2_REPRO, E2_ADI, E2_NAIVE),
    titles = c("Penalized Wald (this paper)", "Repro (Awan and Wang)",
               "PB-ADI (Wang, Chang and Awan)", "PB-naive (naive plug-in + F)"),
    greyfn = grey_reject,
    ticks  = c(ALPHA_12, 0.5, 1),
    key_label = "Rejection probability",
    file   = "fig_exp2_grid.pdf",
    xlab   = "Sample size n", ylab = expression("True " * beta[1]^"*"))
  
  write_figure_wrapper("fig_exp2_grid.pdf",
                       sprintf(paste0("Experiment 2. Rejection probability of $H_0:\\beta_1=0$ at ",
                                      "the $%.2f$ level. Rows are the true slope, columns the ",
                                      "sample size. The grey scale runs from $0$ to $1$ and breaks ",
                                      "at $%.2f$, so cells at or below the nominal level are light ",
                                      "and cells above it dark: the bottom row of each panel is ",
                                      "empirical size and the rest is power."), ALPHA_12, ALPHA_12),
                       "fig:exp2-grid", width = 0.98)
}

fig_exp2_size <- function() {
  
  path <- file.path(DIR_FIG, "fig_exp2_size.pdf")
  pdf(path, width = 6.4, height = 3.7)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  par(mar = c(4.0, 4.2, 1.0, 1.0), mgp = c(2.5, 0.7, 0), tcl = -0.25, las = 1)
  ymax <- max(c(T1_OURS, T1_REPRO, T1_ADI, T1_NAIVE), na.rm = TRUE) * 1.10
  
  plot(NA, xlim = range(E2_N8), ylim = c(0, ymax), log = "x",
       xlab = "Sample size n", ylab = expression("Empirical size at " * beta[1]^"*" == 0),
       axes = FALSE, cex.lab = 0.95)
  axis(1, at = E2_N8, labels = E2_N8, cex.axis = 0.8)
  axis(2, cex.axis = 0.8)
  box(col = "grey60", lwd = 0.6)
  abline(h = pretty(c(0, ymax)), col = "grey92", lwd = 0.6)
  abline(h = ALPHA_12, lty = 2, col = "grey25", lwd = 1.1)
  abline(h = SZB_12,  lty = 3, col = "grey45", lwd = 1.1)
  
  series <- list(
    list(v = T1_OURS,  col = "black",   pch = 19, lty = 1,
         lab = "Penalized Wald (this paper)"),
    list(v = T1_REPRO, col = "grey40",  pch = 15, lty = 1,
         lab = "Repro (Awan and Wang)"),
    list(v = T1_ADI,   col = "grey55",  pch = 17, lty = 2,
         lab = "PB-ADI (Wang, Chang and Awan)"),
    list(v = T1_NAIVE, col = "grey70",  pch = 18, lty = 4,
         lab = "PB-naive (naive plug-in)"))
  
  for (s in series) {
    keep <- !is.na(s$v)
    lines(E2_N8[keep], s$v[keep], col = s$col, lty = s$lty, lwd = 1.6)
    points(E2_N8[keep], s$v[keep], col = s$col, pch = s$pch, cex = 0.85)
  }
  legend("topleft", bty = "n", cex = 0.68,
         legend = c(vapply(series, function(s) s$lab, character(1)),
                    "nominal 0.05", sprintf("exact bound %.5f", SZB_12)),
         col = c(vapply(series, function(s) s$col, character(1)), "grey25", "grey45"),
         lty = c(vapply(series, function(s) s$lty, numeric(1)), 2, 3),
         pch = c(vapply(series, function(s) s$pch, numeric(1)), NA, NA),
         lwd = 1.4)
  
  write_figure_wrapper("fig_exp2_size.pdf",
                       sprintf(paste0("Experiment 2. Empirical size at $\\beta_1^\\ast=0$ against ",
                                      "sample size. The dashed line is the nominal level $%.2f$ ",
                                      "and the dotted line the exact finite-$R$ bound $%.5f$. Our ",
                                      "design stops at $n=1000$, so the two largest sample sizes ",
                                      "are shown for the published methods only."),
                               ALPHA_12, SZB_12),
                       "fig:exp2-size", width = 0.80)
}

## ---------------------------------------------------------------------------
## 4.6  Experiment 3: coverage grid and width curves
## ---------------------------------------------------------------------------

fig_exp3_coverage <- function() {
  gc3 <- function(v) grey_cov(v, nominal = 1 - ALPHA_3)
  fig_heat_2x2(
    panels = list(E3_COV_OURS, E3_COV_REPRO, E3_COV_ADI, E3_COV_NAIVE),
    titles = c("Penalized Wald", "Repro", "PB-ADI", "PB-naive"),
    greyfn = gc3,
    ticks  = c(0.2, 0.6, 1 - ALPHA_3),
    key_label = "Coverage",
    file   = "fig_exp3_coverage.pdf",
    xlab   = "Sample size n", ylab = expression(epsilon),
    dec    = c(3, 2, 2, 2),
    width  = 6.8, height = 4.4)
  
  write_figure_wrapper("fig_exp3_coverage.pdf",
                       sprintf(paste0("Experiment 3. Empirical coverage of the $90\\%%$ interval ",
                                      "for $\\beta_1^\\ast$. Rows are the privacy parameter, ",
                                      "columns the sample size. The grey scale breaks at the ",
                                      "nominal level $%.2f$, so cells below it are light."),
                               1 - ALPHA_3),
                       "fig:exp3-coverage", width = 0.92)
}

fig_exp3_width <- function() {
  
  path <- file.path(DIR_FIG, "fig_exp3_width.pdf")
  pdf(path, width = 8.4, height = 3.1)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  panels <- list("Penalized Wald (this paper)"      = E3_WID_OURS,
                 "Repro (Awan and Wang)"            = E3_WID_REPRO,
                 "PB-ADI (Wang, Chang and Awan)"    = E3_WID_ADI)
  cols <- c("black", "grey40", "grey55", "grey70")
  ltys <- c(1, 2, 3, 4)
  pchs <- c(19, 15, 17, 18)
  yt   <- c(0.5, 1, 2, 5, 10, 20)
  
  par(mfrow = c(1, 3), mar = c(3.9, 3.9, 2.0, 0.8),
      mgp = c(2.4, 0.6, 0), tcl = -0.22, las = 1)
  
  for (p in seq_along(panels)) {
    W <- panels[[p]]
    plot(NA, xlim = range(E3_N), ylim = c(0.3, 30), log = "xy",
         xlab = "n", ylab = "CI width", main = names(panels)[p],
         cex.main = 0.88, font.main = 1, cex.lab = 0.9, axes = FALSE)
    axis(1, at = E3_N, labels = E3_N, cex.axis = 0.78)
    axis(2, at = yt, labels = yt, cex.axis = 0.78)
    box(col = "grey60", lwd = 0.6)
    abline(h = yt, col = "grey93", lwd = 0.5)
    abline(h = E3_BOXW, lty = 2, col = "grey60", lwd = 1.0)
    for (i in seq_along(E3_EPS)) {
      lines(E3_N, W[i, ], col = cols[i], lty = ltys[i], lwd = 1.6)
      points(E3_N, W[i, ], col = cols[i], pch = pchs[i], cex = 0.8)
    }
    if (p == 1L)
      legend("bottomleft", bty = "n", cex = 0.7,
             legend = c(parse(text = paste0("epsilon == ", E3_EPS)),
                        expression("search box")),
             col = c(cols[seq_along(E3_EPS)], "grey60"),
             lty = c(ltys[seq_along(E3_EPS)], 2),
             pch = c(pchs[seq_along(E3_EPS)], NA), lwd = 1.4)
  }
  
  write_figure_wrapper("fig_exp3_width.pdf",
                       sprintf(paste0("Experiment 3. Average width of the $90\\%%$ interval for ",
                                      "$\\beta_1^\\ast$, on logarithmic axes, one curve per ",
                                      "privacy parameter. The dashed horizontal line is the total ",
                                      "width $%g$ of the search box $[%g,%g]$: a curve touching it ",
                                      "is reporting the box rather than an interval, which is the ",
                                      "case at $\\varepsilon=0.1$ for both our procedure and ",
                                      "Repro. \\textsc{pb-naive} is omitted because its widths are ",
                                      "the smallest of the four, which follows from its ",
                                      "under-coverage rather than from efficiency."),
                               E3_BOXW, E3_BOX[1], E3_BOX[2]),
                       "fig:exp3-width", width = 0.98)
}

## ---------------------------------------------------------------------------
## 4.7  Sensitivity: lambda
## ---------------------------------------------------------------------------

fig_lambda <- function() {
  
  path <- file.path(DIR_FIG, "fig_sens_lambda.pdf")
  pdf(path, width = 6.2, height = 3.6)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  x <- seq_along(LAMBDA$lambda) - 1          # geometric grid: ordinal spacing
  par(mar = c(4.0, 4.2, 1.0, 1.0), mgp = c(2.5, 0.7, 0), tcl = -0.25, las = 1)
  yl <- range(c(LAMBDA$w_mu, LAMBDA$w_sig)) * c(0.92, 1.08)
  
  plot(NA, xlim = range(x), ylim = yl, axes = FALSE,
       xlab = expression(lambda[n] * " (ordinal spacing)"),
       ylab = "Average interval width", cex.lab = 0.95)
  axis(1, at = x, labels = formatC(LAMBDA$lambda, format = "g"), cex.axis = 0.8)
  axis(2, cex.axis = 0.8)
  box(col = "grey60", lwd = 0.6)
  abline(h = pretty(yl), col = "grey93", lwd = 0.6)
  abline(v = x[LAMBDA$selected], lty = 3, col = "grey45", lwd = 1.2)
  
  lines(x, LAMBDA$w_mu,  col = "black",  lwd = 1.7)
  points(x, LAMBDA$w_mu, col = "black",  pch = 19, cex = 0.85)
  lines(x, LAMBDA$w_sig, col = "grey45", lwd = 1.7, lty = 2)
  points(x, LAMBDA$w_sig, col = "grey45", pch = 15, cex = 0.85)
  
  legend("topright", bty = "n", cex = 0.72,
         legend = c(expression(mu * " target"), expression(sigma * " target"),
                    "selected"),
         col = c("black", "grey45", "grey45"), lty = c(1, 2, 3),
         pch = c(19, 15, NA), lwd = 1.5)
  
  write_figure_wrapper("fig_sens_lambda.pdf",
                       sprintf(paste0("Sensitivity to the penalty weight $\\lambda_n$: average ",
                                      "interval width for each target. Grid points are drawn at ",
                                      "equal ordinal spacing because the design grid is ",
                                      "geometric; the dotted line marks the selected value ",
                                      "$\\lambda_n=%g$. The $\\mu$-target width falls from $%.3f$ ",
                                      "at $\\lambda_n=0$ to $%.3f$ at the selected value and then ",
                                      "rises slowly, so the penalty is not cosmetic: with no ",
                                      "Mahalanobis term the profile objective is flat in the ",
                                      "nuisance direction. Coverage clears the ",
                                      "Monte-Carlo-adjusted floor at every grid point, so the ",
                                      "figure is about width alone."),
                               LAMBDA$lambda[LAMBDA$selected][1],
                               LAMBDA$w_mu[1], LAMBDA$w_mu[LAMBDA$selected][1]),
                       "fig:sens-lambda", width = 0.78)
}

## ---------------------------------------------------------------------------
## 4.8  Sensitivity: clamping range, size and power
## ---------------------------------------------------------------------------

fig_delta <- function() {
  
  fig_heat_1x3(
    panels = list(DEL_OURS_H0, DEL_REPRO_H0, DEL_PB_H0),
    titles = c("Penalized Wald", "Repro", "Parametric bootstrap"),
    greyfn = grey_reject, ticks = c(ALPHA_12, 0.5, 1),
    key_label = "Type I error",
    file = "fig_sens_delta_size.pdf",
    xlab = "Sample size n", ylab = "Clamping bound")
  
  write_figure_wrapper("fig_sens_delta_size.pdf",
                       sprintf(paste0("Empirical size at $\\beta_1^\\ast=0$ as a function of the ",
                                      "clamping range. The grey scale breaks at $%.2f$, so any ",
                                      "dark cell is a violation of the nominal level. Published ",
                                      "panels are from Awan and Wang (2025, Fig.~6), restricted ",
                                      "to the cells of our design."), ALPHA_12),
                       "fig:sens-delta-size", width = 0.98)
  
  fig_heat_1x3(
    panels = list(DEL_OURS_H1, DEL_REPRO_H1, DEL_PB_H1),
    titles = c("Penalized Wald", "Repro", "Parametric bootstrap"),
    greyfn = grey_reject, ticks = c(ALPHA_12, 0.5, 1),
    key_label = "Power",
    file = "fig_sens_delta_power.pdf",
    xlab = "Sample size n", ylab = "Clamping bound")
  
  write_figure_wrapper("fig_sens_delta_power.pdf",
                       paste0("Power at $\\beta_1^\\ast=1$ as a function of the clamping range. ",
                              "The three panels agree that $\\Delta\\in[0.5,2]$ is the usable ",
                              "range and that power collapses by $\\Delta=10$, where the privacy ",
                              "noise scales as $\\Delta^2$. The bootstrap's apparent advantage at ",
                              "$\\Delta=2$ cannot be read as power, since its size at that ",
                              "$\\Delta$ is inflated."),
                       "fig:sens-delta-power", width = 0.98)
}

fig_delta_full <- function() {
  
  path <- file.path(DIR_FIG, "fig_sens_delta_full.pdf")
  pdf(path, width = 8.8, height = 3.4)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.26))
  heat_panel(DEL_REPRO_H0_FULL, grey_reject, main = "Repro (Awan and Wang)",
             xlab = "Sample size n", ylab = "Clamping bound",
             cex_cell = 0.5, cex_ax = 0.62)
  heat_panel(DEL_PB_H0_FULL, grey_reject,
             main = "Parametric bootstrap (Awan and Wang)",
             xlab = "Sample size n", ylab = "Clamping bound",
             cex_cell = 0.5, cex_ax = 0.62)
  colour_key(grey_reject, c(ALPHA_12, 0.5, 1), "Type I error")
  
  write_figure_wrapper("fig_sens_delta_full.pdf",
                       sprintf(paste0("The published clamping grid in full, for context (Awan and ",
                                      "Wang 2025, Fig.~6). The parametric bootstrap's type-I error ",
                                      "grows with the sample size rather than shrinking, reaching ",
                                      "$%.3f$ at $\\Delta=1$, $n=5000$, while the repro column ",
                                      "stays at or below $%.3f$ everywhere. Our design covers the ",
                                      "left half of these grids; the two largest sample sizes and ",
                                      "the intermediate clamping values $\\Delta\\in\\{0.8,1.5\\}$ ",
                                      "were not run, and those are the cells in which the ",
                                      "bootstrap fails hardest."),
                               max(DEL_PB_H0_FULL), max(DEL_REPRO_H0_FULL)),
                       "fig:pbfull", width = 0.98)
}

## ---------------------------------------------------------------------------
## 4.9  Clamping boundary: joint coverage grid, and the soft-clamp comparison
## ---------------------------------------------------------------------------

fig_pbadi_joint <- function() {
  
  fig_heat_1x3(
    panels = list(G1_JOINT$mah, G1_JOINT$eff, G1_JOINT$pb),
    ## titles pass through c(), so they must be plain character: mixing a
    ## string with an expression() would deparse the expression to its source
    titles = c("Mahalanobis (Repro)", "Penalized Wald (mu target)", "PB-ADI"),
    greyfn = grey_cov, ticks = c(0.2, 0.6, 0.95),
    key_label = "Joint coverage",
    file = "fig_pbadi_joint.pdf",
    xlab = expression("Privacy " * epsilon), ylab = expression(mu^"*"),
    dec = rep(2, 3))
  
  write_figure_wrapper("fig_pbadi_joint.pdf",
                       sprintf(paste0("Study 1. Joint coverage of the nominal $95\\%%$ region for ",
                                      "$(\\mu,\\sigma)$. Rows are the true mean, columns the ",
                                      "privacy budget. The grey scale breaks at $0.95$, so light ",
                                      "cells are under-coverage. \\textsc{pb-adi} degrades in both ",
                                      "directions at once, reaching $%.2f$ at $\\mu^\\ast=3$, ",
                                      "$\\varepsilon=0.1$. Mahalanobis depth is flat at ",
                                      "$%.2f$--$%.2f$ over the whole grid. The Penalized Wald ",
                                      "statistic is flat until the last row, where it drops to ",
                                      "$%.2f$."),
                               min(G1_JOINT$pb), min(G1_JOINT$mah), max(G1_JOINT$mah),
                               min(G1_JOINT$eff)),
                       "fig:pbadi-joint", width = 0.98)
}

fig_pbadi_soft <- function() {
  
  path <- file.path(DIR_FIG, "fig_pbadi_soft.pdf")
  pdf(path, width = 6.4, height = 3.7)
  on.exit({ dev.off(); message("  [pdf]   ", path) }, add = TRUE)
  
  eps <- as.numeric(PB_EPS)
  par(mar = c(4.0, 4.2, 1.0, 1.0), mgp = c(2.5, 0.7, 0), tcl = -0.25, las = 1)
  plot(NA, xlim = range(eps), ylim = c(0.10, 1.03), log = "x", axes = FALSE,
       xlab = expression("Privacy budget " * epsilon),
       ylab = expression("Joint coverage at " * mu^"*" == 3), cex.lab = 0.95)
  axis(1, at = eps, labels = PB_EPS, cex.axis = 0.8)
  axis(2, cex.axis = 0.8)
  box(col = "grey60", lwd = 0.6)
  abline(h = pretty(c(0.1, 1)), col = "grey93", lwd = 0.6)
  abline(h = BND_12, lty = 2, col = "grey40", lwd = 1.1)
  
  mm   <- c("mah", "eff", "pb")
  cols <- c(mah = "black", eff = "grey40", pb = "grey60")
  pchs <- c(mah = 19, eff = 15, pb = 17)
  for (m in mm) {
    H <- S2[[paste0(m, ".hard")]]$joint
    S <- S2[[paste0(m, ".soft")]]$joint
    lines(eps, H, col = cols[m], lwd = 1.7, lty = 1)
    points(eps, H, col = cols[m], pch = pchs[m], cex = 0.85)
    lines(eps, S, col = cols[m], lwd = 1.7, lty = 2)
    points(eps, S, col = cols[m], pch = pchs[m], cex = 0.85, bg = "white")
  }
  legend("bottomright", bty = "n", cex = 0.66, ncol = 2,
         legend = c("Mahalanobis, hard", "Mahalanobis, soft",
                    "Penalized Wald, hard", "Penalized Wald, soft",
                    "PB-ADI, hard", "PB-ADI, soft"),
         col = rep(cols[mm], each = 2), lty = rep(c(1, 2), 3),
         pch = rep(pchs[mm], each = 2), lwd = 1.5)
  
  write_figure_wrapper("fig_pbadi_soft.pdf",
                       sprintf(paste0("Study 2. Joint coverage at $\\mu^\\ast=3$ under hard ",
                                      "(solid) and soft (dashed) clamping, at identical privacy ",
                                      "cost. The dashed horizontal line is the finite-$R$ bound ",
                                      "$%.5f$. Softening the clamp restores the Penalized Wald ",
                                      "statistic to $%.2f$--$%.2f$ from $%.2f$--$%.2f$, with the ",
                                      "largest repair exactly where the hard-clamp loss was ",
                                      "largest. \\textsc{pb-adi} improves but is not repaired, and ",
                                      "Mahalanobis depth is unchanged, as it must be, since it ",
                                      "never differentiates the mechanism."),
                               BND_12,
                               min(S2$eff.soft$joint), max(S2$eff.soft$joint),
                               min(S2$eff.hard$joint), max(S2$eff.hard$joint)),
                       "fig:pbadi-soft", width = 0.80)
}


## ===========================================================================
## 5.  DRIVER
## ===========================================================================

TABLES <- c("tab_exp1_coverage", "tab_exp1_targeted",
            "tab_lambda", "tab_raux", "tab_rpc",
            "tab_exp2_power", "tab_delta",
            "tab_exp3",
            "tab_pbadi_grid", "tab_pbadi_soft")

FIGURES <- c("fig_exp2_grid", "fig_exp2_size",
             "fig_exp3_coverage", "fig_exp3_width",
             "fig_lambda", "fig_delta", "fig_delta_full",
             "fig_pbadi_joint", "fig_pbadi_soft")

## Produces everything.  A writer that stops with an error is reported and the
## run continues, so one problem cannot cost the whole set; anything that did
## not complete is repeated at the end.  Graphics devices are closed on the
## way out even after a failure.
make_all <- function(what = c("both", "tables", "figures")) {
  
  what <- match.arg(what)
  todo <- switch(what,
                 both    = c(TABLES, FIGURES),
                 tables  = TABLES,
                 figures = FIGURES)
  failed <- character(0)
  
  for (f in todo) {
    ok <- tryCatch({ do.call(f, list()); TRUE },
                   error = function(e) {
                     message("  [ERROR] ", f, "(): ", conditionMessage(e))
                     FALSE })
    if (!ok) {
      failed <- c(failed, f)
      while (dev.cur() > 1L) dev.off()   # a failed pdf() would stay open
    }
  }
  
  message("\n-----------------------------------------------------------")
  message("  tables  in ", normalizePath(DIR_TAB, mustWork = FALSE), ": ",
          length(list.files(DIR_TAB, pattern = "\\.tex$")))
  message("  figures in ", normalizePath(DIR_FIG, mustWork = FALSE), ": ",
          length(list.files(DIR_FIG, pattern = "\\.pdf$")))
  if (length(failed))
    message("  did not complete: ", paste(failed, collapse = ", "))
  message("-----------------------------------------------------------")
  
  invisible(failed)
}

## ---------------------------------------------------------------------------
## A short console report of the quantities a reader will want to check first,
## computed from the same values that go into the tables.
## ---------------------------------------------------------------------------

report <- function() {
  
  message("\nReference levels")
  message(sprintf("  Exp 1 and 2: alpha = %.2f, R = %d, cut-off = %d, size bound = %.5f, coverage bound = %.5f",
                  ALPHA_12, R_SYN, rank_cutoff(R_SYN, ALPHA_12), SZB_12, BND_12))
  message(sprintf("  Exp 3      : alpha = %.2f, R = %d, cut-off = %d, coverage bound = %.5f",
                  ALPHA_3, R_SYN, rank_cutoff(R_SYN, ALPHA_3), BND_3))
  message(sprintf("  MC s.e. at the bound: %.4f at 1000 reps, %.4f at 100 reps",
                  mc_se(BND_12, 1000), mc_se(BND_12, 100)))
  
  message("\nExperiment 1: width against Mahalanobis")
  for (i in seq_len(nrow(EXP1_TGT)))
    message(sprintf("  %-6s Mahalanobis %.4f -> Penalized Wald %.4f, %.1f%% narrower, coverage %.3f",
                    EXP1_TGT$target[i], EXP1_TGT$mah[i], EXP1_TGT$pw[i],
                    EXP1_TGT$reduce[i], EXP1_TGT$cov[i]))
  message(sprintf("  joint coverage, all four methods: %.3f to %.3f (bound %.5f)",
                  min(EXP1$cov_jt), max(EXP1$cov_jt), BND_12))
  
  message("\nExperiment 2: empirical size")
  sz <- E2_OURS["0.0", ]
  message(sprintf("  ours   %.3f to %.3f  (bound %.5f) -- %s",
                  min(sz), max(sz), SZB_12,
                  if (max(sz) <= SZB_12) "all cells at or below the bound"
                  else "AT LEAST ONE CELL EXCEEDS THE BOUND"))
  message(sprintf("  repro  %.3f to %.3f", min(E2_REPRO["0.0", ]), max(E2_REPRO["0.0", ])))
  message(sprintf("  pb-adi %.3f to %.3f", min(E2_ADI["0.0", ]),   max(E2_ADI["0.0", ])))
  
  message("\nExperiment 3: box saturation")
  sat <- which(E3_WID_OURS >= E3_SAT, arr.ind = TRUE)
  message(sprintf("  %d of %d cells have width >= %.1f against a box of width %g",
                  nrow(sat), length(E3_WID_OURS), E3_SAT, E3_BOXW))
  if (nrow(sat))
    for (k in seq_len(nrow(sat)))
      message(sprintf("    eps = %s, n = %s : width %.2f",
                      rownames(E3_WID_OURS)[sat[k, 1]],
                      colnames(E3_WID_OURS)[sat[k, 2]],
                      E3_WID_OURS[sat[k, 1], sat[k, 2]]))
  message(sprintf("  coverage %.3f to %.3f against a bound of %.5f",
                  min(E3_COV_OURS), max(E3_COV_OURS), BND_3))
  
  message("\nClamping boundary at mu* = 3, hard clamp")
  message(sprintf("  Mahalanobis      joint coverage %.2f to %.2f",
                  min(S2$mah.hard$joint), max(S2$mah.hard$joint)))
  message(sprintf("  Penalized Wald   joint coverage %.2f to %.2f  <- %.1f s.e. below the bound at worst",
                  min(S2$eff.hard$joint), max(S2$eff.hard$joint),
                  (BND_12 - min(S2$eff.hard$joint)) / mc_se(BND_12, 100)))
  message(sprintf("  PB-ADI           joint coverage %.2f to %.2f",
                  min(S2$pb.hard$joint), max(S2$pb.hard$joint)))
  message("  soft clamp, same privacy cost")
  message(sprintf("  Penalized Wald   joint coverage %.2f to %.2f",
                  min(S2$eff.soft$joint), max(S2$eff.soft$joint)))
  message(sprintf("  PB-ADI           joint coverage %.2f to %.2f",
                  min(S2$pb.soft$joint), max(S2$pb.soft$joint)))
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## Autorun.  Sourcing the file produces everything and prints the report; set
##     options(results.autorun = FALSE)
## beforehand to load the functions and data without writing anything.
## ---------------------------------------------------------------------------

if (isTRUE(getOption("results.autorun", default = TRUE))) {
  message("writing tables and figures")
  make_all()
  report()
}
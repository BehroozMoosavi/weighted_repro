# ==============================================================================
# Script: exp1_R_point_cloud_sensitivity.R
# Description: Evaluates the sensitivity of the Efficient Repro / Operator-II 
#              confidence intervals to the primary repro point-cloud size (R) 
#              in Experiment 1. Compares empirical coverage against exact 
#              finite-R reference levels (R + 2 - cutoff) / (R + 1).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment & Thread Controls
# ------------------------------------------------------------------------------

Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ------------------------------------------------------------------------------
# 1. Experimental Parameters & Search Bounds
# ------------------------------------------------------------------------------

upper_clamp <- 3
lower_clamp <- 0
n           <- 100L
ep          <- 1
alpha       <- 0.05
tol         <- 1e-8
bisection_tol <- 0.1

population_mu    <- 1
population_sigma <- 1

h_fd <- 1e-3

mu_search_lower <- -2
mu_search_upper <- 5
sa_search_lower <- 0.1
sa_search_upper <- 5

# ------------------------------------------------------------------------------
# 2. Package Dependencies & Core Engine Utilities
# ------------------------------------------------------------------------------

req <- c("foreach", "doSNOW")
miss <- req[!vapply(req, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "))

suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doSNOW))

parse_first_int <- function(x) {
  if (!nzchar(x)) return(NA_integer_)
  y <- suppressWarnings(as.integer(sub("[^0-9].*$", "", x)))
  if (length(y) == 0L || is.na(y) || y < 1L) return(NA_integer_)
  y
}

slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_PER_TASK", unset = ""))
if (is.na(slurm_cpus)) slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_ON_NODE", unset = ""))
if (is.na(slurm_cpus)) slurm_cpus <- parallel::detectCores(logical = FALSE)
if (is.na(slurm_cpus) || slurm_cpus < 1L) slurm_cpus <- 1L

sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
  raw <- sa * data_randomness + mu
  clamped_data <- pmax(pmin(raw, upper_clamp), lower_clamp)
  row_mean <- rowMeans(clamped_data)
  n_cols <- ncol(clamped_data)
  centered_ss <- pmax(rowSums(clamped_data^2) - n_cols * row_mean^2, 0)
  row_var <- centered_ss / (n_cols - 1)
  mean_noise_scale <- (upper_clamp - lower_clamp) / (n * ep)
  var_noise_scale <- (upper_clamp - lower_clamp)^2 / (n * ep)
  cbind(
    row_mean + mean_noise_scale * privacy_noises[, 1],
    row_var  + var_noise_scale  * privacy_noises[, 2]
  )
}

safe_inv2 <- function(A, ridge = 1e-10) {
  A <- 0.5 * (A + t(A)) + ridge * diag(2L)
  a <- A[1, 1]; b <- A[1, 2]; d <- A[2, 2]
  det_A <- a * d - b^2
  if (!is.finite(det_A) || det_A <= ridge^2) {
    eg <- eigen(A, symmetric = TRUE)
    ev <- pmax(eg$values, ridge)
    return(eg$vectors %*% (t(eg$vectors) / ev))
  }
  matrix(c(d, -b, -b, a), nrow = 2L) / det_A
}

compute_aux_derivatives <- function(mu, sa, aux_dr, aux_pn) {
  h <- h_fd
  synth_center <- sdp_vec(aux_dr, aux_pn, sa, mu)
  
  if (sa - h > sa_search_lower && sa + h < sa_search_upper) {
    plus <- sdp_vec(aux_dr, aux_pn, sa + h, mu)
    minus <- sdp_vec(aux_dr, aux_pn, sa - h, mu)
    d_sa <- (colMeans(plus) - colMeans(minus)) / (2 * h)
  } else if (sa - h <= sa_search_lower) {
    plus <- sdp_vec(aux_dr, aux_pn, sa + h, mu)
    d_sa <- (colMeans(plus) - colMeans(synth_center)) / h
  } else {
    minus <- sdp_vec(aux_dr, aux_pn, sa - h, mu)
    d_sa <- (colMeans(synth_center) - colMeans(minus)) / h
  }
  
  if (mu - h > mu_search_lower && mu + h < mu_search_upper) {
    plus <- sdp_vec(aux_dr, aux_pn, sa, mu + h)
    minus <- sdp_vec(aux_dr, aux_pn, sa, mu - h)
    d_mu <- (colMeans(plus) - colMeans(minus)) / (2 * h)
  } else if (mu - h <= mu_search_lower) {
    plus <- sdp_vec(aux_dr, aux_pn, sa, mu + h)
    d_mu <- (colMeans(plus) - colMeans(synth_center)) / h
  } else {
    minus <- sdp_vec(aux_dr, aux_pn, sa, mu - h)
    d_mu <- (colMeans(synth_center) - colMeans(minus)) / h
  }
  
  S <- cov(synth_center)
  S <- 0.5 * (S + t(S)) + 1e-8 * diag(2L)
  list(
    d_mu = matrix(d_mu, ncol = 1L),
    d_sa = matrix(d_sa, ncol = 1L),
    cov_inv = safe_inv2(S, 1e-8)
  )
}

efficient_direction <- function(mu, sa, interest, aux_dr, aux_pn, cache) {
  key <- sprintf("%s_%.10g_%.10g", interest, mu, sa)
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))
  
  deriv <- compute_aux_derivatives(mu, sa, aux_dr, aux_pn)
  if (interest == "mu") {
    nuisance <- deriv$d_sa
    target_raw <- deriv$d_mu
  } else {
    nuisance <- deriv$d_mu
    target_raw <- deriv$d_sa
  }
  
  nuisance_norm <- as.numeric(t(nuisance) %*% deriv$cov_inv %*% nuisance)
  if (!is.finite(nuisance_norm) || nuisance_norm <= 1e-10) {
    direction <- c(NA_real_, NA_real_)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  coef <- as.numeric(t(nuisance) %*% deriv$cov_inv %*% target_raw) / nuisance_norm
  target <- target_raw - coef * nuisance
  direction <- as.numeric(deriv$cov_inv %*% target)
  if (any(!is.finite(direction)) || sum(direction^2) < 1e-20) direction <- c(NA_real_, NA_real_)
  assign(key, direction, envir = cache)
  direction
}

efficient_depth <- function(synth, mu, sa, interest, aux_dr, aux_pn, cache, lambda_n) {
  direction <- efficient_direction(mu, sa, interest, aux_dr, aux_pn, cache)
  center <- colMeans(synth)
  centered <- sweep(synth, 2L, center, "-")
  covariance <- crossprod(centered) / nrow(synth)
  covariance <- 0.5 * (covariance + t(covariance)) + 1e-10 * diag(2L)
  cov_inv <- tryCatch(safe_inv2(covariance, 1e-10), error = function(e) NULL)
  eff_scale <- if (all(is.finite(direction))) as.numeric(t(direction) %*% covariance %*% direction) else NA_real_
  if (any(!is.finite(direction)) || is.null(cov_inv) || !is.finite(eff_scale) || eff_scale <= 1e-14) {
    return(rep(NA_real_, nrow(synth)))
  }
  projection <- as.numeric(centered %*% direction)
  q_eff <- projection^2 / eff_scale
  q_mah <- rowSums((centered %*% cov_inv) * centered)
  depth <- 1 / (1 + q_eff + lambda_n * q_mah)
  depth[!is.finite(depth)] <- 0
  depth
}

score_mu_sa <- function(optim_par, data_randomness, privacy_noises, dp_statistic,
                        depth_type, aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic) {
  mu <- optim_par[1]; sa <- optim_par[2]
  interest <- if (depth_type == "efficient_mu") "mu" else "sigma"
  synth <- rbind(sdp_vec(data_randomness, privacy_noises, sa, mu), dp_statistic)
  D <- efficient_depth(synth, mu, sa, interest, aux_dr, aux_pn, eff_cache, lambda_n)
  if (any(!is.finite(D))) return(NA_real_)
  obs <- R_synthetic + 1L
  r <- rank(D, ties.method = "max")[obs]
  -(r + D[obs])
}

score_sa_mu <- function(optim_par, data_randomness, privacy_noises, dp_statistic,
                        depth_type, aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic) {
  sa <- optim_par[1]; mu <- optim_par[2]
  score_mu_sa(c(mu, sa), data_randomness, privacy_noises, dp_statistic, depth_type,
              aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic)
}

accept <- function(optim_par, data_randomness, privacy_noises, dp_statistic,
                   search_lower, search_upper, nuisance_lower, nuisance_upper,
                   depth_type, score_func, aux_dr, aux_pn, eff_cache,
                   lambda_n, R_synthetic, acceptance_threshold) {
  val <- score_func(optim_par, data_randomness, privacy_noises, dp_statistic,
                    depth_type, aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic)
  if (is.finite(val) && (-val) >= acceptance_threshold) return(optim_par)
  
  opt <- tryCatch(
    optim(
      par = optim_par, fn = score_func, method = "L-BFGS-B",
      lower = c(search_lower, nuisance_lower), upper = c(search_upper, nuisance_upper),
      data_randomness = data_randomness, privacy_noises = privacy_noises,
      dp_statistic = dp_statistic, depth_type = depth_type, aux_dr = aux_dr,
      aux_pn = aux_pn, eff_cache = eff_cache, lambda_n = lambda_n,
      R_synthetic = R_synthetic
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$value)) return(c(NA_real_, NA_real_))
  if ((-opt$value) >= acceptance_threshold) return(opt$par)
  c(NA_real_, NA_real_)
}

getConfidenceInterval <- function(optim_par, dp_statistic, data_randomness, privacy_noises,
                                  search_lower, search_upper, nuisance_lower, nuisance_upper,
                                  depth_type, score_func, aux_dr, aux_pn, eff_cache,
                                  lambda_n, R_synthetic, acceptance_threshold) {
  optim_par <- accept(
    optim_par, data_randomness, privacy_noises, dp_statistic,
    search_lower, search_upper, nuisance_lower, nuisance_upper,
    depth_type, score_func, aux_dr, aux_pn, eff_cache,
    lambda_n, R_synthetic, acceptance_threshold
  )
  if (is.na(optim_par[1])) return(c(NA_real_, NA_real_))
  t_mid <- optim_par[1]
  
  l_low <- search_lower; l_up <- t_mid - tol
  while (l_up - l_low > bisection_tol) {
    l_mid <- (l_low + l_up) / 2
    a <- accept(
      c(l_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic,
      l_low, l_mid, nuisance_lower, nuisance_upper, depth_type, score_func,
      aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic, acceptance_threshold
    )
    if (!is.na(a[1])) l_up <- l_mid - tol else l_low <- l_mid
  }
  
  r_low <- t_mid + tol; r_up <- search_upper
  while (r_up - r_low > bisection_tol) {
    r_mid <- (r_low + r_up) / 2
    a <- accept(
      c(r_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic,
      r_mid, r_up, nuisance_lower, nuisance_upper, depth_type, score_func,
      aux_dr, aux_pn, eff_cache, lambda_n, R_synthetic, acceptance_threshold
    )
    if (!is.na(a[1])) r_low <- r_mid + tol else r_up <- r_mid
  }
  c(l_low, r_up)
}

run_mc <- function(nSIM, R_synthetic, R_aux, lambda_n, seed_offset = 123L) {
  acceptance_threshold <- floor(alpha * (R_synthetic + 1L)) + 1L
  n_workers <- min(124L, slurm_cpus, nSIM)
  
  cl <- parallel::makePSOCKcluster(n_workers)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = nSIM, style = 3)
  opts <- list(progress = function(k) setTxtProgressBar(pb, k))
  
  export_names <- c(
    "upper_clamp", "lower_clamp", "n", "ep", "alpha", "tol", "bisection_tol",
    "population_mu", "population_sigma", "h_fd", "mu_search_lower", "mu_search_upper",
    "sa_search_lower", "sa_search_upper", "R_synthetic", "R_aux", "lambda_n",
    "acceptance_threshold", "seed_offset", "sdp_vec", "safe_inv2",
    "compute_aux_derivatives", "efficient_direction", "efficient_depth",
    "score_mu_sa", "score_sa_mu", "accept", "getConfidenceInterval"
  )
  parallel::clusterExport(cl, export_names, envir = environment())
  
  ans <- foreach(s = seq_len(nSIM), .combine = "rbind", .options.snow = opts) %dopar% {
    set.seed(s + seed_offset)
    
    raw_dr <- matrix(rnorm(n), 1L, n)
    raw_pn <- matrix(rnorm(2L), 1L, 2L)
    dp_statistic <- as.numeric(sdp_vec(raw_dr, raw_pn, population_sigma, population_mu))
    
    data_randomness <- matrix(rnorm(n * R_synthetic), R_synthetic, n)
    privacy_noises <- matrix(rnorm(2L * R_synthetic), R_synthetic, 2L)
    aux_dr <- matrix(rnorm(n * R_aux), R_aux, n)
    aux_pn <- matrix(rnorm(2L * R_aux), R_aux, 2L)
    
    cache_mu <- new.env(parent = emptyenv())
    mu_ci <- getConfidenceInterval(
      c(1, 1), dp_statistic, data_randomness, privacy_noises,
      mu_search_lower, mu_search_upper, sa_search_lower, sa_search_upper,
      "efficient_mu", score_mu_sa, aux_dr, aux_pn, cache_mu,
      lambda_n, R_synthetic, acceptance_threshold
    )
    if (any(is.na(mu_ci))) {
      cov_mu <- 0; width_mu <- NA_real_; fail_mu <- 1L
    } else {
      cov_mu <- as.numeric(mu_ci[1] <= population_mu && population_mu <= mu_ci[2])
      width_mu <- diff(mu_ci); fail_mu <- 0L
    }
    
    cache_sigma <- new.env(parent = emptyenv())
    sigma_ci <- getConfidenceInterval(
      c(1, 1), dp_statistic, data_randomness, privacy_noises,
      sa_search_lower, sa_search_upper, mu_search_lower, mu_search_upper,
      "efficient_sigma", score_sa_mu, aux_dr, aux_pn, cache_sigma,
      lambda_n, R_synthetic, acceptance_threshold
    )
    if (any(is.na(sigma_ci))) {
      cov_sigma <- 0; width_sigma <- NA_real_; fail_sigma <- 1L
    } else {
      cov_sigma <- as.numeric(sigma_ci[1] <= population_sigma && population_sigma <= sigma_ci[2])
      width_sigma <- diff(sigma_ci); fail_sigma <- 0L
    }
    
    c(cov_mu = cov_mu, width_mu = width_mu, fail_mu = fail_mu,
      cov_sigma = cov_sigma, width_sigma = width_sigma, fail_sigma = fail_sigma)
  }
  
  try(close(pb), silent = TRUE)
  parallel::stopCluster(cl)
  ans
}

summarize_mc <- function(ans) {
  cm <- ans[, "cov_mu"]; wm <- ans[, "width_mu"]; fm <- ans[, "fail_mu"]
  cs <- ans[, "cov_sigma"]; ws <- ans[, "width_sigma"]; fs <- ans[, "fail_sigma"]
  p_mu <- mean(cm); p_sig <- mean(cs)
  nwm <- sum(is.finite(wm)); nws <- sum(is.finite(ws))
  data.frame(
    coverage_mu = p_mu,
    coverage_mu_se = sqrt(p_mu * (1 - p_mu) / length(cm)),
    width_mu = mean(wm, na.rm = TRUE),
    width_mu_se = if (nwm >= 2L) sd(wm, na.rm = TRUE) / sqrt(nwm) else NA_real_,
    failure_rate_mu = mean(fm),
    coverage_sigma = p_sig,
    coverage_sigma_se = sqrt(p_sig * (1 - p_sig) / length(cs)),
    width_sigma = mean(ws, na.rm = TRUE),
    width_sigma_se = if (nws >= 2L) sd(ws, na.rm = TRUE) / sqrt(nws) else NA_real_,
    failure_rate_sigma = mean(fs)
  )
}

# ------------------------------------------------------------------------------
# 3. Graphics & Publication Utilities
# ------------------------------------------------------------------------------

COL_MU     <- "#0072B2"
COL_SIGMA  <- "#D55E00"
COL_REF    <- "#3D3D3D"
COL_FINITE <- "#CC79A7"
COL_SELECT <- "#009E73"
COL_GRID   <- "#E8E8E8"
COL_AXIS   <- "#4A4A4A"

PCH_MU     <- 21L
PCH_SIGMA  <- 24L
PCH_REF    <- 22L

FIG_FAMILY    <- "Helvetica"
FIG_POINTSIZE <- 9
FIG_DPI       <- 600

FIG_CAIRO         <- isTRUE(capabilities("cairo"))
FIG_UNICODE_GREEK <- FIG_CAIRO
FIG_TIKZ          <- requireNamespace("tikzDevice", quietly = TRUE)

MGP_TICK  <- 0.5
YLAB_GAP  <- 0.55
YLAB_TEXT <- 1.15
XLAB_LINE <- 1.95

fig_symbol_report <- function() {
  route <- if (isTRUE(FIG_TIKZ)) "tikz (.tex, exact LaTeX \\mu)"
  else if (isTRUE(FIG_UNICODE_GREEK)) "cairo (italic Unicode mu)"
  else "base (UPRIGHT Symbol mu -- use \\upmu in captions)"
  cat("Figure symbol route:", route, "\n")
  invisible(route)
}

alpha_col <- function(col, a) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(m[1L], m[2L], m[3L], alpha = a)
}

tikz_on <- function() isTRUE(getOption("exp1.tikz", FALSE))

pm <- function(template) {
  reps <- if (isTRUE(FIG_UNICODE_GREEK))
    c(mu = 'italic("\u03bc")', sigma = 'italic("\u03c3")', lambda = 'italic("\u03bb")')
  else
    c(mu = 'mu', sigma = 'sigma', lambda = 'lambda')
  for (k in names(reps))
    template <- gsub(paste0("{", k, "}"), reps[[k]], template, fixed = TRUE)
  parse(text = template)
}

as_plotmath <- function(x) if (is.character(x)) pm(x) else x

pick <- function(math, tex) {
  if (tikz_on() && !is.null(tex)) tex else as_plotmath(math)
}

append_lab <- function(v, math, tex) if (tikz_on()) c(v, tex) else c(v, math)

pct_labels <- function(v, dec = 0L) {
  fmt <- if (tikz_on()) paste0("%.", dec, "f\\%%") else paste0("%.", dec, "f%%")
  sprintf(fmt, 100 * v)
}

label_lines <- function(labels, cex = 0.92) {
  labels <- gsub("\\\\", "", as.character(labels))
  if (!length(labels)) return(0)
  w <- vapply(strsplit(labels, ""), function(ch) {
    sum(ifelse(ch == ".", 0.278,
               ifelse(ch == "%", 0.889,
                      ifelse(ch == "-", 0.333, 0.556))))
  }, numeric(1))
  max(w) * cex / 1.2
}

mar_left_for <- function(...) {
  max(vapply(list(...), function(l) label_lines(l), numeric(1))) +
    MGP_TICK + YLAB_GAP + YLAB_TEXT
}

ylab_line_for <- function(labels) MGP_TICK + label_lines(labels) + YLAB_GAP

panel_titles <- function(xlab, ylab, ylab_line) {
  mtext(xlab, side = 1, line = XLAB_LINE, cex = par("cex.lab"), col = "black")
  mtext(ylab, side = 2, line = ylab_line, las = 0, cex = par("cex.lab"), col = "black")
}

save_figure <- function(draw_fun, pdf_file, png_file, width, height,
                        pointsize = FIG_POINTSIZE, dpi = FIG_DPI,
                        tex_file = sub("\\.pdf$", ".tex", pdf_file)) {
  options(exp1.tikz = FALSE)
  if (isTRUE(FIG_CAIRO)) {
    grDevices::cairo_pdf(pdf_file, width = width, height = height,
                         pointsize = pointsize, family = FIG_FAMILY)
    draw_fun()
    invisible(grDevices::dev.off())
  } else {
    grDevices::pdf(pdf_file, width = width, height = height,
                   family = FIG_FAMILY, pointsize = pointsize, useDingbats = FALSE)
    draw_fun()
    invisible(grDevices::dev.off())
    try(grDevices::embedFonts(pdf_file), silent = TRUE)
  }
  
  png_args <- list(filename = png_file,
                   width  = round(width  * dpi),
                   height = round(height * dpi),
                   res = dpi, pointsize = pointsize, bg = "white")
  if (isTRUE(FIG_CAIRO)) png_args$type <- "cairo"
  do.call(grDevices::png, png_args)
  draw_fun()
  invisible(grDevices::dev.off())
  
  tex_written <- NA_character_
  if (isTRUE(FIG_TIKZ) && !is.null(tex_file)) {
    options(exp1.tikz = TRUE)
    ok <- tryCatch({
      tikzDevice::tikz(tex_file, width = width, height = height,
                       pointsize = pointsize, standAlone = FALSE, sanitize = FALSE)
      draw_fun()
      invisible(grDevices::dev.off())
      TRUE
    }, error = function(e) {
      try(while (grDevices::dev.cur() > 1L) grDevices::dev.off(), silent = TRUE)
      warning("tikz output skipped: ", conditionMessage(e), call. = FALSE)
      FALSE
    })
    options(exp1.tikz = FALSE)
    if (isTRUE(ok)) tex_written <- tex_file
  }
  invisible(c(pdf = pdf_file, png = png_file, tex = tex_written))
}

fig_par <- function(mar) {
  par(mar = mar, mgp = c(3, MGP_TICK, 0), tcl = -0.22, las = 1,
      family = "sans", cex.axis = 0.92, cex.lab = 1.0, cex.main = 1.0,
      xaxs = "i", yaxs = "i", lend = "butt",
      fg = COL_AXIS, col.axis = COL_AXIS, col.lab = "black")
}

errbar <- function(x, lo, hi, col, lwd = 1.2, cap_frac = 0.010) {
  w <- cap_frac * diff(par("usr")[1:2])
  segments(x, lo, x, hi, col = col, lwd = lwd)
  segments(x - w, lo, x + w, lo, col = col, lwd = lwd)
  segments(x - w, hi, x + w, hi, col = col, lwd = lwd)
  invisible(NULL)
}

inner_ticks <- function(lim, n = 5L) {
  at <- pretty(lim, n = n)
  at[at >= lim[1L] & at <= lim[2L]]
}

infer_nsim <- function(p, se) {
  ok <- is.finite(p) & is.finite(se) & se > 0 & p > 0 & p < 1
  if (!any(ok)) return(NA_real_)
  round(stats::median(p[ok] * (1 - p[ok]) / se[ok]^2))
}

grid_positions <- function(x, even_spacing = TRUE) {
  xv <- as.numeric(x)
  list(xv = xv,
       xp = if (isTRUE(even_spacing)) seq_along(xv) else xv,
       labels = format(xv, trim = TRUE, drop0trailing = TRUE))
}

pad_range <- function(xp) {
  xr <- range(xp, finite = TRUE)
  p <- if (diff(xr) > 0) 0.07 * diff(xr) else 0.5
  c(xr[1L] - p, xr[2L] + p)
}

plot_two_panel <- function(results, x, xlab, pdf_file, png_file,
                           nominal = 0.95, reference = NULL,
                           ref_label = "finite-R reference",
                           ref_label_pm = NULL, ref_label_tex = NULL,
                           xlab_tex = NULL, highlight = NULL,
                           even_spacing = TRUE, width_from_zero = TRUE,
                           fig_width = 7.2, fig_height = 4.1) {
  g <- grid_positions(x, even_spacing)
  xv <- g$xv; xp <- g$xp; xlabels <- g$labels
  xlim <- pad_range(xp)
  
  hp <- NULL
  if (!is.null(highlight)) {
    idx <- which(abs(xv - as.numeric(highlight)[1L]) < 1e-12)
    if (length(idx)) hp <- xp[idx[1L]]
  }
  refv <- if (is.null(reference)) NULL else as.numeric(reference)
  
  cm <- results$coverage_mu;    cms <- results$coverage_mu_se
  cs <- results$coverage_sigma; css <- results$coverage_sigma_se
  mu_lo <- pmax(0, cm - 1.96 * cms); mu_hi <- pmin(1, cm + 1.96 * cms)
  sg_lo <- pmax(0, cs - 1.96 * css); sg_hi <- pmin(1, cs + 1.96 * css)
  
  nrep <- infer_nsim(c(cm, cs), c(cms, css))
  mc_half <- if (is.finite(nrep) && nrep > 0)
    1.96 * sqrt(nominal * (1 - nominal) / nrep) else NA_real_
  
  ycand <- c(mu_lo, mu_hi, sg_lo, sg_hi, nominal, refv)
  if (is.finite(mc_half)) ycand <- c(ycand, nominal - mc_half, nominal + mc_half)
  ylo <- min(ycand, na.rm = TRUE); yhi <- max(ycand, na.rm = TRUE)
  ypad <- max(0.005, 0.11 * (yhi - ylo))
  ylim_a <- c(max(0, ylo - ypad), min(1, yhi + ypad))
  
  yat_a <- inner_ticks(ylim_a, 5L)
  dec_a <- if (length(yat_a) >= 2L && min(diff(yat_a)) < 0.01) 1L else 0L
  
  wm <- results$width_mu;    wms <- results$width_mu_se
  ws <- results$width_sigma; wss <- results$width_sigma_se
  wm_lo <- pmax(0, wm - 1.96 * wms); wm_hi <- wm + 1.96 * wms
  ws_lo <- pmax(0, ws - 1.96 * wss); ws_hi <- ws + 1.96 * wss
  
  wtop <- max(c(wm_hi, ws_hi, wm, ws), na.rm = TRUE)
  wbot <- min(c(wm_lo, ws_lo, wm, ws), na.rm = TRUE)
  ylim_b <- if (isTRUE(width_from_zero)) c(0, 1.07 * wtop) else {
    p <- max(1e-6, 0.12 * (wtop - wbot)); c(max(0, wbot - p), wtop + p)
  }
  yat_b <- inner_ticks(ylim_b, 5L)
  
  draw <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    labs_a <- pct_labels(yat_a, dec_a)
    labs_b <- format(yat_b, trim = TRUE)
    ml <- mar_left_for(labs_a, labs_b)
    xl <- pick(xlab, xlab_tex)
    
    layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L, byrow = TRUE), heights = c(1, lcm(1.0)))
    par(cex = 1)
    fig_par(mar = c(3.4, ml, 1.7, 0.9))
    
    # Panel (a): Empirical coverage
    plot.new(); plot.window(xlim = xlim, ylim = ylim_a)
    if (is.finite(mc_half))
      rect(xlim[1L], nominal - mc_half, xlim[2L], nominal + mc_half,
           col = alpha_col(COL_REF, 0.09), border = NA)
    abline(h = yat_a, col = COL_GRID, lwd = 0.6)
    abline(v = xp,    col = COL_GRID, lwd = 0.6)
    if (!is.null(hp))
      abline(v = hp, col = alpha_col(COL_SELECT, 0.75), lwd = 1.2, lty = 3)
    abline(h = nominal, col = COL_REF, lty = 2, lwd = 1.1)
    if (!is.null(refv)) {
      lines(xp, refv, col = COL_FINITE, lwd = 1.2)
      points(xp, refv, pch = PCH_REF, col = COL_FINITE, bg = "white", cex = 0.85, lwd = 0.9)
    }
    
    errbar(xp, mu_lo, mu_hi, COL_MU)
    errbar(xp, sg_lo, sg_hi, COL_SIGMA)
    lines(xp, cm, col = COL_MU,    lwd = 1.8)
    lines(xp, cs, col = COL_SIGMA, lwd = 1.8)
    points(xp, cm, pch = PCH_MU,    bg = COL_MU,    col = "white", cex = 1.15, lwd = 0.7)
    points(xp, cs, pch = PCH_SIGMA, bg = COL_SIGMA, col = "white", cex = 1.15, lwd = 0.7)
    
    axis(1, at = xp,    labels = xlabels, lwd = 0, lwd.ticks = 0.7)
    axis(2, at = yat_a, labels = labs_a,  lwd = 0, lwd.ticks = 0.7)
    box(bty = "l", col = COL_AXIS, lwd = 0.8)
    panel_titles(xl, "Empirical coverage", ylab_line_for(labs_a))
    mtext("(a)", side = 3, adj = 0, line = 0.35, font = 2, cex = 0.98, col = "black")
    
    # Panel (b): Average CI width
    par(cex = 1)
    plot.new(); plot.window(xlim = xlim, ylim = ylim_b)
    abline(h = yat_b, col = COL_GRID, lwd = 0.6)
    abline(v = xp,    col = COL_GRID, lwd = 0.6)
    if (!is.null(hp))
      abline(v = hp, col = alpha_col(COL_SELECT, 0.75), lwd = 1.2, lty = 3)
    
    errbar(xp, wm_lo, wm_hi, COL_MU)
    errbar(xp, ws_lo, ws_hi, COL_SIGMA)
    lines(xp, wm, col = COL_MU,    lwd = 1.8)
    lines(xp, ws, col = COL_SIGMA, lwd = 1.8)
    points(xp, wm, pch = PCH_MU,    bg = COL_MU,    col = "white", cex = 1.15, lwd = 0.7)
    points(xp, ws, pch = PCH_SIGMA, bg = COL_SIGMA, col = "white", cex = 1.15, lwd = 0.7)
    
    axis(1, at = xp,    labels = xlabels, lwd = 0, lwd.ticks = 0.7)
    axis(2, at = yat_b, labels = labs_b,  lwd = 0, lwd.ticks = 0.7)
    box(bty = "l", col = COL_AXIS, lwd = 0.8)
    panel_titles(xl, "Average CI width", ylab_line_for(labs_b))
    mtext("(b)", side = 3, adj = 0, line = 0.35, font = 2, cex = 0.98, col = "black")
    
    # Shared legend row
    leg_txt <- if (tikz_on())
      c("$\\mu$ of interest", "$\\sigma$ of interest", sprintf("nominal $%g\\%%$", 100 * nominal))
    else
      c(pm('{mu}*" of interest"'), pm('{sigma}*" of interest"'),
        parse(text = sprintf('"nominal %g%%"', 100 * nominal)))
    leg_col <- c(COL_MU, COL_SIGMA, COL_REF)
    leg_pch <- c(PCH_MU, PCH_SIGMA, NA)
    leg_bg  <- c(COL_MU, COL_SIGMA, NA)
    leg_lty <- c(1, 1, 2)
    leg_lwd <- c(1.8, 1.8, 1.1)
    
    if (!is.null(refv)) {
      leg_txt <- append_lab(
        leg_txt,
        if (is.null(ref_label_pm)) as.expression(ref_label) else pm(ref_label_pm),
        if (is.null(ref_label_tex)) ref_label else ref_label_tex)
      leg_col <- c(leg_col, COL_FINITE)
      leg_pch <- c(leg_pch, PCH_REF); leg_bg <- c(leg_bg, "white")
      leg_lty <- c(leg_lty, 1);       leg_lwd <- c(leg_lwd, 1.2)
    }
    if (!is.null(hp)) {
      leg_txt <- append_lab(leg_txt, expression("selected value"), "selected value")
      leg_col <- c(leg_col, COL_SELECT)
      leg_pch <- c(leg_pch, NA); leg_bg <- c(leg_bg, NA)
      leg_lty <- c(leg_lty, 3);  leg_lwd <- c(leg_lwd, 1.2)
    }
    
    par(mar = c(0, 0, 0, 0), cex = 1)
    plot.new()
    legend("center", horiz = TRUE, bty = "n",
           cex = if (length(leg_txt) > 3L) 0.80 else 0.85,
           legend = leg_txt, col = leg_col, pt.bg = leg_bg, pch = leg_pch,
           lty = leg_lty, lwd = leg_lwd, pt.cex = 1.05,
           seg.len = 1.6, x.intersp = 0.7, text.col = "black")
  }
  
  save_figure(draw, pdf_file, png_file, width = fig_width, height = fig_height)
}

plot_failure_rate <- function(results, x, xlab, pdf_file, png_file,
                              xlab_tex = NULL, highlight = NULL, even_spacing = TRUE,
                              fig_width = 5.2, fig_height = 3.8) {
  g <- grid_positions(x, even_spacing)
  xv <- g$xv; xp <- g$xp; xlabels <- g$labels
  xlim <- pad_range(xp)
  
  hp <- NULL
  if (!is.null(highlight)) {
    idx <- which(abs(xv - as.numeric(highlight)[1L]) < 1e-12)
    if (length(idx)) hp <- xp[idx[1L]]
  }
  
  fm <- results$failure_rate_mu
  fs <- results$failure_rate_sigma
  ftop <- max(c(fm, fs, 0), na.rm = TRUE)
  ylim <- c(0, max(0.01, 1.10 * ftop))
  yat <- inner_ticks(ylim, 5L)
  
  draw <- function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    labs <- pct_labels(yat, 1L)
    fig_par(mar = c(3.4, mar_left_for(labs), 1.7, 0.9))
    
    plot.new(); plot.window(xlim = xlim, ylim = ylim)
    abline(h = yat, col = COL_GRID, lwd = 0.6)
    abline(v = xp,  col = COL_GRID, lwd = 0.6)
    if (!is.null(hp))
      abline(v = hp, col = alpha_col(COL_SELECT, 0.75), lwd = 1.2, lty = 3)
    
    lines(xp, fm, col = COL_MU,    lwd = 1.8)
    lines(xp, fs, col = COL_SIGMA, lwd = 1.8)
    points(xp, fm, pch = PCH_MU,    bg = COL_MU,    col = "white", cex = 1.15, lwd = 0.7)
    points(xp, fs, pch = PCH_SIGMA, bg = COL_SIGMA, col = "white", cex = 1.15, lwd = 0.7)
    
    axis(1, at = xp,  labels = xlabels, lwd = 0, lwd.ticks = 0.7)
    axis(2, at = yat, labels = labs,    lwd = 0, lwd.ticks = 0.7)
    box(bty = "l", col = COL_AXIS, lwd = 0.8)
    panel_titles(pick(xlab, xlab_tex), "Interval-construction failure rate", ylab_line_for(labs))
    
    legend("topright", bty = "n", cex = 0.82, text.col = "black",
           legend = if (tikz_on())
             c("$\\mu$ of interest", "$\\sigma$ of interest")
           else
             c(pm('{mu}*" of interest"'), pm('{sigma}*" of interest"')),
           col = c("white", "white"), pt.bg = c(COL_MU, COL_SIGMA),
           pch = c(PCH_MU, PCH_SIGMA), lty = c(1, 1), lwd = c(1.8, 1.8),
           pt.cex = 1.1, seg.len = 1.6, x.intersp = 0.75)
  }
  
  save_figure(draw, pdf_file, png_file, width = fig_width, height = fig_height)
}

fig_symbol_report()

# ------------------------------------------------------------------------------
# 4. Point-Cloud R Sensitivity Execution
# ------------------------------------------------------------------------------

PROJECT_DIR <- path.expand("~/R_Simuls/sensitivity_analysis/exp1/R_point_cloud")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

nSIM     <- 1000L
R_aux    <- 400L
lambda_n <- 1 / log(n)
R_values <- c(50L, 100L, 200L, 400L, 800L)

res_list <- vector("list", length(R_values))
for (j in seq_along(R_values)) {
  R_synthetic <- R_values[j]
  cutoff <- floor(alpha * (R_synthetic + 1L)) + 1L
  finite_R_reference <- (R_synthetic + 2L - cutoff) / (R_synthetic + 1L)
  cat("\n### R_synthetic =", R_synthetic,
      " cutoff =", cutoff,
      " finite-R ref =", finite_R_reference, "###\n")
  raw <- run_mc(
    nSIM = nSIM, R_synthetic = R_synthetic, R_aux = R_aux,
    lambda_n = lambda_n, seed_offset = 123L
  )
  write.csv(raw, file.path(PROJECT_DIR, sprintf("R_synthetic_%d_raw.csv", R_synthetic)), row.names = FALSE)
  sm <- summarize_mc(raw)
  sm$R_synthetic <- R_synthetic
  sm$rank_cutoff <- cutoff
  sm$finite_R_reference <- finite_R_reference
  sm$R_aux <- R_aux
  sm$lambda_n <- lambda_n
  res_list[[j]] <- sm
}

ALL_RESULTS <- do.call(rbind, res_list)
rownames(ALL_RESULTS) <- NULL
write.csv(ALL_RESULTS, file.path(PROJECT_DIR, "R_point_cloud_sensitivity_results.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Output Graphics & Captions Export
# ------------------------------------------------------------------------------

plot_two_panel(
  ALL_RESULTS,
  x = ALL_RESULTS$R_synthetic,
  xlab = '"Repro point-cloud size "*italic(R)',
  xlab_tex = "Repro point-cloud size $R$",
  pdf_file = file.path(PROJECT_DIR, "R_point_cloud_sensitivity_publication.pdf"),
  png_file = file.path(PROJECT_DIR, "R_point_cloud_sensitivity_publication.png"),
  nominal = 1 - alpha,
  reference = ALL_RESULTS$finite_R_reference,
  ref_label = "exact finite-R level",
  ref_label_pm = '"exact finite-"*italic(R)*" level"',
  ref_label_tex = "exact finite-$R$ level"
)

plot_failure_rate(
  ALL_RESULTS,
  x = ALL_RESULTS$R_synthetic,
  xlab = '"Repro point-cloud size "*italic(R)',
  xlab_tex = "Repro point-cloud size $R$",
  pdf_file = file.path(PROJECT_DIR, "R_point_cloud_failure_rate_publication.pdf"),
  png_file = file.path(PROJECT_DIR, "R_point_cloud_failure_rate_publication.png")
)

writeLines(c(
  sprintf(paste0("\\caption{Sensitivity of the efficient-depth repro-sample intervals to the ",
                 "repro point-cloud size $R$. Location--scale normal model with $n = %d$, ",
                 "$(\\mu^*, \\sigma^*) = (%g, %g)$, clamping to $[%g, %g]$, ",
                 "$\\varepsilon = %g$, $R_{\\mathrm{aux}} = %d$ auxiliary draws, ",
                 "$\\lambda_n = 1/\\log n = %.4f$, over %d Monte Carlo replicates per grid ",
                 "point. (a) Empirical coverage; error bars are $\\pm 1.96$ binomial standard ",
                 "errors, the shaded band is the Monte Carlo tolerance region around the ",
                 "nominal $%g\\%%$ level, and the purple curve is the exact finite-$R$ level ",
                 "$(R+2-c)/(R+1)$ with rank cutoff ",
                 "$c = \\lfloor \\alpha(R+1) \\rfloor + 1$, which is the benchmark the ",
                 "procedure actually targets at each $R$. (b) Average interval width with ",
                 "$\\pm 1.96$ standard errors. The horizontal axis uses $\\log_2$ spacing. ",
                 "Failed constructions are scored as non-coverage.}"),
          n, population_mu, population_sigma, lower_clamp, upper_clamp, ep,
          R_aux, lambda_n, nSIM, 100 * (1 - alpha)),
  "",
  paste0("% Rank cutoffs and exact finite-R levels: ",
         paste(sprintf("$R = %d$: $c = %d$, level $= %.4f$",
                       ALL_RESULTS$R_synthetic, ALL_RESULTS$rank_cutoff,
                       ALL_RESULTS$finite_R_reference), collapse = "; "), "."),
  "",
  sprintf(paste0("\\caption{Proportion of Monte Carlo replicates in which the interval ",
                 "construction returned no interval, for $\\mu$ and for $\\sigma$ as the ",
                 "parameter of interest, over %d replicates per grid point.}"), nSIM)
), file.path(PROJECT_DIR, "figure_captions.tex"))

cat("\n============================================================\n")
cat("POINT-CLOUD R SENSITIVITY COMPLETE\n")
print(ALL_RESULTS)

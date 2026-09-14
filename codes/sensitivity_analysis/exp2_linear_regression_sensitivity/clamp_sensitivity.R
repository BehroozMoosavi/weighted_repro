# ==============================================================================
# Script: clamp_sensitivity.R
# Description: Evaluates the sensitivity of the Efficient Repro test for 
#              H0: beta1 = 0 in private linear regression across varying clamping 
#              thresholds (Delta in {0.5, 1.0, 2.0, 5.0, 10.0}) and sample 
#              sizes (n in {100, 200, 500, 1000}).
#
# Regimes:
#   beta1_true = 0 : Empirical Type-I error rate
#   beta1_true = 1 : Empirical power
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Directory Initialization
# ------------------------------------------------------------------------------

PROJECT_DIR <- path.expand("~/R_Simuls/sensitivity_analysis/exp2")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

# ------------------------------------------------------------------------------
# 2. Fixed Parameters & Nuisance Bounds
# ------------------------------------------------------------------------------

mu_x     <- 0.5
tau_x    <- 1.0
beta0    <- -0.5
sa_eps   <- 0.5
alpha    <- 0.05

R_synthetic <- 200L
R_aux       <- 400L
h_fd        <- 1e-3
ep_fixed    <- 1.0
reps        <- 1000L

# Nuisance parameter vector: (beta0, mu, tau, sa)
nuisance_start <- c(beta0, mu_x, tau_x, sa_eps)
nuisance_lower <- c(-5.0, -5.0, 0.001, 0.001)
nuisance_upper <- c( 5.0,  5.0, 5.0,   5.0)

rank_threshold <- floor(alpha * (R_synthetic + 1L)) + 1L

cat(sprintf("R = %d\n", R_synthetic))
cat(sprintf("R_aux = %d\n", R_aux))
cat(sprintf("alpha = %.3f\n", alpha))
cat(sprintf("rank threshold = %d\n", rank_threshold))
cat(sprintf("reps = %d\n\n", reps))

# ------------------------------------------------------------------------------
# 3. Environment & Parallel Backend Setup
# ------------------------------------------------------------------------------

packages_needed <- c("foreach", "doSNOW", "parallelly")
packages_missing <- packages_needed[
  !(packages_needed %in% installed.packages()[, "Package"])
]

if (length(packages_missing) > 0L) {
  install.packages(packages_missing, dependencies = TRUE)
}

for (pkg in packages_needed) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)

available_cores <- suppressWarnings(
  as.integer(parallelly::availableCores(omit = 1))
)

if (!is.finite(available_cores) || available_cores < 1L) {
  available_cores <- 1L
}

n.cores <- max(1L, min(124L, available_cores, reps))
cat(sprintf("Using %d PSOCK workers.\n\n", n.cores))

cl <- parallel::makePSOCKcluster(n.cores)
doSNOW::registerDoSNOW(cl)

# ------------------------------------------------------------------------------
# 4. Clamping Helpers & Differentially Private Summaries
# ------------------------------------------------------------------------------

clamp_val <- function(x, a, b) {
  pmin(pmax(x, a), b)
}

clamp_nuisance <- function(x) {
  pmin(pmax(x, nuisance_lower), nuisance_upper)
}

sdp_vec_lr <- function(
    ux,
    uy,
    noise,
    delta,
    theta,
    ep,
    n
) {
  beta1_ <- theta[1]
  beta0_ <- theta[2]
  mu_    <- theta[3]
  tau_   <- theta[4]
  sa_    <- theta[5]
  
  x <- ux * tau_ + mu_
  y <- beta0_ + beta1_ * x + uy * sa_
  
  xc  <- clamp_val(x, -delta, delta)
  yc  <- clamp_val(y, -delta, delta)
  xyc <- clamp_val(x * y, -delta^2, delta^2)
  
  ep_split <- ep / sqrt(5)
  
  xbar  <- rowMeans(xc)        + (2 * delta   / (n * ep_split)) * noise[, 1]
  ybar  <- rowMeans(yc)        + (2 * delta   / (n * ep_split)) * noise[, 2]
  x2bar <- rowMeans(xc^2)      + (delta^2     / (n * ep_split)) * noise[, 3]
  y2bar <- rowMeans(yc^2)      + (delta^2     / (n * ep_split)) * noise[, 4]
  xybar <- rowMeans(xyc)       + (2 * delta^2 / (n * ep_split)) * noise[, 5]
  
  cbind(xbar, ybar, x2bar, y2bar, xybar)
}

# ------------------------------------------------------------------------------
# 5. Matrix Inversion Utilities
# ------------------------------------------------------------------------------

safe_inv_general <- function(A, ridge = 1e-8) {
  if (any(!is.finite(A))) {
    stop("Non-finite matrix.")
  }
  
  A <- 0.5 * (A + t(A))
  eg <- eigen(A, symmetric = TRUE)
  vals <- pmax(eg$values, ridge)
  
  eg$vectors %*% diag(1 / vals, nrow = length(vals)) %*% t(eg$vectors)
}

# ------------------------------------------------------------------------------
# 6. Efficient Direction via Independent Auxiliary Draws
# ------------------------------------------------------------------------------

efficient_direction_lr <- function(
    theta,
    aux_ux,
    aux_uy,
    aux_N,
    delta,
    ep,
    n,
    cache
) {
  key <- paste(sprintf("%.10g", theta), collapse = "_")
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  
  beta1_ <- theta[1]
  beta0_ <- theta[2]
  mu_    <- theta[3]
  tau_   <- theta[4]
  sa_    <- theta[5]
  
  if (tau_ <= 0 || sa_ <= 0 || any(!is.finite(theta))) {
    out <- rep(NA_real_, 5)
    assign(key, out, envir = cache)
    return(out)
  }
  
  h <- h_fd
  center_cloud <- sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, theta, ep, n)
  
  if (any(!is.finite(center_cloud))) {
    out <- rep(NA_real_, 5)
    assign(key, out, envir = cache)
    return(out)
  }
  
  center_mean <- colMeans(center_cloud)
  
  finite_diff <- function(idx) {
    tp <- theta
    tm <- theta
    tp[idx] <- tp[idx] + h
    tm[idx] <- tm[idx] - h
    
    if (idx %in% c(4, 5) && tm[idx] <= 0.001) {
      up <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tp, ep, n))
      return((up - center_mean) / h)
    }
    
    up <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tp, ep, n))
    dn <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tm, ep, n))
    (up - dn) / (2 * h)
  }
  
  d_beta1 <- finite_diff(1)
  D_eta   <- cbind(finite_diff(2), finite_diff(3), finite_diff(4), finite_diff(5))
  
  if (any(!is.finite(d_beta1)) || any(!is.finite(D_eta))) {
    out <- rep(NA_real_, 5)
    assign(key, out, envir = cache)
    return(out)
  }
  
  Sigma <- cov(center_cloud)
  Sigma <- 0.5 * (Sigma + t(Sigma)) + 1e-8 * diag(5)
  Sigma_inv <- tryCatch(safe_inv_general(Sigma, 1e-8), error = function(e) NULL)
  
  if (is.null(Sigma_inv)) {
    out <- rep(NA_real_, 5)
    assign(key, out, envir = cache)
    return(out)
  }
  
  DtSD <- t(D_eta) %*% Sigma_inv %*% D_eta
  DtSd <- t(D_eta) %*% Sigma_inv %*% d_beta1
  DtSD_inv <- tryCatch(safe_inv_general(DtSD, 1e-10), error = function(e) NULL)
  
  if (is.null(DtSD_inv)) {
    out <- rep(NA_real_, 5)
    assign(key, out, envir = cache)
    return(out)
  }
  
  coef <- DtSD_inv %*% DtSd
  target <- d_beta1 - D_eta %*% coef
  direction <- as.numeric(Sigma_inv %*% target)
  
  if (any(!is.finite(direction)) || sum(direction^2) < 1e-20) {
    direction <- rep(NA_real_, 5)
  }
  
  assign(key, direction, envir = cache)
  direction
}

# ------------------------------------------------------------------------------
# 7. Efficient Depth Evaluation
# ------------------------------------------------------------------------------

efficient_depth_lr <- function(
    synth,
    theta,
    aux_ux,
    aux_uy,
    aux_N,
    delta,
    ep,
    n,
    cache
) {
  direction <- efficient_direction_lr(theta, aux_ux, aux_uy, aux_N, delta, ep, n, cache)
  if (any(!is.finite(direction))) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  center <- colMeans(synth)
  centered <- sweep(synth, 2, center, "-")
  
  Sigma <- crossprod(centered) / nrow(synth)
  Sigma <- 0.5 * (Sigma + t(Sigma)) + 1e-10 * diag(5)
  Sigma_inv <- tryCatch(safe_inv_general(Sigma, 1e-10), error = function(e) NULL)
  if (is.null(Sigma_inv)) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  eff_scale <- as.numeric(t(direction) %*% Sigma %*% direction)
  if (!is.finite(eff_scale) || eff_scale <= 1e-14) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  projection <- as.numeric(centered %*% direction)
  Q_eff <- projection^2 / eff_scale
  Q_full <- rowSums((centered %*% Sigma_inv) * centered)
  
  lambda_current <- 1 / log(n)
  penalty <- Q_eff + lambda_current * Q_full
  depth <- 1 / (1 + penalty)
  
  if (any(!is.finite(depth))) {
    return(rep(NA_real_, nrow(synth)))
  }
  depth
}

# ------------------------------------------------------------------------------
# 8. Test Statistic Scoring & Nuisance Profiling
# ------------------------------------------------------------------------------

score_Eff_theta <- function(
    theta,
    ux,
    uy,
    noise,
    delta,
    s_dp,
    ep,
    n,
    aux_ux,
    aux_uy,
    aux_N,
    cache
) {
  synth <- tryCatch(
    sdp_vec_lr(ux, uy, noise, delta, theta, ep, n),
    error = function(e) NULL
  )
  if (is.null(synth) || any(!is.finite(synth))) {
    return(NA_real_)
  }
  
  synth_full <- rbind(synth, s_dp)
  depth <- efficient_depth_lr(synth_full, theta, aux_ux, aux_uy, aux_N, delta, ep, n, cache)
  if (any(!is.finite(depth))) {
    return(NA_real_)
  }
  
  obs_idx <- nrow(synth_full)
  r <- rank(depth, ties.method = "max")[obs_idx]
  depth_obs <- depth[obs_idx]
  if (!is.finite(r) || !is.finite(depth_obs)) {
    return(NA_real_)
  }
  
  -(r + depth_obs)
}

score_Eff_nuisance <- function(
    nuisance,
    beta1,
    ux,
    uy,
    noise,
    delta,
    s_dp,
    ep,
    n,
    aux_ux,
    aux_uy,
    aux_N,
    cache
) {
  theta <- c(beta1, nuisance)
  score <- score_Eff_theta(
    theta, ux, uy, noise, delta, s_dp, ep, n, aux_ux, aux_uy, aux_N, cache
  )
  if (!is.finite(score)) return(1e12)
  score
}

make_data_anchor <- function(beta1, s_dp) {
  xbar  <- s_dp[1]
  ybar  <- s_dp[2]
  x2bar <- s_dp[3]
  y2bar <- s_dp[4]
  
  mu_hat <- xbar
  var_x_hat <- x2bar - xbar^2
  if (!is.finite(var_x_hat)) var_x_hat <- 1
  
  tau_hat <- sqrt(max(var_x_hat, 0.001^2))
  beta0_hat <- ybar - beta1 * mu_hat
  
  var_y_hat <- y2bar - ybar^2
  if (!is.finite(var_y_hat)) var_y_hat <- sa_eps^2
  
  sa2_hat <- var_y_hat - beta1^2 * tau_hat^2
  sa_hat  <- sqrt(max(sa2_hat, 0.001^2))
  
  clamp_nuisance(c(beta0_hat, mu_hat, tau_hat, sa_hat))
}

# ------------------------------------------------------------------------------
# 9. Multi-Start Acceptance Search Routine
# ------------------------------------------------------------------------------

accept_Efficient <- function(
    beta1,
    nuisance_start_,
    ux,
    uy,
    noise,
    delta,
    s_dp,
    ep,
    n,
    aux_ux,
    aux_uy,
    aux_N,
    cache
) {
  anchor_start  <- make_data_anchor(beta1, s_dp)
  neutral_start <- c(0, 0, 1, 0.5)
  
  starts <- rbind(
    clamp_nuisance(nuisance_start_),
    anchor_start,
    clamp_nuisance(neutral_start)
  )
  starts <- unique(round(starts, 10))
  if (is.null(dim(starts))) starts <- matrix(starts, nrow = 1)
  
  any_resolved <- FALSE
  any_failure  <- FALSE
  
  # Step 1: Direct evaluations at starting points
  for (k in seq_len(nrow(starts))) {
    nuisance_k <- starts[k, , drop = TRUE]
    theta_k    <- c(beta1, nuisance_k)
    raw_score  <- score_Eff_theta(
      theta_k, ux, uy, noise, delta, s_dp, ep, n, aux_ux, aux_uy, aux_N, cache
    )
    
    if (is.finite(raw_score)) {
      score_val <- -raw_score
      any_resolved <- TRUE
      if (score_val >= rank_threshold) return(TRUE)
    } else {
      any_failure <- TRUE
    }
  }
  
  # Step 2: Multi-start L-BFGS-B local optimization
  for (k in seq_len(nrow(starts))) {
    start_k <- starts[k, , drop = TRUE]
    res <- tryCatch(
      optim(
        par = start_k, fn = score_Eff_nuisance, method = "L-BFGS-B",
        lower = nuisance_lower, upper = nuisance_upper,
        control = list(maxit = 100),
        beta1 = beta1, ux = ux, uy = uy, noise = noise, delta = delta,
        s_dp = s_dp, ep = ep, n = n, aux_ux = aux_ux, aux_uy = aux_uy,
        aux_N = aux_N, cache = cache
      ),
      error = function(e) NULL
    )
    
    if (is.null(res) || !is.finite(res$value) || res$value >= 1e11) {
      any_failure <- TRUE
      next
    }
    
    theta_hat <- c(beta1, res$par)
    raw_final <- score_Eff_theta(
      theta_hat, ux, uy, noise, delta, s_dp, ep, n, aux_ux, aux_uy, aux_N, cache
    )
    if (!is.finite(raw_final)) {
      any_failure <- TRUE
      next
    }
    
    final_score <- -raw_final
    any_resolved <- TRUE
    if (final_score >= rank_threshold) return(TRUE)
  }
  
  if (any_failure)  return(NA)
  if (any_resolved) return(FALSE)
  NA
}

# ------------------------------------------------------------------------------
# 10. Single Monte Carlo Cell Execution
# ------------------------------------------------------------------------------

run_power_sweep <- function(
    n_obs,
    beta1_true,
    delta,
    ep
) {
  theta_truth <- c(beta1_true, beta0, mu_x, tau_x, sa_eps)
  
  export_names <- c(
    "clamp_val", "clamp_nuisance", "sdp_vec_lr", "safe_inv_general",
    "efficient_direction_lr", "efficient_depth_lr", "score_Eff_theta",
    "score_Eff_nuisance", "make_data_anchor", "accept_Efficient",
    "h_fd", "R_synthetic", "R_aux", "nuisance_start", "nuisance_lower",
    "nuisance_upper", "rank_threshold", "sa_eps", "reps"
  )
  
  out <- foreach::foreach(
    rep_idx = seq_len(reps),
    .combine = rbind,
    .export = export_names,
    .inorder = FALSE
  ) %dopar% {
    set.seed(rep_idx + 1000L)
    
    # Observed private statistics
    ux1 <- matrix(rnorm(n_obs), nrow = 1)
    uy1 <- matrix(rnorm(n_obs), nrow = 1)
    N1  <- matrix(rnorm(5), nrow = 1)
    s_dp <- sdp_vec_lr(ux1, uy1, N1, delta, theta_truth, ep, n_obs)[1, ]
    
    # Primary repro cloud
    ux_repro <- matrix(rnorm(R_synthetic * n_obs), nrow = R_synthetic)
    uy_repro <- matrix(rnorm(R_synthetic * n_obs), nrow = R_synthetic)
    N_repro  <- matrix(rnorm(R_synthetic * 5), nrow = R_synthetic)
    
    # Independent auxiliary cloud
    aux_ux <- matrix(rnorm(R_aux * n_obs), nrow = R_aux)
    aux_uy <- matrix(rnorm(R_aux * n_obs), nrow = R_aux)
    aux_N  <- matrix(rnorm(R_aux * 5), nrow = R_aux)
    
    eff_cache <- new.env(parent = emptyenv())
    accepted <- accept_Efficient(
      beta1 = 0, nuisance_start_ = nuisance_start, ux = ux_repro, uy = uy_repro,
      noise = N_repro, delta = delta, s_dp = s_dp, ep = ep, n = n_obs,
      aux_ux = aux_ux, aux_uy = aux_uy, aux_N = aux_N, cache = eff_cache
    )
    
    if (isTRUE(accepted)) {
      c(rejected_conservative = 0, rejected_resolved = 0, failed = 0)
    } else if (identical(accepted, FALSE)) {
      c(rejected_conservative = 1, rejected_resolved = 1, failed = 0)
    } else {
      c(rejected_conservative = 0, rejected_resolved = NA_real_, failed = 1)
    }
  }
  
  rejection_conservative <- mean(out[, "rejected_conservative"])
  resolved_idx <- is.finite(out[, "rejected_resolved"])
  rejection_resolved <- if (any(resolved_idx)) mean(out[resolved_idx, "rejected_resolved"]) else NA_real_
  
  failure_rate <- mean(out[, "failed"])
  n_failed     <- sum(out[, "failed"])
  n_resolved   <- reps - n_failed
  
  list(
    rejection_conservative = rejection_conservative,
    rejection_resolved     = rejection_resolved,
    failure_rate           = failure_rate,
    n_failed               = n_failed,
    n_resolved             = n_resolved
  )
}

# ------------------------------------------------------------------------------
# 11. Sensitivity Grid Loop & Results Persistence
# ------------------------------------------------------------------------------

n_list     <- c(100L, 200L, 500L, 1000L)
beta1_list <- c(0, 1)
Delta_list <- c(0.5, 1.0, 2.0, 5.0, 10.0)

results <- data.frame(
  epsilon                        = numeric(0),
  beta1_true                     = numeric(0),
  Delta                          = numeric(0),
  n                              = integer(0),
  rejection_probability          = numeric(0),
  mc_se                          = numeric(0),
  resolved_rejection_probability = numeric(0),
  failure_rate                   = numeric(0),
  n_failed                       = integer(0),
  n_resolved                     = integer(0),
  elapsed_seconds                = numeric(0),
  stringsAsFactors               = FALSE
)

total_cells <- length(beta1_list) * length(Delta_list) * length(n_list)
cell_counter <- 0L

for (beta1_true in beta1_list) {
  for (delta_current in Delta_list) {
    for (n_obs in n_list) {
      cell_counter <- cell_counter + 1L
      
      cat("\n============================================================\n")
      cat(sprintf("CELL %d / %d\n", cell_counter, total_cells))
      cat(sprintf("beta1_true = %.1f | Delta = %.1f | n = %d | epsilon = %.1f\n",
                  beta1_true, delta_current, n_obs, ep_fixed))
      cat("============================================================\n")
      
      t0 <- Sys.time()
      ans <- run_power_sweep(
        n_obs = n_obs, beta1_true = beta1_true, delta = delta_current, ep = ep_fixed
      )
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      
      p_hat <- ans$rejection_conservative
      mc_se <- sqrt(p_hat * (1 - p_hat) / reps)
      
      new_row <- data.frame(
        epsilon                        = ep_fixed,
        beta1_true                     = beta1_true,
        Delta                          = delta_current,
        n                              = n_obs,
        rejection_probability          = p_hat,
        mc_se                          = mc_se,
        resolved_rejection_probability = ans$rejection_resolved,
        failure_rate                   = ans$failure_rate,
        n_failed                       = ans$n_failed,
        n_resolved                     = ans$n_resolved,
        elapsed_seconds                = elapsed
      )
      results <- rbind(results, new_row)
      
      result_label <- if (beta1_true == 0) "TYPE-I ERROR" else "POWER"
      cat(sprintf("%s = %.3f | MC SE = %.4f | failure = %.3f | elapsed = %.1f sec\n",
                  result_label, p_hat, mc_se, ans$failure_rate, elapsed))
      
      write.csv(
        results,
        file.path(PROJECT_DIR, "exp2_delta_sensitivity_long.csv"),
        row.names = FALSE
      )
    }
  }
}

try(parallel::stopCluster(cl), silent = TRUE)
cat("\nSensitivity simulation completed successfully.\n")

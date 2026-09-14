# ==============================================================================
# Script: linearreg_eff.R
# Description: Private linear regression inference under the Efficient Repro 
#              framework. Conducts multi-start profiling over nuisance parameters 
#              to test H0: beta1 = 0. Conservative power accounting treats 
#              numerical non-resolutions as non-rejections.
#
# Parameter layout:
#   theta = (beta1, beta0, mu, tau, sa)
#   Target:    beta1 (evaluated at H0: beta1 = 0)
#   Nuisance:  (beta0, mu, tau, sa)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Directory Initialization
# ------------------------------------------------------------------------------

PROJECT_DIR <- path.expand("~/R_Simuls/Linearegression")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

RESULTS_DIR <- file.path(PROJECT_DIR, "results_eff_corrected")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Simulation Settings & Nuisance Constraints
# ------------------------------------------------------------------------------

mu_x   <- 0.5
tau_x  <- 1.0
beta0  <- -0.5
sa_eps <- 0.5

alpha <- 0.05

R_synthetic <- 200L
R_aux       <- 400L

h_fd <- 1e-3

delta <- 2.0
ep    <- 1.0

reps <- 1000L

# Optimization search box for nuisance vector (beta0, mu, tau, sa)
nuisance_start <- c(
  beta0,
  mu_x,
  tau_x,
  sa_eps
)

nuisance_lower <- c(
  -5.0,
  -5.0,
  0.001,
  0.001
)

nuisance_upper <- c(
  5.0,
  5.0,
  5.0,
  5.0
)

rank_threshold <- floor(
  alpha * (R_synthetic + 1L)
) + 1L

cat(sprintf("R = %d\n", R_synthetic))
cat(sprintf("R_aux = %d\n", R_aux))
cat(sprintf("alpha = %.3f\n", alpha))
cat(sprintf("rank threshold = %d\n\n", rank_threshold))

# ------------------------------------------------------------------------------
# 3. Environment & Parallel Backend Setup
# ------------------------------------------------------------------------------

list.of.packages <- c("foreach", "doSNOW", "parallelly")
new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages) > 0L) {
  install.packages(new.packages, dependencies = TRUE)
}

for (package.i in list.of.packages) {
  suppressPackageStartupMessages(
    library(package.i, character.only = TRUE)
  )
}

# Restrict worker threads to avoid BLAS/MKL CPU oversubscription
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

n.cores <- 124L
cl <- parallel::makePSOCKcluster(n.cores)
doSNOW::registerDoSNOW(cl)

# ------------------------------------------------------------------------------
# 4. Helper Utilities
# ------------------------------------------------------------------------------

clamp_val <- function(x, a, b) {
  pmin(pmax(x, a), b)
}

clamp_nuisance <- function(x) {
  pmin(pmax(x, nuisance_lower), nuisance_upper)
}

# ------------------------------------------------------------------------------
# 5. Differentially Private Summary Statistic Mechanism
# ------------------------------------------------------------------------------

sdp_vec_lr <- function(
    ux,
    uy,
    noise,
    delta,
    theta,
    ep,
    n
) {
  beta1  <- theta[1]
  beta0_ <- theta[2]
  mu     <- theta[3]
  tau    <- theta[4]
  sa     <- theta[5]
  
  x <- ux * tau + mu
  y <- beta0_ + beta1 * x + uy * sa
  
  xc <- clamp_val(x, -delta, delta)
  yc <- clamp_val(y, -delta, delta)
  xyc <- clamp_val(x * y, -delta^2, delta^2)
  
  ep_split <- ep / sqrt(5)
  
  xbar <- rowMeans(xc) + (2 * delta / (n * ep_split)) * noise[, 1]
  ybar <- rowMeans(yc) + (2 * delta / (n * ep_split)) * noise[, 2]
  x2bar <- rowMeans(xc^2) + (delta^2 / (n * ep_split)) * noise[, 3]
  y2bar <- rowMeans(yc^2) + (delta^2 / (n * ep_split)) * noise[, 4]
  xybar <- rowMeans(xyc) + (2 * delta^2 / (n * ep_split)) * noise[, 5]
  
  cbind(
    xbar,
    ybar,
    x2bar,
    y2bar,
    xybar
  )
}

# ------------------------------------------------------------------------------
# 6. Matrix Inversion with Eigenvalue Regularization
# ------------------------------------------------------------------------------

safe_inv_general <- function(
    A,
    ridge = 1e-8
) {
  if (any(!is.finite(A))) {
    stop("Non-finite matrix in safe_inv_general().")
  }
  
  A_sym <- 0.5 * (A + t(A))
  eg <- eigen(A_sym, symmetric = TRUE)
  
  if (
    any(!is.finite(eg$values)) ||
    any(!is.finite(eg$vectors))
  ) {
    stop("Non-finite eigendecomposition.")
  }
  
  eigvals_floored <- pmax(eg$values, ridge)
  
  eg$vectors %*%
    diag(1 / eigvals_floored, nrow = length(eigvals_floored)) %*%
    t(eg$vectors)
}

# ------------------------------------------------------------------------------
# 7. Efficient Direction Computation via Auxiliary Cloud
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
  
  h <- h_fd
  
  beta1  <- theta[1]
  beta0_ <- theta[2]
  mu     <- theta[3]
  tau    <- theta[4]
  sa     <- theta[5]
  
  if (
    !is.finite(tau) ||
    !is.finite(sa) ||
    tau <= 0 ||
    sa <= 0
  ) {
    direction <- rep(NA_real_, 5)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  # Auxiliary cloud centering
  x_center <- aux_ux * tau + mu
  xc_center <- clamp_val(x_center, -delta, delta)
  xc2_center <- xc_center^2
  
  ep_split <- ep / sqrt(5)
  
  y_center <- beta0_ + beta1 * x_center + aux_uy * sa
  yc_center <- clamp_val(y_center, -delta, delta)
  xyc_center <- clamp_val(x_center * y_center, -delta^2, delta^2)
  
  xbar_c <- rowMeans(xc_center) + (2 * delta / (n * ep_split)) * aux_N[, 1]
  ybar_c <- rowMeans(yc_center) + (2 * delta / (n * ep_split)) * aux_N[, 2]
  x2bar_c <- rowMeans(xc2_center) + (delta^2 / (n * ep_split)) * aux_N[, 3]
  y2bar_c <- rowMeans(yc_center^2) + (delta^2 / (n * ep_split)) * aux_N[, 4]
  xybar_c <- rowMeans(xyc_center) + (2 * delta^2 / (n * ep_split)) * aux_N[, 5]
  
  center_cloud <- cbind(xbar_c, ybar_c, x2bar_c, y2bar_c, xybar_c)
  
  if (any(!is.finite(center_cloud))) {
    direction <- rep(NA_real_, 5)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  center_mean <- colMeans(center_cloud)
  
  # Fast parameter evaluation (beta1, beta0, sa do not alter X)
  eval_fast <- function(theta_pert) {
    beta1p <- theta_pert[1]
    beta0p <- theta_pert[2]
    sap    <- theta_pert[5]
    
    y <- beta0p + beta1p * x_center + aux_uy * sap
    yc <- clamp_val(y, -delta, delta)
    xyc <- clamp_val(x_center * y, -delta^2, delta^2)
    
    ybar <- rowMeans(yc) + (2 * delta / (n * ep_split)) * aux_N[, 2]
    y2bar <- rowMeans(yc^2) + (delta^2 / (n * ep_split)) * aux_N[, 4]
    xybar <- rowMeans(xyc) + (2 * delta^2 / (n * ep_split)) * aux_N[, 5]
    
    colMeans(cbind(xbar_c, ybar, x2bar_c, y2bar, xybar))
  }
  
  grad_wrt_fast <- function(idx, lower_bound_check = NULL) {
    tp <- theta
    tm <- theta
    tp[idx] <- theta[idx] + h
    tm[idx] <- theta[idx] - h
    
    if (
      !is.null(lower_bound_check) &&
      (theta[idx] - h) <= lower_bound_check
    ) {
      up <- eval_fast(tp)
      return((up - center_mean) / h)
    }
    
    up <- eval_fast(tp)
    dn <- eval_fast(tm)
    (up - dn) / (2 * h)
  }
  
  grad_wrt_full <- function(idx, lower_bound_check = NULL) {
    tp <- theta
    tm <- theta
    tp[idx] <- theta[idx] + h
    tm[idx] <- theta[idx] - h
    
    if (
      !is.null(lower_bound_check) &&
      (theta[idx] - h) <= lower_bound_check
    ) {
      up <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tp, ep, n))
      return((up - center_mean) / h)
    }
    
    up <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tp, ep, n))
    dn <- colMeans(sdp_vec_lr(aux_ux, aux_uy, aux_N, delta, tm, ep, n))
    (up - dn) / (2 * h)
  }
  
  # Derivative assembly
  d_beta1 <- grad_wrt_fast(1)
  D_nuisance <- cbind(
    grad_wrt_fast(2),        # beta0
    grad_wrt_full(3),        # mu
    grad_wrt_full(4, 0.001), # tau
    grad_wrt_fast(5, 0.001)  # sa
  )
  
  if (any(!is.finite(d_beta1)) || any(!is.finite(D_nuisance))) {
    direction <- rep(NA_real_, 5)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  aux_cov <- cov(center_cloud)
  aux_cov <- 0.5 * (aux_cov + t(aux_cov)) + 1e-8 * diag(5)
  Sigma_inv <- tryCatch(safe_inv_general(aux_cov, 1e-8), error = function(e) NULL)
  
  if (is.null(Sigma_inv)) {
    direction <- rep(NA_real_, 5)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  # GLS partialling-out
  DtSD <- t(D_nuisance) %*% Sigma_inv %*% D_nuisance
  DtSd <- t(D_nuisance) %*% Sigma_inv %*% d_beta1
  DtSD_inv <- tryCatch(safe_inv_general(DtSD, 1e-10), error = function(e) NULL)
  
  if (is.null(DtSD_inv)) {
    direction <- rep(NA_real_, 5)
    assign(key, direction, envir = cache)
    return(direction)
  }
  
  coef <- DtSD_inv %*% DtSd
  target <- d_beta1 - D_nuisance %*% coef
  direction <- as.numeric(Sigma_inv %*% target)
  
  if (any(!is.finite(direction)) || sum(direction^2) < 1e-20) {
    direction <- rep(NA_real_, 5)
  }
  
  assign(key, direction, envir = cache)
  direction
}

# ------------------------------------------------------------------------------
# 8. Efficient Repro Depth Scoring
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
  direction <- efficient_direction_lr(
    theta,
    aux_ux,
    aux_uy,
    aux_N,
    delta,
    ep,
    n,
    cache
  )
  
  if (any(!is.finite(direction))) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  center <- colMeans(synth)
  centered <- sweep(synth, 2, center, "-")
  
  covariance <- crossprod(centered) / nrow(synth)
  covariance <- 0.5 * (covariance + t(covariance)) + 1e-10 * diag(5)
  
  cov_inv <- tryCatch(safe_inv_general(covariance, 1e-10), error = function(e) NULL)
  if (is.null(cov_inv)) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  eff_scale <- as.numeric(t(direction) %*% covariance %*% direction)
  if (!is.finite(eff_scale) || eff_scale <= 1e-14) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  lambda_current <- 1 / log(n)
  projection <- as.numeric(centered %*% direction)
  
  Q_eff <- projection^2 / eff_scale
  Q_full <- rowSums((centered %*% cov_inv) * centered)
  pen <- Q_eff + lambda_current * Q_full
  depth <- 1 / (1 + pen)
  
  if (any(!is.finite(depth))) {
    return(rep(NA_real_, nrow(synth)))
  }
  
  depth
}

# ------------------------------------------------------------------------------
# 9. Test Criterion Scoring Functions
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
  D_synth <- efficient_depth_lr(
    synth_full,
    theta,
    aux_ux,
    aux_uy,
    aux_N,
    delta,
    ep,
    n,
    cache
  )
  
  if (any(!is.finite(D_synth))) {
    return(NA_real_)
  }
  
  obs_idx <- R_synthetic + 1L
  r <- rank(D_synth, ties.method = "max")[obs_idx]
  depth_obs <- D_synth[obs_idx]
  
  if (!is.finite(r) || !is.finite(depth_obs)) {
    return(NA_real_)
  }
  
  # Return negated rank score for minimizers
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
  )
  
  if (!is.finite(score)) {
    return(1e12)
  }
  
  score
}

# Moment-based nuisance parameter initialization
make_data_anchor <- function(beta1, s_dp) {
  xbar  <- s_dp[1]
  ybar  <- s_dp[2]
  x2bar <- s_dp[3]
  y2bar <- s_dp[4]
  
  mu_hat <- xbar
  var_x_hat <- x2bar - xbar^2
  if (!is.finite(var_x_hat)) {
    var_x_hat <- 1
  }
  
  tau_hat <- sqrt(max(var_x_hat, 0.001^2))
  beta0_hat <- ybar - beta1 * mu_hat
  
  var_y_hat <- y2bar - ybar^2
  if (!is.finite(var_y_hat)) {
    var_y_hat <- sa_eps^2
  }
  
  sa2_hat <- var_y_hat - beta1^2 * tau_hat^2
  sa_hat <- sqrt(max(sa2_hat, 0.001^2))
  
  clamp_nuisance(c(beta0_hat, mu_hat, tau_hat, sa_hat))
}

# ------------------------------------------------------------------------------
# 10. Multi-Start Acceptance Search
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
  anchor_start <- make_data_anchor(beta1, s_dp)
  neutral_start <- c(0.0, 0.0, 1.0, 0.5)
  
  starts <- rbind(
    clamp_nuisance(nuisance_start_),
    anchor_start,
    clamp_nuisance(neutral_start)
  )
  
  starts <- unique(round(starts, digits = 10))
  if (is.null(dim(starts))) {
    starts <- matrix(starts, nrow = 1)
  }
  
  any_resolved <- FALSE
  any_failure <- FALSE
  
  # Step 1: Direct evaluation across start anchors
  for (k in seq_len(nrow(starts))) {
    nuisance_k <- starts[k, , drop = TRUE]
    theta_k <- c(beta1, nuisance_k)
    
    direct <- -score_Eff_theta(
      theta_k,
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
    )
    
    if (is.finite(direct)) {
      any_resolved <- TRUE
      if (direct >= rank_threshold) {
        return(TRUE)
      }
    } else {
      any_failure <- TRUE
    }
  }
  
  # Step 2: Multi-start L-BFGS-B local optimization
  for (k in seq_len(nrow(starts))) {
    start_k <- starts[k, , drop = TRUE]
    
    res <- tryCatch(
      optim(
        par = start_k,
        fn = score_Eff_nuisance,
        method = "L-BFGS-B",
        lower = nuisance_lower,
        upper = nuisance_upper,
        control = list(maxit = 100),
        beta1 = beta1,
        ux = ux,
        uy = uy,
        noise = noise,
        delta = delta,
        s_dp = s_dp,
        ep = ep,
        n = n,
        aux_ux = aux_ux,
        aux_uy = aux_uy,
        aux_N = aux_N,
        cache = cache
      ),
      error = function(e) NULL
    )
    
    if (is.null(res) || !is.finite(res$value) || res$value >= 1e11) {
      any_failure <- TRUE
      next
    }
    
    theta_hat <- c(beta1, res$par)
    final_score <- -score_Eff_theta(
      theta_hat,
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
    )
    
    if (!is.finite(final_score)) {
      any_failure <- TRUE
      next
    }
    
    any_resolved <- TRUE
    if (final_score >= rank_threshold) {
      return(TRUE)
    }
  }
  
  if (any_failure) {
    return(NA)
  }
  
  if (any_resolved) {
    return(FALSE)
  }
  
  NA
}

# ------------------------------------------------------------------------------
# 11. Replication Sweep Routine
# ------------------------------------------------------------------------------

run_power_sweep <- function(
    n_obs,
    beta1_true,
    delta,
    ep
) {
  theta_truth <- c(
    beta1_true,
    beta0,
    mu_x,
    tau_x,
    sa_eps
  )
  
  out <- foreach(
    rep_idx = seq_len(reps),
    .combine = rbind,
    .export = setdiff(
      ls(envir = .GlobalEnv),
      c("cl", "n_obs", "beta1_true", "delta", "ep")
    )
  ) %dopar% {
    set.seed(rep_idx + 1000L)
    
    ux1 <- matrix(rnorm(n_obs), nrow = 1, ncol = n_obs)
    uy1 <- matrix(rnorm(n_obs), nrow = 1, ncol = n_obs)
    N1  <- matrix(rnorm(5), nrow = 1, ncol = 5)
    
    s_dp <- sdp_vec_lr(ux1, uy1, N1, delta, theta_truth, ep, n_obs)[1, ]
    
    ux_repro <- matrix(rnorm(R_synthetic * n_obs), nrow = R_synthetic, ncol = n_obs)
    uy_repro <- matrix(rnorm(R_synthetic * n_obs), nrow = R_synthetic, ncol = n_obs)
    N_repro  <- matrix(rnorm(R_synthetic * 5), nrow = R_synthetic, ncol = 5)
    
    aux_ux <- matrix(rnorm(R_aux * n_obs), nrow = R_aux, ncol = n_obs)
    aux_uy <- matrix(rnorm(R_aux * n_obs), nrow = R_aux, ncol = n_obs)
    aux_N  <- matrix(rnorm(R_aux * 5), nrow = R_aux, ncol = 5)
    
    eff_cache <- new.env(parent = emptyenv())
    
    accepted <- accept_Efficient(
      beta1 = 0.0,
      nuisance_start_ = nuisance_start,
      ux = ux_repro,
      uy = uy_repro,
      noise = N_repro,
      delta = delta,
      s_dp = s_dp,
      ep = ep,
      n = n_obs,
      aux_ux = aux_ux,
      aux_uy = aux_uy,
      aux_N = aux_N,
      cache = eff_cache
    )
    
    if (isTRUE(accepted)) {
      c(rejected_conservative = 0, rejected_resolved = 0, failed = 0)
    } else if (identical(accepted, FALSE)) {
      c(rejected_conservative = 1, rejected_resolved = 1, failed = 0)
    } else {
      c(rejected_conservative = 0, rejected_resolved = NA_real_, failed = 1)
    }
  }
  
  power_conservative <- mean(out[, "rejected_conservative"])
  resolved_idx <- is.finite(out[, "rejected_resolved"])
  power_resolved <- if (any(resolved_idx)) {
    mean(out[resolved_idx, "rejected_resolved"])
  } else {
    NA_real_
  }
  
  failure_rate <- mean(out[, "failed"])
  n_failed <- sum(out[, "failed"])
  n_resolved <- reps - n_failed
  
  list(
    power_conservative = power_conservative,
    power_resolved = power_resolved,
    failure_rate = failure_rate,
    n_failed = n_failed,
    n_resolved = n_resolved
  )
}

# ------------------------------------------------------------------------------
# 12. Grid Execution & Results Persistence
# ------------------------------------------------------------------------------

n_list <- c(100, 200, 300, 400, 500, 1000)
beta1_list <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
ep_list <- c(1)

for (ep in ep_list) {
  power_table <- matrix(
    NA_real_,
    nrow = length(n_list),
    ncol = length(beta1_list),
    dimnames = list(n_list, beta1_list)
  )
  
  power_resolved_table <- power_table
  failure_table <- power_table
  failed_count_table <- power_table
  
  for (j in seq_along(beta1_list)) {
    beta1_true <- beta1_list[j]
    
    for (i in seq_along(n_list)) {
      n_obs <- n_list[i]
      t0 <- Sys.time()
      
      ans <- run_power_sweep(
        n_obs = n_obs,
        beta1_true = beta1_true,
        delta = delta,
        ep = ep
      )
      
      power_table[i, j] <- ans$power_conservative
      power_resolved_table[i, j] <- ans$power_resolved
      failure_table[i, j] <- ans$failure_rate
      failed_count_table[i, j] <- ans$n_failed
      
      elapsed <- as.numeric(Sys.time() - t0, units = "secs")
      cat(
        sprintf(
          paste0(
            "[Efficient corrected] ep=%s beta1=%s n=%d | ",
            "power=%.3f | resolved-power=%.3f | failure=%.3f (%d/%d) | %.1fs\n"
          ),
          ep,
          beta1_true,
          n_obs,
          ans$power_conservative,
          ans$power_resolved,
          ans$failure_rate,
          ans$n_failed,
          reps,
          elapsed
        )
      )
      
      write.csv(t(power_table), file.path(RESULTS_DIR, "Efficient_INT.csv"))
      write.csv(t(power_resolved_table), file.path(RESULTS_DIR, "Efficient_INT_resolved_only.csv"))
      write.csv(t(failure_table), file.path(RESULTS_DIR, "Efficient_failure_rate.csv"))
      write.csv(t(failed_count_table), file.path(RESULTS_DIR, "Efficient_failure_count.csv"))
    }
  }
  
  cat("\n============================================================\n")
  cat("CONSERVATIVE POWER / REJECTION TABLE\n")
  cat("============================================================\n")
  print(t(power_table))
  
  cat("\n============================================================\n")
  cat("POWER AMONG RESOLVED REPLICATIONS ONLY\n")
  cat("============================================================\n")
  print(t(power_resolved_table))
  
  cat("\n============================================================\n")
  cat("NUMERICAL FAILURE RATE\n")
  cat("============================================================\n")
  print(t(failure_table))
  cat("\nSaved in:\n", RESULTS_DIR, "\n\n")
}

try(stopCluster(cl), silent = TRUE)
cat("\nDone.\n")

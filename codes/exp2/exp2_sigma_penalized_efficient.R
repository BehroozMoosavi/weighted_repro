# ==============================================================================
# Script: exp2_sigma_penalized_efficient.R
# Description: Finite-sample Penalized Wald Repro confidence interval evaluation 
#              and level-set geometry for sigma in differentially private linear 
#              regression. Includes Halton profiling, compass search refinement, 
#              Monte Carlo sensitivity audits, and empirical inclusion mapping.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment, Threads & Packages
# ------------------------------------------------------------------------------

Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1"
)

required_packages <- c("foreach", "doSNOW")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doSNOW))

PROJECT_DIR <- path.expand("~/R_Simuls/linear_sigma_penalized_efficient")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

# ------------------------------------------------------------------------------
# 2. Parameters & Simulation Configuration
# ------------------------------------------------------------------------------

n_obs <- 100L

beta1_true <-  1.0
beta0_true <- -0.5
mu_true    <-  0.5
tau_true   <-  1.0
sigma_true <-  0.5

theta_true <- c(
  beta1_true,
  beta0_true,
  mu_true,
  tau_true,
  sigma_true
)

THETA_NAMES <- c("beta1", "beta0", "mu", "tau", "sigma")

Delta <- 2.0
mu_gdp <- 1.0
mu_coordinate <- mu_gdp / sqrt(5)

R_synthetic <- 200L
R_aux <- 400L
alpha <- 0.05
nSIM <- 1000L

p_dim <- 1L
q_dim <- 4L
m_dim <- 5L
k_dim <- 5L

TARGET_IDX   <- 5L
NUISANCE_IDX <- 1L:4L

lambda_n <- 1 / log(n_obs)
penalty_share_reference <- (m_dim * lambda_n) / (p_dim + m_dim * lambda_n)

rank_threshold <- floor(alpha * (R_synthetic + 1L)) + 1L
finite_R_reference_level <- (R_synthetic + 2L - rank_threshold) / (R_synthetic + 1L)

SPD_REL_TOL       <- 1e-10
SPD_ABS_TOL       <- 1e-14
DIRECTION_TOL     <- 1e-20
VAR_TOL           <- 1e-14
SEARCH_TIE_WEIGHT <- 0.999

theta_lower <- c(-3.00, -3.00, -3.00, 0.05, 0.05)
theta_upper <- c( 3.00,  3.00,  3.00, 3.00, 1.50)

sigma_lower <- theta_lower[TARGET_IDX]
sigma_upper <- theta_upper[TARGET_IDX]

SIGMA_SCAN_STEP <- 0.05
CI_TOL          <- 0.01

sigma_scan <- seq(sigma_lower, sigma_upper, by = SIGMA_SCAN_STEP)
if (tail(sigma_scan, 1L) < sigma_upper - 1e-12) {
  sigma_scan <- c(sigma_scan, sigma_upper)
}

N_HALTON_STARTS   <- 24L
N_REFINE_STARTS   <- 6L
PS_STEP0          <- c(0.40, 0.40, 0.40, 0.30)
PS_STEP_MIN       <- c(0.025, 0.025, 0.025, 0.025)
PS_MAX_PASSES     <- 50L
CACHE_MAX_ENTRIES <- 2000L

RUN_SEARCH_AUDIT    <- TRUE
N_AUDIT_REPS        <- 100L
AUDIT_HALTON_STARTS <- 64L
AUDIT_REFINE_STARTS <- 12L

RUN_LEVELSET    <- TRUE
N_LEVELSET_REPS <- 100L
LS_IDX_1        <- 1L
LS_IDX_2        <- 2L
LS_HALF_1       <- 0.80
LS_HALF_2       <- 0.80
LS_N_1          <- 41L
LS_N_2          <- 41L

PUB_FONT       <- "Helvetica"
PUB_CEX_AXIS   <- 0.90
PUB_CEX_LAB    <- 1.00
PUB_CEX_PANEL  <- 0.94
PUB_CEX_LEGEND <- 0.76

COL_BLUE   <- "#0072B2"
COL_ORANGE <- "#D55E00"
COL_GREEN  <- "#009E73"
COL_PURPLE <- "#CC79A7"
COL_GRAY   <- "#666666"
COL_LIGHT  <- "#E6E6E6"
COL_BLACK  <- "#000000"

# ------------------------------------------------------------------------------
# 3. Utility Routines & Safe Linear Algebra
# ------------------------------------------------------------------------------

clamp_val <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

in_box <- function(theta) {
  all(is.finite(theta)) &&
    all(theta >= theta_lower) &&
    all(theta <= theta_upper)
}

chol_inverse <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || nrow(A) < 1L || any(!is.finite(A))) {
    return(NULL)
  }
  A <- 0.5 * (A + t(A))
  RR <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(RR)) {
    return(NULL)
  }
  out <- tryCatch(chol2inv(RR), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) {
    return(NULL)
  }
  out
}

strict_spd_inverse <- function(
    A,
    rel_tol = SPD_REL_TOL,
    abs_tol = SPD_ABS_TOL
) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || nrow(A) < 1L || any(!is.finite(A))) {
    return(NULL)
  }
  A <- 0.5 * (A + t(A))
  ee <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ee) || any(!is.finite(ee$values))) {
    return(NULL)
  }
  max_eig <- max(ee$values)
  min_eig <- min(ee$values)
  if (!is.finite(max_eig) || !is.finite(min_eig) || max_eig <= abs_tol) {
    return(NULL)
  }
  threshold <- max(abs_tol, rel_tol * max_eig)
  if (min_eig <= threshold) {
    return(NULL)
  }
  out <- ee$vectors %*% diag(1 / ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  if (any(!is.finite(out))) {
    return(NULL)
  }
  out
}

# ------------------------------------------------------------------------------
# 4. Cache System & Halton Sequence Generation
# ------------------------------------------------------------------------------

CACHE_COUNTER_NAMES <- c(
  "COUNT_entries",
  "COUNT_eval",
  "COUNT_unresolved",
  "COUNT_unique_geomfail",
  "COUNT_identity_wald",
  "COUNT_identity_mh"
)

cache_new <- function() {
  e <- new.env(parent = emptyenv())
  assign("COUNT_entries", 0L, envir = e)
  assign("COUNT_eval", 0L, envir = e)
  assign("COUNT_unresolved", 0L, envir = e)
  assign("COUNT_unique_geomfail", 0L, envir = e)
  assign("COUNT_identity_wald", 0, envir = e)
  assign("COUNT_identity_mh", 0, envir = e)
  e
}

cache_bump <- function(cache, name, by = 1L) {
  assign(name, get(name, envir = cache, inherits = FALSE) + by, envir = cache)
  invisible(NULL)
}

cache_max_update <- function(cache, name, value) {
  if (!is.finite(value)) return(invisible(NULL))
  old <- get(name, envir = cache, inherits = FALSE)
  if (value > old) assign(name, value, envir = cache)
  invisible(NULL)
}

cache_clear_entries <- function(cache) {
  all_names <- ls(cache, all.names = TRUE)
  drop_names <- setdiff(all_names, CACHE_COUNTER_NAMES)
  if (length(drop_names) > 0L) {
    rm(list = drop_names, envir = cache)
  }
  assign("COUNT_entries", 0L, envir = cache)
  invisible(NULL)
}

cache_counters <- function(cache) {
  list(
    n_eval = get("COUNT_eval", envir = cache, inherits = FALSE),
    n_unresolved = get("COUNT_unresolved", envir = cache, inherits = FALSE),
    n_unique_geomfail = get("COUNT_unique_geomfail", envir = cache, inherits = FALSE),
    identity_wald = get("COUNT_identity_wald", envir = cache, inherits = FALSE),
    identity_mh = get("COUNT_identity_mh", envir = cache, inherits = FALSE)
  )
}

theta_cache_key <- function(theta) {
  paste(formatC(theta, digits = 15, format = "fg", flag = "#"), collapse = "|")
}

radical_inverse <- function(index, base) {
  if (index <= 0) return(0)
  result <- 0
  factor <- 1 / base
  i <- index
  while (i > 0) {
    result <- result + factor * (i %% base)
    i <- floor(i / base)
    factor <- factor / base
  }
  result
}

halton_matrix <- function(n_points, dimension, start_index = 1L) {
  primes <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L)
  if (dimension > length(primes)) stop("Not enough Halton bases.")
  H <- matrix(NA_real_, nrow = n_points, ncol = dimension)
  for (i in seq_len(n_points)) {
    idx <- start_index + i - 1L
    for (j in seq_len(dimension)) {
      H[i, j] <- radical_inverse(idx, primes[j])
    }
  }
  H
}

halton_nuisance_starts <- function(n_points, start_index = 1L) {
  H <- halton_matrix(n_points, q_dim, start_index)
  lower <- theta_lower[NUISANCE_IDX]
  upper <- theta_upper[NUISANCE_IDX]
  sweep(sweep(H, 2L, upper - lower, "*"), 2L, lower, "+")
}

# ------------------------------------------------------------------------------
# 5. DGP & Differentially Private Summaries
# ------------------------------------------------------------------------------

draw_seed_cloud <- function(R, n) {
  list(
    ux    = matrix(rnorm(R * n), nrow = R, ncol = n),
    uy    = matrix(rnorm(R * n), nrow = R, ncol = n),
    noise = matrix(rnorm(R * 5L), nrow = R, ncol = 5L)
  )
}

noise_scale <- c(
  2 * Delta / (n_obs * mu_coordinate),
  2 * Delta / (n_obs * mu_coordinate),
  Delta^2 / (n_obs * mu_coordinate),
  Delta^2 / (n_obs * mu_coordinate),
  2 * Delta^2 / (n_obs * mu_coordinate)
)

release_from_seed <- function(
    seed,
    theta,
    want_latent = FALSE
) {
  beta1 <- theta[1L]
  beta0 <- theta[2L]
  mu    <- theta[3L]
  tau   <- theta[4L]
  sigma <- theta[5L]
  
  R_current <- nrow(seed$ux)
  if (!all(is.finite(theta)) || tau <= 0 || sigma <= 0) {
    return(list(S = matrix(NA_real_, nrow = R_current, ncol = 5L)))
  }
  
  x <- mu + tau * seed$ux
  y <- beta0 + beta1 * x + sigma * seed$uy
  
  xc  <- clamp_val(x, -Delta, Delta)
  yc  <- clamp_val(y, -Delta, Delta)
  x2c <- clamp_val(x^2, 0, Delta^2)
  y2c <- clamp_val(y^2, 0, Delta^2)
  xyc <- clamp_val(x * y, -Delta^2, Delta^2)
  
  S <- cbind(
    rowMeans(xc)  + noise_scale[1L] * seed$noise[, 1L],
    rowMeans(yc)  + noise_scale[2L] * seed$noise[, 2L],
    rowMeans(x2c) + noise_scale[3L] * seed$noise[, 3L],
    rowMeans(y2c) + noise_scale[4L] * seed$noise[, 4L],
    rowMeans(xyc) + noise_scale[5L] * seed$noise[, 5L]
  )
  colnames(S) <- c("xbar", "ybar", "x2bar", "y2bar", "xybar")
  
  if (!want_latent) return(list(S = S))
  list(S = S, x = x, y = y, xc = xc, yc = yc)
}

# ------------------------------------------------------------------------------
# 6. Pathwise Analytic Jacobian
# ------------------------------------------------------------------------------

jacobian_analytic <- function(
    seed,
    theta,
    latent
) {
  beta1 <- theta[1L]
  x  <- latent$x
  y  <- latent$y
  xc <- latent$xc
  yc <- latent$yc
  
  if (is.null(x) || is.null(y) || is.null(xc) || is.null(yc)) {
    return(NULL)
  }
  
  I_x  <- abs(x) < Delta
  I_y  <- abs(y) < Delta
  I_xy <- abs(x * y) < Delta^2
  
  dx_list <- list(0, 0, 1, seed$ux, 0)
  dy_list <- list(x, 1, beta1, beta1 * seed$ux, seed$uy)
  
  J <- matrix(0, nrow = 5L, ncol = 5L)
  for (j in seq_len(5L)) {
    dx <- dx_list[[j]]
    dy <- dy_list[[j]]
    dx_zero <- (length(dx) == 1L && isTRUE(dx == 0))
    
    J[1L, j] <- if (dx_zero) 0 else mean(I_x * dx)
    J[2L, j] <- mean(I_y * dy)
    J[3L, j] <- if (dx_zero) 0 else mean(2 * xc * I_x * dx)
    J[4L, j] <- mean(2 * yc * I_y * dy)
    J[5L, j] <- if (dx_zero) mean(I_xy * x * dy) else mean(I_xy * (y * dx + x * dy))
  }
  
  if (any(!is.finite(J))) return(NULL)
  dimnames(J) <- list(c("xbar", "ybar", "x2bar", "y2bar", "xybar"), THETA_NAMES)
  J
}

# ------------------------------------------------------------------------------
# 7. Wald Geometry & Penalized Wald Depth
# ------------------------------------------------------------------------------

wald_geometry <- function(
    theta,
    aux_seed,
    cache
) {
  if (!in_box(theta)) return(NULL)
  key <- theta_cache_key(theta)
  
  if (exists(key, envir = cache, inherits = FALSE)) {
    cached <- get(key, envir = cache, inherits = FALSE)
    if (isTRUE(cached$failed)) return(NULL)
    return(cached)
  }
  
  store_failure <- function() {
    if (get("COUNT_entries", envir = cache, inherits = FALSE) >= CACHE_MAX_ENTRIES) {
      cache_clear_entries(cache)
    }
    assign(key, list(failed = TRUE), envir = cache)
    cache_bump(cache, "COUNT_entries")
    cache_bump(cache, "COUNT_unique_geomfail")
    NULL
  }
  
  aux <- release_from_seed(aux_seed, theta, want_latent = TRUE)
  if (any(!is.finite(aux$S))) return(store_failure())
  
  Sigma_aux <- tryCatch(cov(aux$S), error = function(e) NULL)
  if (is.null(Sigma_aux) || any(!is.finite(Sigma_aux))) return(store_failure())
  
  Sigma_aux <- 0.5 * (Sigma_aux + t(Sigma_aux))
  Sigma_aux_inv <- chol_inverse(Sigma_aux)
  if (is.null(Sigma_aux_inv)) return(store_failure())
  
  J <- jacobian_analytic(aux_seed, theta, aux)
  if (is.null(J)) return(store_failure())
  
  J_target   <- J[, TARGET_IDX, drop = FALSE]
  J_nuisance <- J[, NUISANCE_IDX, drop = FALSE]
  
  nuisance_gram <- t(J_nuisance) %*% Sigma_aux_inv %*% J_nuisance
  nuisance_gram_inv <- strict_spd_inverse(nuisance_gram)
  if (is.null(nuisance_gram_inv)) return(store_failure())
  
  J_bar <- J_target - J_nuisance %*% nuisance_gram_inv %*% t(J_nuisance) %*% Sigma_aux_inv %*% J_target
  direction <- as.numeric(Sigma_aux_inv %*% J_bar)
  
  if (any(!is.finite(direction)) || sum(direction^2) <= DIRECTION_TOL) {
    return(store_failure())
  }
  
  I_eff <- as.numeric(t(J_bar) %*% Sigma_aux_inv %*% J_bar)
  if (!is.finite(I_eff) || I_eff <= VAR_TOL) return(store_failure())
  
  full_information <- t(J) %*% Sigma_aux_inv %*% J
  if (is.null(strict_spd_inverse(full_information))) return(store_failure())
  
  out <- list(
    failed    = FALSE,
    direction = direction,
    J         = J,
    J_bar     = as.numeric(J_bar),
    I_eff     = I_eff,
    Sigma_aux = Sigma_aux
  )
  
  if (get("COUNT_entries", envir = cache, inherits = FALSE) >= CACHE_MAX_ENTRIES) {
    cache_clear_entries(cache)
  }
  assign(key, out, envir = cache)
  cache_bump(cache, "COUNT_entries")
  out
}

penalized_wald_depth <- function(
    cloud,
    theta,
    aux_seed,
    cache
) {
  cloud <- as.matrix(cloud)
  N_cloud <- nrow(cloud)
  if (N_cloud < 2L || ncol(cloud) != k_dim || any(!is.finite(cloud))) return(NULL)
  
  geom <- wald_geometry(theta, aux_seed, cache)
  if (is.null(geom)) return(NULL)
  
  centre <- colMeans(cloud)
  centered <- sweep(cloud, 2L, centre, "-")
  
  Sigma_cloud <- crossprod(centered) / N_cloud
  Sigma_cloud <- 0.5 * (Sigma_cloud + t(Sigma_cloud))
  Sigma_cloud_inv <- chol_inverse(Sigma_cloud)
  if (is.null(Sigma_cloud_inv)) return(NULL)
  
  # Wald component
  d <- geom$direction
  wald_variance <- as.numeric(t(d) %*% Sigma_cloud %*% d)
  if (!is.finite(wald_variance) || wald_variance <= VAR_TOL) return(NULL)
  
  projection <- as.numeric(centered %*% d)
  Q_wald <- projection^2 / wald_variance
  if (any(!is.finite(Q_wald))) return(NULL)
  Q_wald <- pmax(Q_wald, 0)
  
  # Projected Mahalanobis component
  J <- geom$J
  full_gram <- t(J) %*% Sigma_cloud_inv %*% J
  full_gram_inv <- strict_spd_inverse(full_gram)
  if (is.null(full_gram_inv)) return(NULL)
  
  M_proj <- Sigma_cloud_inv %*% J %*% full_gram_inv %*% t(J) %*% Sigma_cloud_inv
  M_proj <- 0.5 * (M_proj + t(M_proj))
  if (any(!is.finite(M_proj))) return(NULL)
  
  Q_mh <- rowSums((centered %*% M_proj) * centered)
  if (any(!is.finite(Q_mh))) return(NULL)
  Q_mh <- pmax(Q_mh, 0)
  
  # Total penalized Wald statistic
  Q_pw <- Q_wald + lambda_n * Q_mh
  if (any(!is.finite(Q_pw))) return(NULL)
  
  depth <- 1 / (1 + Q_pw)
  if (any(!is.finite(depth)) || any(depth <= 0) || any(depth > 1)) return(NULL)
  
  cache_max_update(cache, "COUNT_identity_wald", abs(mean(Q_wald) - p_dim))
  cache_max_update(cache, "COUNT_identity_mh", abs(mean(Q_mh) - m_dim))
  
  list(depth = depth, Q_wald = Q_wald, Q_mh = Q_mh, Q_pw = Q_pw)
}

# ------------------------------------------------------------------------------
# 8. Candidate Compatibility & Nuisance Profiling
# ------------------------------------------------------------------------------

unresolved_answer <- list(
  resolved   = FALSE,
  accepted   = FALSE,
  rank_exact = NA_integer_,
  score      = -Inf,
  depth_obs  = NA_real_,
  Q_pw       = NA_real_,
  Q_wald     = NA_real_,
  Q_mh       = NA_real_
)

candidate_compatibility <- function(
    theta,
    syn_seed,
    s_obs,
    aux_seed,
    cache
) {
  cache_bump(cache, "COUNT_eval")
  if (!in_box(theta)) {
    cache_bump(cache, "COUNT_unresolved")
    return(unresolved_answer)
  }
  
  sim <- release_from_seed(syn_seed, theta)$S
  if (any(!is.finite(sim))) {
    cache_bump(cache, "COUNT_unresolved")
    return(unresolved_answer)
  }
  
  cloud <- rbind(sim, s_obs)
  dep <- penalized_wald_depth(cloud, theta, aux_seed, cache)
  if (is.null(dep)) {
    cache_bump(cache, "COUNT_unresolved")
    return(unresolved_answer)
  }
  
  D <- dep$depth
  if (length(D) != (R_synthetic + 1L)) {
    cache_bump(cache, "COUNT_unresolved")
    return(unresolved_answer)
  }
  
  obs_idx <- R_synthetic + 1L
  depth_obs <- D[obs_idx]
  rank_exact <- as.integer(sum(D <= depth_obs))
  
  list(
    resolved   = TRUE,
    accepted   = (rank_exact >= rank_threshold),
    rank_exact = rank_exact,
    score      = rank_exact + SEARCH_TIE_WEIGHT * depth_obs,
    depth_obs  = depth_obs,
    Q_pw       = dep$Q_pw[obs_idx],
    Q_wald     = dep$Q_wald[obs_idx],
    Q_mh       = dep$Q_mh[obs_idx]
  )
}

data_nuisance_anchor <- function(s_obs) {
  xbar  <- s_obs[1L]
  ybar  <- s_obs[2L]
  x2bar <- s_obs[3L]
  xybar <- s_obs[5L]
  
  var_x <- x2bar - xbar^2
  if (!is.finite(var_x) || var_x <= 1e-4) var_x <- 1
  
  beta1_hat <- (xybar - xbar * ybar) / var_x
  if (!is.finite(beta1_hat)) beta1_hat <- 0
  
  beta0_hat <- ybar - beta1_hat * xbar
  if (!is.finite(beta0_hat)) beta0_hat <- 0
  
  tau_hat <- sqrt(max(var_x, theta_lower[4L]^2))
  out <- c(beta1_hat, beta0_hat, xbar, tau_hat)
  clamp_val(out, theta_lower[NUISANCE_IDX], theta_upper[NUISANCE_IDX])
}

nuisance_start_matrix <- function(s_obs, n_halton, extra_start = NULL) {
  anchor   <- data_nuisance_anchor(s_obs)
  neutral  <- clamp_val(c(0, 0, 0, 1), theta_lower[NUISANCE_IDX], theta_upper[NUISANCE_IDX])
  midpoint <- (theta_lower[NUISANCE_IDX] + theta_upper[NUISANCE_IDX]) / 2
  
  starts <- rbind(anchor, neutral, midpoint, halton_nuisance_starts(n_halton))
  if (!is.null(extra_start) && length(extra_start) == q_dim && all(is.finite(extra_start))) {
    starts <- rbind(clamp_val(extra_start, theta_lower[NUISANCE_IDX], theta_upper[NUISANCE_IDX]), starts)
  }
  
  starts <- unique(round(starts, 12))
  if (is.null(dim(starts))) starts <- matrix(starts, nrow = 1L)
  starts
}

compass_from_start <- function(
    eta_start,
    sigma,
    syn_seed,
    s_obs,
    aux_seed,
    cache
) {
  theta_of <- function(eta) {
    theta <- numeric(5L)
    theta[NUISANCE_IDX] <- eta
    theta[TARGET_IDX]   <- sigma
    theta
  }
  
  ans <- candidate_compatibility(theta_of(eta_start), syn_seed, s_obs, aux_seed, cache)
  if (isTRUE(ans$accepted)) {
    return(list(accepted = TRUE, resolved = TRUE, eta = eta_start, score = ans$score))
  }
  if (!isTRUE(ans$resolved)) {
    return(list(accepted = FALSE, resolved = FALSE, eta = eta_start, score = -Inf))
  }
  
  eta <- as.numeric(eta_start)
  best_score <- ans$score
  step <- PS_STEP0
  pass <- 0L
  
  while (any(step >= PS_STEP_MIN) && pass < PS_MAX_PASSES) {
    pass <- pass + 1L
    improved <- FALSE
    
    for (j in seq_len(q_dim)) {
      if (step[j] < PS_STEP_MIN[j]) next
      
      for (direction_sign in c(-1, 1)) {
        candidate <- eta
        candidate[j] <- candidate[j] + direction_sign * step[j]
        lo <- theta_lower[NUISANCE_IDX[j]]
        hi <- theta_upper[NUISANCE_IDX[j]]
        if (candidate[j] < lo || candidate[j] > hi) next
        
        candidate_answer <- candidate_compatibility(theta_of(candidate), syn_seed, s_obs, aux_seed, cache)
        if (isTRUE(candidate_answer$accepted)) {
          return(list(accepted = TRUE, resolved = TRUE, eta = candidate, score = candidate_answer$score))
        }
        if (!isTRUE(candidate_answer$resolved)) next
        
        if (candidate_answer$score > best_score) {
          best_score <- candidate_answer$score
          eta <- candidate
          improved <- TRUE
        }
      }
    }
    if (!improved) step <- step / 2
  }
  
  list(accepted = FALSE, resolved = TRUE, eta = eta, score = best_score)
}

profile_at_sigma <- function(
    sigma,
    syn_seed,
    s_obs,
    aux_seed,
    cache,
    n_halton = N_HALTON_STARTS,
    n_refine = N_REFINE_STARTS,
    warm_start = NULL
) {
  starts <- nuisance_start_matrix(s_obs, n_halton, warm_start)
  n_starts <- nrow(starts)
  
  theta_of <- function(eta) {
    theta <- numeric(5L)
    theta[NUISANCE_IDX] <- eta
    theta[TARGET_IDX]   <- sigma
    theta
  }
  
  scores <- rep(-Inf, n_starts)
  resolved_flag <- logical(n_starts)
  
  for (k in seq_len(n_starts)) {
    eta <- as.numeric(starts[k, ])
    ans <- candidate_compatibility(theta_of(eta), syn_seed, s_obs, aux_seed, cache)
    
    if (isTRUE(ans$resolved)) {
      resolved_flag[k] <- TRUE
      scores[k] <- ans$score
    }
    if (isTRUE(ans$accepted)) {
      return(list(resolved = TRUE, accepted = TRUE, eta = eta, best_score = ans$score))
    }
  }
  
  valid_idx <- which(resolved_flag)
  if (length(valid_idx) == 0L) {
    return(list(resolved = FALSE, accepted = FALSE, eta = rep(NA_real_, q_dim), best_score = -Inf))
  }
  
  order_idx  <- valid_idx[order(scores[valid_idx], decreasing = TRUE)]
  refine_idx <- head(order_idx, min(n_refine, length(order_idx)))
  best_eta   <- as.numeric(starts[refine_idx[1L], ])
  best_score <- scores[refine_idx[1L]]
  
  for (k in refine_idx) {
    out <- compass_from_start(as.numeric(starts[k, ]), sigma, syn_seed, s_obs, aux_seed, cache)
    if (isTRUE(out$accepted)) {
      return(list(resolved = TRUE, accepted = TRUE, eta = out$eta, best_score = out$score))
    }
    if (isTRUE(out$resolved) && out$score > best_score) {
      best_score <- out$score
      best_eta   <- out$eta
    }
  }
  
  list(resolved = TRUE, accepted = FALSE, eta = best_eta, best_score = best_score)
}

# ------------------------------------------------------------------------------
# 9. Confidence Interval Inversion & Boundary Search
# ------------------------------------------------------------------------------

count_true_components <- function(x) {
  x <- as.logical(x)
  if (length(x) == 0L || !any(x)) return(0L)
  component_starts <- which(x & c(TRUE, !head(x, -1L)))
  length(component_starts)
}

refine_boundary <- function(
    sigma_rejected,
    sigma_accepted,
    syn_seed,
    s_obs,
    aux_seed,
    cache,
    warm_start = NULL
) {
  lower <- min(sigma_rejected, sigma_accepted)
  upper <- max(sigma_rejected, sigma_accepted)
  acceptance_is_on_right <- (sigma_rejected < sigma_accepted)
  
  while (upper - lower > CI_TOL) {
    midpoint <- 0.5 * (lower + upper)
    fit <- profile_at_sigma(midpoint, syn_seed, s_obs, aux_seed, cache, warm_start = warm_start)
    if (!isTRUE(fit$resolved)) return(NA_real_)
    
    if (isTRUE(fit$accepted)) {
      warm_start <- fit$eta
      if (acceptance_is_on_right) upper <- midpoint else lower <- midpoint
    } else {
      if (acceptance_is_on_right) lower <- midpoint else upper <- midpoint
    }
  }
  
  if (acceptance_is_on_right) lower else upper
}

construct_sigma_CI <- function(
    syn_seed,
    s_obs,
    aux_seed,
    cache
) {
  n_grid <- length(sigma_scan)
  resolved <- logical(n_grid)
  accepted <- logical(n_grid)
  eta_at_grid <- matrix(NA_real_, nrow = n_grid, ncol = q_dim)
  warm <- NULL
  
  for (j in seq_len(n_grid)) {
    fit <- profile_at_sigma(sigma_scan[j], syn_seed, s_obs, aux_seed, cache, warm_start = warm)
    resolved[j] <- isTRUE(fit$resolved)
    accepted[j] <- isTRUE(fit$accepted)
    if (length(fit$eta) == q_dim && all(is.finite(fit$eta))) {
      eta_at_grid[j, ] <- fit$eta
    }
    if (isTRUE(fit$accepted)) {
      warm <- fit$eta
    }
  }
  
  accepted_resolved <- resolved & accepted
  acc_idx <- which(accepted_resolved)
  n_accepted_grid <- length(acc_idx)
  n_components <- count_true_components(accepted_resolved)
  
  failed_result <- list(
    CI              = c(NA_real_, NA_real_),
    failed          = TRUE,
    n_accepted_grid = n_accepted_grid,
    n_components    = n_components,
    disconnected    = as.numeric(n_components > 1L)
  )
  
  if (n_accepted_grid == 0L || n_components > 1L) {
    return(failed_result)
  }
  
  left_idx  <- min(acc_idx)
  right_idx <- max(acc_idx)
  
  # Lower bound resolution
  if (left_idx == 1L) {
    lower <- sigma_lower
  } else {
    nb <- left_idx - 1L
    if (!resolved[nb]) return(failed_result)
    if (accepted[nb]) {
      lower <- sigma_scan[nb]
    } else {
      lower <- refine_boundary(
        sigma_scan[nb],
        sigma_scan[left_idx],
        syn_seed,
        s_obs,
        aux_seed,
        cache,
        warm_start = eta_at_grid[left_idx, ]
      )
    }
  }
  
  # Upper bound resolution
  if (right_idx == n_grid) {
    upper <- sigma_upper
  } else {
    nb <- right_idx + 1L
    if (!resolved[nb]) return(failed_result)
    if (accepted[nb]) {
      upper <- sigma_scan[nb]
    } else {
      upper <- refine_boundary(
        sigma_scan[nb],
        sigma_scan[right_idx],
        syn_seed,
        s_obs,
        aux_seed,
        cache,
        warm_start = eta_at_grid[right_idx, ]
      )
    }
  }
  
  if (!is.finite(lower) || !is.finite(upper) || lower > upper) {
    return(failed_result)
  }
  
  lower <- max(lower, sigma_lower)
  upper <- min(upper, sigma_upper)
  
  list(
    CI              = c(lower, upper),
    failed          = FALSE,
    n_accepted_grid = n_accepted_grid,
    n_components    = n_components,
    disconnected    = 0
  )
}

# ------------------------------------------------------------------------------
# 10. Monte Carlo Replication Routine
# ------------------------------------------------------------------------------

one_replication <- function(rep_idx) {
  set.seed(rep_idx + 1000L)
  
  bad <- c(
    rep_idx           = rep_idx,
    lower             = NA_real_,
    upper             = NA_real_,
    width             = NA_real_,
    covered           = 0,
    failed            = 1,
    disconnected      = NA_real_,
    n_components      = NA_real_,
    n_accepted_grid   = NA_real_,
    n_eval            = NA_real_,
    n_unresolved      = NA_real_,
    n_unique_geomfail = NA_real_,
    identity_wald     = NA_real_,
    identity_mh       = NA_real_
  )
  
  obs_seed <- draw_seed_cloud(1L, n_obs)
  s_obs <- as.numeric(release_from_seed(obs_seed, theta_true)$S)
  if (any(!is.finite(s_obs))) return(bad)
  
  syn_seed <- draw_seed_cloud(R_synthetic, n_obs)
  aux_seed <- draw_seed_cloud(R_aux, n_obs)
  cache    <- cache_new()
  
  result <- tryCatch(
    construct_sigma_CI(syn_seed, s_obs, aux_seed, cache),
    error = function(e) {
      list(
        CI              = c(NA_real_, NA_real_),
        failed          = TRUE,
        n_accepted_grid = NA_real_,
        n_components    = NA_real_,
        disconnected    = NA_real_
      )
    }
  )
  
  cnt <- cache_counters(cache)
  CI  <- result$CI
  
  if (isTRUE(result$failed) || length(CI) != 2L || any(!is.finite(CI))) {
    bad["n_accepted_grid"]   <- result$n_accepted_grid
    bad["n_components"]      <- result$n_components
    bad["disconnected"]      <- result$disconnected
    bad["n_eval"]            <- cnt$n_eval
    bad["n_unresolved"]      <- cnt$n_unresolved
    bad["n_unique_geomfail"] <- cnt$n_unique_geomfail
    bad["identity_wald"]     <- cnt$identity_wald
    bad["identity_mh"]       <- cnt$identity_mh
    return(bad)
  }
  
  c(
    rep_idx           = rep_idx,
    lower             = CI[1L],
    upper             = CI[2L],
    width             = CI[2L] - CI[1L],
    covered           = as.numeric(CI[1L] <= sigma_true && sigma_true <= CI[2L]),
    failed            = 0,
    disconnected      = result$disconnected,
    n_components      = result$n_components,
    n_accepted_grid   = result$n_accepted_grid,
    n_eval            = cnt$n_eval,
    n_unresolved      = cnt$n_unresolved,
    n_unique_geomfail = cnt$n_unique_geomfail,
    identity_wald     = cnt$identity_wald,
    identity_mh       = cnt$identity_mh
  )
}

# ------------------------------------------------------------------------------
# 11. Main Monte Carlo Execution
# ------------------------------------------------------------------------------

parse_first_int <- function(x) {
  if (!nzchar(x)) return(NA_integer_)
  out <- suppressWarnings(as.integer(sub("[^0-9].*$", "", x)))
  if (is.na(out) || out < 1L) return(NA_integer_)
  out
}

slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_PER_TASK", unset = ""))
if (is.na(slurm_cpus)) slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_ON_NODE", unset = ""))
if (is.na(slurm_cpus)) slurm_cpus <- parallel::detectCores(logical = FALSE)
if (is.na(slurm_cpus) || slurm_cpus < 1L) slurm_cpus <- 1L

n.cores <- min(124L, slurm_cpus, nSIM)

EXPORT_NAMES <- c(
  "n_obs", "beta1_true", "beta0_true", "mu_true", "tau_true", "sigma_true",
  "theta_true", "THETA_NAMES", "Delta", "mu_gdp", "mu_coordinate", "noise_scale",
  "R_synthetic", "R_aux", "alpha", "nSIM", "p_dim", "q_dim", "m_dim", "k_dim",
  "lambda_n", "penalty_share_reference", "rank_threshold", "finite_R_reference_level",
  "SPD_REL_TOL", "SPD_ABS_TOL", "DIRECTION_TOL", "VAR_TOL", "SEARCH_TIE_WEIGHT",
  "theta_lower", "theta_upper", "TARGET_IDX", "NUISANCE_IDX", "sigma_lower",
  "sigma_upper", "SIGMA_SCAN_STEP", "CI_TOL", "sigma_scan", "N_HALTON_STARTS",
  "N_REFINE_STARTS", "PS_STEP0", "PS_STEP_MIN", "PS_MAX_PASSES", "CACHE_MAX_ENTRIES",
  "CACHE_COUNTER_NAMES", "AUDIT_HALTON_STARTS", "AUDIT_REFINE_STARTS", "clamp_val",
  "in_box", "chol_inverse", "strict_spd_inverse", "cache_new", "cache_bump",
  "cache_max_update", "cache_clear_entries", "cache_counters", "theta_cache_key",
  "radical_inverse", "halton_matrix", "halton_nuisance_starts", "draw_seed_cloud",
  "release_from_seed", "jacobian_analytic", "wald_geometry", "penalized_wald_depth",
  "unresolved_answer", "candidate_compatibility", "data_nuisance_anchor",
  "nuisance_start_matrix", "compass_from_start", "profile_at_sigma",
  "count_true_components", "refine_boundary", "construct_sigma_CI", "one_replication"
)

cat("\n============================================================\n")
cat("EXPERIMENT 2: PENALIZED WALD -- SIGMA TARGET\n")
cat("============================================================\n")
cat(sprintf("n                            %d\n", n_obs))
cat(sprintf("true sigma                   %.2f\n", sigma_true))
cat(sprintf("GDP privacy                  %.2f\n", mu_gdp))
cat(sprintf("Delta                        %.2f\n", Delta))
cat(sprintf("R / R_aux                    %d / %d\n", R_synthetic, R_aux))
cat(sprintf("lambda_n                     %.6f\n", lambda_n))
cat(sprintf("exact rank cutoff            %d\n", rank_threshold))
cat(sprintf("finite-R reference           %.6f\n", finite_R_reference_level))
cat(sprintf("Monte Carlo reps             %d\n", nSIM))
cat(sprintf("workers                      %d\n", n.cores))
cat("============================================================\n\n")

start_time <- Sys.time()
cl <- parallel::makePSOCKcluster(n.cores)
doSNOW::registerDoSNOW(cl)
parallel::clusterExport(cl, EXPORT_NAMES, envir = environment())

pb <- txtProgressBar(max = nSIM, style = 3)
progress <- function(k) setTxtProgressBar(pb, k)

SIM_RESULTS <- foreach(
  rep_idx = seq_len(nSIM),
  .combine = "rbind",
  .options.snow = list(progress = progress)
) %dopar% {
  one_replication(rep_idx)
}

try(close(pb), silent = TRUE)
parallel::stopCluster(cl)

SIM_RESULTS <- as.data.frame(SIM_RESULTS)

# ------------------------------------------------------------------------------
# 12. Main Results Aggregation & Table Output
# ------------------------------------------------------------------------------

coverage    <- mean(SIM_RESULTS$covered)
coverage_se <- sqrt(coverage * (1 - coverage) / nSIM)
n_covered   <- sum(SIM_RESULTS$covered)

cp_lower <- if (n_covered == 0L) 0 else qbeta(0.025, n_covered, nSIM - n_covered + 1L)
cp_upper <- if (n_covered == nSIM) 1 else qbeta(0.975, n_covered + 1L, nSIM - n_covered)

failure_rate      <- mean(SIM_RESULTS$failed)
disconnected_rate <- mean(SIM_RESULTS$disconnected == 1, na.rm = TRUE)

ok_width      <- is.finite(SIM_RESULTS$width)
n_valid_width <- sum(ok_width)

mean_width    <- if (n_valid_width > 0L) mean(SIM_RESULTS$width[ok_width]) else NA_real_
median_width  <- if (n_valid_width > 0L) median(SIM_RESULTS$width[ok_width]) else NA_real_
width_q25     <- if (n_valid_width > 0L) unname(quantile(SIM_RESULTS$width[ok_width], 0.25)) else NA_real_
width_q75     <- if (n_valid_width > 0L) unname(quantile(SIM_RESULTS$width[ok_width], 0.75)) else NA_real_
mean_width_se <- if (n_valid_width >= 2L) sd(SIM_RESULTS$width[ok_width]) / sqrt(n_valid_width) else NA_real_

median_evals     <- median(SIM_RESULTS$n_eval, na.rm = TRUE)
total_unresolved <- sum(SIM_RESULTS$n_unresolved, na.rm = TRUE)
total_evals      <- sum(SIM_RESULTS$n_eval, na.rm = TRUE)
unresolved_fraction <- if (total_evals > 0) total_unresolved / total_evals else NA_real_

total_unique_geomfail <- sum(SIM_RESULTS$n_unique_geomfail, na.rm = TRUE)
worst_identity_wald   <- suppressWarnings(max(SIM_RESULTS$identity_wald, na.rm = TRUE))
worst_identity_mh     <- suppressWarnings(max(SIM_RESULTS$identity_mh, na.rm = TRUE))

if (!is.finite(worst_identity_wald)) worst_identity_wald <- NA_real_
if (!is.finite(worst_identity_mh)) worst_identity_mh <- NA_real_

SUMMARY <- data.frame(
  Method                   = "Penalized Wald",
  Target                   = "sigma",
  n                        = n_obs,
  GDP                      = mu_gdp,
  Delta                    = Delta,
  R                        = R_synthetic,
  R_aux                    = R_aux,
  Lambda                   = lambda_n,
  Coverage                 = coverage,
  MC_SE_Coverage           = coverage_se,
  CP_Lower                 = cp_lower,
  CP_Upper                 = cp_upper,
  Mean_Width               = mean_width,
  MC_SE_Width              = mean_width_se,
  Median_Width             = median_width,
  Width_Q25                = width_q25,
  Width_Q75                = width_q75,
  Failure_Rate             = failure_rate,
  Disconnected_Rate        = disconnected_rate,
  Unresolved_Eval_Fraction = unresolved_fraction,
  N                        = nSIM,
  stringsAsFactors         = FALSE
)

write.csv(SIM_RESULTS, file.path(PROJECT_DIR, "exp2_sigma_raw.csv"), row.names = FALSE)
write.csv(SUMMARY, file.path(PROJECT_DIR, "Table_Exp2_sigma_summary.csv"), row.names = FALSE)

latex_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Finite-sample performance of the Penalized Wald Repro confidence interval for $\\sigma$.}",
  "\\label{tab:exp2-sigma-penalized-wald}",
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "Method & Coverage & 95\\% MC CI & Mean width & Median width & Failure & Disconnected & $N$ \\\\",
  "\\midrule",
  sprintf(
    paste0("Penalized Wald & %.3f & [%.3f, %.3f] & %.3f & %.3f & %.3f & %.3f & %d \\\\"),
    coverage, cp_lower, cp_upper, mean_width, median_width, failure_rate, disconnected_rate, nSIM
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex_table, con = file.path(PROJECT_DIR, "Table_Exp2_sigma_summary.tex"))

# ------------------------------------------------------------------------------
# 13. Performance Diagnostic Plot
# ------------------------------------------------------------------------------

draw_performance_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(
    mfrow = c(1, 2),
    family = "sans",
    mar = c(4.2, 4.4, 2.0, 0.8),
    oma = c(0, 0, 0, 0),
    mgp = c(2.45, 0.65, 0),
    tcl = -0.22,
    las = 1,
    cex.axis = PUB_CEX_AXIS,
    cex.lab = PUB_CEX_LAB
  )
  
  # Panel (a): Empirical coverage
  y_lower <- max(0, min(coverage, cp_lower, 1 - alpha, finite_R_reference_level) - 0.025)
  plot(
    1, coverage, type = "n", xlim = c(0.72, 1.28), ylim = c(y_lower, 1.002),
    xaxt = "n", xlab = "", ylab = "Empirical coverage", bty = "l"
  )
  abline(h = 1 - alpha, lty = 3, lwd = 1.25, col = COL_GRAY)
  abline(h = finite_R_reference_level, lty = 2, lwd = 1.45, col = COL_BLACK)
  arrows(x0 = 1, y0 = cp_lower, x1 = 1, y1 = cp_upper, code = 3, angle = 90, length = 0.05, lwd = 1.6, col = COL_BLUE)
  points(1, coverage, pch = 19, cex = 1.1, col = COL_BLUE)
  axis(1, at = 1, labels = "Penalized Wald", tick = FALSE)
  
  legend(
    "bottomleft",
    legend = c(
      sprintf("Coverage = %.3f", coverage),
      sprintf("95%% MC interval = [%.3f, %.3f]", cp_lower, cp_upper),
      sprintf("Finite-R reference = %.4f", finite_R_reference_level)
    ),
    pch = c(19, NA, NA),
    lty = c(NA, 1, 2),
    col = c(COL_BLUE, COL_BLUE, COL_BLACK),
    lwd = c(NA, 1.5, 1.45),
    bty = "n",
    cex = PUB_CEX_LEGEND,
    x.intersp = 0.7,
    y.intersp = 1.0
  )
  mtext("(a) Coverage", side = 3, adj = 0, line = 0.3, font = 2, cex = PUB_CEX_PANEL)
  
  # Panel (b): Confidence interval width
  widths <- SIM_RESULTS$width[is.finite(SIM_RESULTS$width)]
  if (length(widths) >= 2L && diff(range(widths)) > 0) {
    width_breaks <- if (IQR(widths) > 0) "FD" else 12L
    hist(
      widths, breaks = width_breaks, probability = TRUE, col = "gray92", border = "white",
      xlab = expression("Confidence interval width for " * sigma), ylab = "Density", main = "", bty = "l"
    )
    if (length(unique(widths)) >= 3L) {
      lines(density(widths), lwd = 1.7, col = COL_BLUE)
    }
    abline(v = mean_width, col = COL_ORANGE, lwd = 1.7)
    abline(v = median_width, col = COL_GREEN, lwd = 1.7, lty = 2)
    legend(
      "topright",
      legend = c(sprintf("Mean = %.3f", mean_width), sprintf("Median = %.3f", median_width)),
      col = c(COL_ORANGE, COL_GREEN),
      lty = c(1, 2),
      lwd = 1.7,
      bty = "n",
      cex = PUB_CEX_LEGEND
    )
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient variation in finite interval widths", cex = 0.88)
  }
  mtext("(b) Confidence interval width", side = 3, adj = 0, line = 0.3, font = 2, cex = PUB_CEX_PANEL)
}

pdf(file.path(PROJECT_DIR, "Figure_Exp2_sigma_performance.pdf"), width = 8.6, height = 3.8, family = PUB_FONT, useDingbats = FALSE)
draw_performance_figure()
invisible(dev.off())

png(file.path(PROJECT_DIR, "Figure_Exp2_sigma_performance.png"), width = 2580, height = 1140, res = 300, bg = "white")
draw_performance_figure()
invisible(dev.off())

# ------------------------------------------------------------------------------
# 14. Empirical Repro Level-Set Study
# ------------------------------------------------------------------------------

grid_1 <- seq(max(theta_lower[LS_IDX_1], theta_true[LS_IDX_1] - LS_HALF_1),
              min(theta_upper[LS_IDX_1], theta_true[LS_IDX_1] + LS_HALF_1), length.out = LS_N_1)
grid_2 <- seq(max(theta_lower[LS_IDX_2], theta_true[LS_IDX_2] - LS_HALF_2),
              min(theta_upper[LS_IDX_2], theta_true[LS_IDX_2] + LS_HALF_2), length.out = LS_N_2)

one_levelset_replication <- function(rep_idx) {
  set.seed(100000L + rep_idx)
  ng <- LS_N_1 * LS_N_2
  
  obs_seed <- draw_seed_cloud(1L, n_obs)
  s_obs <- as.numeric(release_from_seed(obs_seed, theta_true)$S)
  if (any(!is.finite(s_obs))) {
    return(list(
      Q_wald      = rep(NA_real_, ng),
      lambda_Q_mh = rep(NA_real_, ng),
      Q_pw        = rep(NA_real_, ng),
      accepted    = rep(0, ng),
      resolved    = rep(0, ng)
    ))
  }
  
  syn_seed <- draw_seed_cloud(R_synthetic, n_obs)
  aux_seed <- draw_seed_cloud(R_aux, n_obs)
  cache    <- cache_new()
  
  Q_wald_vec      <- rep(NA_real_, ng)
  lambda_Q_mh_vec <- rep(NA_real_, ng)
  Q_pw_vec        <- rep(NA_real_, ng)
  accepted_vec    <- rep(0, ng)
  resolved_vec    <- rep(0, ng)
  
  index <- 0L
  for (j in seq_along(grid_2)) {
    for (i in seq_along(grid_1)) {
      index <- index + 1L
      theta_ij <- theta_true
      theta_ij[LS_IDX_1] <- grid_1[i]
      theta_ij[LS_IDX_2] <- grid_2[j]
      
      ans <- candidate_compatibility(theta_ij, syn_seed, s_obs, aux_seed, cache)
      if (isTRUE(ans$resolved)) {
        resolved_vec[index]    <- 1
        accepted_vec[index]    <- as.numeric(ans$accepted)
        Q_wald_vec[index]      <- ans$Q_wald
        lambda_Q_mh_vec[index] <- lambda_n * ans$Q_mh
        Q_pw_vec[index]        <- ans$Q_pw
      }
    }
  }
  
  list(
    Q_wald      = Q_wald_vec,
    lambda_Q_mh = lambda_Q_mh_vec,
    Q_pw        = Q_pw_vec,
    accepted    = accepted_vec,
    resolved    = resolved_vec
  )
}

if (isTRUE(RUN_LEVELSET)) {
  cat("\n============================================================\n")
  cat("EMPIRICAL REPRO LEVEL-SET STUDY\n")
  cat("============================================================\n")
  cat(sprintf("replications                 %d\n", N_LEVELSET_REPS))
  cat(sprintf("grid                         %d x %d\n", LS_N_1, LS_N_2))
  cat(sprintf("total evaluations            %d\n", N_LEVELSET_REPS * LS_N_1 * LS_N_2))
  cat("============================================================\n\n")
  
  levelset_cores <- min(n.cores, N_LEVELSET_REPS)
  cl_ls <- parallel::makePSOCKcluster(levelset_cores)
  doSNOW::registerDoSNOW(cl_ls)
  
  LEVELSET_EXPORTS <- unique(c(
    EXPORT_NAMES, "N_LEVELSET_REPS", "LS_IDX_1", "LS_IDX_2",
    "LS_N_1", "LS_N_2", "grid_1", "grid_2", "one_levelset_replication"
  ))
  parallel::clusterExport(cl_ls, LEVELSET_EXPORTS, envir = environment())
  
  pb_ls <- txtProgressBar(max = N_LEVELSET_REPS, style = 3)
  progress_ls <- function(k) setTxtProgressBar(pb_ls, k)
  
  LEVELSET_RESULTS <- foreach(
    rep_idx = seq_len(N_LEVELSET_REPS),
    .combine = "c",
    .multicombine = TRUE,
    .options.snow = list(progress = progress_ls)
  ) %dopar% {
    list(one_levelset_replication(rep_idx))
  }
  
  try(close(pb_ls), silent = TRUE)
  parallel::stopCluster(cl_ls)
  
  ng <- LS_N_1 * LS_N_2
  QW_SUM       <- numeric(ng)
  QMH_SUM      <- numeric(ng)
  QPW_SUM      <- numeric(ng)
  ACCEPT_SUM   <- numeric(ng)
  RESOLVED_SUM <- numeric(ng)
  
  for (z in LEVELSET_RESULTS) {
    good <- as.logical(z$resolved)
    if (any(good)) {
      QW_SUM[good]  <- QW_SUM[good]  + z$Q_wald[good]
      QMH_SUM[good] <- QMH_SUM[good] + z$lambda_Q_mh[good]
      QPW_SUM[good] <- QPW_SUM[good] + z$Q_pw[good]
    }
    ACCEPT_SUM   <- ACCEPT_SUM   + z$accepted
    RESOLVED_SUM <- RESOLVED_SUM + z$resolved
  }
  
  MEAN_QW          <- ifelse(RESOLVED_SUM > 0, QW_SUM / RESOLVED_SUM, NA_real_)
  MEAN_LAMBDA_QMH  <- ifelse(RESOLVED_SUM > 0, QMH_SUM / RESOLVED_SUM, NA_real_)
  MEAN_QPW         <- ifelse(RESOLVED_SUM > 0, QPW_SUM / RESOLVED_SUM, NA_real_)
  EMPIRICAL_INCLUSION <- ACCEPT_SUM / N_LEVELSET_REPS
  RESOLVED_FRACTION   <- RESOLVED_SUM / N_LEVELSET_REPS
  FAILURE_FRACTION    <- 1 - RESOLVED_FRACTION
  
  Q_W_MEAN_GRID         <- matrix(MEAN_QW, nrow = LS_N_1, ncol = LS_N_2)
  LAMBDA_Q_MH_MEAN_GRID <- matrix(MEAN_LAMBDA_QMH, nrow = LS_N_1, ncol = LS_N_2)
  Q_PW_MEAN_GRID        <- matrix(MEAN_QPW, nrow = LS_N_1, ncol = LS_N_2)
  INCLUSION_GRID        <- matrix(EMPIRICAL_INCLUSION, nrow = LS_N_1, ncol = LS_N_2)
  RESOLVED_GRID         <- matrix(RESOLVED_FRACTION, nrow = LS_N_1, ncol = LS_N_2)
  
  EMPIRICAL_LEVELSET <- expand.grid(beta1 = grid_1, beta0 = grid_2)
  names(EMPIRICAL_LEVELSET) <- c(THETA_NAMES[LS_IDX_1], THETA_NAMES[LS_IDX_2])
  EMPIRICAL_LEVELSET$mean_Q_wald                   <- MEAN_QW
  EMPIRICAL_LEVELSET$mean_lambda_Q_mh              <- MEAN_LAMBDA_QMH
  EMPIRICAL_LEVELSET$mean_Q_pw                     <- MEAN_QPW
  EMPIRICAL_LEVELSET$empirical_inclusion_probability <- EMPIRICAL_INCLUSION
  EMPIRICAL_LEVELSET$resolved_fraction             <- RESOLVED_FRACTION
  EMPIRICAL_LEVELSET$failure_fraction              <- FAILURE_FRACTION
  
  write.csv(EMPIRICAL_LEVELSET, file.path(PROJECT_DIR, "exp2_sigma_empirical_levelset.csv"), row.names = FALSE)
  
  true_i <- which.min(abs(grid_1 - theta_true[LS_IDX_1]))
  true_j <- which.min(abs(grid_2 - theta_true[LS_IDX_2]))
  true_inclusion <- INCLUSION_GRID[true_i, true_j]
  
  LEVELSET_SUMMARY <- data.frame(
    Method                 = "Penalized Wald",
    Type                   = "Empirical Repro inclusion probability",
    Axis_1                 = THETA_NAMES[LS_IDX_1],
    Axis_2                 = THETA_NAMES[LS_IDX_2],
    Sigma_Fixed            = sigma_true,
    Mu_Fixed               = mu_true,
    Tau_Fixed              = tau_true,
    Lambda                 = lambda_n,
    Levelset_Replications  = N_LEVELSET_REPS,
    Grid_N1                = LS_N_1,
    Grid_N2                = LS_N_2,
    Inclusion_At_Truth     = true_inclusion,
    Min_Resolved_Fraction  = min(RESOLVED_FRACTION, na.rm = TRUE),
    Mean_Resolved_Fraction = mean(RESOLVED_FRACTION, na.rm = TRUE),
    stringsAsFactors       = FALSE
  )
  write.csv(LEVELSET_SUMMARY, file.path(PROJECT_DIR, "exp2_sigma_empirical_levelset_summary.csv"), row.names = FALSE)
  
  axis_label <- function(name) {
    switch(name,
           beta1 = expression(beta[1]),
           beta0 = expression(beta[0]),
           mu    = expression(mu),
           tau   = expression(tau),
           sigma = expression(sigma),
           name)
  }
  
  lab_1 <- axis_label(THETA_NAMES[LS_IDX_1])
  lab_2 <- axis_label(THETA_NAMES[LS_IDX_2])
  truth_x <- theta_true[LS_IDX_1]
  truth_y <- theta_true[LS_IDX_2]
  
  publication_contour_levels <- function(Z) {
    good <- Z[is.finite(Z)]
    if (length(good) < 20L || diff(range(good)) <= 0) return(NULL)
    levels <- unname(quantile(good, probs = c(0.25, 0.45, 0.65, 0.82), na.rm = TRUE))
    levels <- sort(unique(levels))
    if (length(levels) < 2L) return(NULL)
    levels
  }
  
  LEV_W  <- publication_contour_levels(Q_W_MEAN_GRID)
  LEV_MH <- publication_contour_levels(LAMBDA_Q_MH_MEAN_GRID)
  LEV_PW <- publication_contour_levels(Q_PW_MEAN_GRID)
  
  setup_geometry_panel <- function() {
    plot(NA_real_, NA_real_, xlim = range(grid_1), ylim = range(grid_2),
         xaxs = "i", yaxs = "i", xlab = lab_1, ylab = lab_2, bty = "l", las = 1)
  }
  
  draw_truth <- function() {
    points(truth_x, truth_y, pch = 4, cex = 1.18, lwd = 2.2, col = COL_BLACK)
  }
  
  # Three-panel geometry figure
  draw_geometry_figure <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(
      mfrow = c(1, 3), family = "sans", mar = c(4.0, 4.0, 2.0, 0.45),
      oma = c(0, 0, 0, 0), mgp = c(2.35, 0.62, 0), tcl = -0.20, las = 1,
      cex.axis = PUB_CEX_AXIS, cex.lab = PUB_CEX_LAB
    )
    
    # Panel (a): Wald component
    setup_geometry_panel()
    if (!is.null(LEV_W)) {
      contour(x = grid_1, y = grid_2, z = Q_W_MEAN_GRID, levels = LEV_W,
              add = TRUE, drawlabels = FALSE, col = COL_BLUE, lwd = 1.55)
    }
    draw_truth()
    mtext("(a) Wald component", side = 3, adj = 0, line = 0.30, font = 2, cex = PUB_CEX_PANEL)
    
    # Panel (b): Penalization component
    setup_geometry_panel()
    if (!is.null(LEV_MH)) {
      contour(x = grid_1, y = grid_2, z = LAMBDA_Q_MH_MEAN_GRID, levels = LEV_MH,
              add = TRUE, drawlabels = FALSE, col = COL_ORANGE, lwd = 1.55)
    }
    draw_truth()
    mtext("(b) Penalization component", side = 3, adj = 0, line = 0.30, font = 2, cex = PUB_CEX_PANEL)
    
    # Panel (c): Penalized Wald and Repro set
    setup_geometry_panel()
    if (!is.null(LEV_PW)) {
      contour(x = grid_1, y = grid_2, z = Q_PW_MEAN_GRID, levels = LEV_PW,
              add = TRUE, drawlabels = FALSE, col = COL_PURPLE, lwd = 1.30)
    }
    
    inclusion_levels <- c(0.50, 0.80, 0.90, 0.95)
    inclusion_lty    <- c(3, 2, 5, 1)
    inclusion_lwd    <- c(1.25, 1.45, 1.65, 2.15)
    finite_inc       <- INCLUSION_GRID[is.finite(INCLUSION_GRID)]
    shown            <- logical(length(inclusion_levels))
    
    for (k in seq_along(inclusion_levels)) {
      lev <- inclusion_levels[k]
      shown[k] <- (length(finite_inc) > 0L && min(finite_inc) < lev && max(finite_inc) > lev)
      if (shown[k]) {
        contour(x = grid_1, y = grid_2, z = INCLUSION_GRID, levels = lev,
                add = TRUE, drawlabels = FALSE, col = COL_BLACK,
                lty = inclusion_lty[k], lwd = inclusion_lwd[k])
      }
    }
    draw_truth()
    
    shown_idx <- which(shown)
    legend_text <- c(expression(Q^{PW}), sprintf("Inclusion %.2f", inclusion_levels[shown_idx]), "True parameter")
    legend(
      "topright",
      legend = legend_text,
      col = c(COL_PURPLE, rep(COL_BLACK, length(shown_idx)), COL_BLACK),
      lty = c(1, inclusion_lty[shown_idx], NA),
      lwd = c(1.30, inclusion_lwd[shown_idx], 2),
      pch = c(NA, rep(NA, length(shown_idx)), 4),
      pt.cex = 0.9,
      bty = "n",
      cex = 0.65,
      x.intersp = 0.65,
      y.intersp = 0.92
    )
    mtext("(c) Penalized Wald and Repro set", side = 3, adj = 0, line = 0.30, font = 2, cex = PUB_CEX_PANEL)
  }
  
  pdf(file.path(PROJECT_DIR, "Figure_Exp2_sigma_empirical_geometry.pdf"), width = 11.2, height = 3.7, family = PUB_FONT, useDingbats = FALSE)
  draw_geometry_figure()
  invisible(dev.off())
  
  png(file.path(PROJECT_DIR, "Figure_Exp2_sigma_empirical_geometry.png"), width = 3360, height = 1110, res = 300, bg = "white")
  draw_geometry_figure()
  invisible(dev.off())
  
  # Dedicated empirical Repro inclusion figure
  prob_levels <- seq(0, 1, by = 0.05)
  prob_cols   <- hcl.colors(length(prob_levels) - 1L, palette = "YlOrRd")
  inclusion_levels <- c(0.50, 0.80, 0.90, 0.95)
  inclusion_lty    <- c(3, 2, 5, 1)
  inclusion_lwd    <- c(1.25, 1.45, 1.70, 2.25)
  
  draw_probability_panel <- function() {
    par(family = "sans", mar = c(4.2, 4.4, 0.8, 0.35), mgp = c(2.45, 0.65, 0),
        tcl = -0.22, las = 1, cex.axis = PUB_CEX_AXIS, cex.lab = PUB_CEX_LAB)
    
    image(x = grid_1, y = grid_2, z = INCLUSION_GRID, breaks = prob_levels,
          col = prob_cols, xlim = range(grid_1), ylim = range(grid_2),
          xaxs = "i", yaxs = "i", xlab = lab_1, ylab = lab_2, axes = FALSE, useRaster = TRUE)
    axis(1)
    axis(2, las = 1)
    box()
    
    finite_inc <- INCLUSION_GRID[is.finite(INCLUSION_GRID)]
    shown <- logical(length(inclusion_levels))
    for (k in seq_along(inclusion_levels)) {
      lev <- inclusion_levels[k]
      shown[k] <- (length(finite_inc) > 0L && min(finite_inc) < lev && max(finite_inc) > lev)
      if (shown[k]) {
        contour(x = grid_1, y = grid_2, z = INCLUSION_GRID, levels = lev,
                add = TRUE, drawlabels = FALSE, col = COL_BLACK,
                lty = inclusion_lty[k], lwd = inclusion_lwd[k])
      }
    }
    points(truth_x, truth_y, pch = 4, cex = 1.25, lwd = 2.5, col = COL_BLACK)
    
    shown_idx <- which(shown)
    legend(
      "bottomleft", inset = 0.02,
      legend = c(sprintf("Inclusion %.2f", inclusion_levels[shown_idx]), "True parameter"),
      col = c(rep(COL_BLACK, length(shown_idx)), COL_BLACK),
      lty = c(inclusion_lty[shown_idx], NA),
      lwd = c(inclusion_lwd[shown_idx], 2),
      pch = c(rep(NA, length(shown_idx)), 4),
      pt.cex = 0.9, bty = "n", cex = 0.73, x.intersp = 0.7, y.intersp = 0.93
    )
  }
  
  draw_probability_colorbar <- function() {
    par(family = "sans", mar = c(4.2, 0.7, 0.8, 3.2), tcl = -0.22)
    yvals <- seq(0, 1, length.out = 201L)
    image(x = c(0, 1), y = yvals, z = matrix(rep(yvals, each = 2L), nrow = 2L),
          col = prob_cols, breaks = prob_levels, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
    axis(4, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)), las = 1, cex.axis = PUB_CEX_AXIS)
    box()
    mtext(expression(hat(pi)(beta[1], beta[0])), side = 3, line = 0.25, cex = 0.88)
  }
  
  draw_empirical_repro_figure <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit({ layout(1); par(old_par) }, add = TRUE)
    layout(matrix(c(1, 2), nrow = 1L), widths = c(6.0, 0.82))
    draw_probability_panel()
    draw_probability_colorbar()
  }
  
  pdf(file.path(PROJECT_DIR, "Figure_Exp2_sigma_empirical_Repro_region.pdf"), width = 6.8, height = 5.0, family = PUB_FONT, useDingbats = FALSE)
  draw_empirical_repro_figure()
  invisible(dev.off())
  
  png(file.path(PROJECT_DIR, "Figure_Exp2_sigma_empirical_Repro_region.png"), width = 2040, height = 1500, res = 300, bg = "white")
  draw_empirical_repro_figure()
  invisible(dev.off())
  
  draw_empirical_repro_bw <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(family = "sans", mar = c(4.2, 4.4, 0.8, 0.8), mgp = c(2.45, 0.65, 0),
        tcl = -0.22, las = 1, cex.axis = PUB_CEX_AXIS, cex.lab = PUB_CEX_LAB)
    plot(NA, xlim = range(grid_1), ylim = range(grid_2), xaxs = "i", yaxs = "i", xlab = lab_1, ylab = lab_2, bty = "l")
    
    finite_inc <- INCLUSION_GRID[is.finite(INCLUSION_GRID)]
    shown <- logical(length(inclusion_levels))
    for (k in seq_along(inclusion_levels)) {
      lev <- inclusion_levels[k]
      shown[k] <- (length(finite_inc) > 0L && min(finite_inc) < lev && max(finite_inc) > lev)
      if (shown[k]) {
        contour(x = grid_1, y = grid_2, z = INCLUSION_GRID, levels = lev,
                add = TRUE, drawlabels = FALSE, col = COL_BLACK,
                lty = inclusion_lty[k], lwd = inclusion_lwd[k])
      }
    }
    points(truth_x, truth_y, pch = 4, cex = 1.25, lwd = 2.5, col = COL_BLACK)
    
    shown_idx <- which(shown)
    legend(
      "topright",
      legend = c(sprintf("Inclusion %.2f", inclusion_levels[shown_idx]), "True parameter"),
      col = COL_BLACK,
      lty = c(inclusion_lty[shown_idx], NA),
      lwd = c(inclusion_lwd[shown_idx], 2),
      pch = c(rep(NA, length(shown_idx)), 4),
      bty = "n",
      cex = PUB_CEX_LEGEND
    )
  }
  
  pdf(file.path(PROJECT_DIR, "Figure_Exp2_sigma_empirical_Repro_region_BW.pdf"), width = 5.8, height = 4.7, family = PUB_FONT, useDingbats = FALSE)
  draw_empirical_repro_bw()
  invisible(dev.off())
}

# ------------------------------------------------------------------------------
# 15. Nuisance-Search Profiling Audit
# ------------------------------------------------------------------------------

if (isTRUE(RUN_SEARCH_AUDIT)) {
  cat("\n============================================================\n")
  cat("RUNNING NUISANCE-SEARCH AUDIT\n")
  cat("============================================================\n")
  
  audit_one_rep <- function(rep_idx) {
    set.seed(rep_idx + 1000L)
    obs_seed <- draw_seed_cloud(1L, n_obs)
    s_obs <- as.numeric(release_from_seed(obs_seed, theta_true)$S)
    bad <- c(rep_idx = rep_idx, n_checked = 0, n_overturned = 0, failed = 1)
    if (any(!is.finite(s_obs))) return(bad)
    
    syn_seed <- draw_seed_cloud(R_synthetic, n_obs)
    aux_seed <- draw_seed_cloud(R_aux, n_obs)
    
    tryCatch({
      cache_prod   <- cache_new()
      cache_strong <- cache_new()
      n_checked    <- 0L
      n_overturned <- 0L
      
      for (sigma_candidate in sigma_scan) {
        prod_fit <- profile_at_sigma(sigma_candidate, syn_seed, s_obs, aux_seed, cache_prod)
        if (!isTRUE(prod_fit$resolved) || isTRUE(prod_fit$accepted)) next
        
        n_checked <- n_checked + 1L
        strong_fit <- profile_at_sigma(
          sigma_candidate, syn_seed, s_obs, aux_seed, cache_strong,
          n_halton = AUDIT_HALTON_STARTS, n_refine = AUDIT_REFINE_STARTS
        )
        if (isTRUE(strong_fit$resolved) && isTRUE(strong_fit$accepted)) {
          n_overturned <- n_overturned + 1L
        }
      }
      c(rep_idx = rep_idx, n_checked = n_checked, n_overturned = n_overturned, failed = 0)
    }, error = function(e) bad)
  }
  
  n_audit_cores <- min(n.cores, N_AUDIT_REPS)
  cl_audit <- parallel::makePSOCKcluster(n_audit_cores)
  doSNOW::registerDoSNOW(cl_audit)
  parallel::clusterExport(cl_audit, EXPORT_NAMES, envir = environment())
  parallel::clusterExport(cl_audit, "audit_one_rep", envir = environment())
  
  pb_audit <- txtProgressBar(max = N_AUDIT_REPS, style = 3)
  progress_audit <- function(k) setTxtProgressBar(pb_audit, k)
  
  AUDIT <- foreach(
    rep_idx = seq_len(N_AUDIT_REPS),
    .combine = "rbind",
    .options.snow = list(progress = progress_audit)
  ) %dopar% {
    audit_one_rep(rep_idx)
  }
  
  try(close(pb_audit), silent = TRUE)
  parallel::stopCluster(cl_audit)
  
  AUDIT <- as.data.frame(AUDIT)
  audit_checked    <- sum(AUDIT$n_checked, na.rm = TRUE)
  audit_overturned <- sum(AUDIT$n_overturned, na.rm = TRUE)
  overturn_rate    <- if (audit_checked > 0L) audit_overturned / audit_checked else NA_real_
  
  AUDIT_SUMMARY <- data.frame(
    Audit_Replications            = N_AUDIT_REPS,
    Production_Halton             = N_HALTON_STARTS,
    Production_Refinements        = N_REFINE_STARTS,
    Audit_Halton                  = AUDIT_HALTON_STARTS,
    Audit_Refinements             = AUDIT_REFINE_STARTS,
    Rejected_Sigma_Values_Checked = audit_checked,
    Rejections_Overturned         = audit_overturned,
    Overturn_Rate                 = overturn_rate,
    Audit_Failures                = sum(AUDIT$failed, na.rm = TRUE),
    stringsAsFactors              = FALSE
  )
  
  write.csv(AUDIT, file.path(PROJECT_DIR, "exp2_sigma_search_audit_raw.csv"), row.names = FALSE)
  write.csv(AUDIT_SUMMARY, file.path(PROJECT_DIR, "exp2_sigma_search_audit_summary.csv"), row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("NUMERICAL PROFILING AUDIT\n")
  cat("============================================================\n")
  print(t(AUDIT_SUMMARY))
  cat("============================================================\n")
}

# ------------------------------------------------------------------------------
# 16. Metadata Export & Final Reporting
# ------------------------------------------------------------------------------

SETTINGS <- data.frame(
  Setting = c(
    "n_obs", "beta1_true", "beta0_true", "mu_true", "tau_true", "sigma_true",
    "Delta", "mu_GDP", "R_synthetic", "R_aux", "alpha", "lambda_n", "nSIM",
    "sigma_scan_step", "CI_tol", "N_halton_starts", "N_refine_starts",
    "N_levelset_reps", "levelset_grid_1", "levelset_grid_2", "N_audit_reps"
  ),
  Value = c(
    n_obs, beta1_true, beta0_true, mu_true, tau_true, sigma_true,
    Delta, mu_gdp, R_synthetic, R_aux, alpha, lambda_n, nSIM,
    SIGMA_SCAN_STEP, CI_TOL, N_HALTON_STARTS, N_REFINE_STARTS,
    N_LEVELSET_REPS, LS_N_1, LS_N_2, N_AUDIT_REPS
  ),
  stringsAsFactors = FALSE
)
write.csv(SETTINGS, file.path(PROJECT_DIR, "exp2_sigma_settings.csv"), row.names = FALSE)

caption_performance <- paste0(
  "\\caption{",
  "\\textbf{Finite-sample performance of the $\\sigma$-target Penalized Wald Repro confidence interval.} ",
  "Panel (a) reports empirical coverage of the true regression-error standard deviation $\\sigma^*=0.5$, ",
  "together with an exact 95\\% Clopper--Pearson Monte Carlo interval. The horizontal reference lines ",
  "indicate the nominal 95\\% level and the corresponding finite-$R$ reference level induced by the ",
  "inclusive-rank rule. Panel (b) displays the empirical distribution of confidence-interval widths among ",
  "replications yielding finite intervals. The solid and dashed vertical lines indicate the mean and median ",
  "widths, respectively. Numerical failures and disconnected profiled acceptance sets are counted as noncoverage.",
  "}"
)

caption_geometry <- paste0(
  "\\caption{",
  "\\textbf{Geometry of the $\\sigma$-target Penalized Wald Repro procedure.} ",
  "The figure displays a two-dimensional nuisance-parameter section in $(\\beta_1,\\beta_0)$, with ",
  "$(\\mu,\\tau,\\sigma)=(0.5,1,0.5)$ fixed at their data-generating values. Panel (a) displays level sets of ",
  "the Monte Carlo mean nuisance-orthogonal Wald component $Q^W$. Panel (b) displays level sets of the Monte ",
  "Carlo mean penalization component $\\lambda_nQ^{MH}$. Panel (c) displays level sets of the Monte Carlo mean ",
  "Penalized Wald statistic $Q^{PW}=Q^W+\\lambda_nQ^{MH}$ together with level sets of the empirical Repro ",
  "inclusion probability. The black cross denotes the data-generating nuisance parameter $(\\beta_1^*,\\beta_0^*)=(1,-0.5)$.",
  "}"
)

caption_empirical <- paste0(
  "\\caption{",
  "\\textbf{Empirical Repro inclusion probability for the $\\sigma$-target Penalized Wald procedure.} ",
  "The figure shows a two-dimensional nuisance-parameter section in $(\\beta_1,\\beta_0)$ with ",
  "$(\\mu,\\tau,\\sigma)=(0.5,1,0.5)$ fixed at their data-generating values. At each grid point, the shading represents ",
  "$\\widehat\\pi(\\beta_1,\\beta_0)=B^{-1}\\sum_{b=1}^{B}\\mathbf{1}\\{(\\beta_1,\\beta_0,\\mu^*,\\tau^*,\\sigma^*)\\in C_b\\}$, ",
  "where $C_b$ is the finite-sample Repro acceptance set in Monte Carlo replication $b$. Contours correspond to empirical ",
  "inclusion probabilities $0.50$, $0.80$, $0.90$, and $0.95$, and the black cross denotes the data-generating nuisance ",
  "parameter $(\\beta_1^*,\\beta_0^*)=(1,-0.5)$. The $0.95$ contour is a level set of repeated-sampling inclusion probability ",
  "and is not itself a 95\\% confidence region.",
  "}"
)

writeLines(
  c("% Figure 1", caption_performance, "", "% Figure 2", caption_geometry, "", "% Figure 3", caption_empirical),
  con = file.path(PROJECT_DIR, "Figure_Exp2_sigma_captions.tex")
)

capture.output(sessionInfo(), file = file.path(PROJECT_DIR, "sessionInfo.txt"))

cat("\n\n============================================================\n")
cat("DONE -- EXPERIMENT 2: SIGMA TARGET\n")
cat("============================================================\n")
cat(sprintf("Output directory             %s\n", PROJECT_DIR))
cat(sprintf("n                            %d\n", n_obs))
cat(sprintf("True sigma                   %.4f\n", sigma_true))
cat(sprintf("Coverage                     %.4f\n", coverage))
cat(sprintf("Coverage MC SE               %.4f\n", coverage_se))
cat(sprintf("95%% Clopper-Pearson          [%.4f, %.4f]\n", cp_lower, cp_upper))
cat(sprintf("Finite-R reference           %.6f\n", finite_R_reference_level))
cat(sprintf("Mean CI width                %.4f\n", mean_width))
cat(sprintf("Median CI width              %.4f\n", median_width))
cat(sprintf("Width IQR                    [%.4f, %.4f]\n", width_q25, width_q75))
cat(sprintf("Failure rate                 %.4f\n", failure_rate))
cat(sprintf("Disconnected-set rate        %.4f\n", disconnected_rate))
cat(sprintf("Unresolved eval fraction     %.6f\n", unresolved_fraction))
cat(sprintf("Unique geometry failures     %d\n", as.integer(total_unique_geomfail)))
cat(sprintf("Worst Wald identity error    %.3e\n", worst_identity_wald))
cat(sprintf("Worst MH identity error      %.3e\n", worst_identity_mh))
cat(sprintf("Exact rank cutoff            %d\n", rank_threshold))
cat(sprintf("lambda_n                     %.6f\n", lambda_n))
cat(sprintf("Penalty share                %.3f\n", penalty_share_reference))
cat(sprintf("Median candidate evals       %.0f\n", median_evals))

if (isTRUE(RUN_LEVELSET)) {
  cat(sprintf("Empirical level-set reps     %d\n", N_LEVELSET_REPS))
  cat(sprintf("Empirical level-set grid     %d x %d\n", LS_N_1, LS_N_2))
  cat(sprintf("Inclusion at true point      %.4f\n", true_inclusion))
}

cat(sprintf("Elapsed                      %s\n", format(Sys.time() - start_time)))
cat("============================================================\n")

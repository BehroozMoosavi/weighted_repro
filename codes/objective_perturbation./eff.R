# ==============================================================================
# Script: eff.R
# Description: Evaluates finite-sample confidence intervals for the target 
#              parameter beta1 in logistic regression under differentially private 
#              objective perturbation (Operator II / Efficient Repro method). 
#              Searches parameter boxes for confidence interval inversion rather 
#              than fixed midpoints and records numerical non-resolutions.
#
# Simulation grid:
#   n       = 100, 200, 500, 1000
#   epsilon = 0.1, 0.3, 1, 3
#   reps    = 1000
#   R       = 200 (inference cloud)
#   R_aux   = 100 (auxiliary direction cloud)
#   alpha   = 0.10
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Thread Control & Global Setup
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
# 1. Experimental Settings & Tolerances
# ------------------------------------------------------------------------------

n_values   <- c(100L, 200L, 500L, 1000L)
eps_values <- c(0.1, 0.3, 1, 3)

reps  <- 1000L
R     <- 200L
R_aux <- 100L
alpha <- 0.10

shape1 <- 1.0
shape2 <- 1.0

beta0      <- 0.5
beta1_true <- 2.0

M_DIM <- 4L

# Efficient direction approximation
FD_STEP            <- 0.05
USE_DIRECTION_GRID <- TRUE

AUX_COV_FLOOR <- 1e-10
NUIS_FLOOR    <- 1e-10
INF_COV_FLOOR <- 1e-10
DENOM_TOL     <- 1e-14

NEWTON_MAXIT  <- 100L
NEWTON_TOL    <- 1e-7
NEWTON_ACCEPT <- 1e-5
MAX_STEP      <- 2.0
DET_FLOOR     <- 1e-14

ROBUST_BFGS_MAXIT   <- 300L
ROBUST_NLMINB_MAXIT <- 500L
ROBUST_GRAD_TOL     <- 1e-6

BOX_NM_MAXIT <- 100L

BETA1_LO <- -10
BETA1_HI <-  10

# Search bounds for nuisance vector: (beta0, log(shape1), log(shape2))
NUIS_LO <- c(-10, -5, -5)
NUIS_HI <- c( 10,  5,  5)

CI_TOL <- 1e-4

rank_threshold <- floor(alpha * (R + 1L)) + 1L
rank_reference <- (R + 2L - rank_threshold) / (R + 1L)

PROJECT_DIR <- path.expand("~/R_Simuls/exp3_objectiveperturb/eff")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

FRESH_START <- FALSE

# ------------------------------------------------------------------------------
# 2. Package Dependencies
# ------------------------------------------------------------------------------

required_packages <- c("foreach", "doSNOW", "tictoc")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall once with install.packages()."
  )
}

suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(tictoc))

# ------------------------------------------------------------------------------
# 3. CPU Core Detection
# ------------------------------------------------------------------------------

parse_first_int <- function(x) {
  if (!nzchar(x)) return(NA_integer_)
  
  y <- suppressWarnings(
    as.integer(sub("[^0-9].*$", "", x))
  )
  
  if (length(y) == 0L || is.na(y) || y < 1L) {
    return(NA_integer_)
  }
  
  y
}

slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_PER_TASK", unset = ""))
if (is.na(slurm_cpus)) {
  slurm_cpus <- parse_first_int(Sys.getenv("SLURM_CPUS_ON_NODE", unset = ""))
}
if (is.na(slurm_cpus)) {
  slurm_cpus <- parallel::detectCores(logical = FALSE)
}
if (is.na(slurm_cpus) || slurm_cpus < 1L) {
  slurm_cpus <- 1L
}

n_workers <- min(125L, slurm_cpus, reps)

cat("\n============================================================\n")
cat(" OBJECTIVE PERTURBATION -- CORRECTED OPERATOR II\n")
cat("============================================================\n")
cat(sprintf("Host                 : %s\n", Sys.info()[["nodename"]]))
cat(sprintf("Allocated CPUs       : %d\n", slurm_cpus))
cat(sprintf("Workers              : %d\n", n_workers))
cat(sprintf("Replications / config: %d\n", reps))
cat(sprintf("Number configs       : %d\n", length(n_values) * length(eps_values)))
cat(sprintf("Results root         : %s\n", PROJECT_DIR))
cat("============================================================\n\n")

# ------------------------------------------------------------------------------
# 4. Error Condition Types
# ------------------------------------------------------------------------------

stat_fail <- function(msg = "statistic not evaluable") {
  stop(
    structure(
      list(message = msg, call = NULL),
      class = c("stat_fail", "error", "condition")
    )
  )
}

cert_hit <- function(par, M) {
  stop(
    structure(
      list(
        message = "acceptance certificate",
        call = NULL,
        par = par,
        M = M
      ),
      class = c("cert_hit", "error", "condition")
    )
  )
}

# ------------------------------------------------------------------------------
# 5. K-Norm Geometry & Rejection Sampling
# ------------------------------------------------------------------------------

lower <- function(x) {
  bottom <- (x + 1)^2 - 1
  middle <- x - 1 / 4
  top    <- x^2
  
  bottom * (x <= -1 / 2) +
    middle * (x > -1 / 2 & x < 1 / 2) +
    top * (x >= 1 / 2)
}

upper <- function(x) {
  bottom <- -x^2
  middle <- x + 1 / 4
  top    <- -(x - 1)^2 + 1
  
  bottom * (x <= -1 / 2) +
    middle * (x > -1 / 2 & x < 1 / 2) +
    top * (x >= 1 / 2)
}

draw_U2 <- function(k) {
  acc <- matrix(numeric(0), nrow = 0L, ncol = 2L)
  
  while (nrow(acc) < k) {
    need <- k - nrow(acc)
    bs <- max(100L, 20L * need)
    
    W <- matrix(
      runif(bs * 2L, min = -1, max = 1),
      nrow = bs,
      ncol = 2L
    )
    
    keep <- (W[, 2] >= lower(W[, 1])) & (W[, 2] <= upper(W[, 1]))
    
    if (any(keep)) {
      acc <- rbind(acc, W[keep, , drop = FALSE])
    }
  }
  
  acc[seq_len(k), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# 6. Eigendecomposition-Based Matrix Inversion
# ------------------------------------------------------------------------------

safe_inv <- function(A, floor_value) {
  if (!all(is.finite(A))) {
    return(NULL)
  }
  
  A <- 0.5 * (A + t(A))
  
  eg <- tryCatch(
    eigen(A, symmetric = TRUE),
    error = function(e) NULL
  )
  
  if (
    is.null(eg) ||
    !all(is.finite(eg$values)) ||
    !all(is.finite(eg$vectors))
  ) {
    return(NULL)
  }
  
  ev <- pmax(eg$values, floor_value)
  
  eg$vectors %*%
    diag(1 / ev, nrow = length(ev)) %*%
    t(eg$vectors)
}

# ------------------------------------------------------------------------------
# 7. Objective-Perturbation Optimization Solvers
# ------------------------------------------------------------------------------

softplus <- function(z) {
  out <- numeric(length(z))
  positive <- z > 0
  
  out[positive] <- z[positive] + log1p(exp(-z[positive]))
  out[!positive] <- log1p(exp(z[!positive]))
  
  out
}

robust_objperturb_row <- function(
    x2,
    y,
    Dn,
    noise0,
    noise1,
    start = c(0, 0)
) {
  if (
    !all(is.finite(x2)) ||
    !all(is.finite(y)) ||
    !is.finite(Dn) ||
    !is.finite(noise0) ||
    !is.finite(noise1) ||
    Dn <= 0
  ) {
    return(c(NA_real_, NA_real_))
  }
  
  if (length(start) != 2L || !all(is.finite(start))) {
    start <- c(0, 0)
  }
  
  objective <- function(par) {
    if (!all(is.finite(par))) return(Inf)
    
    b0 <- par[1]
    b1 <- par[2]
    eta <- b0 + b1 * x2
    
    val <- mean(softplus(eta) - y * eta) +
      0.5 * Dn * (b0^2 + b1^2) +
      noise0 * b0 +
      noise1 * b1
    
    if (is.finite(val)) val else Inf
  }
  
  gradient <- function(par) {
    if (!all(is.finite(par))) {
      return(c(NA_real_, NA_real_))
    }
    
    b0 <- par[1]
    b1 <- par[2]
    
    eta <- b0 + b1 * x2
    p <- stats::plogis(eta)
    residual <- p - y
    
    c(
      mean(residual) + Dn * b0 + noise0,
      mean(x2 * residual) + Dn * b1 + noise1
    )
  }
  
  certify <- function(par) {
    if (length(par) != 2L || !all(is.finite(par))) {
      return(FALSE)
    }
    
    g <- gradient(par)
    all(is.finite(g)) && max(abs(g)) <= ROBUST_GRAD_TOL
  }
  
  fit1 <- tryCatch(
    optim(
      par = start,
      fn = objective,
      gr = gradient,
      method = "BFGS",
      control = list(maxit = ROBUST_BFGS_MAXIT, reltol = 1e-12)
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit1) && certify(fit1$par)) {
    return(as.numeric(fit1$par))
  }
  
  start2 <- if (!is.null(fit1) && length(fit1$par) == 2L && all(is.finite(fit1$par))) {
    fit1$par
  } else {
    start
  }
  
  fit2 <- tryCatch(
    nlminb(
      start = start2,
      objective = objective,
      gradient = gradient,
      control = list(
        iter.max = ROBUST_NLMINB_MAXIT,
        eval.max = 2L * ROBUST_NLMINB_MAXIT,
        rel.tol = 1e-12,
        x.tol = 1e-12
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit2) && certify(fit2$par)) {
    return(as.numeric(fit2$par))
  }
  
  c(NA_real_, NA_real_)
}

# ------------------------------------------------------------------------------
# 8. Private Summary Release: s = (beta0~, beta1~, mean(z)~, mean(z^2)~)
# ------------------------------------------------------------------------------

sdp_matrix <- function(
    seed_mat,
    U_mat,
    N1_mat,
    N2_mat,
    ep,
    theta,
    n_val,
    x_pre = NULL
) {
  nR <- nrow(seed_mat)
  bad_output <- matrix(NA_real_, nrow = nR, ncol = M_DIM)
  
  if (length(theta) != 4L || !all(is.finite(theta))) {
    return(bad_output)
  }
  
  beta1_val  <- theta[1]
  beta0_val  <- theta[2]
  shape1_val <- exp(theta[3])
  shape2_val <- exp(theta[4])
  
  if (
    !is.finite(shape1_val) ||
    !is.finite(shape2_val) ||
    shape1_val <= 0 ||
    shape2_val <= 0
  ) {
    return(bad_output)
  }
  
  ep1 <- 0.9 * ep
  ep2 <- ep - ep1
  
  if (
    !is.finite(ep1) ||
    !is.finite(ep2) ||
    ep1 <= 0 ||
    ep2 <= 0
  ) {
    return(bad_output)
  }
  
  x <- if (is.null(x_pre)) {
    matrix(
      qbeta(seed_mat, shape1_val, shape2_val),
      nrow = nR,
      ncol = n_val
    )
  } else {
    x_pre
  }
  
  if (
    !is.matrix(x) ||
    nrow(x) != nR ||
    ncol(x) != n_val ||
    !all(is.finite(x))
  ) {
    return(bad_output)
  }
  
  x2 <- 2 * x - 1
  p_true <- stats::plogis(beta0_val + beta1_val * x2)
  Y <- U_mat <= p_true
  
  # Awan-Slavkovic l_infinity objective perturbation parameters
  lambda_value <- 0.5
  Delta <- lambda_value / (exp(0.15 * ep1) - 1)
  b_scale <- 2 / (0.85 * ep1)
  
  inv_n <- 1 / n_val
  Dn <- Delta * inv_n
  
  b0_noise <- b_scale * N1_mat[, 1] * inv_n
  b1_noise <- b_scale * N1_mat[, 2] * inv_n
  
  b0 <- numeric(nR)
  b1 <- numeric(nR)
  converged <- logical(nR)
  failed <- logical(nR)
  
  # Vectorized Newton-Raphson iterations
  for (iteration in seq_len(NEWTON_MAXIT)) {
    active <- which(!converged & !failed)
    if (length(active) == 0L) break
    
    xa <- x2[active, , drop = FALSE]
    Ya <- Y[active, , drop = FALSE]
    
    b0a <- b0[active]
    b1a <- b1[active]
    
    eta <- xa * b1a + b0a
    p <- stats::plogis(eta)
    residual <- Ya - p
    
    g0 <- -inv_n * rowSums(residual) + Dn * b0a + b0_noise[active]
    g1 <- -inv_n * rowSums(xa * residual) + Dn * b1a + b1_noise[active]
    
    gn <- pmax(abs(g0), abs(g1))
    conv <- is.finite(gn) & gn <= NEWTON_TOL
    
    if (any(conv)) {
      converged[active[conv]] <- TRUE
    }
    
    keep <- which(!conv)
    if (length(keep) == 0L) next
    
    idx <- active[keep]
    xa <- xa[keep, , drop = FALSE]
    p  <- p[keep, , drop = FALSE]
    g0 <- g0[keep]
    g1 <- g1[keep]
    
    W <- p * (1 - p)
    h11 <- inv_n * rowSums(W) + Dn
    h12 <- inv_n * rowSums(W * xa)
    h22 <- inv_n * rowSums(W * xa * xa) + Dn
    
    detH <- h11 * h22 - h12 * h12
    
    good <- (
      is.finite(g0) & is.finite(g1) &
        is.finite(h11) & is.finite(h12) & is.finite(h22) &
        is.finite(detH) & detH > DET_FLOOR
    )
    
    if (any(!good)) {
      failed[idx[!good]] <- TRUE
    }
    
    if (!any(good)) next
    
    idx_good <- idx[good]
    step0 <- (h22[good] * g0[good] - h12[good] * g1[good]) / detH[good]
    step1 <- (-h12[good] * g0[good] + h11[good] * g1[good]) / detH[good]
    
    finite_step <- is.finite(step0) & is.finite(step1)
    if (any(!finite_step)) {
      failed[idx_good[!finite_step]] <- TRUE
    }
    
    if (!any(finite_step)) next
    
    idx_good <- idx_good[finite_step]
    step0 <- step0[finite_step]
    step1 <- step1[finite_step]
    
    step_norm <- sqrt(step0^2 + step1^2)
    damp <- pmin(1, MAX_STEP / pmax(step_norm, 1e-300))
    
    b0[idx_good] <- b0[idx_good] - damp * step0
    b1[idx_good] <- b1[idx_good] - damp * step1
  }
  
  # Check first-order condition
  candidates <- which(!failed)
  if (length(candidates) > 0L) {
    xa <- x2[candidates, , drop = FALSE]
    Ya <- Y[candidates, , drop = FALSE]
    
    eta <- xa * b1[candidates] + b0[candidates]
    p <- stats::plogis(eta)
    residual <- Ya - p
    
    g0 <- -inv_n * rowSums(residual) + Dn * b0[candidates] + b0_noise[candidates]
    g1 <- -inv_n * rowSums(xa * residual) + Dn * b1[candidates] + b1_noise[candidates]
    
    gf <- pmax(abs(g0), abs(g1))
    bad_final <- (
      !is.finite(gf) | gf > NEWTON_ACCEPT |
        !is.finite(b0[candidates]) | !is.finite(b1[candidates])
    )
    
    if (any(bad_final)) {
      failed[candidates[bad_final]] <- TRUE
    }
  }
  
  # Numerical fallback
  fallback_rows <- which(failed)
  if (length(fallback_rows) > 0L) {
    for (ii in fallback_rows) {
      start_i <- c(b0[ii], b1[ii])
      if (!all(is.finite(start_i))) {
        start_i <- c(0, 0)
      }
      
      sol <- robust_objperturb_row(
        x2 = x2[ii, , drop = TRUE],
        y = as.numeric(Y[ii, , drop = TRUE]),
        Dn = Dn,
        noise0 = b0_noise[ii],
        noise1 = b1_noise[ii],
        start = start_i
      )
      
      if (length(sol) == 2L && all(is.finite(sol))) {
        b0[ii] <- sol[1]
        b1[ii] <- sol[2]
        failed[ii] <- FALSE
      }
    }
  }
  
  # Final convergence check
  final_rows <- which(!failed)
  if (length(final_rows) > 0L) {
    xa <- x2[final_rows, , drop = FALSE]
    Ya <- Y[final_rows, , drop = FALSE]
    
    eta <- xa * b1[final_rows] + b0[final_rows]
    p <- stats::plogis(eta)
    residual <- Ya - p
    
    g0 <- -inv_n * rowSums(residual) + Dn * b0[final_rows] + b0_noise[final_rows]
    g1 <- -inv_n * rowSums(xa * residual) + Dn * b1[final_rows] + b1_noise[final_rows]
    
    final_gradient <- pmax(abs(g0), abs(g1))
    final_bad <- (
      !is.finite(final_gradient) | final_gradient > ROBUST_GRAD_TOL |
        !is.finite(b0[final_rows]) | !is.finite(b1[final_rows])
    )
    
    if (any(final_bad)) {
      failed[final_rows[final_bad]] <- TRUE
    }
  }
  
  if (any(failed)) {
    b0[failed] <- NA_real_
    b1[failed] <- NA_real_
  }
  
  s3 <- rowMeans(x) + N2_mat[, 1] / (ep2 * n_val)
  s4 <- rowMeans(x * x) + N2_mat[, 2] / (ep2 * n_val)
  
  cbind(b0, b1, s3, s4)
}

# ------------------------------------------------------------------------------
# 9. Efficient Direction Construction
# ------------------------------------------------------------------------------

efficient_direction <- function(
    theta,
    seed_a,
    U_a,
    N1_a,
    N2_a,
    ep,
    n_val,
    cache
) {
  theta_dir <- if (USE_DIRECTION_GRID) {
    round(theta / DIR_GRID) * DIR_GRID
  } else {
    theta
  }
  
  key <- paste(sprintf("%.15g", theta_dir), collapse = "_")
  
  if (exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  
  bad_direction <- function() {
    ans <- rep(NA_real_, M_DIM)
    assign(key, ans, envir = cache)
    ans
  }
  
  if (length(theta_dir) != 4L || !all(is.finite(theta_dir))) {
    return(bad_direction())
  }
  
  a <- exp(theta_dir[3])
  b <- exp(theta_dir[4])
  
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    return(bad_direction())
  }
  
  x_base <- matrix(
    qbeta(seed_a, a, b),
    nrow = nrow(seed_a),
    ncol = n_val
  )
  
  if (!all(is.finite(x_base))) {
    return(bad_direction())
  }
  
  aux <- sdp_matrix(
    seed_a, U_a, N1_a, N2_a, ep, theta_dir, n_val, x_pre = x_base
  )
  
  if (!all(is.finite(aux))) {
    return(bad_direction())
  }
  
  m0 <- colMeans(aux)
  J  <- matrix(NA_real_, nrow = M_DIM, ncol = 4L)
  
  for (j in 1:4) {
    tp <- theta_dir
    tp[j] <- tp[j] + FD_STEP
    xp <- if (j <= 2L) x_base else NULL
    
    Sp <- sdp_matrix(
      seed_a, U_a, N1_a, N2_a, ep, tp, n_val, x_pre = xp
    )
    
    if (!all(is.finite(Sp))) {
      return(bad_direction())
    }
    
    J[, j] <- (colMeans(Sp) - m0) / FD_STEP
  }
  
  Za <- sweep(aux, 2L, m0, "-")
  Sa <- crossprod(Za) / nrow(aux)
  Sa <- 0.5 * (Sa + t(Sa))
  
  Sai <- safe_inv(Sa, AUX_COV_FLOOR)
  if (is.null(Sai)) {
    return(bad_direction())
  }
  
  Jb <- J[, 1, drop = FALSE]
  Je <- J[, 2:4, drop = FALSE]
  
  Me <- t(Je) %*% Sai %*% Je
  Mei <- safe_inv(Me, NUIS_FLOOR)
  if (is.null(Mei)) {
    return(bad_direction())
  }
  
  nuisance_coef <- Mei %*% t(Je) %*% Sai %*% Jb
  efficient_derivative <- Jb - Je %*% nuisance_coef
  direction <- as.numeric(Sai %*% efficient_derivative)
  
  if (any(!is.finite(direction)) || sum(direction^2) < 1e-20) {
    return(bad_direction())
  }
  
  assign(key, direction, envir = cache)
  direction
}

# ------------------------------------------------------------------------------
# 10. Operator II Criterion Scoring
# ------------------------------------------------------------------------------

score_theta <- function(
    theta,
    seed,
    U,
    N1,
    N2,
    s_dp,
    ep,
    n_val,
    seed_a,
    U_a,
    N1_a,
    N2_a,
    cache,
    trace = NULL
) {
  synthetic <- sdp_matrix(seed, U, N1, N2, ep, theta, n_val)
  if (!all(is.finite(synthetic))) stat_fail("synthetic cloud")
  if (!all(is.finite(s_dp)))      stat_fail("observed summary")
  
  cloud <- rbind(synthetic, s_dp)
  
  direction <- efficient_direction(
    theta, seed_a, U_a, N1_a, N2_a, ep, n_val, cache
  )
  if (!all(is.finite(direction))) stat_fail("efficient direction")
  
  mu <- colMeans(cloud)
  Z  <- sweep(cloud, 2L, mu, "-")
  
  Sigma_hat <- crossprod(Z) / nrow(cloud)
  Sigma_hat <- 0.5 * (Sigma_hat + t(Sigma_hat))
  
  Sigma_inv <- safe_inv(Sigma_hat, INF_COV_FLOOR)
  if (is.null(Sigma_inv)) stat_fail("inference covariance")
  
  Q_full <- rowSums((Z %*% Sigma_inv) * Z)
  if (!all(is.finite(Q_full))) stat_fail("full Mahalanobis")
  
  denominator <- as.numeric(t(direction) %*% Sigma_hat %*% direction)
  if (!is.finite(denominator) || denominator <= DENOM_TOL) {
    stat_fail("efficient denominator")
  }
  
  projection <- as.numeric(Z %*% direction)
  Q_eff <- projection^2 / denominator
  Q_pen <- Q_eff + lambda_n * Q_full
  if (!all(is.finite(Q_pen))) stat_fail("penalized statistic")
  
  depth <- 1 / (1 + Q_pen)
  if (!all(is.finite(depth))) stat_fail("depth")
  
  observed_index <- R + 1L
  observed_rank  <- rank(depth, ties.method = "max")[observed_index]
  M <- observed_rank + depth[observed_index]
  
  if (!is.finite(M)) stat_fail("M score")
  
  if (!is.null(trace)) {
    if (!is.finite(trace$M) || M > trace$M) {
      trace$M <- M
      trace$par <- theta
    }
    if (M >= rank_threshold) {
      cert_hit(theta, M)
    }
  }
  
  -M
}

# Parameter initialization anchor based on moments
anchor_theta <- function(s_dp) {
  beta0_anchor <- max(min(s_dp[1], 10), -10)
  beta1_anchor <- max(min(s_dp[2], 10), -10)
  
  log_shape1 <- 0
  log_shape2 <- 0
  
  m1 <- s_dp[3]
  m2 <- s_dp[4]
  v  <- m2 - m1^2
  
  if (
    is.finite(m1) && is.finite(v) &&
    m1 > 1e-6 && m1 < 1 - 1e-6 && v > 1e-9
  ) {
    total_shape <- (m1 * (1 - m1) / v - 1)
    
    if (is.finite(total_shape) && total_shape > 1e-6) {
      a <- m1 * total_shape
      b <- (1 - m1) * total_shape
      
      if (a > 1e-6 && b > 1e-6) {
        log_shape1 <- max(min(log(a), 5), -5)
        log_shape2 <- max(min(log(b), 5), -5)
      }
    }
  }
  
  c(beta1_anchor, beta0_anchor, log_shape1, log_shape2)
}

# ------------------------------------------------------------------------------
# 11. Parameter Box Acceptance Search
# ------------------------------------------------------------------------------

accept_parameter_box <- function(
    beta_left,
    beta_right,
    start_theta,
    anchor,
    seed,
    U,
    N1,
    N2,
    s_dp,
    ep,
    n_val,
    seed_a,
    U_a,
    N1_a,
    N2_a,
    cache
) {
  if (
    !is.finite(beta_left) ||
    !is.finite(beta_right) ||
    beta_left > beta_right
  ) {
    return(
      list(status = "unresolved", par = rep(NA_real_, 4L), M = NA_real_)
    )
  }
  
  lower_box <- c(beta_left, NUIS_LO)
  upper_box <- c(beta_right, NUIS_HI)
  span_box  <- upper_box - lower_box
  span_safe <- pmax(span_box, 1e-12)
  
  clip_theta <- function(th) {
    pmin(pmax(th, lower_box), upper_box)
  }
  
  beta_mid <- (beta_left + beta_right) / 2
  starts <- rbind(
    clip_theta(start_theta),
    clip_theta(anchor),
    c(beta_mid, anchor[2:4]),
    c(beta_mid, 0, 0, 0)
  )
  
  starts <- unique(round(starts, 12))
  if (is.null(dim(starts))) starts <- matrix(starts, nrow = 1L)
  
  trace     <- new.env(parent = emptyenv())
  trace$M   <- NA_real_
  trace$par <- rep(NA_real_, 4L)
  
  n_valid   <- 0L
  n_failed  <- 0L
  n_optimizer_completed <- 0L
  
  eval_score <- function(theta_candidate) {
    ans <- tryCatch(
      score_theta(
        theta_candidate, seed, U, N1, N2, s_dp, ep, n_val,
        seed_a, U_a, N1_a, N2_a, cache, trace
      ),
      stat_fail = function(e) {
        n_failed <<- n_failed + 1L
        NA_real_
      }
    )
    if (is.finite(ans)) n_valid <<- n_valid + 1L
    ans
  }
  
  # Step 1: Direct evaluations
  for (ii in seq_len(nrow(starts))) {
    direct <- tryCatch(
      eval_score(starts[ii, , drop = TRUE]),
      cert_hit = function(e) e
    )
    
    if (inherits(direct, "cert_hit")) {
      return(list(status = "accept", par = direct$par, M = direct$M))
    }
    
    if (is.finite(direct) && (-direct) >= rank_threshold) {
      return(list(status = "accept", par = starts[ii, , drop = TRUE], M = -direct))
    }
  }
  
  # Step 2: Multi-start Nelder-Mead optimization
  OUT_OF_BOX <- 1e6
  objective_unit_cube <- function(z) {
    if (!all(is.finite(z)) || any(z < 0) || any(z > 1)) return(OUT_OF_BOX)
    
    theta_candidate <- lower_box + z * span_safe
    deg <- span_box <= 0
    if (any(deg)) theta_candidate[deg] <- lower_box[deg]
    
    val <- eval_score(theta_candidate)
    if (!is.finite(val)) return(OUT_OF_BOX)
    val
  }
  
  for (ii in seq_len(nrow(starts))) {
    z_start <- (starts[ii, ] - lower_box) / span_safe
    z_start <- pmin(pmax(z_start, 0), 1)
    
    result <- tryCatch(
      list(
        type = "optimizer",
        opt = optim(
          par = z_start,
          fn = objective_unit_cube,
          method = "Nelder-Mead",
          control = list(maxit = BOX_NM_MAXIT, reltol = 1e-7)
        )
      ),
      cert_hit = function(e) list(type = "certificate", par = e$par, M = e$M),
      error = function(e) list(type = "error", opt = NULL)
    )
    
    if (identical(result$type, "certificate")) {
      return(list(status = "accept", par = result$par, M = result$M))
    }
    
    if (
      identical(result$type, "optimizer") &&
      !is.null(result$opt) &&
      is.finite(result$opt$value) &&
      result$opt$value < OUT_OF_BOX
    ) {
      n_optimizer_completed <- n_optimizer_completed + 1L
      z_best <- pmin(pmax(result$opt$par, 0), 1)
      theta_best <- lower_box + z_best * span_safe
      deg <- span_box <= 0
      if (any(deg)) theta_best[deg] <- lower_box[deg]
      
      final <- tryCatch(eval_score(theta_best), cert_hit = function(e) e)
      if (inherits(final, "cert_hit")) {
        return(list(status = "accept", par = final$par, M = final$M))
      }
      if (is.finite(final) && (-final) >= rank_threshold) {
        return(list(status = "accept", par = theta_best, M = -final))
      }
    }
  }
  
  if (is.finite(trace$M) && trace$M >= rank_threshold) {
    return(list(status = "accept", par = trace$par, M = trace$M))
  }
  
  if (n_valid > 0L && n_optimizer_completed > 0L) {
    return(
      list(
        status = "reject",
        par = rep(NA_real_, 4L),
        M = if (is.finite(trace$M)) trace$M else NA_real_
      )
    )
  }
  
  list(
    status = "unresolved",
    par = rep(NA_real_, 4L),
    M = if (is.finite(trace$M)) trace$M else NA_real_
  )
}

# ------------------------------------------------------------------------------
# 12. Confidence Interval Construction by Bisection Search
# ------------------------------------------------------------------------------

getCI <- function(
    seed,
    U,
    N1,
    N2,
    s_dp,
    ep,
    n_val,
    seed_a,
    U_a,
    N1_a,
    N2_a,
    cache
) {
  anchor <- anchor_theta(s_dp)
  if (!all(is.finite(anchor))) {
    return(list(ci = c(NA_real_, NA_real_), status = "unresolved"))
  }
  
  initial <- accept_parameter_box(
    beta_left = BETA1_LO, beta_right = BETA1_HI, start_theta = anchor,
    anchor = anchor, seed = seed, U = U, N1 = N1, N2 = N2, s_dp = s_dp,
    ep = ep, n_val = n_val, seed_a = seed_a, U_a = U_a, N1_a = N1_a,
    N2_a = N2_a, cache = cache
  )
  
  if (initial$status == "unresolved") {
    return(list(ci = c(NA_real_, NA_real_), status = "unresolved"))
  }
  if (initial$status != "accept") {
    return(list(ci = c(NA_real_, NA_real_), status = "no_accept"))
  }
  
  accepted_theta <- initial$par
  centre <- accepted_theta[1]
  
  # Search left boundary
  left_lower <- BETA1_LO
  left_upper <- centre
  left_start <- accepted_theta
  
  while (left_upper - left_lower > CI_TOL) {
    left_middle <- (left_lower + left_upper) / 2
    res <- accept_parameter_box(
      beta_left = left_lower, beta_right = left_middle, start_theta = left_start,
      anchor = anchor, seed = seed, U = U, N1 = N1, N2 = N2, s_dp = s_dp,
      ep = ep, n_val = n_val, seed_a = seed_a, U_a = U_a, N1_a = N1_a,
      N2_a = N2_a, cache = cache
    )
    
    if (res$status == "unresolved") {
      return(list(ci = c(NA_real_, NA_real_), status = "unresolved"))
    }
    
    if (res$status == "accept") {
      left_upper <- left_middle
      left_start <- res$par
    } else {
      left_lower <- left_middle
    }
  }
  
  # Search right boundary
  right_lower <- centre
  right_upper <- BETA1_HI
  right_start <- accepted_theta
  
  while (right_upper - right_lower > CI_TOL) {
    right_middle <- (right_lower + right_upper) / 2
    res <- accept_parameter_box(
      beta_left = right_middle, beta_right = right_upper, start_theta = right_start,
      anchor = anchor, seed = seed, U = U, N1 = N1, N2 = N2, s_dp = s_dp,
      ep = ep, n_val = n_val, seed_a = seed_a, U_a = U_a, N1_a = N1_a,
      N2_a = N2_a, cache = cache
    )
    
    if (res$status == "unresolved") {
      return(list(ci = c(NA_real_, NA_real_), status = "unresolved"))
    }
    
    if (res$status == "accept") {
      right_lower <- right_middle
      right_start <- res$par
    } else {
      right_upper <- right_middle
    }
  }
  
  list(ci = c(left_lower, right_upper), status = "ok")
}

# ------------------------------------------------------------------------------
# 13. Monte Carlo Grid Execution
# ------------------------------------------------------------------------------

ALL_SUMMARIES <- list()
config_counter <- 0L

for (n_cur in n_values) {
  for (eps_cur in eps_values) {
    config_counter <- config_counter + 1L
    
    n  <- as.integer(n_cur)
    ep <- as.numeric(eps_cur)
    
    theta_true <- c(beta1_true, beta0, log(shape1), log(shape2))
    lambda_n   <- 1 / log(n)
    DIR_GRID   <- sqrt(lambda_n) / 10
    
    ep_label   <- format(ep, scientific = FALSE, trim = TRUE)
    RUN_TAG    <- sprintf("n_%d_eps_%s_eff", n, ep_label)
    
    RESULTS_ROOT <- file.path(PROJECT_DIR, RUN_TAG)
    DIR_REPS     <- file.path(RESULTS_ROOT, "logistic")
    DIR_SUMMARY  <- file.path(RESULTS_ROOT, "summary")
    
    dir.create(DIR_REPS, recursive = TRUE, showWarnings = FALSE)
    dir.create(DIR_SUMMARY, recursive = TRUE, showWarnings = FALSE)
    
    rep_file    <- function(k) file.path(DIR_REPS, sprintf("rep_%04d.csv", k))
    FILE_CIS    <- file.path(DIR_SUMMARY, "summary_CIs.csv")
    FILE_COVER  <- file.path(DIR_SUMMARY, "Repro_logistic.csv")
    FILE_CONFIG <- file.path(RESULTS_ROOT, "config.txt")
    
    existing_reps <- list.files(DIR_REPS, pattern = "^rep_[0-9]+\\.csv$", full.names = TRUE)
    if (FRESH_START) {
      if (length(existing_reps) > 0L) invisible(file.remove(existing_reps))
      if (file.exists(FILE_CIS))   invisible(file.remove(FILE_CIS))
      if (file.exists(FILE_COVER)) invisible(file.remove(FILE_COVER))
    }
    
    writeLines(
      c(
        "objperturb -- corrected Efficient Operator II",
        "",
        sprintf("n                 : %d", n),
        sprintf("epsilon           : %s", ep_label),
        sprintf("reps              : %d", reps),
        sprintf("R                 : %d", R),
        sprintf("R_aux             : %d", R_aux),
        sprintf("alpha             : %.4f", alpha),
        sprintf("lambda_n          : %.10f", lambda_n),
        sprintf("DIR_GRID          : %.10f", DIR_GRID),
        sprintf("USE_DIRECTION_GRID: %s", USE_DIRECTION_GRID),
        sprintf("FD_STEP           : %.10f", FD_STEP),
        sprintf("CI_TOL            : %.10f", CI_TOL),
        sprintf("BOX_NM_MAXIT      : %d", BOX_NM_MAXIT),
        sprintf("seed rule         : rep + 1000"),
        sprintf("workers           : %d", n_workers),
        sprintf("generated         : %s", format(Sys.time()))
      ),
      FILE_CONFIG
    )
    
    cat("\n\n################################################################\n")
    cat(sprintf("CONFIGURATION %d / 16: n = %d, epsilon = %s\n", config_counter, n, ep_label))
    cat(sprintf("Folder: %s\n", RESULTS_ROOT))
    cat("################################################################\n\n")
    
    cl <- parallel::makePSOCKcluster(n_workers)
    doSNOW::registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = reps, style = 3)
    snow_options <- list(progress = function(k) setTxtProgressBar(pb, k))
    
    parallel::clusterExport(
      cl,
      c(
        "lower", "upper", "draw_U2", "safe_inv", "softplus", "robust_objperturb_row",
        "stat_fail", "cert_hit", "sdp_matrix", "efficient_direction", "score_theta",
        "anchor_theta", "accept_parameter_box", "getCI", "R", "R_aux", "alpha", "ep",
        "n", "reps", "theta_true", "beta1_true", "M_DIM", "lambda_n", "DIR_GRID",
        "USE_DIRECTION_GRID", "FD_STEP", "AUX_COV_FLOOR", "NUIS_FLOOR", "INF_COV_FLOOR",
        "DENOM_TOL", "NEWTON_MAXIT", "NEWTON_TOL", "NEWTON_ACCEPT", "MAX_STEP", "DET_FLOOR",
        "ROBUST_BFGS_MAXIT", "ROBUST_NLMINB_MAXIT", "ROBUST_GRAD_TOL", "BOX_NM_MAXIT",
        "BETA1_LO", "BETA1_HI", "NUIS_LO", "NUIS_HI", "CI_TOL", "rank_threshold",
        "rank_reference", "rep_file"
      ),
      envir = environment()
    )
    
    CIs <- foreach(
      rep = seq_len(reps),
      .combine = "rbind",
      .packages = "tictoc",
      .options.snow = snow_options
    ) %dopar% {
      set.seed(rep + 1000L)
      result_file <- rep_file(rep)
      tic()
      
      # Observed data generation
      seed_obs <- matrix(runif(n), nrow = 1L, ncol = n)
      U_obs    <- matrix(runif(n), nrow = 1L, ncol = n)
      U1_obs   <- matrix(runif(2L, min = -1, max = 1), nrow = 1L, ncol = 2L)
      G1_obs   <- rgamma(1L, shape = 3, rate = 1)
      N1_obs   <- U1_obs * G1_obs
      U2_obs   <- draw_U2(1L)
      G2_obs   <- rgamma(1L, shape = 3, rate = 1)
      N2_obs   <- U2_obs * G2_obs
      
      observed_matrix <- sdp_matrix(
        seed_obs, U_obs, N1_obs, N2_obs, ep, theta_true, n
      )
      s_dp <- as.numeric(observed_matrix[1L, ])
      
      # Primary synthetic cloud
      seed_mat <- matrix(runif(n * R), nrow = R, ncol = n)
      U_mat    <- matrix(runif(n * R), nrow = R, ncol = n)
      U1_mat   <- matrix(runif(2L * R, min = -1, max = 1), nrow = R, ncol = 2L)
      G1_mat   <- rgamma(R, shape = 3, rate = 1)
      N1_mat   <- U1_mat * G1_mat
      U2_mat   <- draw_U2(R)
      G2_mat   <- rgamma(R, shape = 3, rate = 1)
      N2_mat   <- U2_mat * G2_mat
      
      # Auxiliary direction cloud
      seed_aux <- matrix(runif(n * R_aux), nrow = R_aux, ncol = n)
      U_aux    <- matrix(runif(n * R_aux), nrow = R_aux, ncol = n)
      U1_aux   <- matrix(runif(2L * R_aux, min = -1, max = 1), nrow = R_aux, ncol = 2L)
      G1_aux   <- rgamma(R_aux, shape = 3, rate = 1)
      N1_aux   <- U1_aux * G1_aux
      U2_aux   <- draw_U2(R_aux)
      G2_aux   <- rgamma(R_aux, shape = 3, rate = 1)
      N2_aux   <- U2_aux * G2_aux
      
      direction_cache <- new.env(parent = emptyenv())
      
      ci_result <- if (!all(is.finite(s_dp))) {
        list(ci = c(NA_real_, NA_real_), status = "unresolved")
      } else {
        tryCatch(
          getCI(
            seed_mat, U_mat, N1_mat, N2_mat, s_dp, ep, n,
            seed_aux, U_aux, N1_aux, N2_aux, direction_cache
          ),
          error = function(e) list(ci = c(NA_real_, NA_real_), status = "unresolved")
        )
      }
      
      timing <- toc(quiet = TRUE)
      elapsed_seconds <- as.numeric(timing$toc - timing$tic)
      
      out_row <- data.frame(
        rep = rep,
        lower = ci_result$ci[1],
        upper = ci_result$ci[2],
        seconds = elapsed_seconds,
        status = ci_result$status,
        stringsAsFactors = FALSE
      )
      
      write.csv(out_row, result_file, row.names = FALSE)
      out_row
    }
    
    parallel::stopCluster(cl)
    try(close(pb), silent = TRUE)
    
    CIs_df <- as.data.frame(CIs, stringsAsFactors = FALSE)
    CIs_df$rep     <- as.integer(CIs_df$rep)
    CIs_df$lower   <- as.numeric(CIs_df$lower)
    CIs_df$upper   <- as.numeric(CIs_df$upper)
    CIs_df$seconds <- as.numeric(CIs_df$seconds)
    CIs_df$status  <- as.character(CIs_df$status)
    
    write.csv(CIs_df, FILE_CIS, row.names = FALSE)
    
    # Coverage calculation
    resolved   <- (CIs_df$status == "ok" & is.finite(CIs_df$lower) & is.finite(CIs_df$upper))
    unresolved <- CIs_df$status == "unresolved"
    no_accept  <- CIs_df$status == "no_accept"
    
    n_resolved   <- sum(resolved)
    n_unresolved <- sum(unresolved)
    n_no_accept  <- sum(no_accept)
    
    covered_resolved   <- (resolved & CIs_df$lower <= beta1_true & CIs_df$upper >= beta1_true)
    n_covered_resolved <- sum(covered_resolved, na.rm = TRUE)
    coverage_resolved  <- if (n_resolved > 0L) n_covered_resolved / n_resolved else NA_real_
    
    coverage_lower_all <- n_covered_resolved / reps
    coverage_upper_all <- (n_covered_resolved + n_unresolved + n_no_accept) / reps
    
    widths <- ifelse(resolved, CIs_df$upper - CIs_df$lower, NA_real_)
    average_width <- if (any(is.finite(widths))) mean(widths, na.rm = TRUE) else NA_real_
    width_se <- if (sum(is.finite(widths)) >= 2L) {
      sd(widths, na.rm = TRUE) / sqrt(sum(is.finite(widths)))
    } else {
      NA_real_
    }
    
    mean_seconds <- mean(CIs_df$seconds, na.rm = TRUE)
    
    if (n_resolved > 0L) {
      coverage_lo <- if (n_covered_resolved == 0L) 0 else qbeta(0.025, n_covered_resolved, n_resolved - n_covered_resolved + 1L)
      coverage_hi <- if (n_covered_resolved == n_resolved) 1 else qbeta(0.975, n_covered_resolved + 1L, n_resolved - n_covered_resolved)
      p_under <- pbinom(n_covered_resolved, n_resolved, rank_reference)
    } else {
      coverage_lo <- NA_real_
      coverage_hi <- NA_real_
      p_under     <- NA_real_
    }
    
    summary_table <- data.frame(
      n                   = n,
      epsilon             = ep,
      reps_requested      = reps,
      reps_produced       = n_resolved,
      R                   = R,
      R_aux               = R_aux,
      alpha               = alpha,
      rank_cutoff         = rank_threshold,
      finite_R_reference  = rank_reference,
      covered             = n_covered_resolved,
      coverage            = coverage_resolved,
      coverage_resolved   = coverage_resolved,
      coverage_lo         = coverage_lo,
      coverage_hi         = coverage_hi,
      coverage_lower_all  = coverage_lower_all,
      coverage_upper_all  = coverage_upper_all,
      p_under             = p_under,
      width_all           = average_width,
      n_width_all         = n_resolved,
      width_real          = average_width,
      n_width_real        = n_resolved,
      width_se            = width_se,
      n_no_accepted_point = n_no_accept,
      n_full_range        = 0L,
      n_unresolved        = n_unresolved,
      failure_rate        = (n_unresolved + n_no_accept) / reps,
      mean_seconds        = mean_seconds,
      DIR_GRID            = DIR_GRID,
      USE_DIRECTION_GRID  = USE_DIRECTION_GRID,
      FD_STEP             = FD_STEP,
      CI_TOL              = CI_TOL,
      BOX_NM_MAXIT        = BOX_NM_MAXIT,
      robust_grad_tol     = ROBUST_GRAD_TOL,
      slurm_cpus          = slurm_cpus,
      workers             = n_workers
    )
    
    write.csv(summary_table, FILE_COVER, row.names = FALSE)
    ALL_SUMMARIES[[length(ALL_SUMMARIES) + 1L]] <- summary_table
    
    cat("\n============================================================\n")
    cat(" CORRECTED EFFICIENT OPERATOR II\n")
    cat("============================================================\n")
    cat(sprintf("n                    : %d\n", n))
    cat(sprintf("epsilon              : %.3f\n", ep))
    cat(sprintf("Requested reps       : %d\n", reps))
    cat(sprintf("Resolved CIs         : %d\n", n_resolved))
    cat(sprintf("Unresolved numerical : %d\n", n_unresolved))
    cat(sprintf("No accepted point    : %d\n", n_no_accept))
    cat(sprintf("Failure rate         : %.4f\n", (n_unresolved + n_no_accept) / reps))
    cat(sprintf("Coverage (resolved)  : %.4f (%d/%d)\n", coverage_resolved, n_covered_resolved, n_resolved))
    cat(sprintf("95%% exact CI         : (%.4f, %.4f)\n", coverage_lo, coverage_hi))
    cat(sprintf("Mean resolved width  : %.4f\n", average_width))
    cat(sprintf("Results folder       : %s\n", RESULTS_ROOT))
    cat("============================================================\n")
  }
}

# ------------------------------------------------------------------------------
# 14. Save Master Summary File
# ------------------------------------------------------------------------------

MASTER_SUMMARY <- do.call(rbind, ALL_SUMMARIES)
rownames(MASTER_SUMMARY) <- NULL
MASTER_SUMMARY_FILE <- file.path(PROJECT_DIR, "all_16_results.csv")
write.csv(MASTER_SUMMARY, MASTER_SUMMARY_FILE, row.names = FALSE)

cat("\n\n################################################################\n")
cat(" ALL 16 CONFIGURATIONS FINISHED\n")
cat("################################################################\n")
cat(sprintf("Master CSV: %s\n\n", MASTER_SUMMARY_FILE))

print(
  MASTER_SUMMARY[
    ,
    c(
      "n", "epsilon", "coverage_resolved", "coverage_lower_all",
      "coverage_upper_all", "width_real", "n_unresolved",
      "n_no_accepted_point", "failure_rate", "mean_seconds"
    )
  ]
)

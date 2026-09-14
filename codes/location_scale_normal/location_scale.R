# ==============================================================================
# Script: location_scale.R
# Description: Single-draw confidence-region construction and marginal CI
#              evaluation for the clamped location-scale normal model.
#              Compares Mahalanobis depth, Efficient Repro (mu-interest),
#              Efficient Repro (sigma-interest), and PB-ADI methods.
# ==============================================================================

PROJECT_DIR <- path.expand("~/R_Simuls/Location_Scale_nomal")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

RESULTS_DIR <- file.path(PROJECT_DIR, "results_apply")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Data Generating Process (DGP) & Differentially Private Summaries
# ------------------------------------------------------------------------------

upper_clamp <- 3
lower_clamp <- 0
n <- 100
ep <- 1

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
    row_var + var_noise_scale * privacy_noises[, 2]
  )
}

# ------------------------------------------------------------------------------
# 2. Simulation & Inference Settings
# ------------------------------------------------------------------------------

value_r <- 100
R_synthetic <- 200
alpha <- 0.05
tol <- 1e-8

population_mu <- 1
population_sigma <- 1

R_aux <- 400
R_indirect_est <- 50
B_paramboot <- 200
h_fd <- 1e-3
lambda_n <- 1 / log(n)
pb_minimum_valid <- max(20L, floor(0.8 * B_paramboot))

acceptance_threshold <- floor(alpha * (R_synthetic + 1)) + 1

# ------------------------------------------------------------------------------
# 3. Environment & Parallel Backend Setup
# ------------------------------------------------------------------------------

list.of.packages <- c("foreach", "doSNOW", "parallelly")
new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]

if (length(new.packages) > 0) {
  install.packages(new.packages, dep = TRUE)
}

for (package.i in list.of.packages) {
  suppressPackageStartupMessages(
    library(package.i, character.only = TRUE)
  )
}

# Prevent BLAS/MKL thread oversubscription during parallel execution
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

n.cores <- min(124, parallel::detectCores() - 1)
cl <- parallel::makePSOCKcluster(n.cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = value_r^2, style = 3)
progress <- function(k) setTxtProgressBar(pb, k)
opts <- list(progress = progress)

# ------------------------------------------------------------------------------
# 4. Matrix Utilities
# ------------------------------------------------------------------------------

safe_inv2 <- function(A, ridge = 1e-10) {
  if (
    is.null(A) ||
    !all(dim(A) == c(2L, 2L)) ||
    any(!is.finite(A))
  ) {
    stop("Invalid matrix in safe_inv2().")
  }
  
  A <- 0.5 * (A + t(A)) + ridge * diag(2L)
  
  a <- A[1L, 1L]
  b <- A[1L, 2L]
  d <- A[2L, 2L]
  
  det_A <- a * d - b^2
  
  if (!is.finite(det_A) || det_A <= ridge^2) {
    eg <- eigen(A, symmetric = TRUE)
    
    if (
      any(!is.finite(eg$values)) ||
      any(!is.finite(eg$vectors))
    ) {
      stop("Invalid eigendecomposition.")
    }
    
    ev <- pmax(eg$values, ridge)
    
    return(
      eg$vectors %*%
        diag(1 / ev, nrow = 2L) %*%
        t(eg$vectors)
    )
  }
  
  matrix(
    c(d, -b, -b, a),
    nrow = 2L
  ) / det_A
}

stabilize_cov2 <- function(A, ridge = 1e-8) {
  if (
    is.null(A) ||
    !all(dim(A) == c(2L, 2L)) ||
    any(!is.finite(A))
  ) {
    return(ridge * diag(2L))
  }
  
  0.5 * (A + t(A)) + ridge * diag(2L)
}

# ------------------------------------------------------------------------------
# 5. Fixed Auxiliary Draws & Caches
# ------------------------------------------------------------------------------

set.seed(300000)

aux_dr <- matrix(
  rnorm(n * R_aux),
  nrow = R_aux,
  ncol = n
)

aux_pn <- matrix(
  rnorm(2 * R_aux),
  nrow = R_aux,
  ncol = 2
)

eff_cache_mu <- new.env(parent = emptyenv())
eff_cache_sigma <- new.env(parent = emptyenv())

# ------------------------------------------------------------------------------
# 6. Efficient Score & Direction Calculations
# ------------------------------------------------------------------------------

compute_aux_derivatives <- function(mu, sa) {
  if (!is.finite(mu) || !is.finite(sa) || sa <= tol) {
    return(NULL)
  }
  
  h <- h_fd
  
  synth_center <- tryCatch(
    sdp_vec(aux_dr, aux_pn, sa, mu),
    error = function(e) NULL
  )
  
  if (
    is.null(synth_center) ||
    any(!is.finite(synth_center))
  ) {
    return(NULL)
  }
  
  center_mean <- colMeans(synth_center)
  
  if (sa - h > tol) {
    synth_sa_plus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa + h,
      mu
    )
    
    synth_sa_minus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa - h,
      mu
    )
    
    d_sa <- (
      colMeans(synth_sa_plus) -
        colMeans(synth_sa_minus)
    ) / (2 * h)
    
  } else {
    synth_sa_plus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa + h,
      mu
    )
    
    d_sa <- (
      colMeans(synth_sa_plus) -
        center_mean
    ) / h
  }
  
  synth_mu_plus <- sdp_vec(
    aux_dr,
    aux_pn,
    sa,
    mu + h
  )
  
  synth_mu_minus <- sdp_vec(
    aux_dr,
    aux_pn,
    sa,
    mu - h
  )
  
  d_mu <- (
    colMeans(synth_mu_plus) -
      colMeans(synth_mu_minus)
  ) / (2 * h)
  
  if (any(!is.finite(c(d_mu, d_sa)))) {
    return(NULL)
  }
  
  aux_cov <- cov(synth_center)
  aux_cov <- 0.5 * (aux_cov + t(aux_cov)) +
    1e-8 * diag(2L)
  
  cov_inv <- tryCatch(
    safe_inv2(aux_cov, 1e-8),
    error = function(e) NULL
  )
  
  if (is.null(cov_inv)) {
    return(NULL)
  }
  
  list(
    d_mu = matrix(d_mu, ncol = 1L),
    d_sa = matrix(d_sa, ncol = 1L),
    cov_inv = cov_inv
  )
}

efficient_direction <- function(mu, sa, interest) {
  cache <- if (interest == "mu") {
    eff_cache_mu
  } else {
    eff_cache_sigma
  }
  
  key <- sprintf(
    "%.10g_%.10g",
    mu,
    sa
  )
  
  if (
    exists(
      key,
      envir = cache,
      inherits = FALSE
    )
  ) {
    return(
      get(
        key,
        envir = cache,
        inherits = FALSE
      )
    )
  }
  
  deriv <- compute_aux_derivatives(
    mu,
    sa
  )
  
  if (is.null(deriv)) {
    direction <- c(
      NA_real_,
      NA_real_
    )
    
    assign(
      key,
      direction,
      envir = cache
    )
    
    return(direction)
  }
  
  if (interest == "mu") {
    nuisance <- deriv$d_sa
    target_raw <- deriv$d_mu
  } else {
    nuisance <- deriv$d_mu
    target_raw <- deriv$d_sa
  }
  
  nuisance_norm <- as.numeric(
    t(nuisance) %*%
      deriv$cov_inv %*%
      nuisance
  )
  
  if (
    !is.finite(nuisance_norm) ||
    nuisance_norm <= 1e-10
  ) {
    direction <- c(
      NA_real_,
      NA_real_
    )
    
    assign(
      key,
      direction,
      envir = cache
    )
    
    return(direction)
  }
  
  projection_coefficient <- as.numeric(
    t(nuisance) %*%
      deriv$cov_inv %*%
      target_raw
  ) / nuisance_norm
  
  target <- target_raw -
    projection_coefficient * nuisance
  
  direction <- as.numeric(
    deriv$cov_inv %*%
      target
  )
  
  if (
    any(!is.finite(direction)) ||
    sum(direction^2) < 1e-20
  ) {
    direction <- c(
      NA_real_,
      NA_real_
    )
  }
  
  assign(
    key,
    direction,
    envir = cache
  )
  
  direction
}

efficient_depth <- function(
    synth,
    mu,
    sa,
    interest
) {
  direction <- efficient_direction(
    mu,
    sa,
    interest
  )
  
  if (any(!is.finite(direction))) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  center <- colMeans(synth)
  
  centered <- sweep(
    synth,
    2L,
    center,
    "-"
  )
  
  covariance <- crossprod(centered) /
    nrow(synth)
  
  covariance <- 0.5 * (
    covariance +
      t(covariance)
  ) +
    1e-10 * diag(2L)
  
  cov_inv <- tryCatch(
    safe_inv2(
      covariance,
      1e-10
    ),
    error = function(e) NULL
  )
  
  if (is.null(cov_inv)) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  eff_scale <- as.numeric(
    t(direction) %*%
      covariance %*%
      direction
  )
  
  if (
    !is.finite(eff_scale) ||
    eff_scale <= 1e-14
  ) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  projection <- as.numeric(
    centered %*%
      direction
  )
  
  efficient_quadratic <- projection^2 /
    eff_scale
  
  mahalanobis_quadratic <- rowSums(
    (
      centered %*%
        cov_inv
    ) *
      centered
  )
  
  penalty <- efficient_quadratic +
    lambda_n *
    mahalanobis_quadratic
  
  depth <- 1 / (
    1 +
      penalty
  )
  
  if (any(!is.finite(depth))) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  depth
}

# ------------------------------------------------------------------------------
# 7. Mahalanobis Depth Functions
# ------------------------------------------------------------------------------

depth_mahalanobis_selfref <- function(synth) {
  m <- colMeans(synth)
  S <- cov(synth)
  
  S_inv <- tryCatch(
    solve(S),
    error = function(e) {
      tryCatch(
        safe_inv2(
          S,
          1e-10
        ),
        error = function(e2) NULL
      )
    }
  )
  
  if (is.null(S_inv)) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  centered <- sweep(
    synth,
    2,
    m,
    "-"
  )
  
  md2 <- rowSums(
    (
      centered %*%
        S_inv
    ) *
      centered
  )
  
  depth <- 1 / (
    1 +
      md2
  )
  
  if (any(!is.finite(depth))) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  depth
}

depth_mahalanobis_cross <- function(
    x,
    synth
) {
  m <- colMeans(synth)
  S <- cov(synth)
  
  S_inv <- tryCatch(
    solve(S),
    error = function(e) {
      tryCatch(
        safe_inv2(
          S,
          1e-10
        ),
        error = function(e2) NULL
      )
    }
  )
  
  if (is.null(S_inv)) {
    return(NA_real_)
  }
  
  centered <- sweep(
    matrix(
      x,
      nrow = 1
    ),
    2,
    m,
    "-"
  )
  
  md2 <- rowSums(
    (
      centered %*%
        S_inv
    ) *
      centered
  )
  
  depth <- 1 / (
    1 +
      md2
  )
  
  if (!is.finite(depth)) {
    return(NA_real_)
  }
  
  depth
}

# ------------------------------------------------------------------------------
# 8. Test Statistic Scoring & Confidence Set Inversion
# ------------------------------------------------------------------------------

score_mu_sa <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type
) {
  mu <- optim_par[1]
  sa <- optim_par[2]
  
  if (
    !is.finite(mu) ||
    !is.finite(sa) ||
    sa <= tol
  ) {
    return(NA_real_)
  }
  
  synth <- tryCatch(
    rbind(
      sdp_vec(
        data_randomness,
        privacy_noises,
        sa,
        mu
      ),
      dp_statistic
    ),
    error = function(e) NULL
  )
  
  if (is.null(synth)) {
    return(NA_real_)
  }
  
  if (
    depth_type %in%
    c(
      "efficient_mu",
      "efficient_sigma"
    )
  ) {
    interest <- if (
      depth_type ==
      "efficient_mu"
    ) {
      "mu"
    } else {
      "sigma"
    }
    
    D_synth <- efficient_depth(
      synth,
      mu,
      sa,
      interest
    )
    
  } else {
    D_synth <- depth_mahalanobis_selfref(
      synth
    )
  }
  
  if (any(!is.finite(D_synth))) {
    return(NA_real_)
  }
  
  obs_index <- R_synthetic + 1L
  
  r <- rank(
    D_synth,
    ties.method = "max"
  )[obs_index]
  
  depth_obs <- D_synth[
    obs_index
  ]
  
  if (
    !is.finite(r) ||
    !is.finite(depth_obs)
  ) {
    return(NA_real_)
  }
  
  -(
    r +
      depth_obs
  )
}

score_sa_mu <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type
) {
  sa <- optim_par[1]
  mu <- optim_par[2]
  
  if (
    !is.finite(mu) ||
    !is.finite(sa) ||
    sa <= tol
  ) {
    return(NA_real_)
  }
  
  score_mu_sa(
    c(mu, sa),
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type
  )
}

score_for_optim <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    score_func
) {
  ans <- score_func(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type
  )
  
  if (!is.finite(ans)) {
    return(1e12)
  }
  
  ans
}

make_box_starts <- function(
    optim_par,
    search_lower,
    search_upper,
    nuisance_lower,
    nuisance_upper
) {
  search_mid <- (
    search_lower +
      search_upper
  ) / 2
  
  nuisance_mid <- (
    nuisance_lower +
      nuisance_upper
  ) / 2
  
  s25 <- search_lower +
    0.25 * (
      search_upper -
        search_lower
    )
  
  s75 <- search_lower +
    0.75 * (
      search_upper -
        search_lower
    )
  
  n25 <- nuisance_lower +
    0.25 * (
      nuisance_upper -
        nuisance_lower
    )
  
  n75 <- nuisance_lower +
    0.75 * (
      nuisance_upper -
        nuisance_lower
    )
  
  supplied <- c(
    min(
      max(
        optim_par[1],
        search_lower
      ),
      search_upper
    ),
    min(
      max(
        optim_par[2],
        nuisance_lower
      ),
      nuisance_upper
    )
  )
  
  starts <- rbind(
    supplied,
    c(search_mid, nuisance_mid),
    c(s25, n25),
    c(s25, n75),
    c(s75, n25),
    c(s75, n75)
  )
  
  starts <- unique(
    round(
      starts,
      12
    )
  )
  
  if (is.null(dim(starts))) {
    starts <- matrix(
      starts,
      nrow = 1L
    )
  }
  
  starts
}

accept_box <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    search_lower,
    search_upper,
    nuisance_lower,
    nuisance_upper,
    depth_type,
    score_func
) {
  starts <- make_box_starts(
    optim_par,
    search_lower,
    search_upper,
    nuisance_lower,
    nuisance_upper
  )
  
  any_resolved <- FALSE
  any_failure <- FALSE
  
  # Step 1: Direct evaluations at candidate grid points
  for (k in seq_len(nrow(starts))) {
    start_k <- starts[
      k,
      ,
      drop = TRUE
    ]
    
    proposed_result <- score_func(
      start_k,
      data_randomness,
      privacy_noises,
      dp_statistic,
      depth_type
    )
    
    if (is.finite(proposed_result)) {
      any_resolved <- TRUE
      
      if (
        (-proposed_result) >=
        acceptance_threshold
      ) {
        return(
          list(
            status = "accept",
            par = start_k
          )
        )
      }
      
    } else {
      any_failure <- TRUE
    }
  }
  
  # Step 2: Multi-start local optimization via L-BFGS-B
  for (k in seq_len(nrow(starts))) {
    start_k <- starts[
      k,
      ,
      drop = TRUE
    ]
    
    opt <- tryCatch(
      optim(
        par = start_k,
        fn = score_for_optim,
        method = "L-BFGS-B",
        lower = c(
          search_lower,
          nuisance_lower
        ),
        upper = c(
          search_upper,
          nuisance_upper
        ),
        control = list(
          maxit = 100
        ),
        data_randomness = data_randomness,
        privacy_noises = privacy_noises,
        dp_statistic = dp_statistic,
        depth_type = depth_type,
        score_func = score_func
      ),
      error = function(e) NULL
    )
    
    if (
      is.null(opt) ||
      !is.finite(opt$value) ||
      opt$value >= 1e11
    ) {
      any_failure <- TRUE
      next
    }
    
    final_score <- score_func(
      opt$par,
      data_randomness,
      privacy_noises,
      dp_statistic,
      depth_type
    )
    
    if (!is.finite(final_score)) {
      any_failure <- TRUE
      next
    }
    
    any_resolved <- TRUE
    
    if (
      (-final_score) >=
      acceptance_threshold
    ) {
      return(
        list(
          status = "accept",
          par = opt$par
        )
      )
    }
  }
  
  if (any_failure) {
    return(
      list(
        status = "unresolved",
        par = c(
          NA_real_,
          NA_real_
        )
      )
    )
  }
  
  if (any_resolved) {
    return(
      list(
        status = "reject",
        par = c(
          NA_real_,
          NA_real_
        )
      )
    )
  }
  
  list(
    status = "unresolved",
    par = c(
      NA_real_,
      NA_real_
    )
  )
}

getConfidenceInterval <- function(
    optim_par,
    dp_statistic,
    data_randomness,
    privacy_noises,
    search_lower,
    search_upper,
    nuisance_lower,
    nuisance_upper,
    depth_type,
    score_func
) {
  initial <- accept_box(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    search_lower,
    search_upper,
    nuisance_lower,
    nuisance_upper,
    depth_type,
    score_func
  )
  
  if (
    initial$status ==
    "unresolved"
  ) {
    cat(
      "Initial CI search unresolved; returning full search range.\n"
    )
    
    return(
      c(
        search_lower,
        search_upper
      )
    )
  }
  
  if (
    initial$status !=
    "accept"
  ) {
    cat(
      "Failed to find an accepted point in CI search.\n"
    )
    
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  initial_par <- initial$par
  theta_accepted_middle_val <- initial_par[1]
  
  # Search left boundary via bisection
  left_lower <- search_lower
  left_upper <- theta_accepted_middle_val - tol
  
  while (
    left_upper -
    left_lower >
    tol
  ) {
    left_middle <- (
      left_lower +
        left_upper
    ) / 2
    
    res <- accept_box(
      initial_par,
      data_randomness,
      privacy_noises,
      dp_statistic,
      left_lower,
      left_middle,
      nuisance_lower,
      nuisance_upper,
      depth_type,
      score_func
    )
    
    if (
      res$status ==
      "accept"
    ) {
      left_upper <- res$par[1] - tol
      
    } else if (
      res$status ==
      "reject"
    ) {
      left_lower <- left_middle
      
    } else {
      left_lower <- search_lower
      break
    }
  }
  
  # Search right boundary via bisection
  right_lower <- theta_accepted_middle_val + tol
  right_upper <- search_upper
  
  while (
    right_upper -
    right_lower >
    tol
  ) {
    right_middle <- (
      right_lower +
        right_upper
    ) / 2
    
    res <- accept_box(
      initial_par,
      data_randomness,
      privacy_noises,
      dp_statistic,
      right_middle,
      right_upper,
      nuisance_lower,
      nuisance_upper,
      depth_type,
      score_func
    )
    
    if (
      res$status ==
      "accept"
    ) {
      right_lower <- res$par[1] + tol
      
    } else if (
      res$status ==
      "reject"
    ) {
      right_upper <- right_middle
      
    } else {
      right_upper <- search_upper
      break
    }
  }
  
  c(
    left_lower,
    right_upper
  )
}

# ------------------------------------------------------------------------------
# 9. Fixed Draw Simulation & Boundary Sweep (Repro Methods)
# ------------------------------------------------------------------------------

dp_statistic <- c(1, 0.75)

d1_lower <- -10
d1_upper <- 10
d2_lower <- tol
d2_upper <- 10

set.seed(1)

data_randomness <- matrix(
  rnorm(n * R_synthetic),
  ncol = n,
  nrow = R_synthetic
)

privacy_noises <- matrix(
  rnorm(R_synthetic * 2),
  ncol = 2,
  nrow = R_synthetic
)

method_colors <- c(
  mahalanobis = "#1b98e0",
  efficient_mu = "#D55E00",
  efficient_sigma = "#CC79A7",
  PB_ADI = "#009E73"
)

all_boundaries <- list()

for (
  depth_type in
  c(
    "mahalanobis",
    "efficient_mu",
    "efficient_sigma"
  )
) {
  cat(depth_type, "\n")
  cat(
    d1_lower,
    d1_upper,
    d2_lower,
    d2_upper,
    "\n"
  )
  
  optim_par0 <- c(1, 1)
  
  d1_new_range <- getConfidenceInterval(
    optim_par0,
    dp_statistic,
    data_randomness,
    privacy_noises,
    d1_lower,
    d1_upper,
    d2_lower,
    d2_upper,
    depth_type,
    score_mu_sa
  )
  
  d2_new_range <- getConfidenceInterval(
    optim_par0,
    dp_statistic,
    data_randomness,
    privacy_noises,
    d2_lower,
    d2_upper,
    d1_lower,
    d1_upper,
    depth_type,
    score_sa_mu
  )
  
  cat(
    d1_new_range,
    d2_new_range,
    "\n"
  )
  
  if (
    any(
      !is.finite(
        c(
          d1_new_range,
          d2_new_range
        )
      )
    )
  ) {
    stop(
      paste(
        "Could not construct finite search ranges for",
        depth_type
      )
    )
  }
  
  boundary <- foreach(
    i = 1:(value_r^2),
    .combine = "rbind",
    .options.snow = opts,
    .export = setdiff(
      ls(envir = .GlobalEnv),
      c("cl", "pb")
    )
  ) %dopar% {
    d1_idx <- (i - 1) %% value_r
    d2_idx <- as.integer(
      (i - 1) / value_r
    )
    
    search_lower <- d1_new_range[1] +
      d1_idx / value_r *
      (
        d1_new_range[2] -
          d1_new_range[1]
      )
    
    search_upper <- d1_new_range[1] +
      (d1_idx + 1) / value_r *
      (
        d1_new_range[2] -
          d1_new_range[1]
      )
    
    nuisance_lower <- d2_new_range[1] +
      d2_idx / value_r *
      (
        d2_new_range[2] -
          d2_new_range[1]
      )
    
    nuisance_upper <- d2_new_range[1] +
      (d2_idx + 1) / value_r *
      (
        d2_new_range[2] -
          d2_new_range[1]
      )
    
    mid_pt <- c(
      (
        search_lower +
          search_upper
      ) / 2,
      (
        nuisance_lower +
          nuisance_upper
      ) / 2
    )
    
    result <- accept_box(
      mid_pt,
      data_randomness,
      privacy_noises,
      dp_statistic,
      search_lower,
      search_upper,
      nuisance_lower,
      nuisance_upper,
      depth_type,
      score_mu_sa
    )
    
    if (
      result$status ==
      "accept"
    ) {
      c(
        search_lower,
        nuisance_lower,
        search_upper,
        nuisance_upper,
        1
      )
      
    } else if (
      result$status ==
      "unresolved"
    ) {
      c(
        search_lower,
        nuisance_lower,
        search_upper,
        nuisance_upper,
        2
      )
      
    } else {
      c(
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        0
      )
    }
  }
  
  boundary_valid <- boundary[
    is.finite(
      boundary[, 1]
    ),
    ,
    drop = FALSE
  ]
  
  accepted_cells <- sum(
    boundary_valid[, 5] == 1
  )
  
  unresolved_cells <- sum(
    boundary_valid[, 5] == 2
  )
  
  write.csv(
    boundary_valid,
    file.path(
      RESULTS_DIR,
      paste(
        depth_type,
        value_r,
        ep,
        "boundary.csv",
        sep = "_"
      )
    ),
    row.names = FALSE
  )
  
  area <- sum(
    (
      boundary_valid[, 3] -
        boundary_valid[, 1]
    ) *
      (
        boundary_valid[, 4] -
          boundary_valid[, 2]
      )
  )
  
  cat(
    depth_type,
    area,
    "\n"
  )
  
  cat(
    "Accepted cells:",
    accepted_cells,
    "\n"
  )
  
  cat(
    "Unresolved cells:",
    unresolved_cells,
    "\n"
  )
  
  all_boundaries[[depth_type]] <- list(
    boundary = boundary_valid[
      ,
      1:4,
      drop = FALSE
    ],
    status = boundary_valid[, 5],
    area = area,
    mu_interval = d1_new_range,
    sigma_interval = d2_new_range,
    n_accepted = accepted_cells,
    n_unresolved = unresolved_cells
  )
  
  pdf(
    file.path(
      RESULTS_DIR,
      paste(
        depth_type,
        value_r,
        ep,
        "region.pdf",
        sep = "_"
      )
    ),
    width = 5,
    height = 5
  )
  
  plot(
    c(0, 4),
    c(0, 4),
    type = "n",
    xlab = "",
    ylab = "",
    main = area
  )
  
  if (
    nrow(
      boundary_valid
    ) >
    0L
  ) {
    idx_accept <- boundary_valid[, 5] == 1
    
    if (any(idx_accept)) {
      rect(
        boundary_valid[
          idx_accept,
          1
        ],
        boundary_valid[
          idx_accept,
          2
        ],
        boundary_valid[
          idx_accept,
          3
        ],
        boundary_valid[
          idx_accept,
          4
        ],
        border = NA,
        col = method_colors[[depth_type]]
      )
    }
    
    idx_unresolved <- boundary_valid[, 5] == 2
    
    if (any(idx_unresolved)) {
      rect(
        boundary_valid[
          idx_unresolved,
          1
        ],
        boundary_valid[
          idx_unresolved,
          2
        ],
        boundary_valid[
          idx_unresolved,
          3
        ],
        boundary_valid[
          idx_unresolved,
          4
        ],
        border = method_colors[[depth_type]],
        lty = 2,
        col = adjustcolor(
          method_colors[[depth_type]],
          alpha.f = 0.20
        )
      )
    }
  }
  
  invisible(
    dev.off()
  )
}

# ------------------------------------------------------------------------------
# 10. Parametric Bootstrap Adaptive Indirect Inference (PB-ADI)
# ------------------------------------------------------------------------------

cat("PB_ADI\n")

adi_score <- function(
    optim_par,
    idr,
    ipn,
    observed
) {
  mu <- optim_par[1]
  sa <- optim_par[2]
  
  if (
    !is.finite(mu) ||
    !is.finite(sa) ||
    sa <= d2_lower
  ) {
    return(1e12)
  }
  
  synth <- sdp_vec(
    idr,
    ipn,
    sa,
    mu
  )
  
  dep <- depth_mahalanobis_cross(
    observed,
    synth
  )
  
  if (!is.finite(dep)) {
    return(1e12)
  }
  
  -dep
}

adi_estimate <- function(
    observed,
    seed
) {
  set.seed(seed)
  
  idr <- matrix(
    rnorm(
      R_indirect_est *
        n
    ),
    nrow = R_indirect_est,
    ncol = n
  )
  
  ipn <- matrix(
    rnorm(
      2 *
        R_indirect_est
    ),
    nrow = R_indirect_est,
    ncol = 2
  )
  
  initial <- c(
    observed[1],
    sqrt(
      max(
        1e-12,
        observed[2]
      )
    )
  )
  
  initial <- pmin(
    pmax(
      initial,
      c(
        d1_lower,
        d2_lower
      )
    ),
    c(
      d1_upper,
      d2_upper
    )
  )
  
  opt <- optim(
    par = initial,
    fn = adi_score,
    lower = c(
      d1_lower,
      d2_lower
    ),
    upper = c(
      d1_upper,
      d2_upper
    ),
    method = "L-BFGS-B",
    idr = idr,
    ipn = ipn,
    observed = observed
  )
  
  if (
    !is.finite(opt$value) ||
    opt$value >= 1e11 ||
    any(!is.finite(opt$par))
  ) {
    stop(
      sprintf(
        "ADI optimization returned non-finite output (code %d): %s",
        opt$convergence,
        opt$message
      )
    )
  }
  
  opt$par
}

theta_hat <- adi_estimate(
  dp_statistic,
  seed = 400000
)

set.seed(500000)

dr_boot <- matrix(
  rnorm(
    B_paramboot *
      n
  ),
  nrow = B_paramboot,
  ncol = n
)

pn_boot <- matrix(
  rnorm(
    2 *
      B_paramboot
  ),
  nrow = B_paramboot,
  ncol = 2
)

boot_summaries <- sdp_vec(
  dr_boot,
  pn_boot,
  theta_hat[2],
  theta_hat[1]
)

boot_results <- foreach(
  b = 1:B_paramboot,
  .combine = "rbind",
  .export = setdiff(
    ls(envir = .GlobalEnv),
    c("cl", "pb")
  )
) %dopar% {
  tryCatch(
    adi_estimate(
      observed = boot_summaries[
        b,
      ],
      seed = 600000 + b
    ),
    error = function(e) {
      c(
        NA_real_,
        NA_real_
      )
    }
  )
}

boot_est <- boot_results[
  apply(
    boot_results,
    1,
    function(x) {
      all(is.finite(x))
    }
  ),
  ,
  drop = FALSE
]

B_valid <- nrow(
  boot_est
)

if (
  B_valid <
  pb_minimum_valid
) {
  stop(
    sprintf(
      "Only %d of %d PB estimates were valid.",
      B_valid,
      B_paramboot
    )
  )
}

qp <- c(
  alpha / 2,
  1 - alpha / 2
)

basic_draws <- -boot_est +
  2 *
  matrix(
    theta_hat,
    nrow(boot_est),
    2,
    byrow = TRUE
  )

mu_interval <- as.numeric(
  quantile(
    basic_draws[, 1],
    probs = qp,
    type = 7
  )
)

mu_interval <- c(
  max(
    d1_lower,
    mu_interval[1]
  ),
  min(
    d1_upper,
    mu_interval[2]
  )
)

sigma_interval <- as.numeric(
  quantile(
    pmax(
      d2_lower,
      basic_draws[, 2]
    ),
    probs = qp,
    type = 7
  )
)

sigma_interval <- c(
  max(
    d2_lower,
    sigma_interval[1]
  ),
  min(
    d2_upper,
    sigma_interval[2]
  )
)

cat(
  "PB_ADI marginal",
  mu_interval,
  sigma_interval,
  "\n"
)

boot_cov <- stabilize_cov2(
  cov(
    basic_draws
  ),
  1e-8
)

boot_cov_inv <- safe_inv2(
  boot_cov,
  1e-8
)

centered_reflected <- sweep(
  basic_draws,
  2,
  theta_hat,
  "-"
)

quadratic <- rowSums(
  (
    centered_reflected %*%
      boot_cov_inv
  ) *
    centered_reflected
)

cutoff <- unname(
  quantile(
    quadratic,
    1 - alpha,
    type = 7
  )
)

pb_area <- pi *
  cutoff *
  sqrt(
    max(
      det(
        boot_cov
      ),
      0
    )
  )

cat(
  "PB_ADI joint area (ellipse)",
  pb_area,
  "\n"
)

eg <- eigen(
  boot_cov,
  symmetric = TRUE
)

ev <- pmax(
  eg$values,
  0
)

cov_sqrt <- eg$vectors %*%
  diag(
    sqrt(ev),
    nrow = 2
  ) %*%
  t(eg$vectors)

theta_grid_pts <- seq(
  0,
  2 * pi,
  length.out = 400
)

ellipse_pts <- t(
  theta_hat +
    sqrt(cutoff) *
    (
      cov_sqrt %*%
        rbind(
          cos(theta_grid_pts),
          sin(theta_grid_pts)
        )
    )
)

write.csv(
  ellipse_pts,
  file.path(
    RESULTS_DIR,
    paste(
      "PB_ADI",
      value_r,
      ep,
      "ellipse_boundary.csv",
      sep = "_"
    )
  ),
  row.names = FALSE
)

write.csv(
  boot_est,
  file.path(
    RESULTS_DIR,
    paste(
      "PB_ADI",
      value_r,
      ep,
      "boot_estimates.csv",
      sep = "_"
    )
  ),
  row.names = FALSE
)

all_boundaries[["PB_ADI"]] <- list(
  boundary = matrix(
    numeric(0),
    ncol = 4
  ),
  area = pb_area,
  mu_interval = mu_interval,
  sigma_interval = sigma_interval
)

pdf(
  file.path(
    RESULTS_DIR,
    paste(
      "PB_ADI",
      value_r,
      ep,
      "region.pdf",
      sep = "_"
    )
  ),
  width = 5,
  height = 5
)

plot(
  c(0, 4),
  c(0, 4),
  type = "n",
  xlab = "",
  ylab = "",
  main = pb_area
)

polygon(
  ellipse_pts[, 1],
  ellipse_pts[, 2],
  border = NA,
  col = method_colors[["PB_ADI"]]
)

points(
  theta_hat[1],
  theta_hat[2],
  pch = 4,
  lwd = 2
)

invisible(
  dev.off()
)

try(
  close(pb),
  silent = TRUE
)

try(
  stopCluster(cl),
  silent = TRUE
)

# ------------------------------------------------------------------------------
# 11. Joint Confidence Region Comparison Plot
# ------------------------------------------------------------------------------

truth_color <- "#111111"

all_x <- c(
  population_mu,
  dp_statistic[1],
  ellipse_pts[, 1]
)

all_y <- c(
  population_sigma,
  dp_statistic[2],
  ellipse_pts[, 2]
)

for (
  m in names(
    all_boundaries
  )
) {
  b <- all_boundaries[[m]]$boundary
  
  if (
    nrow(b) >
    0
  ) {
    all_x <- c(
      all_x,
      b[, 1],
      b[, 3]
    )
    
    all_y <- c(
      all_y,
      b[, 2],
      b[, 4]
    )
  }
}

all_x <- all_x[
  is.finite(
    all_x
  )
]

all_y <- all_y[
  is.finite(
    all_y
  )
]

x_pad <- max(
  0.08 *
    diff(
      range(
        all_x
      )
    ),
  0.05
)

y_pad <- max(
  0.08 *
    diff(
      range(
        all_y
      )
    ),
  0.05
)

xlim_combined <- range(
  all_x
) +
  c(
    -x_pad,
    x_pad
  )

ylim_combined <- range(
  all_y
) +
  c(
    -y_pad,
    y_pad
  )

pdf(
  file.path(
    RESULTS_DIR,
    paste(
      "comparison",
      value_r,
      ep,
      "region.pdf",
      sep = "_"
    )
  ),
  width = 6.4,
  height = 5.6
)

plot(
  NA,
  xlim = xlim_combined,
  ylim = ylim_combined,
  xlab = expression(mu),
  ylab = expression(sigma),
  main = "Confidence region comparison"
)

for (
  m in names(
    all_boundaries
  )
) {
  b <- all_boundaries[[m]]$boundary
  
  if (
    nrow(b) >
    0
  ) {
    rect(
      b[, 1],
      b[, 2],
      b[, 3],
      b[, 4],
      border = NA,
      col = adjustcolor(
        method_colors[[m]],
        alpha.f = 0.35
      )
    )
  }
}

polygon(
  ellipse_pts[, 1],
  ellipse_pts[, 2],
  border = method_colors[["PB_ADI"]],
  col = adjustcolor(
    method_colors[["PB_ADI"]],
    alpha.f = 0.30
  ),
  lwd = 1.8
)

points(
  population_mu,
  population_sigma,
  pch = 8,
  cex = 1.4,
  lwd = 1.8,
  col = truth_color
)

points(
  theta_hat[1],
  theta_hat[2],
  pch = 4,
  cex = 1.2,
  lwd = 1.6,
  col = truth_color
)

legend(
  "topright",
  inset = 0.02,
  bty = "n",
  cex = 0.82,
  legend = c(
    sprintf(
      "Mahalanobis (area=%.3f)",
      all_boundaries[["mahalanobis"]]$area
    ),
    sprintf(
      "Efficient \u03bc (area=%.3f)",
      all_boundaries[["efficient_mu"]]$area
    ),
    sprintf(
      "Efficient \u03c3 (area=%.3f)",
      all_boundaries[["efficient_sigma"]]$area
    ),
    sprintf(
      "PB-ADI (area=%.3f)",
      pb_area
    ),
    "True (\u03bc*, \u03c3*)",
    "PB-ADI theta_hat"
  ),
  pch = c(
    22,
    22,
    22,
    22,
    8,
    4
  ),
  pt.bg = c(
    adjustcolor(
      method_colors[["mahalanobis"]],
      alpha.f = 0.5
    ),
    adjustcolor(
      method_colors[["efficient_mu"]],
      alpha.f = 0.5
    ),
    adjustcolor(
      method_colors[["efficient_sigma"]],
      alpha.f = 0.5
    ),
    adjustcolor(
      method_colors[["PB_ADI"]],
      alpha.f = 0.5
    ),
    NA,
    NA
  ),
  col = c(
    method_colors[["mahalanobis"]],
    method_colors[["efficient_mu"]],
    method_colors[["efficient_sigma"]],
    method_colors[["PB_ADI"]],
    truth_color,
    truth_color
  )
)

invisible(
  dev.off()
)

# ------------------------------------------------------------------------------
# 12. Marginal Confidence Interval Plots
# ------------------------------------------------------------------------------

open_pdf <- function(
    path,
    width = 6.8,
    height = 4.4
) {
  if (
    capabilities(
      "cairo"
    )
  ) {
    grDevices::cairo_pdf(
      filename = path,
      width = width,
      height = height,
      family = "serif",
      onefile = TRUE
    )
    
  } else {
    grDevices::pdf(
      file = path,
      width = width,
      height = height,
      family = "serif",
      useDingbats = FALSE
    )
  }
}

draw_marginal_ci_plot <- function(
    param_symbol,
    truth,
    intervals,
    method_names,
    colors,
    output_file,
    main_title
) {
  n_methods <- length(
    intervals
  )
  
  y_positions <- seq(
    n_methods,
    1
  )
  
  finite_vals <- c(
    unlist(
      intervals
    ),
    truth
  )
  
  finite_vals <- finite_vals[
    is.finite(
      finite_vals
    )
  ]
  
  x_range <- range(
    finite_vals
  )
  
  x_pad <- max(
    0.08 *
      diff(
        x_range
      ),
    0.03
  )
  
  x_limits <- x_range +
    c(
      -x_pad,
      x_pad
    )
  
  open_pdf(
    output_file,
    width = 6.8,
    height = 4.4
  )
  
  par(
    mar = c(
      4.4,
      8.2,
      3.6,
      1.4
    ),
    mgp = c(
      2.6,
      0.75,
      0
    ),
    tcl = -0.25,
    las = 1,
    xaxs = "i",
    cex.axis = 0.95,
    cex.lab = 1.15
  )
  
  plot(
    NA,
    xlim = x_limits,
    ylim = c(
      0.5,
      n_methods + 0.5
    ),
    xlab = param_symbol,
    ylab = "",
    yaxt = "n",
    axes = FALSE
  )
  
  axis(1)
  
  box(
    bty = "l",
    lwd = 0.9
  )
  
  abline(
    h = y_positions,
    col = "grey92",
    lwd = 0.7
  )
  
  abline(
    v = truth,
    lty = 2,
    lwd = 1.2,
    col = "grey35"
  )
  
  for (
    j in
    seq_len(
      n_methods
    )
  ) {
    lo <- intervals[[j]][1]
    hi <- intervals[[j]][2]
    y <- y_positions[j]
    
    segments(
      lo,
      y,
      hi,
      y,
      col = colors[j],
      lwd = 3.2,
      lend = "butt"
    )
    
    segments(
      lo,
      y - 0.09,
      lo,
      y + 0.09,
      col = colors[j],
      lwd = 1.5
    )
    
    segments(
      hi,
      y - 0.09,
      hi,
      y + 0.09,
      col = colors[j],
      lwd = 1.5
    )
  }
  
  axis(
    2,
    at = y_positions,
    labels = method_names,
    tick = FALSE,
    cex.axis = 0.95
  )
  
  title(
    main = main_title,
    line = 1.6,
    cex.main = 1.08
  )
  
  legend(
    "bottomright",
    inset = 0.02,
    bty = "n",
    cex = 0.82,
    legend = as.expression(
      bquote(
        .(param_symbol)^"*"
      )
    ),
    lty = 2,
    lwd = 1.2,
    col = "grey35"
  )
  
  invisible(
    dev.off()
  )
}

method_labels <- c(
  mahalanobis = "Mahalanobis",
  efficient_mu = "Efficient (\u03bc)",
  efficient_sigma = "Efficient (\u03c3)",
  PB_ADI = "PB-ADI"
)

method_order <- c(
  "mahalanobis",
  "efficient_mu",
  "efficient_sigma",
  "PB_ADI"
)

plot_colors <- unname(
  method_colors[
    method_order
  ]
)

plot_labels <- unname(
  method_labels[
    method_order
  ]
)

mu_intervals <- lapply(
  method_order,
  function(m) {
    all_boundaries[[m]]$mu_interval
  }
)

sigma_intervals <- lapply(
  method_order,
  function(m) {
    all_boundaries[[m]]$sigma_interval
  }
)

draw_marginal_ci_plot(
  param_symbol = quote(mu),
  truth = population_mu,
  intervals = mu_intervals,
  method_names = plot_labels,
  colors = plot_colors,
  output_file = file.path(
    RESULTS_DIR,
    paste0(
      "marginal_CI_mu_",
      value_r,
      "_",
      ep,
      ".pdf"
    )
  ),
  main_title = expression(
    "95% Marginal Confidence Intervals for " *
      mu
  )
)

draw_marginal_ci_plot(
  param_symbol = quote(sigma),
  truth = population_sigma,
  intervals = sigma_intervals,
  method_names = plot_labels,
  colors = plot_colors,
  output_file = file.path(
    RESULTS_DIR,
    paste0(
      "marginal_CI_sigma_",
      value_r,
      "_",
      ep,
      ".pdf"
    )
  ),
  main_title = expression(
    "95% Marginal Confidence Intervals for " *
      sigma
  )
)

cat(
  "\nAll outputs written to:",
  RESULTS_DIR,
  "\n"
)

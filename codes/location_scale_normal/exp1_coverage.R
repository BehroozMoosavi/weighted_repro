# ==============================================================================
# Script: exp1_coverage.R
# Description: Multi-replicate (nSIM) coverage and area simulation study
#              (marginal and joint coverage) across four inferential methods:
#              Mahalanobis Repro, Efficient Repro (mu-interest), Efficient
#              Repro (sigma-interest), and PB-ADI. Produces summary tables
#              and publication-ready comparison plots.
# ==============================================================================

PROJECT_DIR <- path.expand("~/R_Simuls/Location_Scale_nomal")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

RESULTS_DIR <- file.path(PROJECT_DIR, "results_coverage")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Data Generating Process (DGP) & Privacy Mechanism
# ------------------------------------------------------------------------------

upper_clamp <- 3
lower_clamp <- 0
n <- 100
ep <- 1

clean_clamp_meanvar <- function(x) {
  clamp_x <- pmax(lower_clamp, pmin(upper_clamp, x))
  c(mean(clamp_x), var(clamp_x))
}

sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
  raw <- sa * data_randomness + mu
  
  clamped_data <- pmax(
    pmin(raw, upper_clamp),
    lower_clamp
  )
  
  row_mean <- rowMeans(clamped_data)
  n_cols <- ncol(clamped_data)
  
  centered_ss <- pmax(
    rowSums(clamped_data^2) -
      n_cols * row_mean^2,
    0
  )
  
  row_var <- centered_ss / (n_cols - 1)
  
  mean_noise_scale <-
    (upper_clamp - lower_clamp) /
    (n * ep)
  
  var_noise_scale <-
    (upper_clamp - lower_clamp)^2 /
    (n * ep)
  
  cbind(
    row_mean +
      mean_noise_scale * privacy_noises[, 1],
    row_var +
      var_noise_scale * privacy_noises[, 2]
  )
}

# ------------------------------------------------------------------------------
# 2. Simulation Parameters & Optimization Bounds
# ------------------------------------------------------------------------------

nSIM <- 1000
value_r <- 50
R_synthetic <- 200
alpha <- 0.05
tol <- 1e-8
bisection_tol <- 10e-4

population_mu <- 1
population_sigma <- 1

R_aux <- 400
R_indirect_est <- 50
B_paramboot <- 200
h_fd <- 1e-3

lambda_n <- 1 / log(n)

pb_minimum_valid <- max(
  20L,
  floor(0.8 * B_paramboot)
)

mu_search_lower <- -2
mu_search_upper <- 5

sa_search_lower <- 0.1
sa_search_upper <- 5

acceptance_threshold <-
  floor(alpha * (R_synthetic + 1)) + 1

result_digits <- 6

# ------------------------------------------------------------------------------
# 3. Environment & Parallel Backend Setup
# ------------------------------------------------------------------------------

list.of.packages <- c(
  "foreach",
  "doSNOW",
  "parallelly"
)

new.packages <- list.of.packages[
  !(list.of.packages %in%
      installed.packages()[, "Package"])
]

if (length(new.packages) > 0) {
  install.packages(
    new.packages,
    dependencies = TRUE
  )
}

for (package.i in list.of.packages) {
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

# Restrict worker multithreading to avoid BLAS/MKL CPU oversubscription
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

n.cores <- min(
  124L,
  max(
    1L,
    parallel::detectCores() - 1L
  )
)

cl <- parallel::makePSOCKcluster(n.cores)
doSNOW::registerDoSNOW(cl)

pb <- txtProgressBar(
  max = nSIM,
  style = 3
)

progress <- function(k) {
  setTxtProgressBar(pb, k)
}

opts <- list(
  progress = progress
)

start_time <- Sys.time()

# ------------------------------------------------------------------------------
# 4. Linear Algebra Utilities
# ------------------------------------------------------------------------------

safe_inv2 <- function(A, ridge = 1e-10) {
  if (
    is.null(A) ||
    !all(dim(A) == c(2L, 2L)) ||
    any(!is.finite(A))
  ) {
    stop("Invalid matrix in safe_inv2().")
  }
  
  A <- 0.5 * (A + t(A)) +
    ridge * diag(2L)
  
  a <- A[1L, 1L]
  b <- A[1L, 2L]
  d <- A[2L, 2L]
  
  det_A <- a * d - b^2
  
  if (
    !is.finite(det_A) ||
    det_A <= ridge^2
  ) {
    eg <- eigen(
      A,
      symmetric = TRUE
    )
    
    if (
      any(!is.finite(eg$values)) ||
      any(!is.finite(eg$vectors))
    ) {
      stop("Invalid eigendecomposition.")
    }
    
    ev <- pmax(
      eg$values,
      ridge
    )
    
    return(
      eg$vectors %*%
        diag(
          1 / ev,
          nrow = 2L
        ) %*%
        t(eg$vectors)
    )
  }
  
  matrix(
    c(
      d,
      -b,
      -b,
      a
    ),
    nrow = 2L
  ) / det_A
}

stabilize_cov2 <- function(A, ridge = 1e-8) {
  if (
    is.null(A) ||
    !all(dim(A) == c(2L, 2L)) ||
    any(!is.finite(A))
  ) {
    return(
      ridge * diag(2L)
    )
  }
  
  0.5 * (A + t(A)) +
    ridge * diag(2L)
}

# ------------------------------------------------------------------------------
# 5. Efficient Direction & Depth Scoring (GLS Direction)
# ------------------------------------------------------------------------------

compute_aux_derivatives <- function(
    mu,
    sa,
    aux_dr,
    aux_pn
) {
  if (
    !is.finite(mu) ||
    !is.finite(sa) ||
    sa < sa_search_lower ||
    sa > sa_search_upper
  ) {
    return(NULL)
  }
  
  h <- h_fd
  
  synth_center <- tryCatch(
    sdp_vec(
      aux_dr,
      aux_pn,
      sa,
      mu
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(synth_center) ||
    any(!is.finite(synth_center))
  ) {
    return(NULL)
  }
  
  center_mean <- colMeans(
    synth_center
  )
  
  # Numerical derivative with respect to sigma
  if (
    sa - h > sa_search_lower &&
    sa + h < sa_search_upper
  ) {
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
    
  } else if (
    sa - h <= sa_search_lower
  ) {
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
    
  } else {
    synth_sa_minus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa - h,
      mu
    )
    
    d_sa <- (
      center_mean -
        colMeans(synth_sa_minus)
    ) / h
  }
  
  # Numerical derivative with respect to mu
  if (
    mu - h > mu_search_lower &&
    mu + h < mu_search_upper
  ) {
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
    
  } else if (
    mu - h <= mu_search_lower
  ) {
    synth_mu_plus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa,
      mu + h
    )
    
    d_mu <- (
      colMeans(synth_mu_plus) -
        center_mean
    ) / h
    
  } else {
    synth_mu_minus <- sdp_vec(
      aux_dr,
      aux_pn,
      sa,
      mu - h
    )
    
    d_mu <- (
      center_mean -
        colMeans(synth_mu_minus)
    ) / h
  }
  
  if (
    any(
      !is.finite(
        c(
          d_mu,
          d_sa
        )
      )
    )
  ) {
    return(NULL)
  }
  
  aux_cov <- cov(
    synth_center
  )
  
  aux_cov <- 0.5 * (
    aux_cov +
      t(aux_cov)
  ) +
    1e-8 * diag(2L)
  
  cov_inv <- tryCatch(
    safe_inv2(
      aux_cov,
      1e-8
    ),
    error = function(e) NULL
  )
  
  if (is.null(cov_inv)) {
    return(NULL)
  }
  
  list(
    d_mu = matrix(
      d_mu,
      ncol = 1L
    ),
    d_sa = matrix(
      d_sa,
      ncol = 1L
    ),
    cov_inv = cov_inv
  )
}

efficient_direction <- function(
    mu,
    sa,
    interest,
    aux_dr,
    aux_pn,
    cache
) {
  key <- sprintf(
    "%.10g_%.10g",
    mu,
    sa
  )
  
  if (
    !is.null(cache) &&
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
    sa,
    aux_dr,
    aux_pn
  )
  
  if (is.null(deriv)) {
    direction <- c(
      NA_real_,
      NA_real_
    )
    
    if (!is.null(cache)) {
      assign(
        key,
        direction,
        envir = cache
      )
    }
    
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
    
    if (!is.null(cache)) {
      assign(
        key,
        direction,
        envir = cache
      )
    }
    
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
  
  if (!is.null(cache)) {
    assign(
      key,
      direction,
      envir = cache
    )
  }
  
  direction
}

efficient_depth <- function(
    synth,
    mu,
    sa,
    interest,
    aux_dr,
    aux_pn,
    cache
) {
  direction <- efficient_direction(
    mu,
    sa,
    interest,
    aux_dr,
    aux_pn,
    cache
  )
  
  if (
    any(
      !is.finite(
        direction
      )
    )
  ) {
    return(
      rep(
        NA_real_,
        nrow(synth)
      )
    )
  }
  
  center <- colMeans(
    synth
  )
  
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
  
  efficient_quadratic <-
    projection^2 /
    eff_scale
  
  mahalanobis_quadratic <- rowSums(
    (
      centered %*%
        cov_inv
    ) *
      centered
  )
  
  penalty <-
    efficient_quadratic +
    lambda_n *
    mahalanobis_quadratic
  
  depth <- 1 / (
    1 +
      penalty
  )
  
  if (
    any(
      !is.finite(
        depth
      )
    )
  ) {
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
# 6. Mahalanobis Depth Scoring
# ------------------------------------------------------------------------------

depth_mahalanobis_selfref <- function(
    synth
) {
  m <- colMeans(
    synth
  )
  
  S <- cov(
    synth
  )
  
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
  
  if (
    any(
      !is.finite(
        depth
      )
    )
  ) {
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
  m <- colMeans(
    synth
  )
  
  S <- cov(
    synth
  )
  
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
# 7. Criterion Functions & Multi-Start Bisection Routines
# ------------------------------------------------------------------------------

score_mu_sa <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    aux_dr,
    aux_pn,
    eff_cache
) {
  mu <- optim_par[1]
  sa <- optim_par[2]
  
  if (
    !is.finite(mu) ||
    !is.finite(sa) ||
    mu < mu_search_lower ||
    mu > mu_search_upper ||
    sa < sa_search_lower ||
    sa > sa_search_upper
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
      interest,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
  } else {
    D_synth <- depth_mahalanobis_selfref(
      synth
    )
  }
  
  if (
    any(
      !is.finite(
        D_synth
      )
    )
  ) {
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
    depth_type,
    aux_dr,
    aux_pn,
    eff_cache
) {
  sa <- optim_par[1]
  mu <- optim_par[2]
  
  score_mu_sa(
    c(
      mu,
      sa
    ),
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    aux_dr,
    aux_pn,
    eff_cache
  )
}

score_for_optim <- function(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    score_func,
    aux_dr,
    aux_pn,
    eff_cache
) {
  ans <- score_func(
    optim_par,
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    aux_dr,
    aux_pn,
    eff_cache
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
    score_func,
    aux_dr,
    aux_pn,
    eff_cache
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
  for (
    k in seq_len(
      nrow(starts)
    )
  ) {
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
      depth_type,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    if (
      is.finite(
        proposed_result
      )
    ) {
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
  for (
    k in seq_len(
      nrow(starts)
    )
  ) {
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
        data_randomness =
          data_randomness,
        privacy_noises =
          privacy_noises,
        dp_statistic =
          dp_statistic,
        depth_type =
          depth_type,
        score_func =
          score_func,
        aux_dr =
          aux_dr,
        aux_pn =
          aux_pn,
        eff_cache =
          eff_cache
      ),
      error = function(e)
        NULL
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
      depth_type,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    if (
      !is.finite(
        final_score
      )
    ) {
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
    score_func,
    aux_dr,
    aux_pn,
    eff_cache
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
    score_func,
    aux_dr,
    aux_pn,
    eff_cache
  )
  
  if (
    initial$status ==
    "unresolved"
  ) {
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
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  initial_par <- initial$par
  t_mid <- initial_par[1]
  
  # Search left boundary via bisection
  l_low <- search_lower
  l_up <- t_mid - tol
  
  while (
    l_up -
    l_low >
    bisection_tol
  ) {
    l_mid <- (
      l_low +
        l_up
    ) / 2
    
    res <- accept_box(
      c(
        l_mid,
        initial_par[2]
      ),
      data_randomness,
      privacy_noises,
      dp_statistic,
      l_low,
      l_mid,
      nuisance_lower,
      nuisance_upper,
      depth_type,
      score_func,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    if (
      res$status ==
      "accept"
    ) {
      l_up <- min(
        l_up,
        res$par[1] - tol
      )
      
    } else if (
      res$status ==
      "reject"
    ) {
      l_low <- l_mid
      
    } else {
      # Conservative unresolved fallback
      l_low <- search_lower
      break
    }
  }
  
  # Search right boundary via bisection
  r_low <- t_mid + tol
  r_up <- search_upper
  
  while (
    r_up -
    r_low >
    bisection_tol
  ) {
    r_mid <- (
      r_low +
        r_up
    ) / 2
    
    res <- accept_box(
      c(
        r_mid,
        initial_par[2]
      ),
      data_randomness,
      privacy_noises,
      dp_statistic,
      r_mid,
      r_up,
      nuisance_lower,
      nuisance_upper,
      depth_type,
      score_func,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    if (
      res$status ==
      "accept"
    ) {
      r_low <- max(
        r_low,
        res$par[1] + tol
      )
      
    } else if (
      res$status ==
      "reject"
    ) {
      r_up <- r_mid
      
    } else {
      # Conservative unresolved fallback
      r_up <- search_upper
      break
    }
  }
  
  c(
    l_low,
    r_up
  )
}

# ------------------------------------------------------------------------------
# 8. Parametric Bootstrap Adaptive Indirect Inference (PB-ADI) Helpers
# ------------------------------------------------------------------------------

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
    mu < mu_search_lower ||
    mu > mu_search_upper ||
    sa < sa_search_lower ||
    sa > sa_search_upper
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
  
  if (
    !is.finite(
      dep
    )
  ) {
    return(1e12)
  }
  
  -dep
}

adi_estimate <- function(
    observed
) {
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
        mu_search_lower,
        sa_search_lower
      )
    ),
    c(
      mu_search_upper,
      sa_search_upper
    )
  )
  
  opt <- tryCatch(
    optim(
      par = initial,
      fn = adi_score,
      lower = c(
        mu_search_lower,
        sa_search_lower
      ),
      upper = c(
        mu_search_upper,
        sa_search_upper
      ),
      method = "L-BFGS-B",
      idr = idr,
      ipn = ipn,
      observed = observed
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(opt) ||
    !is.finite(opt$value) ||
    opt$value >= 1e11 ||
    any(!is.finite(opt$par))
  ) {
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  opt$par
}

# ------------------------------------------------------------------------------
# 9. Main Parallel Simulation Loop
# ------------------------------------------------------------------------------

cat(
  "Starting coverage study: nSIM =",
  nSIM,
  ", value_r =",
  value_r,
  "\n"
)

sim_results <- foreach(
  s = 1:nSIM,
  .combine = "rbind",
  .options.snow = opts,
  .export = setdiff(
    ls(envir = .GlobalEnv),
    c(
      "cl",
      "pb"
    )
  )
) %dopar% {
  set.seed(
    s +
      123
  )
  
  # Generate observed private sample statistics
  raw_dr <- matrix(
    rnorm(n),
    nrow = 1,
    ncol = n
  )
  
  raw_pn <- matrix(
    rnorm(2),
    nrow = 1,
    ncol = 2
  )
  
  dp_statistic <- as.numeric(
    sdp_vec(
      raw_dr,
      raw_pn,
      population_sigma,
      population_mu
    )
  )
  
  # Inference randomness clouds
  data_randomness <- matrix(
    rnorm(
      n *
        R_synthetic
    ),
    ncol = n,
    nrow = R_synthetic
  )
  
  privacy_noises <- matrix(
    rnorm(
      R_synthetic *
        2
    ),
    ncol = 2,
    nrow = R_synthetic
  )
  
  aux_dr <- matrix(
    rnorm(
      n *
        R_aux
    ),
    nrow = R_aux,
    ncol = n
  )
  
  aux_pn <- matrix(
    rnorm(
      2 *
        R_aux
    ),
    nrow = R_aux,
    ncol = 2
  )
  
  row_out <- c()
  
  # Evaluation: Repro methods (Mahalanobis, Efficient-mu, Efficient-sigma)
  for (
    depth_type in
    c(
      "mahalanobis",
      "efficient_mu",
      "efficient_sigma"
    )
  ) {
    eff_cache <- if (
      depth_type %in%
      c(
        "efficient_mu",
        "efficient_sigma"
      )
    ) {
      new.env(
        parent = emptyenv()
      )
    } else {
      NULL
    }
    
    mu_range <- getConfidenceInterval(
      c(
        1,
        1
      ),
      dp_statistic,
      data_randomness,
      privacy_noises,
      mu_search_lower,
      mu_search_upper,
      sa_search_lower,
      sa_search_upper,
      depth_type,
      score_mu_sa,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    sigma_range <- getConfidenceInterval(
      c(
        1,
        1
      ),
      dp_statistic,
      data_randomness,
      privacy_noises,
      sa_search_lower,
      sa_search_upper,
      mu_search_lower,
      mu_search_upper,
      depth_type,
      score_sa_mu,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    if (
      any(
        !is.finite(
          c(
            mu_range,
            sigma_range
          )
        )
      )
    ) {
      row_out <- c(
        row_out,
        cov_mu = NA_real_,
        cov_sigma = NA_real_,
        cov_joint = NA_real_,
        width_mu = NA_real_,
        width_sigma = NA_real_,
        area = NA_real_,
        failure = 1
      )
      
      next
    }
    
    cov_mu <- as.numeric(
      mu_range[1] <= population_mu &&
        population_mu <= mu_range[2]
    )
    
    cov_sigma <- as.numeric(
      sigma_range[1] <= population_sigma &&
        population_sigma <= sigma_range[2]
    )
    
    true_score <- score_mu_sa(
      c(
        population_mu,
        population_sigma
      ),
      data_randomness,
      privacy_noises,
      dp_statistic,
      depth_type,
      aux_dr,
      aux_pn,
      eff_cache
    )
    
    cov_joint <- if (
      is.finite(
        true_score
      )
    ) {
      as.numeric(
        (-true_score) >=
          acceptance_threshold
      )
    } else {
      NA_real_
    }
    
    # Numerical area computation across grid cells
    mu_vals <- seq(
      mu_range[1],
      mu_range[2],
      length.out = value_r
    )
    
    sa_vals <- seq(
      sigma_range[1],
      sigma_range[2],
      length.out = value_r
    )
    
    step_mu <- mu_vals[2] -
      mu_vals[1]
    
    step_sa <- sa_vals[2] -
      sa_vals[1]
    
    total_area <- 0
    unresolved_area_cells <- 0L
    
    for (
      m_i in
      1:(value_r - 1)
    ) {
      for (
        s_i in
        1:(value_r - 1)
      ) {
        mid_pt <- c(
          (
            mu_vals[m_i] +
              mu_vals[m_i + 1]
          ) / 2,
          (
            sa_vals[s_i] +
              sa_vals[s_i + 1]
          ) / 2
        )
        
        score_val <- score_mu_sa(
          mid_pt,
          data_randomness,
          privacy_noises,
          dp_statistic,
          depth_type,
          aux_dr,
          aux_pn,
          eff_cache
        )
        
        if (
          is.finite(
            score_val
          )
        ) {
          if (
            (-score_val) >=
            acceptance_threshold
          ) {
            total_area <-
              total_area +
              step_mu *
              step_sa
          }
          
        } else {
          total_area <-
            total_area +
            step_mu *
            step_sa
          
          unresolved_area_cells <-
            unresolved_area_cells +
            1L
        }
      }
    }
    
    failure_flag <- as.numeric(
      !is.finite(
        cov_joint
      ) ||
        unresolved_area_cells >
        0L
    )
    
    row_out <- c(
      row_out,
      cov_mu = cov_mu,
      cov_sigma = cov_sigma,
      cov_joint = cov_joint,
      width_mu =
        mu_range[2] -
        mu_range[1],
      width_sigma =
        sigma_range[2] -
        sigma_range[1],
      area = total_area,
      failure = failure_flag
    )
  }
  
  # Evaluation: PB-ADI
  theta_hat <- adi_estimate(
    dp_statistic
  )
  
  if (
    any(
      !is.finite(
        theta_hat
      )
    )
  ) {
    row_out <- c(
      row_out,
      pb_cov_mu = NA_real_,
      pb_cov_sigma = NA_real_,
      pb_cov_joint = NA_real_,
      pb_width_mu = NA_real_,
      pb_width_sigma = NA_real_,
      pb_area = NA_real_,
      pb_failure = 1
    )
    
  } else {
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
    
    boot_est <- t(
      sapply(
        seq_len(
          B_paramboot
        ),
        function(b) {
          adi_estimate(
            boot_summaries[
              b,
            ]
          )
        }
      )
    )
    
    boot_est <- boot_est[
      apply(
        boot_est,
        1,
        function(x) {
          all(
            is.finite(
              x
            )
          )
        }
      ),
      ,
      drop = FALSE
    ]
    
    if (
      nrow(
        boot_est
      ) <
      pb_minimum_valid
    ) {
      row_out <- c(
        row_out,
        pb_cov_mu = NA_real_,
        pb_cov_sigma = NA_real_,
        pb_cov_joint = NA_real_,
        pb_width_mu = NA_real_,
        pb_width_sigma = NA_real_,
        pb_area = NA_real_,
        pb_failure = 1
      )
      
    } else {
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
      
      pb_mu_interval <- as.numeric(
        quantile(
          basic_draws[, 1],
          probs = qp,
          type = 7
        )
      )
      
      pb_mu_interval <- c(
        max(
          mu_search_lower,
          pb_mu_interval[1]
        ),
        min(
          mu_search_upper,
          pb_mu_interval[2]
        )
      )
      
      pb_sigma_interval <- as.numeric(
        quantile(
          pmax(
            sa_search_lower,
            basic_draws[, 2]
          ),
          probs = qp,
          type = 7
        )
      )
      
      pb_sigma_interval <- c(
        max(
          sa_search_lower,
          pb_sigma_interval[1]
        ),
        min(
          sa_search_upper,
          pb_sigma_interval[2]
        )
      )
      
      pb_cov_mu <- as.numeric(
        pb_mu_interval[1] <= population_mu &&
          population_mu <= pb_mu_interval[2]
      )
      
      pb_cov_sigma <- as.numeric(
        pb_sigma_interval[1] <= population_sigma &&
          population_sigma <= pb_sigma_interval[2]
      )
      
      boot_cov <- stabilize_cov2(
        cov(
          basic_draws
        ),
        1e-8
      )
      
      boot_cov_inv <- tryCatch(
        safe_inv2(
          boot_cov,
          1e-8
        ),
        error = function(e) NULL
      )
      
      if (
        is.null(
          boot_cov_inv
        )
      ) {
        row_out <- c(
          row_out,
          pb_cov_mu = NA_real_,
          pb_cov_sigma = NA_real_,
          pb_cov_joint = NA_real_,
          pb_width_mu = NA_real_,
          pb_width_sigma = NA_real_,
          pb_area = NA_real_,
          pb_failure = 1
        )
        
      } else {
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
        
        truth_centered <- c(
          population_mu,
          population_sigma
        ) -
          theta_hat
        
        truth_quadratic <- as.numeric(
          t(truth_centered) %*%
            boot_cov_inv %*%
            truth_centered
        )
        
        pb_cov_joint <- as.numeric(
          truth_quadratic <=
            cutoff
        )
        
        row_out <- c(
          row_out,
          pb_cov_mu = pb_cov_mu,
          pb_cov_sigma = pb_cov_sigma,
          pb_cov_joint = pb_cov_joint,
          pb_width_mu =
            pb_mu_interval[2] -
            pb_mu_interval[1],
          pb_width_sigma =
            pb_sigma_interval[2] -
            pb_sigma_interval[1],
          pb_area = pb_area,
          pb_failure = 0
        )
      }
    }
  }
  
  row_out
}

try(
  close(pb),
  silent = TRUE
)

try(
  stopCluster(cl),
  silent = TRUE
)

# ------------------------------------------------------------------------------
# 10. Aggregation & Output Export
# ------------------------------------------------------------------------------

methods_univariate <- c(
  "mahalanobis",
  "efficient_mu",
  "efficient_sigma"
)

method_labels <- c(
  mahalanobis = "Mahalanobis",
  efficient_mu = "Efficient (\u03bc-interest)",
  efficient_sigma = "Efficient (\u03c3-interest)"
)

univariate_rows <- list()

for (
  i in seq_along(
    methods_univariate
  )
) {
  m <- methods_univariate[i]
  prefix <- (
    i -
      1
  ) *
    7
  
  cov_mu_vec <- sim_results[
    ,
    prefix + 1
  ]
  
  cov_sigma_vec <- sim_results[
    ,
    prefix + 2
  ]
  
  cov_joint_vec <- sim_results[
    ,
    prefix + 3
  ]
  
  width_mu_vec <- sim_results[
    ,
    prefix + 4
  ]
  
  width_sigma_vec <- sim_results[
    ,
    prefix + 5
  ]
  
  area_vec <- sim_results[
    ,
    prefix + 6
  ]
  
  failure_vec <- sim_results[
    ,
    prefix + 7
  ]
  
  n_mu_valid <- sum(
    is.finite(
      cov_mu_vec
    )
  )
  
  n_sigma_valid <- sum(
    is.finite(
      cov_sigma_vec
    )
  )
  
  n_joint_valid <- sum(
    is.finite(
      cov_joint_vec
    )
  )
  
  coverage_mu_mean <- mean(
    cov_mu_vec,
    na.rm = TRUE
  )
  
  coverage_sigma_mean <- mean(
    cov_sigma_vec,
    na.rm = TRUE
  )
  
  coverage_joint_mean <- mean(
    cov_joint_vec,
    na.rm = TRUE
  )
  
  univariate_rows[[m]] <- data.frame(
    method =
      method_labels[[m]],
    coverage_mu =
      coverage_mu_mean,
    coverage_mu_se =
      sqrt(
        coverage_mu_mean *
          (
            1 -
              coverage_mu_mean
          ) /
          max(
            n_mu_valid,
            1
          )
      ),
    coverage_sigma =
      coverage_sigma_mean,
    coverage_sigma_se =
      sqrt(
        coverage_sigma_mean *
          (
            1 -
              coverage_sigma_mean
          ) /
          max(
            n_sigma_valid,
            1
          )
      ),
    coverage_joint =
      coverage_joint_mean,
    coverage_joint_se =
      sqrt(
        coverage_joint_mean *
          (
            1 -
              coverage_joint_mean
          ) /
          max(
            n_joint_valid,
            1
          )
      ),
    width_mu =
      mean(
        width_mu_vec,
        na.rm = TRUE
      ),
    width_mu_se =
      sd(
        width_mu_vec,
        na.rm = TRUE
      ) /
      sqrt(
        max(
          sum(
            is.finite(
              width_mu_vec
            )
          ),
          1
        )
      ),
    width_sigma =
      mean(
        width_sigma_vec,
        na.rm = TRUE
      ),
    width_sigma_se =
      sd(
        width_sigma_vec,
        na.rm = TRUE
      ) /
      sqrt(
        max(
          sum(
            is.finite(
              width_sigma_vec
            )
          ),
          1
        )
      ),
    area =
      mean(
        area_vec,
        na.rm = TRUE
      ),
    area_se =
      sd(
        area_vec,
        na.rm = TRUE
      ) /
      sqrt(
        max(
          sum(
            is.finite(
              area_vec
            )
          ),
          1
        )
      ),
    failure_rate =
      mean(
        failure_vec,
        na.rm = TRUE
      ),
    n_joint_valid =
      n_joint_valid
  )
}

pb_prefix <- length(
  methods_univariate
) *
  7

pb_cov_mu_vec <- sim_results[
  ,
  pb_prefix + 1
]

pb_cov_sigma_vec <- sim_results[
  ,
  pb_prefix + 2
]

pb_cov_joint_vec <- sim_results[
  ,
  pb_prefix + 3
]

pb_width_mu_vec <- sim_results[
  ,
  pb_prefix + 4
]

pb_width_sigma_vec <- sim_results[
  ,
  pb_prefix + 5
]

pb_area_vec <- sim_results[
  ,
  pb_prefix + 6
]

pb_failure_vec <- sim_results[
  ,
  pb_prefix + 7
]

pb_n_mu_valid <- sum(
  is.finite(
    pb_cov_mu_vec
  )
)

pb_n_sigma_valid <- sum(
  is.finite(
    pb_cov_sigma_vec
  )
)

pb_n_joint_valid <- sum(
  is.finite(
    pb_cov_joint_vec
  )
)

pb_coverage_mu <- mean(
  pb_cov_mu_vec,
  na.rm = TRUE
)

pb_coverage_sigma <- mean(
  pb_cov_sigma_vec,
  na.rm = TRUE
)

pb_coverage_joint <- mean(
  pb_cov_joint_vec,
  na.rm = TRUE
)

pb_row <- data.frame(
  method =
    "PB-ADI",
  coverage_mu =
    pb_coverage_mu,
  coverage_mu_se =
    sqrt(
      pb_coverage_mu *
        (
          1 -
            pb_coverage_mu
        ) /
        max(
          pb_n_mu_valid,
          1
        )
    ),
  coverage_sigma =
    pb_coverage_sigma,
  coverage_sigma_se =
    sqrt(
      pb_coverage_sigma *
        (
          1 -
            pb_coverage_sigma
        ) /
        max(
          pb_n_sigma_valid,
          1
        )
    ),
  coverage_joint =
    pb_coverage_joint,
  coverage_joint_se =
    sqrt(
      pb_coverage_joint *
        (
          1 -
            pb_coverage_joint
        ) /
        max(
          pb_n_joint_valid,
          1
        )
    ),
  width_mu =
    mean(
      pb_width_mu_vec,
      na.rm = TRUE
    ),
  width_mu_se =
    sd(
      pb_width_mu_vec,
      na.rm = TRUE
    ) /
    sqrt(
      max(
        sum(
          is.finite(
            pb_width_mu_vec
          )
        ),
        1
      )
    ),
  width_sigma =
    mean(
      pb_width_sigma_vec,
      na.rm = TRUE
    ),
  width_sigma_se =
    sd(
      pb_width_sigma_vec,
      na.rm = TRUE
    ) /
    sqrt(
      max(
        sum(
          is.finite(
            pb_width_sigma_vec
          )
        ),
        1
      )
    ),
  area =
    mean(
      pb_area_vec,
      na.rm = TRUE
    ),
  area_se =
    sd(
      pb_area_vec,
      na.rm = TRUE
    ) /
    sqrt(
      max(
        sum(
          is.finite(
            pb_area_vec
          )
        ),
        1
      )
    ),
  failure_rate =
    mean(
      pb_failure_vec,
      na.rm = TRUE
    ),
  n_joint_valid =
    pb_n_joint_valid
)

summary_table <- do.call(
  rbind,
  c(
    univariate_rows,
    list(
      pb_row
    )
  )
)

summary_table[
  ,
  -1
] <- round(
  summary_table[
    ,
    -1
  ],
  result_digits
)

print(
  summary_table
)

write.csv(
  summary_table,
  file.path(
    RESULTS_DIR,
    paste0(
      "coverage_study_n",
      n,
      "_ep",
      ep,
      "_nSIM",
      nSIM,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  sim_results,
  file.path(
    RESULTS_DIR,
    paste0(
      "coverage_study_raw_n",
      n,
      "_ep",
      ep,
      "_nSIM",
      nSIM,
      ".csv"
    )
  ),
  row.names = FALSE
)

cat(
  "\n============================================================\n"
)

cat(
  "Coverage study complete. Elapsed:",
  format(
    Sys.time() -
      start_time
  ),
  "\n"
)

cat(
  "For comparison, the paper's Table 2 (n=100, ep=1) reports:\n"
)

cat(
  "  PB (adaptive indirect): coverage(\u03bc,\u03c3)=0.943(0.007), area=0.339(0.004)\n"
)

cat(
  "  Repro (Mahalanobis):    coverage(\u03bc,\u03c3)=0.971(0.005), area=0.3619(0.004)\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nNumerical failure rates:\n"
)

print(
  summary_table[
    ,
    c(
      "method",
      "failure_rate",
      "n_joint_valid"
    )
  ]
)

# ------------------------------------------------------------------------------
# 11. Plot Generation
# ------------------------------------------------------------------------------

method_colors_named <- c(
  "Mahalanobis" = "#1b98e0",
  "Efficient (\u03bc-interest)" = "#D55E00",
  "Efficient (\u03c3-interest)" = "#CC79A7",
  "PB-ADI" = "#009E73"
)

plot_colors <- unname(
  method_colors_named[
    summary_table$method
  ]
)

plot_colors[
  is.na(
    plot_colors
  )
] <- "#7F7F7F"

alpha_used <- 0.05

CONTENT_SINGLE_W <- 6.8
CONTENT_SINGLE_H <- 5.4

CONTENT_PANEL_W <- 9.6
CONTENT_PANEL_H <- 4.6

draw_topright_guide <- function(
    label,
    cex = 0.78
) {
  usr <- par("usr")
  
  x_pos <- usr[2] -
    0.02 *
    diff(
      usr[1:2]
    )
  
  y_pos <- usr[4] -
    0.04 *
    diff(
      usr[3:4]
    )
  
  text(
    x_pos,
    y_pos,
    label,
    adj = c(
      1,
      1
    ),
    cex = cex,
    col = "grey30",
    xpd = TRUE
  )
}

# Plot 1: Joint area comparison across methods
pdf(
  file.path(
    RESULTS_DIR,
    "joint_area_comparison.pdf"
  ),
  width = CONTENT_SINGLE_W,
  height = CONTENT_SINGLE_H
)

par(
  mai = c(
    1.3,
    1.0,
    0.7,
    0.3
  ),
  mgp = c(
    2.4,
    0.7,
    0
  ),
  tcl = -0.25,
  las = 2,
  cex.axis = 0.85,
  cex.lab = 0.95,
  cex.main = 1.0
)

bar_x <- barplot(
  summary_table$area,
  names.arg = summary_table$method,
  col = plot_colors,
  border = NA,
  ylim = c(
    0,
    max(
      summary_table$area +
        summary_table$area_se,
      na.rm = TRUE
    ) *
      1.2
  ),
  ylab = "Average joint area",
  main = "Average confidence-region area by method"
)

arrows(
  bar_x,
  summary_table$area -
    summary_table$area_se,
  bar_x,
  summary_table$area +
    summary_table$area_se,
  angle = 90,
  code = 3,
  length = 0.05,
  lwd = 1.2
)

box(
  bty = "l"
)

invisible(
  dev.off()
)

# Plot 2: Empirical joint coverage comparison
pdf(
  file.path(
    RESULTS_DIR,
    "joint_coverage_comparison.pdf"
  ),
  width = CONTENT_SINGLE_W,
  height = CONTENT_SINGLE_H
)

par(
  mai = c(
    0.9,
    1.9,
    0.7,
    0.3
  ),
  mgp = c(
    2.4,
    0.7,
    0
  ),
  tcl = -0.25,
  las = 1,
  cex.axis = 0.88,
  cex.lab = 0.95,
  cex.main = 1.0
)

n_methods <- nrow(
  summary_table
)

y_positions <- seq(
  n_methods,
  1
)

x_range <- range(
  c(
    summary_table$coverage_joint -
      summary_table$coverage_joint_se,
    summary_table$coverage_joint +
      summary_table$coverage_joint_se,
    1 -
      alpha_used
  ),
  na.rm = TRUE
)

x_pad <- max(
  0.15 *
    diff(
      x_range
    ),
  0.003
)

plot(
  NA,
  xlim = x_range +
    c(
      -x_pad,
      x_pad
    ),
  ylim = c(
    0.4,
    n_methods +
      0.6
  ),
  xlab = "Joint coverage",
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
  v = 1 -
    alpha_used,
  lty = 2,
  lwd = 1.2,
  col = "grey35"
)

for (
  j in seq_len(
    n_methods
  )
) {
  arrows(
    summary_table$coverage_joint[j] -
      summary_table$coverage_joint_se[j],
    y_positions[j],
    summary_table$coverage_joint[j] +
      summary_table$coverage_joint_se[j],
    y_positions[j],
    angle = 90,
    code = 3,
    length = 0.05,
    col = plot_colors[j],
    lwd = 1.4
  )
  
  points(
    summary_table$coverage_joint[j],
    y_positions[j],
    pch = 19,
    cex = 1.5,
    col = plot_colors[j]
  )
  
  text(
    summary_table$coverage_joint[j],
    y_positions[j] +
      0.26,
    sprintf(
      "%.3f",
      summary_table$coverage_joint[j]
    ),
    cex = 0.78,
    xpd = TRUE
  )
}

axis(
  2,
  at = y_positions,
  labels = summary_table$method,
  tick = FALSE,
  cex.axis = 0.85
)

title(
  main = "Joint coverage by method",
  line = 1.5
)

draw_topright_guide(
  sprintf(
    "nominal %.0f%%",
    100 *
      (
        1 -
          alpha_used
      )
  )
)

invisible(
  dev.off()
)

# Plot 3: Scatter plot of joint coverage vs. region area
pdf(
  file.path(
    RESULTS_DIR,
    "coverage_vs_area_scatter.pdf"
  ),
  width = CONTENT_SINGLE_W,
  height = CONTENT_SINGLE_H
)

par(
  mai = c(
    0.9,
    1.0,
    0.7,
    0.3
  ),
  mgp = c(
    2.4,
    0.7,
    0
  ),
  tcl = -0.25,
  cex.axis = 0.88,
  cex.lab = 0.95,
  cex.main = 1.0
)

x_range <- range(
  c(
    summary_table$area -
      summary_table$area_se,
    summary_table$area +
      summary_table$area_se
  ),
  na.rm = TRUE
)

y_range <- range(
  c(
    summary_table$coverage_joint -
      summary_table$coverage_joint_se,
    summary_table$coverage_joint +
      summary_table$coverage_joint_se,
    1 -
      alpha_used
  ),
  na.rm = TRUE
)

x_pad <- 0.08 *
  diff(
    x_range
  )

y_pad <- 0.10 *
  diff(
    y_range
  )

plot(
  NA,
  xlim = x_range +
    c(
      -x_pad,
      x_pad
    ),
  ylim = y_range +
    c(
      -y_pad,
      1.6 *
        y_pad
    ),
  xlab = "Average area",
  ylab = "Joint coverage",
  main = "Coverage vs. area (nSIM-averaged)"
)

abline(
  h = 1 -
    alpha_used,
  lty = 2,
  lwd = 1.2,
  col = "grey35"
)

for (
  i in seq_len(
    nrow(
      summary_table
    )
  )
) {
  arrows(
    summary_table$area[i] -
      summary_table$area_se[i],
    summary_table$coverage_joint[i],
    summary_table$area[i] +
      summary_table$area_se[i],
    summary_table$coverage_joint[i],
    angle = 90,
    code = 3,
    length = 0.04,
    col = plot_colors[i],
    lwd = 1.3
  )
  
  arrows(
    summary_table$area[i],
    summary_table$coverage_joint[i] -
      summary_table$coverage_joint_se[i],
    summary_table$area[i],
    summary_table$coverage_joint[i] +
      summary_table$coverage_joint_se[i],
    angle = 90,
    code = 3,
    length = 0.04,
    col = plot_colors[i],
    lwd = 1.3
  )
  
  points(
    summary_table$area[i],
    summary_table$coverage_joint[i],
    pch = 19,
    cex = 1.4,
    col = plot_colors[i]
  )
}

draw_topright_guide(
  sprintf(
    "nominal %.0f%%",
    100 *
      (
        1 -
          alpha_used
      )
  )
)

legend(
  "topright",
  inset = c(
    0.02,
    0.09
  ),
  bty = "n",
  cex = 0.78,
  legend = summary_table$method,
  pch = 19,
  col = plot_colors
)

box(
  bty = "l"
)

invisible(
  dev.off()
)

# Plot 4: Marginal coverage comparison (mu and sigma)
draw_coverage_panel <- function(
    values,
    se,
    xlab,
    main_title,
    ref_line
) {
  n_m <- length(
    values
  )
  
  y_pos <- seq(
    n_m,
    1
  )
  
  x_range <- range(
    c(
      values -
        se,
      values +
        se,
      ref_line
    ),
    na.rm = TRUE
  )
  
  x_pad <- max(
    0.15 *
      diff(
        x_range
      ),
    0.003
  )
  
  par(
    mai = c(
      0.9,
      1.8,
      0.6,
      0.25
    ),
    mgp = c(
      2.3,
      0.65,
      0
    ),
    tcl = -0.25,
    las = 1,
    cex.axis = 0.82,
    cex.main = 0.95
  )
  
  plot(
    NA,
    xlim = x_range +
      c(
        -x_pad,
        x_pad
      ),
    ylim = c(
      0.4,
      n_m +
        0.6
    ),
    xlab = xlab,
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
    h = y_pos,
    col = "grey92",
    lwd = 0.7
  )
  
  abline(
    v = ref_line,
    lty = 2,
    lwd = 1.2,
    col = "grey35"
  )
  
  for (
    j in seq_len(
      n_m
    )
  ) {
    arrows(
      values[j] -
        se[j],
      y_pos[j],
      values[j] +
        se[j],
      y_pos[j],
      angle = 90,
      code = 3,
      length = 0.045,
      col = plot_colors[j],
      lwd = 1.3
    )
    
    points(
      values[j],
      y_pos[j],
      pch = 19,
      cex = 1.35,
      col = plot_colors[j]
    )
    
    text(
      values[j],
      y_pos[j] +
        0.26,
      sprintf(
        "%.3f",
        values[j]
      ),
      cex = 0.72,
      xpd = TRUE
    )
  }
  
  axis(
    2,
    at = y_pos,
    labels = summary_table$method,
    tick = FALSE,
    cex.axis = 0.82
  )
  
  title(
    main = main_title,
    line = 1.3,
    cex.main = 0.95
  )
  
  draw_topright_guide(
    sprintf(
      "nominal %.0f%%",
      100 *
        ref_line
    ),
    cex = 0.68
  )
}

pdf(
  file.path(
    RESULTS_DIR,
    "marginal_coverage_comparison.pdf"
  ),
  width = CONTENT_PANEL_W,
  height = CONTENT_PANEL_H
)

par(
  mfrow = c(
    1,
    2
  )
)

draw_coverage_panel(
  summary_table$coverage_mu,
  summary_table$coverage_mu_se,
  expression(
    "Coverage of " *
      mu
  ),
  expression(
    "Marginal coverage: " *
      mu
  ),
  1 -
    alpha_used
)

draw_coverage_panel(
  summary_table$coverage_sigma,
  summary_table$coverage_sigma_se,
  expression(
    "Coverage of " *
      sigma
  ),
  expression(
    "Marginal coverage: " *
      sigma
  ),
  1 -
    alpha_used
)

invisible(
  dev.off()
)

# Plot 5: Marginal interval width comparison (mu and sigma)
draw_marginal_bar <- function(
    values,
    se,
    ylab,
    main_title
) {
  y_top <- max(
    values +
      se,
    na.rm = TRUE
  ) *
    1.2
  
  par(
    mai = c(
      1.3,
      0.9,
      0.6,
      0.25
    ),
    mgp = c(
      2.3,
      0.65,
      0
    ),
    tcl = -0.25,
    las = 2,
    cex.axis = 0.78,
    cex.main = 0.95
  )
  
  bx <- barplot(
    values,
    names.arg = summary_table$method,
    col = plot_colors,
    border = NA,
    ylim = c(
      0,
      y_top
    ),
    ylab = ylab,
    main = main_title
  )
  
  arrows(
    bx,
    values -
      se,
    bx,
    values +
      se,
    angle = 90,
    code = 3,
    length = 0.04,
    lwd = 1.1
  )
  
  box(
    bty = "l"
  )
}

pdf(
  file.path(
    RESULTS_DIR,
    "marginal_width_comparison.pdf"
  ),
  width = CONTENT_PANEL_W,
  height = CONTENT_PANEL_H
)

par(
  mfrow = c(
    1,
    2
  )
)

draw_marginal_bar(
  summary_table$width_mu,
  summary_table$width_mu_se,
  expression(
    "Average width of " *
      mu *
      " CI"
  ),
  expression(
    "Marginal width: " *
      mu
  )
)

draw_marginal_bar(
  summary_table$width_sigma,
  summary_table$width_sigma_se,
  expression(
    "Average width of " *
      sigma *
      " CI"
  ),
  expression(
    "Marginal width: " *
      sigma
  )
)

invisible(
  dev.off()
)

cat(
  "\nPlots written to:",
  RESULTS_DIR,
  "\n"
)

cat(
  "  joint_area_comparison.pdf\n"
)

cat(
  "  joint_coverage_comparison.pdf\n"
)

cat(
  "  coverage_vs_area_scatter.pdf   <- start here\n"
)

cat(
  "  marginal_coverage_comparison.pdf\n"
)

cat(
  "  marginal_width_comparison.pdf\n"
)

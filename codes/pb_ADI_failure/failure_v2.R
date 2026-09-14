# ==============================================================================
# Script: failure_v2.R
# Description: Comparative evaluation of Hard Clamp vs. Soft Clamp (logistic 
#              transformation) under a boundary parameter regime (mu* = 3.0) 
#              across privacy levels epsilon in [0.1, 1.0]. Evaluates stability 
#              and coverage performance for Mahalanobis Repro, Efficient Repro, 
#              and PB-ADI.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Directory Initialization
# ------------------------------------------------------------------------------

PROJECT_DIR <- path.expand("~/R_Simuls/PB_adi_failure")
dir.create(PROJECT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PROJECT_DIR)

RESULTS_DIR <- file.path(PROJECT_DIR, "results_softclamp_study_v2")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Experimental Configurations & Search Bounds
# ------------------------------------------------------------------------------

n <- 100
upper_clamp <- 3
lower_clamp <- 0
population_mu <- 3.0
population_sigma <- 1.0

nSIM <- 100
R_synthetic <- 200
value_r <- 50
alpha <- 0.05
tol <- 1e-8

R_indirect_est <- 50
B_paramboot <- 200
result_digits <- 6

h_fd <- 1e-3
R_aux <- 400
lambda_n <- 1 / log(n)

acceptance_threshold <- floor(alpha * (R_synthetic + 1)) + 1

ep_grid <- c(0.1, 0.2, 0.5, 1.0)
clamp_types <- c("hard", "soft")

# Search bounds: Repro and Efficient
mu_search_lower_repro <- -2
mu_search_upper_repro <- 5
sa_search_lower_repro <- 0.1
sa_search_upper_repro <- 5

# Search bounds: PB-ADI
mu_lower_adi <- -2
mu_upper_adi <- 10
sa_lower_adi <- 1e-6
sa_upper_adi <- 10

# ------------------------------------------------------------------------------
# 3. Package Dependencies & Parallel Setup
# ------------------------------------------------------------------------------

list.of.packages <- c("foreach", "doSNOW", "ddalpha", "parallelly")
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

n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
doSNOW::registerDoSNOW(cl)

# ------------------------------------------------------------------------------
# 4. Clamping Transformations & Private DGP Functions
# ------------------------------------------------------------------------------

# Logistic Soft Clamp Transformation mapping R strictly into (lower_clamp, upper_clamp)
soft_clamp_fn <- function(x, k = 1.5) {
  mid_pt <- (lower_clamp + upper_clamp) / 2
  lower_clamp + (upper_clamp - lower_clamp) / (1 + exp(-k * (x - mid_pt)))
}

sdp_vec_repro <- function(data_randomness, privacy_noises, sa, mu, ep, clamp_type = "hard") {
  n_inner <- dim(data_randomness)[2]
  data <- sa * data_randomness + mu
  
  if (clamp_type == "soft") {
    data_clamp <- soft_clamp_fn(data)
  } else {
    data_clamp <- pmax(pmin(data, upper_clamp), lower_clamp)
  }
  
  s1 <- rowMeans(data_clamp) + (upper_clamp - lower_clamp) / (n_inner * ep) * privacy_noises[, 1]
  n_cols <- ncol(data_clamp)
  centered_ss <- pmax(rowSums(data_clamp^2) - n_cols * (rowMeans(data_clamp))^2, 0)
  s2 <- centered_ss / (n_cols - 1) + (upper_clamp - lower_clamp)^2 / (n_inner * ep) * privacy_noises[, 2]
  
  cbind(s1, s2)
}

# ------------------------------------------------------------------------------
# 5. Matrix Utilities
# ------------------------------------------------------------------------------

safe_inv2 <- function(A, ridge = 1e-10) {
  A <- 0.5 * (A + t(A)) + ridge * diag(2L)
  a <- A[1L, 1L]
  b <- A[1L, 2L]
  d <- A[2L, 2L]
  
  det_A <- a * d - b^2
  
  if (!is.finite(det_A) || det_A <= ridge^2) {
    eg <- eigen(A, symmetric = TRUE)
    ev <- pmax(eg$values, ridge)
    return(eg$vectors %*% (t(eg$vectors) / ev))
  }
  
  matrix(c(d, -b, -b, a), nrow = 2L) / det_A
}

# ------------------------------------------------------------------------------
# 6. Auxiliary Derivatives & Efficient Direction
# ------------------------------------------------------------------------------

compute_aux_derivatives <- function(mu, sa, aux_dr, aux_pn, ep, clamp_type) {
  h <- h_fd
  synth_center <- sdp_vec_repro(aux_dr, aux_pn, sa, mu, ep, clamp_type)
  
  # Numerical derivative with respect to sigma
  if (sa - h > sa_search_lower_repro && sa + h < sa_search_upper_repro) {
    synth_sa_plus <- sdp_vec_repro(aux_dr, aux_pn, sa + h, mu, ep, clamp_type)
    synth_sa_minus <- sdp_vec_repro(aux_dr, aux_pn, sa - h, mu, ep, clamp_type)
    d_sa <- (colMeans(synth_sa_plus) - colMeans(synth_sa_minus)) / (2 * h)
  } else if (sa - h <= sa_search_lower_repro) {
    synth_sa_plus <- sdp_vec_repro(aux_dr, aux_pn, sa + h, mu, ep, clamp_type)
    d_sa <- (colMeans(synth_sa_plus) - colMeans(synth_center)) / h
  } else {
    synth_sa_minus <- sdp_vec_repro(aux_dr, aux_pn, sa - h, mu, ep, clamp_type)
    d_sa <- (colMeans(synth_center) - colMeans(synth_sa_minus)) / h
  }
  
  # Numerical derivative with respect to mu
  if (mu - h > mu_search_lower_repro && mu + h < mu_search_upper_repro) {
    synth_mu_plus <- sdp_vec_repro(aux_dr, aux_pn, sa, mu + h, ep, clamp_type)
    synth_mu_minus <- sdp_vec_repro(aux_dr, aux_pn, sa, mu - h, ep, clamp_type)
    d_mu <- (colMeans(synth_mu_plus) - colMeans(synth_mu_minus)) / (2 * h)
  } else if (mu - h <= mu_search_lower_repro) {
    synth_mu_plus <- sdp_vec_repro(aux_dr, aux_pn, sa, mu + h, ep, clamp_type)
    d_mu <- (colMeans(synth_mu_plus) - colMeans(synth_center)) / h
  } else {
    synth_mu_minus <- sdp_vec_repro(aux_dr, aux_pn, sa, mu - h, ep, clamp_type)
    d_mu <- (colMeans(synth_center) - colMeans(synth_mu_minus)) / h
  }
  
  aux_cov <- cov(synth_center)
  aux_cov <- 0.5 * (aux_cov + t(aux_cov)) + 1e-8 * diag(2L)
  
  list(
    d_mu    = matrix(d_mu, ncol = 1L),
    d_sa    = matrix(d_sa, ncol = 1L),
    cov_inv = safe_inv2(aux_cov, 1e-8)
  )
}

efficient_direction <- function(mu, sa, interest, aux_dr, aux_pn, cache, ep, clamp_type) {
  key <- sprintf("%.10g_%.10g", mu, sa)
  
  if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
    return(get(key, envir = cache, inherits = FALSE))
  }
  
  deriv <- compute_aux_derivatives(mu, sa, aux_dr, aux_pn, ep, clamp_type)
  
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
    if (!is.null(cache)) assign(key, direction, envir = cache)
    return(direction)
  }
  
  projection_coefficient <- as.numeric(t(nuisance) %*% deriv$cov_inv %*% target_raw) / nuisance_norm
  target <- target_raw - projection_coefficient * nuisance
  direction <- as.numeric(deriv$cov_inv %*% target)
  
  if (any(!is.finite(direction)) || sum(direction^2) < 1e-20) {
    direction <- c(NA_real_, NA_real_)
  }
  
  if (!is.null(cache)) assign(key, direction, envir = cache)
  direction
}

efficient_depth <- function(synth, mu, sa, interest, aux_dr, aux_pn, cache, ep, clamp_type) {
  direction <- efficient_direction(mu, sa, interest, aux_dr, aux_pn, cache, ep, clamp_type)
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
  efficient_quadratic <- projection^2 / eff_scale
  mahalanobis_quadratic <- rowSums((centered %*% cov_inv) * centered)
  penalty <- efficient_quadratic + lambda_n * mahalanobis_quadratic
  depth <- 1 / (1 + penalty)
  depth[!is.finite(depth)] <- 0
  depth
}

# ------------------------------------------------------------------------------
# 7. Criterion Scoring & Bisection (Repro Methods)
# ------------------------------------------------------------------------------

score_mu_sa <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type,
                        aux_dr = NULL, aux_pn = NULL, eff_cache = NULL, ep, clamp_type) {
  mu <- optim_par[1]
  sa <- optim_par[2]
  synth <- sdp_vec_repro(data_randomness, privacy_noises, sa, mu, ep, clamp_type)
  synth <- rbind(synth, dp_statistic)
  
  if (depth_type == "efficient_mu") {
    D_synth <- efficient_depth(synth, mu, sa, "mu", aux_dr, aux_pn, eff_cache, ep, clamp_type)
    if (any(!is.finite(D_synth))) return(NA_real_)
  } else {
    D_synth <- depth.Mahalanobis(synth, synth)
  }
  
  r <- rank(D_synth, ties.method = "max")[R_synthetic + 1]
  s <- r + D_synth[R_synthetic + 1]
  -s
}

score_sa_mu <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type,
                        aux_dr = NULL, aux_pn = NULL, eff_cache = NULL, ep, clamp_type) {
  score_mu_sa(
    c(optim_par[2], optim_par[1]),
    data_randomness,
    privacy_noises,
    dp_statistic,
    depth_type,
    aux_dr,
    aux_pn,
    eff_cache,
    ep,
    clamp_type
  )
}

accept <- function(optim_par, data_randomness, privacy_noises, dp_statistic,
                   search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func,
                   aux_dr = NULL, aux_pn = NULL, eff_cache = NULL, ep, clamp_type) {
  proposed_result <- score_func(
    optim_par, data_randomness, privacy_noises, dp_statistic, depth_type,
    aux_dr, aux_pn, eff_cache, ep, clamp_type
  )
  if (is.finite(proposed_result) && (-proposed_result) >= acceptance_threshold) {
    return(optim_par)
  }
  
  opt <- tryCatch(
    optim(
      par = optim_par,
      fn = score_func,
      method = "L-BFGS-B",
      lower = c(search_lower, nuisance_lower),
      upper = c(search_upper, nuisance_upper),
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = dp_statistic,
      depth_type = depth_type,
      aux_dr = aux_dr,
      aux_pn = aux_pn,
      eff_cache = eff_cache,
      ep = ep,
      clamp_type = clamp_type
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || !is.finite(opt$value)) return(c(NA, NA))
  if ((-opt$value) >= acceptance_threshold) return(opt$par)
  c(NA, NA)
}

getConfidenceInterval <- function(optim_par, dp_statistic, data_randomness, privacy_noises,
                                  search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func,
                                  aux_dr = NULL, aux_pn = NULL, eff_cache = NULL, ep, clamp_type) {
  optim_par <- accept(
    optim_par, data_randomness, privacy_noises, dp_statistic,
    search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func,
    aux_dr, aux_pn, eff_cache, ep, clamp_type
  )
  if (is.na(optim_par[1])) return(c(NA, NA))
  
  t_mid_val <- optim_par[1]
  
  # Left boundary bisection
  l_low <- search_lower
  l_up  <- t_mid_val - tol
  while (l_up - l_low > 0.1) {
    l_mid <- (l_low + l_up) / 2
    res_l <- accept(
      c(l_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic,
      l_low, l_mid, nuisance_lower, nuisance_upper, depth_type, score_func,
      aux_dr, aux_pn, eff_cache, ep, clamp_type
    )
    if (!is.na(res_l[1])) l_up <- l_mid - tol else l_low <- l_mid
  }
  
  # Right boundary bisection
  r_low <- t_mid_val + tol
  r_up  <- search_upper
  while (r_up - r_low > 0.1) {
    r_mid <- (r_low + r_up) / 2
    res_r <- accept(
      c(r_mid, optim_par[2]), data_randomness, privacy_noises, dp_statistic,
      r_mid, r_up, nuisance_lower, nuisance_upper, depth_type, score_func,
      aux_dr, aux_pn, eff_cache, ep, clamp_type
    )
    if (!is.na(res_r[1])) r_low <- r_mid + tol else r_up <- r_mid
  }
  
  c(l_low, r_up)
}

# ------------------------------------------------------------------------------
# 8. Parametric Bootstrap Adaptive Indirect Inference (PB-ADI) Routines
# ------------------------------------------------------------------------------

clean_clamp_meanvar <- function(x, clamp_type = "hard") {
  if (clamp_type == "soft") {
    clamp_x <- soft_clamp_fn(x)
  } else {
    clamp_x <- pmax(lower_clamp, pmin(upper_clamp, x))
  }
  c(mean(clamp_x), var(clamp_x))
}

sdp_vec_adi <- function(data_randomness, privacy_noises, sa, mu, clamp_type) {
  data <- sa * data_randomness + mu
  t(apply(data, 1, function(row) clean_clamp_meanvar(row, clamp_type))) + privacy_noises
}

score_adi <- function(optim_par, data_randomness, privacy_noises, dp_statistic, clamp_type) {
  mu <- optim_par[1]
  sa <- optim_par[2]
  synth <- sdp_vec_adi(data_randomness, privacy_noises, sa, mu, clamp_type)
  D_synth <- depth.Mahalanobis(dp_statistic, synth)
  -D_synth
}

solve_meanstd_from_clamp <- function(clamped_meanvar, ep, clamp_type) {
  initialized_value <- clamped_meanvar
  initialized_value[2] <- sqrt(max(1e-12, initialized_value[2]))
  
  sensitivity_adi      <- (upper_clamp - lower_clamp) / n
  sd_of_noise_mean_adi <- sensitivity_adi / ep
  sensitivity_var_adi  <- (upper_clamp - lower_clamp)^2 / n
  sd_of_noise_var_adi  <- sensitivity_var_adi / ep
  
  data_randomness <- matrix(rnorm(n * R_indirect_est), ncol = n, nrow = R_indirect_est)
  privacy_noises  <- matrix(rnorm(2 * R_indirect_est), ncol = 2, nrow = R_indirect_est)
  privacy_noises  <- t(t(privacy_noises) * c(sd_of_noise_mean_adi, sd_of_noise_var_adi))
  
  opt <- tryCatch(
    optim(
      par = initialized_value,
      fn = score_adi,
      lower = c(mu_lower_adi, sa_lower_adi),
      upper = c(mu_upper_adi, sa_upper_adi),
      method = "L-BFGS-B",
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = clamped_meanvar,
      clamp_type = clamp_type
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || !is.finite(opt$value) || any(!is.finite(opt$par))) {
    return(c(NA_real_, NA_real_))
  }
  opt$par
}

# ------------------------------------------------------------------------------
# 9. Main Study Loop: Hard vs. Soft Clamping Comparison
# ------------------------------------------------------------------------------

all_master_rows <- list()
start_time_total <- Sys.time()

for (c_type in clamp_types) {
  for (current_ep in ep_grid) {
    
    cat("\n============================================================\n")
    cat("RUNNING REGIME: Clamp =", toupper(c_type), "| mu* =", population_mu, "| epsilon =", current_ep, "| nSIM =", nSIM, "\n")
    cat("============================================================\n")
    
    pb <- txtProgressBar(max = nSIM, style = 3)
    progress <- function(k) setTxtProgressBar(pb, k)
    opts <- list(progress = progress)
    
    sim_results <- foreach(
      s = 1:nSIM,
      .combine = "rbind",
      .packages = c("ddalpha"),
      .options.snow = opts,
      .export = setdiff(ls(envir = .GlobalEnv), c("cl", "pb"))
    ) %dopar% {
      
      set.seed(s + 123)
      
      # Sample private observations
      raw_data <- rnorm(n, mean = population_mu, sd = population_sigma)
      if (c_type == "soft") {
        data_c <- soft_clamp_fn(raw_data)
      } else {
        data_c <- pmax(pmin(raw_data, upper_clamp), lower_clamp)
      }
      
      obs_s1 <- mean(data_c) + (upper_clamp - lower_clamp) / (n * current_ep) * rnorm(1)
      obs_s2 <- var(data_c) + (upper_clamp - lower_clamp)^2 / (n * current_ep) * rnorm(1)
      obs_stats <- c(obs_s1, obs_s2)
      
      row_out <- c()
      
      # Evaluate Repro (Mahalanobis) and Efficient Repro (mu-interest)
      data_rand  <- matrix(rnorm(n * R_synthetic), ncol = n, nrow = R_synthetic)
      priv_noise <- matrix(rnorm(R_synthetic * 2), ncol = 2, nrow = R_synthetic)
      aux_dr     <- matrix(rnorm(n * R_aux), nrow = R_aux, ncol = n)
      aux_pn     <- matrix(rnorm(2 * R_aux), nrow = R_aux, ncol = 2)
      
      mu_guess <- min(max(obs_stats[1], mu_search_lower_repro), mu_search_upper_repro)
      sa_guess <- min(max(sqrt(max(obs_stats[2], 1e-6)), sa_search_lower_repro), sa_search_upper_repro)
      
      for (depth_type in c("mahalanobis", "efficient_mu")) {
        this_cache <- if (depth_type == "efficient_mu") new.env(parent = emptyenv()) else NULL
        
        d1_rng <- getConfidenceInterval(
          c(mu_guess, sa_guess), obs_stats, data_rand, priv_noise,
          mu_search_lower_repro, mu_search_upper_repro, sa_search_lower_repro, sa_search_upper_repro,
          depth_type, score_mu_sa, aux_dr, aux_pn, this_cache, current_ep, c_type
        )
        d2_rng <- getConfidenceInterval(
          c(sa_guess, mu_guess), obs_stats, data_rand, priv_noise,
          sa_search_lower_repro, sa_search_upper_repro, mu_search_lower_repro, mu_search_upper_repro,
          depth_type, score_sa_mu, aux_dr, aux_pn, this_cache, current_ep, c_type
        )
        
        if (any(is.na(d1_rng)) || any(is.na(d2_rng))) {
          row_out <- c(
            row_out,
            setNames(
              c(0, 0, 0, NA, NA, 0),
              paste0(depth_type, c("_cov_joint", "_cov_mu", "_cov_sigma", "_width_mu", "_width_sigma", "_area"))
            )
          )
          next
        }
        
        cov_mu <- as.numeric(population_mu >= d1_rng[1] && population_mu <= d1_rng[2])
        cov_sigma <- as.numeric(population_sigma >= d2_rng[1] && population_sigma <= d2_rng[2])
        width_mu <- d1_rng[2] - d1_rng[1]
        width_sigma <- d2_rng[2] - d2_rng[1]
        
        # Grid sweep for joint confidence region area
        mu_vals <- seq(d1_rng[1], d1_rng[2], length.out = value_r)
        sa_vals <- seq(d2_rng[1], d2_rng[2], length.out = value_r)
        step_mu <- mu_vals[2] - mu_vals[1]
        step_sa <- sa_vals[2] - sa_vals[1]
        
        total_area <- 0
        covered <- 0
        
        for (m_i in 1:(value_r - 1)) {
          for (s_i in 1:(value_r - 1)) {
            mid_pt <- c((mu_vals[m_i] + mu_vals[m_i + 1]) / 2, (sa_vals[s_i] + sa_vals[s_i + 1]) / 2)
            score_val <- score_mu_sa(
              mid_pt, data_rand, priv_noise, obs_stats, depth_type,
              aux_dr, aux_pn, this_cache, current_ep, c_type
            )
            
            if (is.finite(score_val) && (-score_val) >= acceptance_threshold) {
              total_area <- total_area + (step_mu * step_sa)
              if (
                population_mu >= mu_vals[m_i] && population_mu <= mu_vals[m_i + 1] &&
                population_sigma >= sa_vals[s_i] && population_sigma <= sa_vals[s_i + 1]
              ) {
                covered <- 1
              }
            }
          }
        }
        
        row_out <- c(
          row_out,
          setNames(
            c(covered, cov_mu, cov_sigma, width_mu, width_sigma, total_area),
            paste0(depth_type, c("_cov_joint", "_cov_mu", "_cov_sigma", "_width_mu", "_width_sigma", "_area"))
          )
        )
      }
      
      # Evaluate PB-ADI
      ADI_meansd <- solve_meanstd_from_clamp(obs_stats, current_ep, c_type)
      
      if (any(is.na(ADI_meansd))) {
        row_out <- c(
          row_out,
          setNames(
            rep(NA_real_, 6),
            c("PB_ADI_cov_mu", "PB_ADI_cov_sigma", "PB_ADI_width_mu", "PB_ADI_width_sigma",
              "PB_ADI_cov_joint_box", "PB_ADI_area_box")
          )
        )
      } else {
        clean_means_vars_new <- sapply(
          seq_len(B_paramboot),
          function(x) clean_clamp_meanvar(rnorm(n = n, mean = ADI_meansd[1], sd = ADI_meansd[2]), c_type)
        )
        
        sensitivity_adi      <- (upper_clamp - lower_clamp) / n
        sd_of_noise_mean_adi <- sensitivity_adi / current_ep
        sensitivity_var_adi  <- (upper_clamp - lower_clamp)^2 / n
        sd_of_noise_var_adi  <- sensitivity_var_adi / current_ep
        
        privacy_noises_new <- matrix(rnorm(2 * B_paramboot), ncol = 2, nrow = B_paramboot)
        privacy_noises_new <- t(t(privacy_noises_new) * c(sd_of_noise_mean_adi, sd_of_noise_var_adi))
        noisy_means_vars_new <- t(clean_means_vars_new) + privacy_noises_new
        
        ADI_meansd_pb <- t(
          sapply(seq_len(B_paramboot), function(i) solve_meanstd_from_clamp(noisy_means_vars_new[i, ], current_ep, c_type))
        )
        ADI_meansd_pb <- ADI_meansd_pb[apply(ADI_meansd_pb, 1, function(x) all(is.finite(x))), , drop = FALSE]
        
        if (nrow(ADI_meansd_pb) < floor(0.8 * B_paramboot)) {
          row_out <- c(
            row_out,
            setNames(
              rep(NA_real_, 6),
              c("PB_ADI_cov_mu", "PB_ADI_cov_sigma", "PB_ADI_width_mu", "PB_ADI_width_sigma",
                "PB_ADI_cov_joint_box", "PB_ADI_area_box")
            )
          )
        } else {
          CI_mean_ends <- quantile(2 * ADI_meansd[1] - ADI_meansd_pb[, 1], probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
          CI_std_ends  <- quantile(pmax(0, 2 * ADI_meansd[2] - ADI_meansd_pb[, 2]), probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
          
          pb_cov_mean   <- as.numeric(CI_mean_ends[1] <= population_mu && CI_mean_ends[2] >= population_mu)
          pb_cov_sd     <- as.numeric(CI_std_ends[1] <= population_sigma^2 && CI_std_ends[2] >= population_sigma^2)
          pb_width_mean <- CI_mean_ends[2] - CI_mean_ends[1]
          pb_width_sd   <- CI_std_ends[2] - CI_std_ends[1]
          
          pb_cov_joint_box <- as.numeric(pb_cov_mean == 1 && pb_cov_sd == 1)
          pb_area_box <- pb_width_mean * pb_width_sd
          
          row_out <- c(
            row_out,
            setNames(
              c(pb_cov_mean, pb_cov_sd, pb_width_mean, pb_width_sd, pb_cov_joint_box, pb_area_box),
              c("PB_ADI_cov_mu", "PB_ADI_cov_sigma", "PB_ADI_width_mu", "PB_ADI_width_sigma",
                "PB_ADI_cov_joint_box", "PB_ADI_area_box")
            )
          )
        }
      }
      
      row_out
    }
    
    try(close(pb), silent = TRUE)
    
    # Summary statistics
    col_mean_se <- function(x) {
      m <- mean(x, na.rm = TRUE)
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      c(mean = m, se = se)
    }
    
    summary_rows <- list()
    for (m in c("mahalanobis", "efficient_mu")) {
      cj   <- col_mean_se(sim_results[, paste0(m, "_cov_joint")])
      cmu  <- col_mean_se(sim_results[, paste0(m, "_cov_mu")])
      csig <- col_mean_se(sim_results[, paste0(m, "_cov_sigma")])
      wmu  <- col_mean_se(sim_results[, paste0(m, "_width_mu")])
      wsig <- col_mean_se(sim_results[, paste0(m, "_width_sigma")])
      ar   <- col_mean_se(sim_results[, paste0(m, "_area")])
      
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        clamp_type   = c_type,
        mu_star      = population_mu,
        epsilon      = current_ep,
        method       = ifelse(m == "mahalanobis", "Mahalanobis (Repro)", "Efficient (mu-interest)"),
        cov_joint    = cj["mean"],
        cov_joint_se = cj["se"],
        cov_mu       = cmu["mean"],
        cov_mu_se    = cmu["se"],
        cov_sigma    = csig["mean"],
        cov_sigma_se = csig["se"],
        width_mu     = wmu["mean"],
        width_mu_se  = wmu["se"],
        width_sigma  = wsig["mean"],
        width_sigma_se = wsig["se"],
        area         = ar["mean"],
        area_se      = ar["se"]
      )
    }
    
    pb_cmu    <- col_mean_se(sim_results[, "PB_ADI_cov_mu"])
    pb_csig   <- col_mean_se(sim_results[, "PB_ADI_cov_sigma"])
    pb_cjoint <- col_mean_se(sim_results[, "PB_ADI_cov_joint_box"])
    pb_wmu    <- col_mean_se(sim_results[, "PB_ADI_width_mu"])
    pb_wsig   <- col_mean_se(sim_results[, "PB_ADI_width_sigma"])
    pb_area   <- col_mean_se(sim_results[, "PB_ADI_area_box"])
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      clamp_type   = c_type,
      mu_star      = population_mu,
      epsilon      = current_ep,
      method       = "PB-ADI",
      cov_joint    = pb_cjoint["mean"],
      cov_joint_se = pb_cjoint["se"],
      cov_mu       = pb_cmu["mean"],
      cov_mu_se    = pb_cmu["se"],
      cov_sigma    = pb_csig["mean"],
      cov_sigma_se = pb_csig["se"],
      width_mu     = pb_wmu["mean"],
      width_mu_se  = pb_wmu["se"],
      width_sigma  = pb_wsig["mean"],
      width_sigma_se = pb_wsig["se"],
      area         = pb_area["mean"],
      area_se      = pb_area["se"]
    )
    
    cell_table <- do.call(rbind, summary_rows)
    cell_table[, -(1:4)] <- round(cell_table[, -(1:4)], result_digits)
    print(cell_table, row.names = FALSE)
    
    all_master_rows[[length(all_master_rows) + 1]] <- cell_table
    write.csv(cell_table, file.path(RESULTS_DIR, paste0("softclamp_mu3.0_ep", current_ep, "_", c_type, ".csv")), row.names = FALSE)
  }
}

try(stopCluster(cl), silent = TRUE)

# ------------------------------------------------------------------------------
# 10. Master Summary Table Export
# ------------------------------------------------------------------------------

master_summary_table <- do.call(rbind, all_master_rows)
write.csv(master_summary_table, file.path(RESULTS_DIR, "master_softclamp_summary.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("Soft Clamp Study Complete! Results saved to:\n", file.path(RESULTS_DIR, "master_softclamp_summary.csv"), "\n")
cat("Total Elapsed Time:", format(Sys.time() - start_time_total), "\n")
cat("============================================================\n")

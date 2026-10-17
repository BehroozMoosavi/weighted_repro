# MH, PW and PB-ADI for the two-target private linear regression
# beta = (beta1, beta0), eta = (mu, tau, sigma); n = 100, 200, 500, 1000
rm(list = ls())

N_VALUES  <- c(100L, 200L, 500L, 1000L)
METHODS   <- c("MH", "PW", "PB_ADI")
REPS      <- 1000L
PILOT_RUN <- FALSE
PILOT_REPS <- 2L
MAX_CORES <- 124L
BASE_DIR  <- path.expand("~/R_Simuls/linear_regression_two_target/MH_PW_PBADI")

theta_truth <- c(1, -0.5, 0.5, 1, 0.5)
beta_true   <- theta_truth[1:2]
alpha <- 0.05; delta <- 2; D2 <- delta^2; ep_split <- 1 / sqrt(5)
R_synthetic <- 200L; R_aux <- 400L; h_fd <- 1e-3
rank_threshold <- floor(alpha * (R_synthetic + 1L)) + 1L
nuis_lo <- c(-5, 0.001, 0.001); nuis_hi <- c(5, 5, 5); neutral <- c(0, delta / 2, delta / 2)
MAXIT <- 100L; SIG_ROWS <- 50L; CAL_MAXIT <- 60L; SIG_TOL <- 0.02
SIG_MULT <- c(0.85, 1.18); ETA_PERT <- c(-1, 1); X_CACHE_MAX <- 32L
RIDGE_COV <- 1e-8; RIDGE_INFO <- 1e-10
R_ADI <- 50L; B_BOOT <- 200L; MIN_BOOT_OK <- floor(0.8 * B_BOOT); ADI_MAXIT <- 100L
BOOT_WARM_START <- FALSE   # TRUE = start bootstrap fits at theta_hat (original PB-ADI code)
TH_LO <- c(-5, -5, -5, 0.001, 0.001); TH_HI <- c(5, 5, 5, 5, 5)

# box half-width around beta_true, coarse step, fine step
GRID <- list(`100` = c(2.0, 0.20, 0.05), `200` = c(1.2, 0.10, 0.025),
             `500` = c(0.8, 0.10, 0.02), `1000` = c(0.4, 0.05, 0.01))
MARGIN <- 2L

tag_of <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, compress = FALSE, version = 2)
  substr(unname(tools::md5sum(f)), 1L, 8L)
}
COMMON <- list(theta_truth, alpha, delta, ep_split)
REPRO  <- list(R_synthetic, nuis_lo, nuis_hi, neutral, MAXIT, SIG_ROWS, CAL_MAXIT, SIG_TOL,
               SIG_MULT, ETA_PERT, GRID, MARGIN)
TAG <- list(MH = tag_of(list(COMMON, REPRO)),
            PW = tag_of(list(COMMON, REPRO, R_aux, h_fd, RIDGE_COV, RIDGE_INFO)),
            PB_ADI = paste0(if (BOOT_WARM_START) "warm_" else "",
                            tag_of(list(COMMON, R_ADI, B_BOOT, MIN_BOOT_OK, ADI_MAXIT,
                                        TH_LO, TH_HI, BOOT_WARM_START))))

REPRO_FIELDS <- c("rep", "cov_oracle", "cov_matched", "cov_matched_strict", "cov_refined",
                  "unres_matched", "unres_refined", "area_lower", "area_upper",
                  "n_accepted", "n_unresolved", "boundary_hit", "window_hit",
                  "cal_ok", "tau_hat", "seconds")
PB_FIELDS <- c("rep", "adi_ok", "adi_q", "pb_ok", "covered", "area", "n_boot",
               "boot_q_frac", "seconds")
Q_TOL <- 0.01

for (p in c("foreach", "doSNOW", "parallelly"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(foreach); library(doSNOW) })
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")

# ---------------------------------------------------------------------- shared
clamp <- function(x, a, b) pmin(pmax(x, a), b)
clamp_nuis <- function(x) pmin(pmax(x, nuis_lo), nuis_hi)
safe_inv <- function(A, ridge = 1e-8) {
  if (any(!is.finite(A))) return(NULL)
  e <- tryCatch(eigen(0.5 * (A + t(A)), symmetric = TRUE), error = function(e) NULL)
  if (is.null(e) || any(!is.finite(e$values)) || any(!is.finite(e$vectors))) return(NULL)
  e$vectors %*% (t(e$vectors) / pmax(e$values, ridge))
}
cp_int <- function(k, m) {
  if (m <= 0) return(c(NA_real_, NA_real_))
  c(if (k == 0) 0 else qbeta(alpha / 2, k, m - k + 1),
    if (k == m) 1 else qbeta(1 - alpha / 2, k + 1, m - k))
}
cells <- function(lo, hi, h) {
  z <- (hi - lo) / h
  if (abs(z - round(z)) > 1e-9) stop("grid step does not tile the box")
  lo + (seq_len(round(z)) - 0.5) * h
}
rep_seed <- function(rep, n) as.integer(10000000L + 10000L * n + rep)
noise_sd <- function(n) c(2 * delta, 2 * delta, D2, D2, 2 * D2) / (n * ep_split)

release <- function(ux, uy, noise, th, nsd) {
  x <- ux * th[4] + th[3]; y <- th[2] + th[1] * x + uy * th[5]
  xc <- clamp(x, -delta, delta); yc <- clamp(y, -delta, delta)
  cbind(rowMeans(xc), rowMeans(yc), rowMeans(xc * xc), rowMeans(yc * yc),
        rowMeans(clamp(x * y, -D2, D2))) + sweep(noise, 2L, nsd, "*")
}
one_release <- function(th, n, nsd) {
  ux <- matrix(rnorm(n), 1L); uy <- matrix(rnorm(n), 1L); nz <- matrix(rnorm(5L), 1L)
  as.numeric(release(ux, uy, nz, th, nsd))
}

# ------------------------------------------------------------------ MH and PW
draw_rep <- function(rep, n, aux) {
  set.seed(rep_seed(rep, n)); nsd <- noise_sd(n)
  d <- list(s = one_release(theta_truth, n, nsd))
  d$ux <- matrix(rnorm(R_synthetic * n), R_synthetic)
  d$uy <- matrix(rnorm(R_synthetic * n), R_synthetic)
  d$noise <- matrix(rnorm(R_synthetic * 5L), R_synthetic)
  if (aux) {
    d$ux_a <- matrix(rnorm(R_aux * n), R_aux)
    d$uy_a <- matrix(rnorm(R_aux * n), R_aux)
    d$noise_a <- matrix(rnorm(R_aux * 5L), R_aux)
  }
  m <- seq_len(SIG_ROWS)
  d$ux_s <- d$ux[m, , drop = FALSE]; d$uy_s <- d$uy[m, , drop = FALSE]
  d$noise_s <- d$noise[m, , drop = FALSE]
  d
}

X_at <- function(ctx, mu, tau, cache = TRUE) {
  key <- paste(sprintf("%a", c(mu, tau)), collapse = "_"); e <- ctx$xc
  if (cache && !is.null(v <- e[[key]])) return(v)
  x <- ctx$d$ux * tau + mu; xc <- clamp(x, -delta, delta)
  v <- list(x = x,
            c1 = rowMeans(xc) + ctx$nsd[1] * ctx$d$noise[, 1],
            c3 = rowMeans(xc * xc) + ctx$nsd[3] * ctx$d$noise[, 3],
            xs = x[seq_len(SIG_ROWS), , drop = FALSE])
  if (cache && e$.n < X_CACHE_MAX) { assign(key, v, envir = e); e$.n <- e$.n + 1L }
  v
}

synth <- function(ctx, th, cache) {
  X <- X_at(ctx, th[3], th[4], cache); nz <- ctx$d$noise; s <- ctx$nsd
  y <- th[2] + th[1] * X$x + ctx$d$uy * th[5]; yc <- clamp(y, -delta, delta)
  cbind(X$c1, rowMeans(yc) + s[2] * nz[, 2], X$c3, rowMeans(yc * yc) + s[4] * nz[, 4],
        rowMeans(clamp(X$x * y, -D2, D2)) + s[5] * nz[, 5])
}

pw_direction <- function(ctx, th) {
  d <- ctx$d; nsd <- ctx$nsd
  aux_mean <- function(t) {
    Z <- release(d$ux_a, d$uy_a, d$noise_a, t, nsd)
    if (any(!is.finite(Z))) NULL else colMeans(Z)
  }
  Z0 <- release(d$ux_a, d$uy_a, d$noise_a, th, nsd)
  if (any(!is.finite(Z0))) return(NULL)
  m0 <- colMeans(Z0); Zc <- sweep(Z0, 2L, m0)
  S <- crossprod(Zc) / nrow(Z0); S <- 0.5 * (S + t(S)) + RIDGE_COV * diag(5L)
  Si <- safe_inv(S, RIDGE_COV); if (is.null(Si)) return(NULL)
  J <- matrix(NA_real_, 5L, 5L)
  for (j in 1:5) {
    h <- h_fd * max(1, abs(th[j]))
    if (j >= 4L) { room <- th[j] - nuis_lo[j - 2L]; if (room > 0) h <- min(h, 0.45 * room) }
    if (!is.finite(h) || h <= 0) return(NULL)
    tp <- th; tm <- th; tp[j] <- th[j] + h; tm[j] <- th[j] - h
    mp <- aux_mean(tp); if (is.null(mp)) return(NULL)
    if (j >= 4L && tm[j] <= 0) { J[, j] <- (mp - m0) / h; next }
    mm <- aux_mean(tm); if (is.null(mm)) return(NULL)
    J[, j] <- (mp - mm) / (2 * h)
  }
  Jb <- J[, 1:2, drop = FALSE]; Je <- J[, 3:5, drop = FALSE]
  Ii <- safe_inv(crossprod(Je, Si %*% Je), RIDGE_INFO); if (is.null(Ii)) return(NULL)
  A <- Si %*% (Jb - Je %*% (Ii %*% crossprod(Je, Si %*% Jb)))
  if (is.null(safe_inv(crossprod(A, S %*% A), RIDGE_INFO))) return(NULL)
  A
}

score <- function(ctx, th, cache = TRUE) {
  if (any(!is.finite(th)) || th[4] <= 0 || th[5] <= 0) return(NA_real_)
  if (ctx$method == "PW") { A <- pw_direction(ctx, th); if (is.null(A)) return(NA_real_) }
  sy <- synth(ctx, th, cache); if (any(!is.finite(sy))) return(NA_real_)
  full <- rbind(sy, ctx$d$s); N <- nrow(full)
  Cc <- sweep(full, 2L, colMeans(full))
  S <- crossprod(Cc) / N; S <- 0.5 * (S + t(S)) + 1e-10 * diag(5L)
  Si <- safe_inv(S, 1e-10); if (is.null(Si)) return(NA_real_)
  Q <- rowSums((Cc %*% Si) * Cc)
  if (ctx$method == "PW") {
    Ei <- safe_inv(crossprod(A, S %*% A), RIDGE_INFO); if (is.null(Ei)) return(NA_real_)
    P <- Cc %*% A
    Q <- rowSums((P %*% Ei) * P) + ctx$lambda * Q
  }
  if (any(!is.finite(Q))) return(NA_real_)
  D <- 1 / (1 + Q)
  -(sum(D <= D[N]) + 0.5 * D[N])
}

calibrate <- function(ctx) {
  d <- ctx$d; s <- d$s; nsd <- ctx$nsd
  obj <- function(p) {
    if (any(!is.finite(p)) || p[2] <= 0) return(1e12)
    xc <- clamp(d$ux_s * p[2] + p[1], -delta, delta)
    z <- (mean(rowMeans(xc) + nsd[1] * d$noise_s[, 1]) - s[1])^2 +
      (mean(rowMeans(xc * xc) + nsd[3] * d$noise_s[, 3]) - s[3])^2
    if (is.finite(z)) z else 1e12
  }
  vx <- s[3] - s[1]^2; if (!is.finite(vx) || vx <= 0) vx <- 1
  st <- clamp(c(s[1], sqrt(max(vx, 1e-6))), nuis_lo[1:2], nuis_hi[1:2])
  r <- tryCatch(optim(st, obj, method = "L-BFGS-B", lower = nuis_lo[1:2],
                      upper = nuis_hi[1:2], control = list(maxit = CAL_MAXIT)),
                error = function(e) NULL)
  if (is.null(r) || !is.finite(r$value) || r$value >= 1e11)
    return(list(par = st, resid = sqrt(obj(st)), ok = 0))
  list(par = clamp(r$par, nuis_lo[1:2], nuis_hi[1:2]), resid = sqrt(r$value), ok = 1)
}

solve_sigma <- function(ctx, tg, X) {
  d <- ctx$d; s <- d$s; nsd <- ctx$nsd; xs <- X$xs
  obj <- function(sa) {
    if (!is.finite(sa) || sa <= 0) return(1e12)
    y <- tg[2] + tg[1] * xs + d$uy_s * sa; yc <- clamp(y, -delta, delta)
    z <- (mean(rowMeans(yc) + nsd[2] * d$noise_s[, 2]) - s[2])^2 +
      (mean(rowMeans(yc * yc) + nsd[4] * d$noise_s[, 4]) - s[4])^2 +
      (mean(rowMeans(clamp(xs * y, -D2, D2)) + nsd[5] * d$noise_s[, 5]) - s[5])^2
    if (is.finite(z)) z else 1e12
  }
  r <- tryCatch(optimize(obj, c(nuis_lo[3], nuis_hi[3]), tol = SIG_TOL), error = function(e) NULL)
  if (is.null(r) || !is.finite(r$minimum)) NA_real_ else clamp(r$minimum, nuis_lo[3], nuis_hi[3])
}

accept <- function(ctx, tg, refine = FALSE, prior = NULL) {
  if (isTRUE(prior)) return(TRUE)
  mt <- ctx$mu_tau
  sa <- solve_sigma(ctx, tg, X_at(ctx, mt[1], mt[2])); if (!is.finite(sa)) sa <- delta / 2
  dd <- 1 / sqrt(ctx$n)
  st <- rbind(c(mt, sa), c(mt, sa * SIG_MULT[1]), c(mt, sa * SIG_MULT[2]),
              c(mt[1], mt[2] + dd * ETA_PERT[1], sa), c(mt[1], mt[2] + dd * ETA_PERT[2], sa),
              neutral)
  st <- unique(round(t(apply(st, 1L, clamp_nuis)), 10))
  if (is.null(prior)) {
    resolved <- FALSE; failed <- FALSE
    for (k in seq_len(nrow(st))) {
      v <- -score(ctx, c(tg, st[k, ]))
      if (is.finite(v)) { resolved <- TRUE; if (v >= rank_threshold) return(TRUE) } else failed <- TRUE
    }
  } else {
    failed <- is.na(prior); resolved <- !failed
  }
  if (refine) {
    obj <- function(e) { z <- score(ctx, c(tg, e), FALSE); if (is.finite(z)) z else 1e12 }
    for (k in seq_len(min(3L, nrow(st)))) {
      r <- tryCatch(optim(st[k, ], obj, method = "L-BFGS-B", lower = nuis_lo, upper = nuis_hi,
                          control = list(maxit = MAXIT)), error = function(e) NULL)
      if (is.null(r) || !is.finite(r$value) || r$value >= 1e11) { failed <- TRUE; next }
      v <- -score(ctx, c(tg, r$par), FALSE)
      if (!is.finite(v)) { failed <- TRUE; next }
      resolved <- TRUE
      if (v >= rank_threshold) return(TRUE)
    }
  }
  if (failed) NA else if (resolved) FALSE else NA
}

enc <- function(a) if (isTRUE(a)) 1L else if (identical(a, FALSE)) 0L else -1L

region <- function(ctx) {
  g <- GRID[[as.character(ctx$n)]]; hw <- g[1]; hc <- g[2]; hf <- g[3]
  B1 <- beta_true[1] + c(-hw, hw); B0 <- beta_true[2] + c(-hw, hw)
  C1 <- cells(B1[1], B1[2], hc); C0 <- cells(B0[1], B0[2], hc)
  F1 <- cells(B1[1], B1[2], hf); F0 <- cells(B0[1], B0[2], hf)
  ev <- function(b1, b0) {
    m <- matrix(NA_integer_, length(b0), length(b1))
    for (i in seq_along(b0)) for (j in seq_along(b1)) m[i, j] <- enc(accept(ctx, c(b1[j], b0[i])))
    m
  }
  ins <- which(ev(C1, C0) != 0L, arr.ind = TRUE)
  if (nrow(ins)) {
    mg <- MARGIN * hc
    lo1 <- max(B1[1], min(C1[ins[, 2]]) - hc / 2 - mg); hi1 <- min(B1[2], max(C1[ins[, 2]]) + hc / 2 + mg)
    lo0 <- max(B0[1], min(C0[ins[, 1]]) - hc / 2 - mg); hi0 <- min(B0[2], max(C0[ins[, 1]]) + hc / 2 + mg)
  } else {
    lo1 <- B1[1]; hi1 <- B1[2]; lo0 <- B0[1]; hi0 <- B0[2]
  }
  f1 <- F1[F1 - hf / 2 >= lo1 - 1e-9 & F1 + hf / 2 <= hi1 + 1e-9]
  f0 <- F0[F0 - hf / 2 >= lo0 - 1e-9 & F0 + hf / 2 <= hi0 + 1e-9]
  if (!length(f1)) f1 <- F1
  if (!length(f0)) f0 <- F0
  st <- ev(f1, f0); acc <- st == 1L; inside <- acc | st == -1L
  e1 <- c(min(f1), max(f1)) + c(-1, 1) * hf / 2; e0 <- c(min(f0), max(f0)) + c(-1, 1) * hf / 2
  hit <- (abs(e1[1] - B1[1]) < 1e-9 && any(inside[, 1])) ||
    (abs(e1[2] - B1[2]) < 1e-9 && any(inside[, ncol(inside)])) ||
    (abs(e0[1] - B0[1]) < 1e-9 && any(inside[1, ])) ||
    (abs(e0[2] - B0[2]) < 1e-9 && any(inside[nrow(inside), ]))
  win <- (any(inside[, 1]) && e1[1] > B1[1] + 1e-9) ||
    (any(inside[, ncol(inside)]) && e1[2] < B1[2] - 1e-9) ||
    (any(inside[1, ]) && e0[1] > B0[1] + 1e-9) ||
    (any(inside[nrow(inside), ]) && e0[2] < B0[2] - 1e-9)
  c(area_lower = hf^2 * sum(acc), area_upper = hf^2 * sum(inside),
    n_accepted = sum(acc), n_unresolved = sum(st == -1L),
    boundary_hit = as.numeric(hit), window_hit = as.numeric(win))
}

run_repro <- function(rep, n, method) {
  t0 <- proc.time()[[3]]
  ctx <- new.env()
  ctx$n <- n; ctx$method <- method; ctx$nsd <- noise_sd(n); ctx$lambda <- 1 / log(n)
  ctx$d <- draw_rep(rep, n, method == "PW")
  ctx$xc <- new.env(); ctx$xc$.n <- 0L
  vo <- -score(ctx, theta_truth)
  cal <- calibrate(ctx); ctx$mu_tau <- cal$par
  tm <- accept(ctx, beta_true)
  tr <- accept(ctx, beta_true, refine = TRUE, prior = tm)
  out <- c(rep = rep, cov_oracle = if (is.finite(vo)) as.numeric(vo >= rank_threshold) else NA_real_,
           cov_matched = as.numeric(!identical(tm, FALSE)),
           cov_matched_strict = as.numeric(isTRUE(tm)),
           cov_refined = as.numeric(!identical(tr, FALSE)),
           unres_matched = as.numeric(is.na(tm)), unres_refined = as.numeric(is.na(tr)),
           region(ctx), cal_ok = cal$ok, tau_hat = unname(cal$par[2]),
           seconds = proc.time()[[3]] - t0)
  out[REPRO_FIELDS]
}

# -------------------------------------------------------------------- PB-ADI
safe_cov <- function(X, ridge) { S <- cov(X); 0.5 * (S + t(S)) + ridge * diag(ncol(S)) }
inv_floor <- function(A, ridge) {
  e <- eigen(0.5 * (A + t(A)), symmetric = TRUE)
  e$vectors %*% (t(e$vectors) / pmax(e$values, ridge))
}
adi_obj <- function(th, ux, uy, nz, s, nsd) {
  Q <- tryCatch({
    sy <- release(ux, uy, nz, th, nsd); v <- s - colMeans(sy)
    as.numeric(t(v) %*% inv_floor(safe_cov(sy, 1e-8), 1e-8) %*% v)
  }, error = function(e) Inf)
  if (is.finite(Q)) Q else 1e20
}
adi_init <- function(s, n) {
  t2 <- max(1e-12, s[3] - s[1]^2)
  b1 <- (s[5] - s[1] * s[2]) / t2
  b0 <- (s[2] * s[3] - s[1] * s[5]) / t2
  s2 <- max(1e-12, (s[4] - 2 * b0 * s[2] - 2 * b1 * s[5] + b0^2 + 2 * b0 * b1 * s[1] +
                      b1^2 * s[3]) * n / (n - 2))
  clamp(c(b1, b0, s[1], sqrt(t2 * n / (n - 1)), sqrt(s2)), TH_LO, TH_HI)
}
theta_adi <- function(s, n, nsd, init) {
  ux <- matrix(rnorm(R_ADI * n), R_ADI); uy <- matrix(rnorm(R_ADI * n), R_ADI)
  nz <- matrix(rnorm(R_ADI * 5L), R_ADI)
  best <- Inf; par <- rep(NA_real_, 5L); code <- NA_integer_
  for (st in list(init, c(0, 0, 0, 1, 1))) {
    st <- clamp(st, TH_LO, TH_HI)
    v <- adi_obj(st, ux, uy, nz, s, nsd)
    if (v < best) { best <- v; par <- st; code <- NA_integer_ }
    o <- tryCatch(optim(st, adi_obj, method = "L-BFGS-B", lower = TH_LO, upper = TH_HI,
                        control = list(maxit = ADI_MAXIT),
                        ux = ux, uy = uy, nz = nz, s = s, nsd = nsd),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && all(is.finite(o$par)) && o$value < best) {
      best <- o$value; par <- o$par; code <- o$convergence
    }
  }
  list(theta = par, ok = all(is.finite(par)), code = code, value = best)
}

run_pb <- function(rep, n) {
  t0 <- proc.time()[[3]]
  set.seed(rep_seed(rep, n)); nsd <- noise_sd(n)
  s <- one_release(theta_truth, n, nsd)
  out <- c(rep = rep, adi_ok = 0, adi_q = NA, pb_ok = 0, covered = NA, area = NA,
           n_boot = 0, boot_q_frac = NA, seconds = NA)
  fit <- tryCatch(theta_adi(s, n, nsd, adi_init(s, n)), error = function(e) NULL)
  if (!is.null(fit) && fit$ok) {
    out["adi_ok"] <- 1; th <- fit$theta
    out["adi_q"] <- fit$value
    bt <- matrix(NA_real_, B_BOOT, 5L); qv <- rep(NA_real_, B_BOOT)
    for (b in seq_len(B_BOOT)) {
      sb <- one_release(th, n, nsd)
      f <- tryCatch(theta_adi(sb, n, nsd, if (BOOT_WARM_START) th else adi_init(sb, n)),
                    error = function(e) NULL)
      if (!is.null(f) && f$ok) bt[b, ] <- f$theta
      if (!is.null(f)) qv[b] <- f$value
    }
    out["boot_q_frac"] <- mean(qv > Q_TOL, na.rm = TRUE)
    bt <- bt[stats::complete.cases(bt), , drop = FALSE]
    out["n_boot"] <- nrow(bt)
    if (nrow(bt) >= MIN_BOOT_OK) {
      E <- sweep(bt[, 1:2, drop = FALSE], 2L, th[1:2])
      S <- safe_cov(E, 1e-10); Si <- inv_floor(S, 1e-10)
      q <- rowSums((E %*% Si) * E); q <- q[is.finite(q)]
      if (length(q) >= MIN_BOOT_OK) {
        cut <- sort(q)[min(length(q), ceiling((length(q) + 1) * (1 - alpha)))]
        dv <- beta_true - th[1:2]; qt <- as.numeric(t(dv) %*% Si %*% dv)
        ds <- det(S)
        area <- if (is.finite(ds) && ds > 0 && is.finite(cut)) pi * cut * sqrt(ds) else NA_real_
        out["covered"] <- as.numeric(is.finite(qt) && qt <= cut)
        out["area"] <- area
        out["pb_ok"] <- as.numeric(is.finite(area))
      }
    }
  }
  out["seconds"] <- proc.time()[[3]] - t0
  out[PB_FIELDS]
}

# ------------------------------------------------------------------ summaries
summarize <- function(res, m, n, mins) {
  if (m == "PB_ADI") {
    ok <- is.finite(res$covered) & is.finite(res$area); u <- sum(ok)
    k <- sum(res$covered[ok]); ci <- cp_int(k, u); a <- res$area[ok]
    data.frame(Method = m, n = n, Reps_Used = u, Coverage = k / u, CP_Lo = ci[1], CP_Hi = ci[2],
               Cov_Strict = k / nrow(res),
               Mean_Area_L = mean(a), Mean_Area = mean(a), SE_Area = sd(a) / sqrt(length(a)),
               ADI_Q_Bad = mean(res$adi_q > Q_TOL, na.rm = TRUE),
               Boot_Q_Frac = mean(res$boot_q_frac, na.rm = TRUE),
               PB_Fail = mean(res$pb_ok == 0, na.rm = TRUE),
               Task_Errors = sum(is.na(res$seconds)),
               Sec_Per_Rep = mean(res$seconds, na.rm = TRUE), Wall_Min = mins)
  } else {
    ok <- is.finite(res$cov_matched) & is.finite(res$area_upper); u <- sum(ok)
    km <- sum(res$cov_matched[ok]); ci <- cp_int(km, u)
    o <- res$cov_oracle[is.finite(res$cov_oracle)]; oc <- cp_int(sum(o), length(o))
    aL <- res$area_lower[ok]; aU <- res$area_upper[ok]
    data.frame(Method = m, n = n, Reps_Used = u, Coverage = km / u, CP_Lo = ci[1], CP_Hi = ci[2],
               Cov_Strict = sum(res$cov_matched_strict == 1, na.rm = TRUE) / nrow(res),
               Cov_Refined = mean(res$cov_refined[ok]),
               Cov_Oracle = mean(o), Oracle_CP_Lo = oc[1], Oracle_CP_Hi = oc[2],
               Unres_Truth = mean(res$unres_matched[ok]),
               Mean_Area_L = mean(aL), Mean_Area = mean(aU), SE_Area = sd(aU) / sqrt(length(aU)),
               Boundary_Rate = mean(res$boundary_hit, na.rm = TRUE),
               Window_Rate = mean(res$window_hit, na.rm = TRUE),
               Unresolved_Cells = mean(res$n_unresolved, na.rm = TRUE),
               Cal_Fail = mean(res$cal_ok == 0, na.rm = TRUE),
               Task_Errors = sum(is.na(res$seconds)),
               Sec_Per_Rep = mean(res$seconds, na.rm = TRUE), Wall_Min = mins)
  }
}

# ------------------------------------------------------------------------ run
dir.create(BASE_DIR, recursive = TRUE, showWarnings = FALSE)
reps <- if (PILOT_RUN) PILOT_REPS else REPS
ncores <- max(1L, min(MAX_CORES, parallelly::availableCores(), reps))
EXPORTS <- ls()

cl <- parallel::makePSOCKcluster(ncores)
registerDoSNOW(cl)
parallel::clusterExport(cl, EXPORTS)

SUMMARY <- list()
for (m in METHODS) for (n in N_VALUES) {
  cat(sprintf("\n%s | n = %d | reps = %d | cores = %d\n", m, n, reps, ncores))
  raw_file <- file.path(BASE_DIR, sprintf("%s_raw_n%d_R%d_%s.csv", m, n, reps, TAG[[m]]))
  need <- if (m == "PB_ADI") PB_FIELDS else REPRO_FIELDS
  old <- if (file.exists(raw_file)) tryCatch(read.csv(raw_file), error = function(e) NULL)
  if (!is.null(old) && nrow(old) == reps && all(need %in% names(old))) {
    s <- summarize(old, m, n, NA_real_)
    SUMMARY[[length(SUMMARY) + 1L]] <- s
    cat("already done, read from", raw_file, "\n"); print(s, row.names = FALSE); next
  }
  parallel::clusterExport(cl, c("m", "n"))
  bar <- txtProgressBar(max = reps, style = 3)
  opts <- list(progress = function(k) setTxtProgressBar(bar, k))
  t0 <- Sys.time()
  res <- foreach(r = seq_len(reps), .combine = rbind, .options.snow = opts,
                 .noexport = EXPORTS) %dopar% {
                   tryCatch(if (m == "PB_ADI") run_pb(r, n) else run_repro(r, n, m),
                            error = function(e) {
                              cat(sprintf("[%s n=%d rep=%d] %s\n", m, n, r, conditionMessage(e)),
                                  file = file.path(BASE_DIR, "errors.log"), append = TRUE)
                              f <- if (m == "PB_ADI") PB_FIELDS else REPRO_FIELDS
                              setNames(c(r, rep(NA_real_, length(f) - 1L)), f)
                            })
                 }
  close(bar)
  res <- as.data.frame(rbind(res), row.names = NULL)
  write.csv(res, raw_file, row.names = FALSE)
  s <- summarize(res, m, n, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  SUMMARY[[length(SUMMARY) + 1L]] <- s
  print(s, row.names = FALSE)
}
parallel::stopCluster(cl)

for (m in METHODS) {
  tab <- do.call(rbind, Filter(function(x) x$Method == m, SUMMARY))
  write.csv(tab, file.path(BASE_DIR, sprintf("%s_summary_%s.csv", m, TAG[[m]])), row.names = FALSE)
  cat(sprintf("\n============================== %s ==============================\n", m))
  print(tab, row.names = FALSE, digits = 4)
}

cols <- c("Method", "n", "Reps_Used", "Coverage", "CP_Lo", "CP_Hi", "Mean_Area", "SE_Area",
          "Cov_Strict", "Mean_Area_L")
COMPARE <- do.call(rbind, lapply(SUMMARY, function(x) x[, cols]))
write.csv(COMPARE, file.path(BASE_DIR, sprintf("comparison_%s.csv", tag_of(TAG))),
          row.names = FALSE)
cat("\n============================== COMPARISON ==============================\n")
print(COMPARE, row.names = FALSE, digits = 4)
cat("\nResults: ", BASE_DIR, "\n", sep = "")

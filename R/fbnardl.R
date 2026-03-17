#' Fourier Bootstrap Nonlinear ARDL Estimation
#'
#' Estimates a nonlinear ARDL (NARDL) model augmented with flexible Fourier
#' terms and tests for cointegration using either PSS bounds-testing critical
#' values (FNARDL) or residual-based bootstrap critical values (FBNARDL).
#'
#' @param y Numeric vector. Dependent variable.
#' @param x_decomp Numeric matrix or vector. Variable(s) to decompose into
#'   positive and negative partial sums.
#' @param x_ctrl Numeric matrix or vector, or \code{NULL}. Additional
#'   (non-decomposed) control regressors. Default is \code{NULL}.
#' @param type Character. \code{"fnardl"} (default) for Fourier NARDL with
#'   PSS/Kripfganz-Schneider critical values, or \code{"fbnardl"} for Fourier
#'   Bootstrap NARDL.
#' @param max_lag Integer. Maximum lag order to search over. Default is
#'   \code{4}.
#' @param max_k Numeric. Maximum Fourier frequency. Default is \code{3}.
#' @param ic Character. Information criterion: \code{"aic"} (default) or
#'   \code{"bic"}.
#' @param reps Integer. Number of bootstrap replications for
#'   \code{type = "fbnardl"}. Default is \code{499}.
#' @param no_fourier Logical. If \code{TRUE}, suppress Fourier terms (pure
#'   NARDL). Default is \code{FALSE}.
#'
#' @return A list of class \code{"fbnardl"} containing:
#' \describe{
#'   \item{type}{Character. Estimator type used.}
#'   \item{coefficients}{Numeric vector. Estimated coefficients.}
#'   \item{se}{Numeric vector. Standard errors.}
#'   \item{tstat}{Numeric vector. t-statistics.}
#'   \item{pval}{Numeric vector. p-values.}
#'   \item{best_p}{Integer. Selected lag order for the dependent variable.}
#'   \item{best_q}{Integer vector. Selected lag orders for decomposed variables.}
#'   \item{best_kstar}{Numeric. Selected Fourier frequency.}
#'   \item{ic_val}{Numeric. Best information criterion value.}
#'   \item{aic}{Numeric. AIC of final model.}
#'   \item{bic}{Numeric. BIC of final model.}
#'   \item{r2}{Numeric. R-squared.}
#'   \item{r2_adj}{Numeric. Adjusted R-squared.}
#'   \item{nobs}{Integer. Observations used.}
#'   \item{Fov}{Numeric. Overall F-statistic for cointegration.}
#'   \item{t_dep}{Numeric. t-statistic on lagged dependent variable (ECM).}
#'   \item{Find}{Numeric. F-statistic on lagged independent level terms.}
#'   \item{lr_pos}{Numeric vector. Long-run multipliers for positive components.}
#'   \item{lr_neg}{Numeric vector. Long-run multipliers for negative components.}
#'   \item{bootstrap}{List. Bootstrap critical values (if \code{type = "fbnardl"}).}
#'   \item{residuals}{Numeric vector. Residuals from the final model.}
#' }
#'
#' @details
#' \strong{Decomposition:}
#' Each variable in \code{x_decomp} is split into positive and negative
#' partial sums: \eqn{x^+ = \sum \max(\Delta x, 0)} and
#' \eqn{x^- = \sum \min(\Delta x, 0)}.
#'
#' \strong{Model selection:}
#' A two-step procedure selects the optimal Fourier frequency \eqn{k^*} in
#' Step 1 (minimum SSR from a maximal model, Yilanci et al. 2020), then
#' searches over all \eqn{(p, q, r)} lag combinations using the specified
#' information criterion in Step 2.
#'
#' \strong{Cointegration test:}
#' Three test statistics are reported — the overall F (Fov), the t-statistic
#' on the lagged dependent variable (t_dep), and the F-statistic on lagged
#' levels of independent variables (Find).  For \code{type = "fbnardl"},
#' bootstrap critical values are computed from residual-based resampling
#' (Bertelli et al. 2022).
#'
#' @references
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In W. C. Horrace & R. C. Sickles (Eds.), \emph{Festschrift in Honor of
#' Peter Schmidt}. Springer. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' Yilanci, V., Bozoklu, S., & Gorus, M. S. (2020). Fourier ARDL approach.
#' \emph{Evaluation Review}, 44(5-6), 431–450.
#' \doi{10.1177/0193841X18823479}
#'
#' Bertelli, S., Vacca, G., & Zoia, M. G. (2022). Bootstrap cointegration
#' tests in ARDL models. \emph{Statistical Methods and Applications}, 31,
#' 1231–1268. \doi{10.1007/s10260-021-00571-7}
#'
#' Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions for
#' critical value bounds and approximate p-values in equilibrium correction
#' models. \emph{Oxford Bulletin of Economics and Statistics}, 82(6),
#' 1456–1481. \doi{10.1111/obes.12376}
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' x <- cumsum(rnorm(n))
#' y <- 0.6 * x + rnorm(n, sd = 0.5)
#' res <- fbnardl(y, x, type = "fnardl", max_lag = 2, max_k = 2)
#' print(res)
#'
#' @export
fbnardl <- function(y,
                    x_decomp,
                    x_ctrl    = NULL,
                    type      = c("fnardl", "fbnardl"),
                    max_lag   = 4L,
                    max_k     = 3,
                    ic        = c("aic", "bic"),
                    reps      = 499L,
                    no_fourier = FALSE) {

  type <- match.arg(type)
  ic   <- match.arg(ic)

  ## ---- Input sanitising ----
  y <- as.numeric(y)
  if (!is.matrix(x_decomp)) x_decomp <- matrix(as.numeric(x_decomp), ncol = 1L)
  if (!is.null(x_ctrl) && !is.matrix(x_ctrl))
    x_ctrl <- matrix(as.numeric(x_ctrl), ncol = 1L)

  max_lag <- as.integer(max_lag)
  reps    <- as.integer(reps)
  n_total <- length(y)

  if (n_total < 30L) stop("Insufficient observations (need >= 30).")
  if (max_lag < 1L)  stop("'max_lag' must be >= 1.")

  n_dec  <- ncol(x_decomp)
  n_ctrl <- if (!is.null(x_ctrl)) ncol(x_ctrl) else 0L

  ## ---- Step 0: partial-sum decomposition ----
  x_pos <- x_neg <- matrix(0, nrow = n_total, ncol = n_dec)
  for (j in seq_len(n_dec)) {
    dx        <- c(0, diff(x_decomp[, j]))
    x_pos[, j] <- cumsum(pmax(dx, 0))
    x_neg[, j] <- cumsum(pmin(dx, 0))
  }

  ## ---- Step 1: select k* by minimum SSR ----
  t_idx <- seq_len(n_total)
  best_kstar <- 0
  best_ssr_k <- Inf

  if (no_fourier) {
    k_grid <- 0
  } else {
    k_grid <- seq(0.1, max_k, by = 0.1)
  }

  for (kv in k_grid) {
    regs <- .fbnardl_build_maximal(y, x_pos, x_neg, x_ctrl,
                                   max_lag, kv, t_idx, n_total)
    if (is.null(regs)) next
    fit  <- tryCatch(lm.fit(regs$X, regs$dY), error = function(e) NULL)
    if (is.null(fit)) next
    ssr_k <- sum(fit$residuals^2)
    if (ssr_k < best_ssr_k) { best_ssr_k <- ssr_k; best_kstar <- kv }
  }

  ## ---- Step 2: select (p, q, r) by IC with fixed k* ----
  best_ic  <- Inf
  best_p   <- 1L
  best_q   <- rep(0L, n_dec)
  best_r   <- if (n_ctrl > 0L) rep(0L, n_ctrl) else integer(0L)
  best_fit <- NULL
  best_X   <- NULL
  best_dY  <- NULL

  for (p in seq_len(max_lag)) {
    q_combos <- expand.grid(rep(list(0:max_lag), n_dec))
    r_combos <- if (n_ctrl > 0L) expand.grid(rep(list(0:max_lag), n_ctrl)) else
      data.frame(dummy = 0L)

    for (qi in seq_len(nrow(q_combos))) {
      q_vec <- as.integer(q_combos[qi, ])
      for (ri in seq_len(nrow(r_combos))) {
        r_vec <- if (n_ctrl > 0L) as.integer(r_combos[ri, ]) else integer(0L)

        regs <- .fbnardl_build_model(y, x_pos, x_neg, x_ctrl,
                                     p, q_vec, r_vec,
                                     best_kstar, t_idx, n_total)
        if (is.null(regs)) next
        fit <- tryCatch(lm.fit(regs$X, regs$dY), error = function(e) NULL)
        if (is.null(fit)) next

        n_eff <- length(regs$dY)
        k_p   <- ncol(regs$X)
        if (n_eff - k_p < 5L) next

        ssr  <- sum(fit$residuals^2)
        ic_v <- if (ic == "aic")
          n_eff * log(ssr / n_eff) + 2 * k_p
        else
          n_eff * log(ssr / n_eff) + k_p * log(n_eff)

        if (ic_v < best_ic) {
          best_ic  <- ic_v
          best_p   <- p
          best_q   <- q_vec
          best_r   <- r_vec
          best_fit <- fit
          best_X   <- regs$X
          best_dY  <- regs$dY
        }
      }
    }
  }

  if (is.null(best_fit))
    stop("No valid model found. Try increasing max_lag or reducing the number of regressors.")

  ## ---- Extract statistics from best model ----
  n_eff  <- length(best_dY)
  k_p    <- ncol(best_X)
  ssr    <- sum(best_fit$residuals^2)
  sigma2 <- ssr / (n_eff - k_p)

  XtX_inv <- tryCatch(solve(crossprod(best_X)), error = function(e) NULL)
  if (is.null(XtX_inv))
    stop("Singular design matrix in best model.")

  coefs  <- as.numeric(XtX_inv %*% crossprod(best_X, best_dY))
  se_vec <- sqrt(diag(XtX_inv) * sigma2)
  tstat  <- coefs / se_vec
  df_r   <- n_eff - k_p
  pval   <- 2 * stats::pt(abs(tstat), df = df_r, lower.tail = FALSE)

  names(coefs) <- names(se_vec) <- names(tstat) <- names(pval) <-
    colnames(best_X)

  aic_val <- n_eff * log(ssr / n_eff) + 2 * k_p
  bic_val <- n_eff * log(ssr / n_eff) + k_p * log(n_eff)
  ss_tot  <- sum((best_dY - mean(best_dY))^2)
  r2      <- 1 - ssr / ss_tot
  r2_adj  <- 1 - (1 - r2) * (n_eff - 1) / (df_r)

  ## ---- Long-run multipliers ----
  ## ECM coefficient: coefficient on L.y (last column block)
  ## Identify column names containing "Ly_lag" and "Lpos_" / "Lneg_"
  cn  <- colnames(best_X)
  ecm_idx <- which(grepl("^Ly_lag$", cn))
  lr_pos <- lr_neg <- numeric(n_dec)
  if (length(ecm_idx) == 1L) {
    alpha <- coefs[ecm_idx]
    for (j in seq_len(n_dec)) {
      pos_idx <- which(grepl(paste0("^Lpos_", j, "$"), cn))
      neg_idx <- which(grepl(paste0("^Lneg_", j, "$"), cn))
      if (length(pos_idx) == 1L && abs(alpha) > .Machine$double.eps)
        lr_pos[j] <- -coefs[pos_idx] / alpha
      if (length(neg_idx) == 1L && abs(alpha) > .Machine$double.eps)
        lr_neg[j] <- -coefs[neg_idx] / alpha
    }
  }

  ## ---- Bounds test (Fov, t_dep, Find) ----
  level_cols <- grep("^(Ly_lag|Lpos_|Lneg_|Lctrl_)", cn)
  indep_cols <- grep("^(Lpos_|Lneg_|Lctrl_)", cn)
  ecm_col    <- grep("^Ly_lag$", cn)

  Fov  <- .fbnardl_F_test(best_X, best_dY, coefs, sigma2, level_cols, df_r)
  Find <- .fbnardl_F_test(best_X, best_dY, coefs, sigma2, indep_cols, df_r)
  t_dep <- if (length(ecm_col) == 1L) tstat[ecm_col] else NA_real_

  ## ---- Bootstrap (type = "fbnardl") ----
  boot_res <- NULL
  if (type == "fbnardl" && reps > 0L) {
    message("Computing bootstrap critical values (", reps, " reps)...")
    boot_res <- .fbnardl_bootstrap(best_X, best_dY, coefs, level_cols,
                                   indep_cols, ecm_col, reps, df_r)
  }

  structure(
    list(
      type        = type,
      coefficients = coefs,
      se          = se_vec,
      tstat       = tstat,
      pval        = pval,
      best_p      = best_p,
      best_q      = best_q,
      best_r      = best_r,
      best_kstar  = best_kstar,
      ic_val      = best_ic,
      aic         = aic_val,
      bic         = bic_val,
      r2          = r2,
      r2_adj      = r2_adj,
      nobs        = n_eff,
      Fov         = Fov,
      t_dep       = t_dep,
      Find        = Find,
      lr_pos      = lr_pos,
      lr_neg      = lr_neg,
      bootstrap   = boot_res,
      residuals   = best_fit$residuals,
      ic          = ic
    ),
    class = "fbnardl"
  )
}


# ============================================================
# Internal: build maximal ARDL(p_max, q_max) for k* selection
# ============================================================
#' @keywords internal
.fbnardl_build_maximal <- function(y, x_pos, x_neg, x_ctrl,
                                   max_lag, kv, t_idx, n_total) {
  n_dec  <- ncol(x_pos)
  n_ctrl <- if (!is.null(x_ctrl)) ncol(x_ctrl) else 0L

  dy <- diff(y)
  T1 <- length(dy)

  ## We need max_lag rows of lead to start
  start <- max_lag + 1L
  if (start > T1) return(NULL)
  idx   <- start:T1    # index in dy (length = T1 - max_lag)
  n_eff <- length(idx)
  if (n_eff < 10L) return(NULL)

  dY <- dy[idx]

  cols <- list()

  ## Lagged differences of y: L1.dy ... Lp.dy
  for (j in seq_len(max_lag))
    cols[[paste0("Ldy_", j)]] <- dy[idx - j]

  ## Decomposed pos/neg: D.x+ D.x- L1.D.x+ ... Lq.D.x+
  for (d in seq_len(n_dec)) {
    dx_pos <- diff(x_pos[, d])
    dx_neg <- diff(x_neg[, d])
    for (j in 0:max_lag) {
      lbl <- if (j == 0) "" else paste0("L", j, ".")
      cols[[paste0(lbl, "Dpos_", d)]] <- dx_pos[idx - j]
      cols[[paste0(lbl, "Dneg_", d)]] <- dx_neg[idx - j]
    }
  }

  ## Controls
  if (n_ctrl > 0L) {
    for (d in seq_len(n_ctrl)) {
      dx_c <- diff(x_ctrl[, d])
      for (j in 0:max_lag) {
        lbl <- if (j == 0) "" else paste0("L", j, ".")
        cols[[paste0(lbl, "Dctrl_", d)]] <- dx_c[idx - j]
      }
    }
  }

  ## Lagged levels (ECM)
  cols[["Ly_lag"]] <- y[idx]   # L.y in ECM context (idx maps to t-1 in original)
  for (d in seq_len(n_dec)) {
    cols[[paste0("Lpos_", d)]] <- x_pos[idx, d]
    cols[[paste0("Lneg_", d)]] <- x_neg[idx, d]
  }
  if (n_ctrl > 0L) {
    for (d in seq_len(n_ctrl))
      cols[[paste0("Lctrl_", d)]] <- x_ctrl[idx, d]
  }

  ## Fourier terms
  if (kv > 0) {
    cols[["sin_k"]] <- sin(2 * pi * kv * t_idx[idx + max_lag] / n_total)
    cols[["cos_k"]] <- cos(2 * pi * kv * t_idx[idx + max_lag] / n_total)
  }

  ## Constant
  cols[["const"]] <- rep(1, n_eff)

  X <- do.call(cbind, cols)
  ok <- complete.cases(cbind(dY, X))
  if (sum(ok) < 10L) return(NULL)
  list(X = X[ok, ], dY = dY[ok])
}


# ============================================================
# Internal: build specific ARDL(p, q, r) model
# ============================================================
#' @keywords internal
.fbnardl_build_model <- function(y, x_pos, x_neg, x_ctrl,
                                  p, q_vec, r_vec,
                                  kstar, t_idx, n_total) {
  n_dec  <- ncol(x_pos)
  n_ctrl <- if (!is.null(x_ctrl)) ncol(x_ctrl) else 0L
  max_lag <- max(p, max(q_vec), if (length(r_vec) > 0) max(r_vec) else 0L)

  dy <- diff(y)
  T1 <- length(dy)
  start <- max_lag + 1L
  if (start > T1) return(NULL)
  idx   <- start:T1
  n_eff <- length(idx)
  if (n_eff < 10L) return(NULL)

  dY <- dy[idx]
  cols <- list()

  for (j in seq_len(p))
    cols[[paste0("Ldy_", j)]] <- dy[idx - j]

  for (d in seq_len(n_dec)) {
    dx_pos <- diff(x_pos[, d])
    dx_neg <- diff(x_neg[, d])
    qi <- q_vec[d]
    for (j in 0:qi) {
      lbl <- if (j == 0) "" else paste0("L", j, ".")
      cols[[paste0(lbl, "Dpos_", d)]] <- dx_pos[idx - j]
      cols[[paste0(lbl, "Dneg_", d)]] <- dx_neg[idx - j]
    }
  }

  if (n_ctrl > 0L) {
    for (d in seq_len(n_ctrl)) {
      dx_c <- diff(x_ctrl[, d])
      rj <- r_vec[d]
      for (j in 0:rj) {
        lbl <- if (j == 0) "" else paste0("L", j, ".")
        cols[[paste0(lbl, "Dctrl_", d)]] <- dx_c[idx - j]
      }
    }
  }

  cols[["Ly_lag"]] <- y[idx]
  for (d in seq_len(n_dec)) {
    cols[[paste0("Lpos_", d)]] <- x_pos[idx, d]
    cols[[paste0("Lneg_", d)]] <- x_neg[idx, d]
  }
  if (n_ctrl > 0L) {
    for (d in seq_len(n_ctrl))
      cols[[paste0("Lctrl_", d)]] <- x_ctrl[idx, d]
  }

  if (kstar > 0) {
    cols[["sin_k"]] <- sin(2 * pi * kstar * t_idx[idx + max_lag] / n_total)
    cols[["cos_k"]] <- cos(2 * pi * kstar * t_idx[idx + max_lag] / n_total)
  }

  cols[["const"]] <- rep(1, n_eff)

  X  <- do.call(cbind, cols)
  ok <- complete.cases(cbind(dY, X))
  if (sum(ok) < 10L) return(NULL)
  list(X = X[ok, ], dY = dY[ok])
}


# ============================================================
# Internal: F-test for a subset of columns (Wald test)
# ============================================================
#' @keywords internal
.fbnardl_F_test <- function(X, dY, coefs, sigma2, col_idx, df_r) {
  if (length(col_idx) == 0L) return(NA_real_)
  R <- matrix(0, nrow = length(col_idx), ncol = ncol(X))
  for (i in seq_along(col_idx)) R[i, col_idx[i]] <- 1
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NA_real_)
  Rb   <- R %*% coefs
  RVRT <- R %*% XtX_inv %*% t(R)
  RVRT_inv <- tryCatch(solve(RVRT), error = function(e) NULL)
  if (is.null(RVRT_inv)) return(NA_real_)
  Fstat <- as.numeric(t(Rb) %*% RVRT_inv %*% Rb) / (length(col_idx) * sigma2)
  Fstat
}


# ============================================================
# Internal: residual-based bootstrap (Bertelli et al. 2022)
# ============================================================
#' @keywords internal
.fbnardl_bootstrap <- function(X, dY, coefs, level_cols,
                                indep_cols, ecm_col, reps, df_r) {
  n_eff  <- nrow(X)
  resids <- as.numeric(dY - X %*% coefs)
  resids <- resids - mean(resids)

  Fov_boot  <- numeric(reps)
  t_boot    <- numeric(reps)
  Find_boot <- numeric(reps)

  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NULL)

  for (b in seq_len(reps)) {
    ## Resample residuals under H0 (unit root: no level terms)
    dY_b   <- X %*% coefs + sample(resids, n_eff, replace = TRUE)
    fit_b  <- tryCatch(lm.fit(X, as.numeric(dY_b)), error = function(e) NULL)
    if (is.null(fit_b)) next

    coefs_b  <- as.numeric(fit_b$coefficients)
    ssr_b    <- sum(fit_b$residuals^2)
    sigma2_b <- ssr_b / df_r

    Fov_boot[b]  <- .fbnardl_F_test(X, dY_b, coefs_b, sigma2_b, level_cols, df_r)
    Find_boot[b] <- .fbnardl_F_test(X, dY_b, coefs_b, sigma2_b, indep_cols, df_r)
    se_b <- sqrt(diag(XtX_inv) * sigma2_b)
    t_boot[b] <- if (length(ecm_col) == 1L && se_b[ecm_col] > 0)
      coefs_b[ecm_col] / se_b[ecm_col] else NA_real_
  }

  ## Remove NAs
  Fov_boot  <- Fov_boot[is.finite(Fov_boot)]
  t_boot    <- t_boot[is.finite(t_boot)]
  Find_boot <- Find_boot[is.finite(Find_boot)]

  list(
    Fov_cv10  = stats::quantile(Fov_boot,  0.90, names = FALSE),
    Fov_cv05  = stats::quantile(Fov_boot,  0.95, names = FALSE),
    Fov_cv01  = stats::quantile(Fov_boot,  0.99, names = FALSE),
    t_cv10    = stats::quantile(t_boot,    0.10, names = FALSE),
    t_cv05    = stats::quantile(t_boot,    0.05, names = FALSE),
    t_cv01    = stats::quantile(t_boot,    0.01, names = FALSE),
    Find_cv10 = stats::quantile(Find_boot, 0.90, names = FALSE),
    Find_cv05 = stats::quantile(Find_boot, 0.95, names = FALSE),
    Find_cv01 = stats::quantile(Find_boot, 0.99, names = FALSE)
  )
}


#' Print Method for fbnardl Objects
#'
#' @param x An object of class \code{"fbnardl"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.fbnardl <- function(x, ...) {
  cat("Fourier Bootstrap Nonlinear ARDL\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Type         : %s\n", toupper(x$type)))
  cat(sprintf("Best p       : %d\n", x$best_p))
  cat(sprintf("Best k*      : %.2f\n", x$best_kstar))
  cat(sprintf("IC (%s)     : %.4f\n", x$ic, x$ic_val))
  cat(sprintf("R-squared    : %.4f\n", x$r2))
  cat(sprintf("Adj R-sq     : %.4f\n", x$r2_adj))
  cat(sprintf("Observations : %d\n", x$nobs))
  cat("\nCointegration Test Statistics:\n")
  cat(sprintf("  Fov  : %.4f\n", x$Fov))
  cat(sprintf("  t_dep: %.4f\n", x$t_dep))
  cat(sprintf("  Find : %.4f\n", x$Find))
  if (!is.null(x$bootstrap)) {
    b <- x$bootstrap
    cat("\nBootstrap Critical Values (10%/5%/1%):\n")
    cat(sprintf("  Fov  : %.3f / %.3f / %.3f\n", b$Fov_cv10, b$Fov_cv05, b$Fov_cv01))
    cat(sprintf("  t    : %.3f / %.3f / %.3f\n", b$t_cv10,   b$t_cv05,   b$t_cv01))
    cat(sprintf("  Find : %.3f / %.3f / %.3f\n", b$Find_cv10, b$Find_cv05, b$Find_cv01))
  }
  if (any(x$lr_pos != 0) || any(x$lr_neg != 0)) {
    cat("\nLong-Run Multipliers:\n")
    for (j in seq_along(x$lr_pos))
      cat(sprintf("  Var %d: LR(+) = %.4f  LR(-) = %.4f\n", j, x$lr_pos[j], x$lr_neg[j]))
  }
  invisible(x)
}

#' Compute Bounds or Bootstrap Cointegration Test
#'
#' @param model Estimated lm object
#' @param data Data frame with decomposed variables
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param control_vars Names of control variables
#' @param type "fnardl" for analytical CV, "fbnardl" for bootstrap
#' @param reps Number of bootstrap replications
#' @param fourier Logical, include Fourier terms
#' @param kstar Optimal Fourier frequency
#' @param p Optimal lag for dependent variable
#' @param q Optimal lags for decomposed variables
#' @param r Optimal lags for control variables
#' @param nobs Number of observations
#' @return List with test statistics and critical values
#' @keywords internal
compute_bounds_test <- function(model, data, depvar, decompose, control_vars,
                                 type, reps, fourier, kstar, p, q, r, nobs) {

  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  df_r <- model$df.residual

  # Identify level terms
  level_vars <- c("L_y")
  for (var in decompose) {
    level_vars <- c(level_vars, paste0("L_", var, "_pos"), paste0("L_", var, "_neg"))
  }
  for (var in control_vars) {
    level_vars <- c(level_vars, paste0("L_", var))
  }

  # Independent level vars (excluding L_y)
  indep_level_vars <- setdiff(level_vars, "L_y")

  # Filter to existing variables
  level_vars <- level_vars[level_vars %in% names(coefs)]
  indep_level_vars <- indep_level_vars[indep_level_vars %in% names(coefs)]

  # Overall F-test (Fov): joint test of all level variables
  if (length(level_vars) > 0) {
    Fov <- compute_joint_f_test(model, level_vars)
  } else {
    Fov <- list(F = NA, df1 = NA, df2 = NA, p = NA)
  }

  # t-test on ECM coefficient (L_y)
  if ("L_y" %in% names(coefs)) {
    t_dep <- coefs["L_y"] / sqrt(vcov_mat["L_y", "L_y"])
    t_dep_p <- 2 * stats::pt(abs(t_dep), df_r, lower.tail = FALSE)
  } else {
    t_dep <- t_dep_p <- NA
  }

  # F-test on independent level variables (Find)
  if (length(indep_level_vars) > 0) {
    Find <- compute_joint_f_test(model, indep_level_vars)
  } else {
    Find <- list(F = NA, df1 = NA, df2 = NA, p = NA)
  }

  # Number of long-run forcing variables (for PSS bounds)
  k <- length(indep_level_vars)

  if (type == "fnardl") {
    # Analytical critical values (Kripfganz & Schneider, 2020)
    cv <- get_pss_critical_values(k)

    result <- list(
      type = "fnardl",
      Fov = Fov$F,
      Fov_df1 = Fov$df1,
      Fov_df2 = Fov$df2,
      Fov_p = Fov$p,
      t_dep = t_dep,
      t_dep_p = t_dep_p,
      Find = Find$F,
      Find_df1 = Find$df1,
      Find_df2 = Find$df2,
      Find_p = Find$p,
      k = k,
      cv = cv,
      method = "PSS Bounds Testing (Pesaran, Shin & Smith, 2001)",
      cv_source = "Kripfganz & Schneider (2020)"
    )
  } else {
    # Bootstrap critical values (Bertelli, Vacca & Zoia, 2022)
    message("\n  Computing bootstrap distributions...")

    bs_result <- bootstrap_cointegration_test(
      model = model,
      data = data,
      depvar = depvar,
      decompose = decompose,
      control_vars = control_vars,
      fourier = fourier,
      kstar = kstar,
      p = p, q = q, r = r,
      nobs = nobs,
      reps = reps,
      level_vars = level_vars,
      indep_level_vars = indep_level_vars
    )

    result <- list(
      type = "fbnardl",
      Fov = Fov$F,
      Fov_cv10 = bs_result$Fov_cv10,
      Fov_cv05 = bs_result$Fov_cv05,
      Fov_cv01 = bs_result$Fov_cv01,
      Fov_pval = bs_result$Fov_pval,
      t_dep = t_dep,
      t_cv10 = bs_result$t_cv10,
      t_cv05 = bs_result$t_cv05,
      t_cv01 = bs_result$t_cv01,
      t_pval = bs_result$t_pval,
      Find = Find$F,
      Find_cv10 = bs_result$Find_cv10,
      Find_cv05 = bs_result$Find_cv05,
      Find_cv01 = bs_result$Find_cv01,
      Find_pval = bs_result$Find_pval,
      k = k,
      reps = reps,
      method = "Bootstrap Cointegration Test",
      cv_source = "Bertelli, Vacca & Zoia (2022)"
    )
  }

  return(result)
}


#' Compute Joint F-Test
#'
#' @param model lm object
#' @param vars Vector of variable names to test
#' @return List with F statistic, df, and p-value
#' @keywords internal
compute_joint_f_test <- function(model, vars) {
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)

  # Keep only variables that exist
  vars <- vars[vars %in% names(coefs)]

  if (length(vars) == 0) {
    return(list(F = NA, df1 = NA, df2 = NA, p = NA))
  }

  # Restriction matrix: R %*% beta = 0
  R <- matrix(0, nrow = length(vars), ncol = length(coefs))
  rownames(R) <- vars
  colnames(R) <- names(coefs)

  for (i in seq_along(vars)) {
    R[i, vars[i]] <- 1
  }

  # Beta vector (subset)
  beta <- coefs[vars]

  # Variance-covariance of restricted coefficients
  vcov_sub <- vcov_mat[vars, vars]

  # F statistic: (R*beta)' * (R*V*R')^-1 * (R*beta) / q
  q <- length(vars)
  df_r <- model$df.residual

  tryCatch({
    vcov_inv <- solve(vcov_sub)
    F_stat <- as.numeric(t(beta) %*% vcov_inv %*% beta) / q
    p_val <- stats::pf(F_stat, q, df_r, lower.tail = FALSE)

    return(list(F = F_stat, df1 = q, df2 = df_r, p = p_val))
  }, error = function(e) {
    return(list(F = NA, df1 = q, df2 = df_r, p = NA))
  })
}


#' Get PSS Critical Values
#'
#' Approximate critical values from Pesaran, Shin & Smith (2001) Table CI(iii)
#' Case III: unrestricted intercept, no trend
#'
#' @param k Number of long-run forcing variables
#' @return Data frame with I(0) and I(1) bounds at 10%, 5%, 1% levels
#' @keywords internal
get_pss_critical_values <- function(k) {

  # I(0) bounds (lower)
  lb <- switch(as.character(k),
               "1" = c(4.04, 4.94, 6.84),
               "2" = c(3.17, 3.79, 5.15),
               "3" = c(2.72, 3.23, 4.29),
               "4" = c(2.45, 2.86, 3.74),
               "5" = c(2.26, 2.62, 3.41),
               "6" = c(2.12, 2.45, 3.15),
               "7" = c(2.03, 2.32, 2.96),
               c(1.95, 2.22, 2.79))  # k >= 8

  # I(1) bounds (upper)
  ub <- switch(as.character(k),
               "1" = c(4.78, 5.73, 7.84),
               "2" = c(4.14, 4.85, 6.36),
               "3" = c(3.77, 4.35, 5.61),
               "4" = c(3.52, 4.01, 5.06),
               "5" = c(3.35, 3.79, 4.68),
               "6" = c(3.22, 3.61, 4.43),
               "7" = c(3.13, 3.50, 4.26),
               c(3.06, 3.39, 4.10))  # k >= 8

  data.frame(
    level = c("10%", "5%", "1%"),
    I0_bound = lb,
    I1_bound = ub
  )
}


#' Bootstrap Cointegration Test
#'
#' Implements the bootstrap procedure of Bertelli, Vacca & Zoia (2022)
#' and McNown, Sam & Goh (2018)
#'
#' @param model Estimated lm object
#' @param data Data frame
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param control_vars Names of control variables
#' @param fourier Logical
#' @param kstar Fourier frequency
#' @param p,q,r Lag orders
#' @param nobs Number of observations
#' @param reps Number of bootstrap replications
#' @param level_vars Names of all level variables
#' @param indep_level_vars Names of independent level variables
#' @return List with bootstrap critical values and p-values
#' @keywords internal
bootstrap_cointegration_test <- function(model, data, depvar, decompose, control_vars,
                                          fourier, kstar, p, q, r, nobs, reps,
                                          level_vars, indep_level_vars) {

  # Extract model components
  y <- stats::model.response(stats::model.frame(model))
  X <- stats::model.matrix(model)
  resid <- stats::residuals(model)
  fitted_vals <- stats::fitted(model)

  n <- length(y)

  # Storage for bootstrap statistics
  Fov_boot <- numeric(reps)
  t_boot <- numeric(reps)
  Find_boot <- numeric(reps)

  # Original statistics
  coefs_orig <- stats::coef(model)
  vcov_orig <- stats::vcov(model)

  # Bootstrap loop
  for (b in seq_len(reps)) {
    # Resample residuals with replacement
    resid_boot <- sample(resid, n, replace = TRUE)

    # Generate bootstrap y under H0 (no cointegration)
    # Under H0: level coefficients are zero
    # y_boot = X %*% beta_restricted + resid_boot

    # Create restricted model (setting level coefficients to zero)
    beta_restricted <- coefs_orig
    beta_restricted[level_vars[level_vars %in% names(beta_restricted)]] <- 0

    y_boot <- X %*% beta_restricted + resid_boot

    # Estimate model with bootstrap data
    model_data_boot <- data.frame(y_boot = y_boot, X[, -1, drop = FALSE])
    names(model_data_boot)[1] <- "dy"

    formula_boot <- stats::as.formula(paste("dy ~",
                                             paste(colnames(X)[-1], collapse = " + ")))

    tryCatch({
      fit_boot <- stats::lm(formula_boot, data = model_data_boot)

      # Compute bootstrap test statistics
      Fov_boot[b] <- compute_joint_f_test(fit_boot, level_vars)$F
      Find_boot[b] <- compute_joint_f_test(fit_boot, indep_level_vars)$F

      if ("L_y" %in% names(stats::coef(fit_boot))) {
        t_boot[b] <- stats::coef(fit_boot)["L_y"] /
          sqrt(stats::vcov(fit_boot)["L_y", "L_y"])
      } else {
        t_boot[b] <- NA
      }
    }, error = function(e) {
      Fov_boot[b] <- NA
      t_boot[b] <- NA
      Find_boot[b] <- NA
    })
  }

  # Remove NA values
  Fov_boot <- Fov_boot[!is.na(Fov_boot)]
  t_boot <- t_boot[!is.na(t_boot)]
  Find_boot <- Find_boot[!is.na(Find_boot)]

  # Compute critical values
  Fov_cv <- if (length(Fov_boot) > 0) {
    stats::quantile(Fov_boot, c(0.90, 0.95, 0.99))
  } else {
    c(NA, NA, NA)
  }

  # For t-test, use lower tail (negative values)
  t_cv <- if (length(t_boot) > 0) {
    stats::quantile(t_boot, c(0.10, 0.05, 0.01))
  } else {
    c(NA, NA, NA)
  }

  Find_cv <- if (length(Find_boot) > 0) {
    stats::quantile(Find_boot, c(0.90, 0.95, 0.99))
  } else {
    c(NA, NA, NA)
  }

  # Compute p-values
  Fov_orig <- compute_joint_f_test(model, level_vars)$F
  t_orig <- coefs_orig["L_y"] / sqrt(vcov_orig["L_y", "L_y"])
  Find_orig <- compute_joint_f_test(model, indep_level_vars)$F

  Fov_pval <- if (length(Fov_boot) > 0 && !is.na(Fov_orig)) {
    mean(Fov_boot >= Fov_orig, na.rm = TRUE)
  } else {
    NA
  }

  t_pval <- if (length(t_boot) > 0 && !is.na(t_orig)) {
    mean(t_boot <= t_orig, na.rm = TRUE)  # Lower tail for t-test
  } else {
    NA
  }

  Find_pval <- if (length(Find_boot) > 0 && !is.na(Find_orig)) {
    mean(Find_boot >= Find_orig, na.rm = TRUE)
  } else {
    NA
  }

  return(list(
    Fov_cv10 = Fov_cv[1],
    Fov_cv05 = Fov_cv[2],
    Fov_cv01 = Fov_cv[3],
    Fov_pval = Fov_pval,
    t_cv10 = t_cv[1],
    t_cv05 = t_cv[2],
    t_cv01 = t_cv[3],
    t_pval = t_pval,
    Find_cv10 = Find_cv[1],
    Find_cv05 = Find_cv[2],
    Find_cv01 = Find_cv[3],
    Find_pval = Find_pval
  ))
}

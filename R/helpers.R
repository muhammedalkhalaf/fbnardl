#' Decompose Variables into Positive and Negative Partial Sums
#'
#' @param data Data frame containing the variables
#' @param decompose Character vector of variable names to decompose
#' @return Data frame with original variables plus _pos and _neg partial sums
#' @importFrom utils head
#' @keywords internal
decompose_partials <- function(data, decompose) {
  result <- data

  for (var in decompose) {
    x <- data[[var]]
    dx <- c(NA, diff(x))

    # Positive partial sum: cumsum of max(dx, 0)
    pos_changes <- pmax(dx, 0)
    pos_changes[is.na(pos_changes)] <- 0
    x_pos <- cumsum(pos_changes)

    # Negative partial sum: cumsum of min(dx, 0)
    neg_changes <- pmin(dx, 0)
    neg_changes[is.na(neg_changes)] <- 0
    x_neg <- cumsum(neg_changes)

    result[[paste0(var, "_pos")]] <- x_pos
    result[[paste0(var, "_neg")]] <- x_neg
  }

  return(result)
}


#' Select Optimal Model via Two-Step Procedure
#'
#' Step 1: Select k* by minimum SSR (Yilanci et al. 2020)
#' Step 2: Select (p, q, r) by AIC/BIC with fixed k*
#'
#' @param data Data frame with decomposed variables
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param control_vars Names of control variables
#' @param maxlag Maximum lag order
#' @param maxk Maximum Fourier frequency
#' @param ic Information criterion ("aic" or "bic")
#' @param fourier Logical, include Fourier terms
#' @param nobs Number of observations
#' @return List with optimal p, q, r, kstar, ic_value, and ssr_by_k
#' @keywords internal
select_optimal_model <- function(data, depvar, decompose, control_vars,
                                  maxlag, maxk, ic, fourier, nobs) {

  # Generate time index
  data$t_index <- seq_len(nrow(data))

  # Generate first difference of dependent variable
  y <- data[[depvar]]
  dy <- c(NA, diff(y))
  data$dy <- dy

  # Build Fourier frequency grid
  if (fourier) {
    k_grid <- seq(0.1, maxk, by = 0.1)
  } else {
    k_grid <- 0
  }

  # ============================================================================
  # Step 1: Select k* by minimum SSR
  # ============================================================================
  message("  Step 1: Selecting k* by minimum SSR (Yilanci et al. 2020)...")

  ssr_results <- data.frame(k = k_grid, ssr = NA_real_)
  best_ssr <- Inf
  best_kstar <- k_grid[1]

  for (ki in seq_along(k_grid)) {
    k <- k_grid[ki]

    # Generate Fourier terms
    if (k > 0) {
      data$fourier_sin <- sin(2 * pi * k * data$t_index / nobs)
      data$fourier_cos <- cos(2 * pi * k * data$t_index / nobs)
    }

    # Build maximal model formula
    formula_parts <- c()

    # Lagged differences of dependent variable
    for (j in 1:maxlag) {
      lag_name <- paste0("dy_L", j)
      data[[lag_name]] <- c(rep(NA, j), head(dy, -j))
      formula_parts <- c(formula_parts, lag_name)
    }

    # Decomposed variables (pos and neg) at maxlag
    for (var in decompose) {
      for (suffix in c("_pos", "_neg")) {
        vname <- paste0(var, suffix)
        x <- data[[vname]]
        dx <- c(NA, diff(x))

        for (j in 0:maxlag) {
          if (j == 0) {
            dvar_name <- paste0("d_", vname)
            data[[dvar_name]] <- dx
          } else {
            dvar_name <- paste0("d_", vname, "_L", j)
            data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
          }
          formula_parts <- c(formula_parts, dvar_name)
        }
      }
    }

    # Control variables at maxlag
    for (var in control_vars) {
      x <- data[[var]]
      dx <- c(NA, diff(x))

      for (j in 0:maxlag) {
        if (j == 0) {
          dvar_name <- paste0("d_", var)
          data[[dvar_name]] <- dx
        } else {
          dvar_name <- paste0("d_", var, "_L", j)
          data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
        }
        formula_parts <- c(formula_parts, dvar_name)
      }
    }

    # Lagged levels (ECM terms)
    # L.depvar
    data$L_y <- c(NA, head(y, -1))
    formula_parts <- c(formula_parts, "L_y")

    # L.decomposed_pos, L.decomposed_neg
    for (var in decompose) {
      for (suffix in c("_pos", "_neg")) {
        vname <- paste0(var, suffix)
        x <- data[[vname]]
        lag_name <- paste0("L_", vname)
        data[[lag_name]] <- c(NA, head(x, -1))
        formula_parts <- c(formula_parts, lag_name)
      }
    }

    # L.controls
    for (var in control_vars) {
      x <- data[[var]]
      lag_name <- paste0("L_", var)
      data[[lag_name]] <- c(NA, head(x, -1))
      formula_parts <- c(formula_parts, lag_name)
    }

    # Fourier terms
    if (k > 0) {
      formula_parts <- c(formula_parts, "fourier_sin", "fourier_cos")
    }

    # Fit model
    formula_str <- paste("dy ~", paste(formula_parts, collapse = " + "))
    tryCatch({
      fit <- stats::lm(stats::as.formula(formula_str), data = data,
                       na.action = stats::na.omit)
      ssr <- sum(stats::residuals(fit)^2)
      ssr_results$ssr[ki] <- ssr

      if (ssr < best_ssr) {
        best_ssr <- ssr
        best_kstar <- k
      }
    }, error = function(e) {
      ssr_results$ssr[ki] <- NA
    })
  }

  message("  Optimal k* = ", sprintf("%.2f", best_kstar),
          " (min SSR = ", sprintf("%.4f", best_ssr), ")")

  # ============================================================================
  # Step 2: Select (p, q, r) by IC with fixed k*
  # ============================================================================
  message("  Step 2: Selecting (p, q, r) by ", toupper(ic), " with fixed k*...")

  # Regenerate Fourier terms for best k*
  if (best_kstar > 0) {
    data$fourier_sin <- sin(2 * pi * best_kstar * data$t_index / nobs)
    data$fourier_cos <- cos(2 * pi * best_kstar * data$t_index / nobs)
  }

  best_ic <- Inf
  best_p <- 1
  best_q <- rep(0, length(decompose))
  best_r <- if (length(control_vars) > 0) rep(0, length(control_vars)) else integer(0)
  total_models <- 0

  # Grid search over all lag combinations
  for (p in 1:maxlag) {
    # Generate all combinations for q and r
    n_decomp <- length(decompose)
    n_ctrl <- length(control_vars)

    if (n_decomp > 0) {
      q_combos <- expand.grid(replicate(n_decomp, 0:maxlag, simplify = FALSE))
    } else {
      q_combos <- data.frame(matrix(nrow = 1, ncol = 0))
    }

    if (n_ctrl > 0) {
      r_combos <- expand.grid(replicate(n_ctrl, 0:maxlag, simplify = FALSE))
    } else {
      r_combos <- data.frame(matrix(nrow = 1, ncol = 0))
    }

    # Combine q and r combinations
    for (qi in seq_len(nrow(q_combos))) {
      for (ri in seq_len(nrow(r_combos))) {
        total_models <- total_models + 1

        q_vals <- if (n_decomp > 0) as.numeric(q_combos[qi, ]) else integer(0)
        r_vals <- if (n_ctrl > 0) as.numeric(r_combos[ri, ]) else integer(0)

        # Build formula
        formula_parts <- c()

        # Lagged differences of dependent variable
        for (j in 1:p) {
          lag_name <- paste0("dy_L", j)
          if (!lag_name %in% names(data)) {
            data[[lag_name]] <- c(rep(NA, j), head(dy, -j))
          }
          formula_parts <- c(formula_parts, lag_name)
        }

        # Decomposed variables with variable-specific q
        for (vi in seq_along(decompose)) {
          var <- decompose[vi]
          qi_val <- q_vals[vi]

          for (suffix in c("_pos", "_neg")) {
            vname <- paste0(var, suffix)
            x <- data[[vname]]
            dx <- c(NA, diff(x))

            for (j in 0:qi_val) {
              if (j == 0) {
                dvar_name <- paste0("d_", vname)
                if (!dvar_name %in% names(data)) {
                  data[[dvar_name]] <- dx
                }
              } else {
                dvar_name <- paste0("d_", vname, "_L", j)
                if (!dvar_name %in% names(data)) {
                  data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
                }
              }
              formula_parts <- c(formula_parts, dvar_name)
            }
          }
        }

        # Control variables with variable-specific r
        for (vi in seq_along(control_vars)) {
          var <- control_vars[vi]
          ri_val <- r_vals[vi]

          x <- data[[var]]
          dx <- c(NA, diff(x))

          for (j in 0:ri_val) {
            if (j == 0) {
              dvar_name <- paste0("d_", var)
              if (!dvar_name %in% names(data)) {
                data[[dvar_name]] <- dx
              }
            } else {
              dvar_name <- paste0("d_", var, "_L", j)
              if (!dvar_name %in% names(data)) {
                data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
              }
            }
            formula_parts <- c(formula_parts, dvar_name)
          }
        }

        # Lagged levels
        formula_parts <- c(formula_parts, "L_y")

        for (var in decompose) {
          formula_parts <- c(formula_parts,
                             paste0("L_", var, "_pos"),
                             paste0("L_", var, "_neg"))
        }

        for (var in control_vars) {
          formula_parts <- c(formula_parts, paste0("L_", var))
        }

        # Fourier terms
        if (best_kstar > 0) {
          formula_parts <- c(formula_parts, "fourier_sin", "fourier_cos")
        }

        # Fit model and compute IC
        formula_str <- paste("dy ~", paste(formula_parts, collapse = " + "))
        tryCatch({
          fit <- stats::lm(stats::as.formula(formula_str), data = data,
                           na.action = stats::na.omit)

          n_used <- stats::nobs(fit)
          k_params <- length(stats::coef(fit))
          ssr <- sum(stats::residuals(fit)^2)

          if (ic == "aic") {
            ic_val <- n_used * log(ssr / n_used) + 2 * k_params
          } else {
            ic_val <- n_used * log(ssr / n_used) + k_params * log(n_used)
          }

          if (ic_val < best_ic) {
            best_ic <- ic_val
            best_p <- p
            best_q <- q_vals
            best_r <- r_vals
          }
        }, error = function(e) {
          # Skip invalid models
        })
      }
    }
  }

  return(list(
    p = best_p,
    q = best_q,
    r = best_r,
    kstar = best_kstar,
    ic_value = best_ic,
    total_models = total_models,
    ssr_by_k = ssr_results
  ))
}


#' Estimate Final Model
#'
#' @param data Data frame with decomposed variables
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param control_vars Names of control variables
#' @param p Optimal lag for dependent variable
#' @param q Optimal lags for decomposed variables
#' @param r Optimal lags for control variables
#' @param kstar Optimal Fourier frequency
#' @param nobs Number of observations
#' @param fourier Logical, include Fourier terms
#' @return lm object with estimated model
#' @keywords internal
estimate_final_model <- function(data, depvar, decompose, control_vars,
                                  p, q, r, kstar, nobs, fourier) {

  # Generate time index
  data$t_index <- seq_len(nrow(data))

  # Generate first difference of dependent variable
  y <- data[[depvar]]
  dy <- c(NA, diff(y))
  data$dy <- dy

  # Generate Fourier terms
  if (fourier && kstar > 0) {
    data$fourier_sin <- sin(2 * pi * kstar * data$t_index / nobs)
    data$fourier_cos <- cos(2 * pi * kstar * data$t_index / nobs)
  }

  # Build formula
  formula_parts <- c()

  # Lagged differences of dependent variable
  for (j in 1:p) {
    lag_name <- paste0("dy_L", j)
    data[[lag_name]] <- c(rep(NA, j), head(dy, -j))
    formula_parts <- c(formula_parts, lag_name)
  }

  # Decomposed variables with variable-specific q
  for (vi in seq_along(decompose)) {
    var <- decompose[vi]
    qi_val <- q[vi]

    for (suffix in c("_pos", "_neg")) {
      vname <- paste0(var, suffix)
      x <- data[[vname]]
      dx <- c(NA, diff(x))

      for (j in 0:qi_val) {
        if (j == 0) {
          dvar_name <- paste0("d_", vname)
          data[[dvar_name]] <- dx
        } else {
          dvar_name <- paste0("d_", vname, "_L", j)
          data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
        }
        formula_parts <- c(formula_parts, dvar_name)
      }
    }
  }

  # Control variables with variable-specific r
  for (vi in seq_along(control_vars)) {
    var <- control_vars[vi]
    ri_val <- r[vi]

    x <- data[[var]]
    dx <- c(NA, diff(x))

    for (j in 0:ri_val) {
      if (j == 0) {
        dvar_name <- paste0("d_", var)
        data[[dvar_name]] <- dx
      } else {
        dvar_name <- paste0("d_", var, "_L", j)
        data[[dvar_name]] <- c(rep(NA, j), head(dx, -j))
      }
      formula_parts <- c(formula_parts, dvar_name)
    }
  }

  # Lagged levels (ECM terms)
  data$L_y <- c(NA, head(y, -1))
  formula_parts <- c(formula_parts, "L_y")

  for (var in decompose) {
    for (suffix in c("_pos", "_neg")) {
      vname <- paste0(var, suffix)
      x <- data[[vname]]
      lag_name <- paste0("L_", vname)
      data[[lag_name]] <- c(NA, head(x, -1))
      formula_parts <- c(formula_parts, lag_name)
    }
  }

  for (var in control_vars) {
    x <- data[[var]]
    lag_name <- paste0("L_", var)
    data[[lag_name]] <- c(NA, head(x, -1))
    formula_parts <- c(formula_parts, lag_name)
  }

  # Fourier terms
  if (fourier && kstar > 0) {
    formula_parts <- c(formula_parts, "fourier_sin", "fourier_cos")
  }

  # Fit model
  formula_str <- paste("dy ~", paste(formula_parts, collapse = " + "))
  fit <- stats::lm(stats::as.formula(formula_str), data = data,
                   na.action = stats::na.omit)

  return(fit)
}


#' Compute Short-Run and Long-Run Multipliers
#'
#' @param model Estimated lm object
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param control_vars Names of control variables
#' @param q Optimal lags for decomposed variables
#' @param r Optimal lags for control variables
#' @return List with SR and LR multipliers for each variable
#' @keywords internal
compute_multipliers <- function(model, depvar, decompose, control_vars, q, r) {

  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  df_r <- model$df.residual

  result <- list()

  # ECM coefficient (speed of adjustment)
  ecm_coef <- coefs["L_y"]

  # Multipliers for decomposed variables
  for (vi in seq_along(decompose)) {
    var <- decompose[vi]
    qi_val <- q[vi]

    # Short-run multipliers (sum of differenced coefficients)
    sr_pos <- 0
    sr_neg <- 0

    for (j in 0:qi_val) {
      if (j == 0) {
        pos_name <- paste0("d_", var, "_pos")
        neg_name <- paste0("d_", var, "_neg")
      } else {
        pos_name <- paste0("d_", var, "_pos_L", j)
        neg_name <- paste0("d_", var, "_neg_L", j)
      }

      if (pos_name %in% names(coefs)) {
        sr_pos <- sr_pos + coefs[pos_name]
      }
      if (neg_name %in% names(coefs)) {
        sr_neg <- sr_neg + coefs[neg_name]
      }
    }

    # Long-run multipliers (using delta method)
    lr_pos_name <- paste0("L_", var, "_pos")
    lr_neg_name <- paste0("L_", var, "_neg")

    if (lr_pos_name %in% names(coefs) && !is.na(ecm_coef) && ecm_coef != 0) {
      lr_pos <- -coefs[lr_pos_name] / ecm_coef
      lr_neg <- -coefs[lr_neg_name] / ecm_coef

      # Delta method for standard errors
      lr_pos_se <- compute_lr_se(coefs, vcov_mat, lr_pos_name, "L_y")
      lr_neg_se <- compute_lr_se(coefs, vcov_mat, lr_neg_name, "L_y")

      lr_pos_t <- lr_pos / lr_pos_se
      lr_neg_t <- lr_neg / lr_neg_se
      lr_pos_p <- 2 * stats::pt(abs(lr_pos_t), df_r, lower.tail = FALSE)
      lr_neg_p <- 2 * stats::pt(abs(lr_neg_t), df_r, lower.tail = FALSE)
    } else {
      lr_pos <- lr_neg <- lr_pos_se <- lr_neg_se <- NA
      lr_pos_t <- lr_neg_t <- lr_pos_p <- lr_neg_p <- NA
    }

    result[[var]] <- list(
      sr_pos = sr_pos,
      sr_neg = sr_neg,
      lr_pos = lr_pos,
      lr_neg = lr_neg,
      lr_pos_se = lr_pos_se,
      lr_neg_se = lr_neg_se,
      lr_pos_t = lr_pos_t,
      lr_neg_t = lr_neg_t,
      lr_pos_p = lr_pos_p,
      lr_neg_p = lr_neg_p
    )
  }

  # Multipliers for control variables
  for (vi in seq_along(control_vars)) {
    var <- control_vars[vi]
    ri_val <- r[vi]

    # Short-run multiplier
    sr <- 0
    for (j in 0:ri_val) {
      if (j == 0) {
        name <- paste0("d_", var)
      } else {
        name <- paste0("d_", var, "_L", j)
      }
      if (name %in% names(coefs)) {
        sr <- sr + coefs[name]
      }
    }

    # Long-run multiplier
    lr_name <- paste0("L_", var)
    if (lr_name %in% names(coefs) && !is.na(ecm_coef) && ecm_coef != 0) {
      lr <- -coefs[lr_name] / ecm_coef
      lr_se <- compute_lr_se(coefs, vcov_mat, lr_name, "L_y")
      lr_t <- lr / lr_se
      lr_p <- 2 * stats::pt(abs(lr_t), df_r, lower.tail = FALSE)
    } else {
      lr <- lr_se <- lr_t <- lr_p <- NA
    }

    result[[var]] <- list(
      sr = sr,
      lr = lr,
      lr_se = lr_se,
      lr_t = lr_t,
      lr_p = lr_p
    )
  }

  return(result)
}


#' Compute Long-Run Standard Error via Delta Method
#'
#' @param coefs Coefficient vector
#' @param vcov_mat Variance-covariance matrix
#' @param beta_name Name of numerator coefficient
#' @param alpha_name Name of denominator coefficient (ECM)
#' @return Standard error of long-run multiplier
#' @keywords internal
compute_lr_se <- function(coefs, vcov_mat, beta_name, alpha_name) {
  beta <- coefs[beta_name]
  alpha <- coefs[alpha_name]

  if (is.na(alpha) || alpha == 0) return(NA)

  # Gradient of -beta/alpha
  grad_beta <- -1 / alpha
  grad_alpha <- beta / (alpha^2)

  # Variance via delta method
  var_beta <- vcov_mat[beta_name, beta_name]
  var_alpha <- vcov_mat[alpha_name, alpha_name]
  cov_ab <- vcov_mat[beta_name, alpha_name]

  var_lr <- grad_beta^2 * var_beta + grad_alpha^2 * var_alpha +
    2 * grad_beta * grad_alpha * cov_ab

  if (var_lr < 0) return(NA)

  return(sqrt(var_lr))
}


#' Compute Asymmetry Tests
#'
#' @param model Estimated lm object
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @return List with Wald test results for each variable
#' @keywords internal
compute_asymmetry_tests <- function(model, depvar, decompose) {

  result <- list()
  coefs <- stats::coef(model)
  vcov_mat <- stats::vcov(model)
  df_r <- model$df.residual

  for (var in decompose) {
    # Short-run asymmetry: test d_var_pos = d_var_neg (contemporaneous)
    pos_name <- paste0("d_", var, "_pos")
    neg_name <- paste0("d_", var, "_neg")

    if (pos_name %in% names(coefs) && neg_name %in% names(coefs)) {
      diff_coef <- coefs[pos_name] - coefs[neg_name]
      var_diff <- vcov_mat[pos_name, pos_name] + vcov_mat[neg_name, neg_name] -
        2 * vcov_mat[pos_name, neg_name]

      if (var_diff > 0) {
        sr_f <- (diff_coef^2) / var_diff
        sr_p <- stats::pf(sr_f, 1, df_r, lower.tail = FALSE)
      } else {
        sr_f <- sr_p <- NA
      }
    } else {
      sr_f <- sr_p <- NA
    }

    # Long-run asymmetry: test LR_pos = LR_neg
    lr_pos_name <- paste0("L_", var, "_pos")
    lr_neg_name <- paste0("L_", var, "_neg")
    alpha_name <- "L_y"

    if (all(c(lr_pos_name, lr_neg_name, alpha_name) %in% names(coefs))) {
      alpha <- coefs[alpha_name]
      beta_pos <- coefs[lr_pos_name]
      beta_neg <- coefs[lr_neg_name]

      if (!is.na(alpha) && alpha != 0) {
        # LR_pos - LR_neg = (-beta_pos/alpha) - (-beta_neg/alpha) = (beta_neg - beta_pos)/alpha
        # Gradient w.r.t. beta_pos: 1/alpha
        # Gradient w.r.t. beta_neg: -1/alpha
        # Gradient w.r.t. alpha: (beta_pos - beta_neg)/alpha^2

        diff_lr <- (beta_neg - beta_pos) / alpha
        grad <- c(1/alpha, -1/alpha, (beta_pos - beta_neg) / (alpha^2))

        var_idx <- c(lr_pos_name, lr_neg_name, alpha_name)
        vcov_sub <- vcov_mat[var_idx, var_idx]
        var_diff_lr <- as.numeric(t(grad) %*% vcov_sub %*% grad)

        if (var_diff_lr > 0) {
          lr_chi2 <- (diff_lr^2) / var_diff_lr
          lr_p <- stats::pchisq(lr_chi2, 1, lower.tail = FALSE)
        } else {
          lr_chi2 <- lr_p <- NA
        }
      } else {
        lr_chi2 <- lr_p <- NA
      }
    } else {
      lr_chi2 <- lr_p <- NA
    }

    result[[var]] <- list(
      sr_f = sr_f,
      sr_p = sr_p,
      lr_chi2 = lr_chi2,
      lr_p = lr_p
    )
  }

  return(result)
}

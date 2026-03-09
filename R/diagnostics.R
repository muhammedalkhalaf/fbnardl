#' Compute Diagnostic Tests
#'
#' @param model Estimated lm object
#' @return List with diagnostic test results
#' @keywords internal
compute_diagnostics <- function(model) {

  resid <- stats::residuals(model)
  n <- length(resid)

  result <- list()

  # ============================================================================
  # 1. Breusch-Godfrey LM Test for Serial Correlation (order 1 and 4)
  # ============================================================================
  result$bg1 <- bg_test(model, 1)
  result$bg4 <- bg_test(model, 4)

  # ============================================================================
  # 2. ARCH LM Test (Engle, 1982)
  # ============================================================================
  result$arch1 <- arch_test(resid, 1)
  result$arch4 <- arch_test(resid, 4)

  # ============================================================================
  # 3. Breusch-Pagan Test for Heteroskedasticity
  # ============================================================================
  result$bp <- bp_test(model)

  # ============================================================================
  # 4. Jarque-Bera Normality Test
  # ============================================================================
  result$jb <- jb_test(resid)

  # ============================================================================
  # 5. Ramsey RESET Test (functional form)
  # ============================================================================
  result$reset <- reset_test(model)

  return(result)
}


#' Breusch-Godfrey LM Test for Serial Correlation
#'
#' @param model lm object
#' @param order Lag order for test
#' @return List with LM statistic, df, and p-value
#' @keywords internal
bg_test <- function(model, order) {
  resid <- stats::residuals(model)
  n <- length(resid)
  X <- stats::model.matrix(model)

  # Create lagged residuals
  resid_lags <- matrix(NA, nrow = n, ncol = order)
  for (j in seq_len(order)) {
    resid_lags[, j] <- c(rep(NA, j), head(resid, -j))
  }

  # Auxiliary regression: e_t on X and lagged residuals
  aux_data <- data.frame(e = resid, X[, -1, drop = FALSE], resid_lags)
  names(aux_data)[(ncol(X) + 1):ncol(aux_data)] <- paste0("e_L", seq_len(order))

  aux_formula <- stats::as.formula(paste("e ~",
                                          paste(c(names(aux_data)[2:(ncol(X))],
                                                  paste0("e_L", seq_len(order))),
                                                collapse = " + ")))

  tryCatch({
    aux_fit <- stats::lm(aux_formula, data = aux_data, na.action = stats::na.omit)
    r2 <- summary(aux_fit)$r.squared
    n_aux <- stats::nobs(aux_fit)
    lm_stat <- n_aux * r2
    p_val <- stats::pchisq(lm_stat, order, lower.tail = FALSE)

    return(list(
      name = paste0("Breusch-Godfrey LM(", order, ")"),
      statistic = lm_stat,
      df = order,
      p_value = p_val
    ))
  }, error = function(e) {
    return(list(
      name = paste0("Breusch-Godfrey LM(", order, ")"),
      statistic = NA,
      df = order,
      p_value = NA
    ))
  })
}


#' ARCH LM Test
#'
#' @param resid Residual vector
#' @param order Lag order for test
#' @return List with LM statistic, df, and p-value
#' @keywords internal
arch_test <- function(resid, order) {
  n <- length(resid)
  e2 <- resid^2

  # Create lagged squared residuals
  e2_lags <- matrix(NA, nrow = n, ncol = order)
  for (j in seq_len(order)) {
    e2_lags[, j] <- c(rep(NA, j), head(e2, -j))
  }

  # Auxiliary regression
  aux_data <- data.frame(e2 = e2, e2_lags)
  names(aux_data)[-1] <- paste0("e2_L", seq_len(order))

  aux_formula <- stats::as.formula(paste("e2 ~",
                                          paste(paste0("e2_L", seq_len(order)),
                                                collapse = " + ")))

  tryCatch({
    aux_fit <- stats::lm(aux_formula, data = aux_data, na.action = stats::na.omit)
    r2 <- summary(aux_fit)$r.squared
    n_aux <- stats::nobs(aux_fit)
    lm_stat <- n_aux * r2
    p_val <- stats::pchisq(lm_stat, order, lower.tail = FALSE)

    return(list(
      name = paste0("ARCH LM(", order, ")"),
      statistic = lm_stat,
      df = order,
      p_value = p_val
    ))
  }, error = function(e) {
    return(list(
      name = paste0("ARCH LM(", order, ")"),
      statistic = NA,
      df = order,
      p_value = NA
    ))
  })
}


#' Breusch-Pagan Test for Heteroskedasticity
#'
#' @param model lm object
#' @return List with LM statistic, df, and p-value
#' @keywords internal
bp_test <- function(model) {
  resid <- stats::residuals(model)
  X <- stats::model.matrix(model)
  n <- length(resid)

  e2 <- resid^2
  e2_scaled <- e2 / mean(e2)

  # Regress scaled squared residuals on X
  aux_data <- data.frame(e2_scaled = e2_scaled, X[, -1, drop = FALSE])

  aux_formula <- stats::as.formula(paste("e2_scaled ~",
                                          paste(names(aux_data)[-1], collapse = " + ")))

  tryCatch({
    aux_fit <- stats::lm(aux_formula, data = aux_data)
    ess <- sum((stats::fitted(aux_fit) - mean(e2_scaled))^2)
    bp_stat <- ess / 2
    df <- ncol(X) - 1
    p_val <- stats::pchisq(bp_stat, df, lower.tail = FALSE)

    return(list(
      name = "Breusch-Pagan",
      statistic = bp_stat,
      df = df,
      p_value = p_val
    ))
  }, error = function(e) {
    return(list(
      name = "Breusch-Pagan",
      statistic = NA,
      df = NA,
      p_value = NA
    ))
  })
}


#' Jarque-Bera Normality Test
#'
#' @param resid Residual vector
#' @return List with JB statistic, df, and p-value
#' @keywords internal
jb_test <- function(resid) {
  n <- length(resid)
  m2 <- mean(resid^2)
  m3 <- mean(resid^3)
  m4 <- mean(resid^4)

  skewness <- m3 / (m2^1.5)
  kurtosis <- m4 / (m2^2)

  jb_stat <- n / 6 * (skewness^2 + (kurtosis - 3)^2 / 4)
  p_val <- stats::pchisq(jb_stat, 2, lower.tail = FALSE)

  return(list(
    name = "Jarque-Bera",
    statistic = jb_stat,
    df = 2,
    p_value = p_val,
    skewness = skewness,
    kurtosis = kurtosis
  ))
}


#' Ramsey RESET Test
#'
#' @param model lm object
#' @param power Powers of fitted values to include (default: 2:3)
#' @return List with F statistic, df, and p-value
#' @keywords internal
reset_test <- function(model, power = 2:3) {
  y <- stats::model.response(stats::model.frame(model))
  X <- stats::model.matrix(model)
  fitted_vals <- stats::fitted(model)

  # Add powers of fitted values
  fitted_powers <- sapply(power, function(p) fitted_vals^p)
  colnames(fitted_powers) <- paste0("fitted_", power)

  # Augmented model
  X_aug <- cbind(X, fitted_powers)

  tryCatch({
    # Original model RSS
    rss0 <- sum(stats::residuals(model)^2)
    df0 <- model$df.residual

    # Augmented model
    aug_fit <- stats::lm.fit(X_aug, y)
    rss1 <- sum(aug_fit$residuals^2)
    df1 <- length(y) - ncol(X_aug)

    # F test
    q <- ncol(fitted_powers)
    f_stat <- ((rss0 - rss1) / q) / (rss1 / df1)
    p_val <- stats::pf(f_stat, q, df1, lower.tail = FALSE)

    return(list(
      name = "Ramsey RESET",
      statistic = f_stat,
      df1 = q,
      df2 = df1,
      p_value = p_val
    ))
  }, error = function(e) {
    return(list(
      name = "Ramsey RESET",
      statistic = NA,
      df1 = NA,
      df2 = NA,
      p_value = NA
    ))
  })
}


#' Compute Dynamic Multipliers
#'
#' Following Shin, Yu & Greenwood-Nimmo (2014), compute cumulative
#' dynamic multipliers for positive and negative shocks
#'
#' @param model Estimated lm object
#' @param depvar Name of dependent variable
#' @param decompose Names of decomposed variables
#' @param p Lag order for dependent variable
#' @param horizon Number of periods for multiplier computation
#' @return List with dynamic multiplier paths
#' @keywords internal
compute_dynamic_multipliers <- function(model, depvar, decompose, p, horizon) {

  coefs <- stats::coef(model)
  result <- list()

  # ECM coefficient
  rho <- coefs["L_y"]
  if (is.na(rho)) rho <- 0

  # Autoregressive coefficients
  phi <- numeric(p)
  for (j in 1:p) {
    name <- paste0("dy_L", j)
    if (name %in% names(coefs)) {
      phi[j] <- coefs[name]
    }
  }

  for (var in decompose) {
    # Level coefficients
    beta_pos <- coefs[paste0("L_", var, "_pos")]
    beta_neg <- coefs[paste0("L_", var, "_neg")]

    if (is.na(beta_pos)) beta_pos <- 0
    if (is.na(beta_neg)) beta_neg <- 0

    # Short-run coefficients (contemporaneous)
    theta_pos <- coefs[paste0("d_", var, "_pos")]
    theta_neg <- coefs[paste0("d_", var, "_neg")]

    if (is.na(theta_pos)) theta_pos <- 0
    if (is.na(theta_neg)) theta_neg <- 0

    # Compute impulse responses via recursive formula
    # h_t = theta_0 + sum_j(phi_j * h_{t-j}) + rho * H_{t-1} + beta
    # where H_t = cumulative sum of h

    h_pos <- numeric(horizon + 1)
    h_neg <- numeric(horizon + 1)
    H_pos <- numeric(horizon + 1)  # Cumulative
    H_neg <- numeric(horizon + 1)

    for (t in 1:(horizon + 1)) {
      if (t == 1) {
        h_pos[t] <- theta_pos
        h_neg[t] <- theta_neg
      } else {
        # Autoregressive component
        ar_pos <- 0
        ar_neg <- 0
        for (j in 1:min(p, t - 1)) {
          ar_pos <- ar_pos + phi[j] * h_pos[t - j]
          ar_neg <- ar_neg + phi[j] * h_neg[t - j]
        }

        # ECM component
        ecm_pos <- rho * H_pos[t - 1] + beta_pos
        ecm_neg <- rho * H_neg[t - 1] + beta_neg

        h_pos[t] <- ar_pos + ecm_pos
        h_neg[t] <- ar_neg + ecm_neg
      }

      # Cumulative sum
      H_pos[t] <- sum(h_pos[1:t])
      H_neg[t] <- sum(h_neg[1:t])
    }

    result[[var]] <- list(
      horizon = 0:horizon,
      h_pos = h_pos,
      h_neg = h_neg,
      H_pos = H_pos,
      H_neg = H_neg,
      H_diff = H_pos - H_neg  # Asymmetry
    )
  }

  return(result)
}

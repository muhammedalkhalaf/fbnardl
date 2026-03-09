#' Fourier Bootstrap Nonlinear ARDL Estimation
#'
#' @description
#' Estimates a Fourier Bootstrap Nonlinear Autoregressive Distributed Lag
#' (FBNARDL) model for analyzing asymmetric cointegration relationships.
#' Combines the NARDL framework of Shin, Yu & Greenwood-Nimmo (2014) with
#' Fourier approximation for smooth structural breaks and optional bootstrap
#' inference.
#'
#' @details
#' The FBNARDL model extends the standard ARDL bounds testing approach by:
#' \enumerate{
#'   \item Decomposing specified regressors into positive and negative partial
#'         sums to capture asymmetric effects
#'   \item Adding Fourier terms (sine and cosine) to approximate smooth
#'         structural breaks without specifying break dates
#'   \item Providing either analytical critical values (Kripfganz & Schneider,
#'         2020) or bootstrap critical values (Bertelli, Vacca & Zoia, 2022)
#' }
#'
#' The general model specification is:
#' \deqn{\Delta y_t = \alpha + \sum_{j=1}^{p} \phi_j \Delta y_{t-j} +
#'       \sum_{j=0}^{q} (\theta^+_j \Delta x^+_{t-j} + \theta^-_j \Delta x^-_{t-j}) +
#'       \rho y_{t-1} + \beta^+ x^+_{t-1} + \beta^- x^-_{t-1} +
#'       \gamma_1 \sin(2\pi k^* t/T) + \gamma_2 \cos(2\pi k^* t/T) + \varepsilon_t}
#'
#' where \eqn{x^+} and \eqn{x^-} are positive and negative partial sums of the
#' decomposed variable, and \eqn{k^*} is the optimal Fourier frequency.
#'
#' The two-step selection procedure follows Yilanci et al. (2020):
#' \enumerate{
#'   \item Select optimal Fourier frequency \eqn{k^*} by minimizing SSR
#'   \item Select optimal lag orders (p, q, r) by AIC/BIC with fixed \eqn{k^*}
#' }
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} where y is the
#'   dependent variable and x1, x2, ... are the regressors.
#' @param data A data frame containing the time series variables.
#' @param decompose Character vector of variable names to decompose into
#'   positive and negative partial sums for asymmetric effects.
#' @param type Character string specifying the model type. Either \code{"fnardl"}
#'   for Fourier NARDL with analytical critical values (Kripfganz & Schneider,
#'   2020) or \code{"fbnardl"} for Fourier Bootstrap NARDL with bootstrap
#'   critical values (Bertelli, Vacca & Zoia, 2022). Default is \code{"fnardl"}.
#' @param maxlag Integer specifying the maximum lag order to search for the
#'   dependent variable and all regressors. Default is 4.
#' @param maxk Numeric specifying the maximum Fourier frequency to search.
#'   Default is 3. Frequencies are searched in increments of 0.1.
#' @param ic Character string specifying the information criterion for lag
#'   selection. Either \code{"aic"} or \code{"bic"}. Default is \code{"aic"}.
#' @param reps Integer specifying the number of bootstrap replications when
#'   \code{type = "fbnardl"}
#' @param level Numeric specifying the confidence level for intervals.
#'   Default is 0.95.
#' @param fourier Logical. If \code{FALSE}, no Fourier terms are included
#'   (pure NARDL). Default is \code{TRUE}.
#' @param horizon Integer specifying the horizon for dynamic multiplier
#'   computation. Default is 20.
#' @param graph Logical. If \code{TRUE}, produce diagnostic and multiplier
#'   plots. Default is \code{TRUE}.
#'
#' @return An object of class \code{"fbnardl"} containing:
#' \describe{
#'   \item{coefficients}{Named vector of estimated coefficients}
#'   \item{vcov}{Variance-covariance matrix of coefficients}
#'   \item{residuals}{Model residuals}
#'   \item{fitted.values}{Fitted values}
#'   \item{lags}{List with optimal lag orders (p, q for each decomposed var, r
#'     for controls)}
#'   \item{kstar}{Optimal Fourier frequency}
#'   \item{ic_value}{Information criterion value of selected model}
#'   \item{multipliers}{List with short-run and long-run multipliers}
#'   \item{asymmetry_tests}{Wald test results for asymmetry}
#'   \item{bounds_test}{Bounds/bootstrap cointegration test results}
#'   \item{diagnostics}{Diagnostic test results}
#'   \item{dynamic_multipliers}{Dynamic multiplier paths}
#'   \item{nobs}{Number of observations used}
#'   \item{df.residual}{Residual degrees of freedom}
#'   \item{r.squared}{R-squared}
#'   \item{adj.r.squared}{Adjusted R-squared}
#'   \item{sigma}{Residual standard error}
#'   \item{call}{The matched call}
#'   \item{formula}{The model formula}
#'   \item{model}{The model frame}
#'   \item{decomposed_vars}{Names of decomposed variables}
#'   \item{control_vars}{Names of control (non-decomposed) variables}
#'   \item{type}{Model type used}
#' }
#'
#' @references
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework.
#' In R. C. Sickles & W. C. Horrace (Eds.), \emph{Festschrift in Honor of
#' Peter Schmidt} (pp. 281-314). Springer. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
#' to the analysis of level relationships. \emph{Journal of Applied
#' Econometrics}, 16(3), 289-326. \doi{10.1002/jae.616}
#'
#' Yilanci, V., Bozoklu, S., & Gorus, M. S. (2020). Are BRICS countries
#' pollution havens? Evidence from a bootstrap ARDL bounds testing approach
#' with a Fourier function. \emph{Sustainable Cities and Society}, 55, 102035.
#' \doi{10.1016/j.scs.2020.102035}
#'
#' Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions
#' for critical value bounds and approximate p-values in equilibrium
#' correction models. \emph{Oxford Bulletin of Economics and Statistics},
#' 82(6), 1456-1481. \doi{10.1111/obes.12377}
#'
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration. \emph{Applied
#' Economics}, 50(13), 1509-1521. \doi{10.1080/00036846.2017.1366643}
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(oil_inflation)
#'
#' # Estimate FNARDL model (analytical critical values)
#' model_fnardl <- fbnardl(
#'   inflation ~ oil_price + gdp_growth,
#'   data = oil_inflation,
#'   decompose = "oil_price",
#'   type = "fnardl",
#'   maxlag = 4,
#'   maxk = 3
#' )
#' summary(model_fnardl)
#'
#' # Estimate FBNARDL model (bootstrap critical values)
#' model_fbnardl <- fbnardl(
#'   inflation ~ oil_price + gdp_growth,
#'   data = oil_inflation,
#'   decompose = "oil_price",
#'   type = "fbnardl",
#'   maxlag = 4,
#'   reps = 999
#' )
#' summary(model_fbnardl)
#'
#' # Extract dynamic multipliers
#' plot(model_fnardl, type = "multipliers")
#' }
#'
#' @export
fbnardl <- function(formula, data, decompose, type = c("fnardl", "fbnardl"),
                    maxlag = 4, maxk = 3, ic = c("aic", "bic"), reps = 999,
                    level = 0.95, fourier = TRUE, horizon = 20,
                    graph = TRUE) {

  # Match arguments
  type <- match.arg(type)
  ic <- match.arg(ic)
  call <- match.call()

 # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (maxlag < 1) {
    stop("'maxlag' must be >= 1")
  }
  if (maxk < 0.1 && fourier) {
    stop("'maxk' must be >= 0.1 when fourier = TRUE")
  }
  if (reps < 100 && type == "fbnardl") {
    stop("'reps' must be >= 100 for bootstrap inference")
  }

  # Parse formula
  formula_vars <- all.vars(formula)
  depvar <- formula_vars[1]
  all_regressors <- formula_vars[-1]

  # Validate decompose variables
  if (missing(decompose) || length(decompose) == 0) {
    stop("'decompose' must specify at least one variable for asymmetric decomposition")
  }
  if (!all(decompose %in% all_regressors)) {
    stop("All variables in 'decompose' must be in the formula")
  }

  # Separate decomposed and control variables
  control_vars <- setdiff(all_regressors, decompose)

  # Extract data
  model_data <- data[, c(depvar, all_regressors), drop = FALSE]
  model_data <- stats::na.omit(model_data)
  nobs <- nrow(model_data)

  if (nobs < 30) {
    stop("Too few observations (", nobs, "). Need at least 30.")
  }

  # ============================================================================
  # Message: Model estimation starting
  # ============================================================================
  message("\n", paste(rep("=", 70), collapse = ""))
  if (type == "fbnardl") {
    message("  Fourier Bootstrap NARDL (FBNARDL) Estimation")
  } else {
    message("  Fourier NARDL (FNARDL) Estimation")
  }
  message(paste(rep("=", 70), collapse = ""))
  message("  Dependent variable : ", depvar)
  message("  Decomposed var(s)  : ", paste(decompose, collapse = ", "))
  if (length(control_vars) > 0) {
    message("  Control var(s)     : ", paste(control_vars, collapse = ", "))
  }
  message("  Max lag (p,q,r)    : ", maxlag)
  if (fourier) {
    message("  Max Fourier freq   : ", maxk)
  }
  message("  Info criterion     : ", toupper(ic))
  message("  Observations       : ", nobs)
  if (type == "fbnardl") {
    message("  Bootstrap reps     : ", reps)
  }
  message(paste(rep("=", 70), collapse = ""))

  # ============================================================================
  # Decompose variables into positive/negative partial sums
  # ============================================================================
  decomposed_data <- decompose_partials(model_data, decompose)

  # ============================================================================
  # Two-step selection: k* by min SSR, then (p,q,r) by IC
  # ============================================================================
  message("\n  Searching for optimal model...")

  selection_result <- select_optimal_model(
    data = decomposed_data,
    depvar = depvar,
    decompose = decompose,
    control_vars = control_vars,
    maxlag = maxlag,
    maxk = maxk,
    ic = ic,
    fourier = fourier,
    nobs = nobs
  )

  message("  Models evaluated   : ", selection_result$total_models)
  message("  Best ", toupper(ic), "          : ",
          sprintf("%.4f", selection_result$ic_value))

  # ============================================================================
  # Estimate final model
  # ============================================================================
  final_model <- estimate_final_model(
    data = decomposed_data,
    depvar = depvar,
    decompose = decompose,
    control_vars = control_vars,
    p = selection_result$p,
    q = selection_result$q,
    r = selection_result$r,
    kstar = selection_result$kstar,
    nobs = nobs,
    fourier = fourier
  )

  # ============================================================================
  # Compute multipliers
  # ============================================================================
  multipliers <- compute_multipliers(
    model = final_model,
    depvar = depvar,
    decompose = decompose,
    control_vars = control_vars,
    q = selection_result$q,
    r = selection_result$r
  )

  # ============================================================================
  # Asymmetry tests
  # ============================================================================
  asymmetry_tests <- compute_asymmetry_tests(
    model = final_model,
    depvar = depvar,
    decompose = decompose
  )

  # ============================================================================
  # Bounds/bootstrap cointegration test
  # ============================================================================
  bounds_test <- compute_bounds_test(
    model = final_model,
    data = decomposed_data,
    depvar = depvar,
    decompose = decompose,
    control_vars = control_vars,
    type = type,
    reps = reps,
    fourier = fourier,
    kstar = selection_result$kstar,
    p = selection_result$p,
    q = selection_result$q,
    r = selection_result$r,
    nobs = nobs
  )

  # ============================================================================
  # Diagnostic tests
  # ============================================================================
  diagnostics <- compute_diagnostics(final_model)

  # ============================================================================
  # Dynamic multipliers
  # ============================================================================
  dynamic_mult <- compute_dynamic_multipliers(
    model = final_model,
    depvar = depvar,
    decompose = decompose,
    p = selection_result$p,
    horizon = horizon
  )

  # ============================================================================
  # Prepare output object
  # ============================================================================
  result <- list(
    coefficients = stats::coef(final_model),
    vcov = stats::vcov(final_model),
    residuals = stats::residuals(final_model),
    fitted.values = stats::fitted(final_model),
    lags = list(
      p = selection_result$p,
      q = selection_result$q,
      r = selection_result$r
    ),
    kstar = selection_result$kstar,
    ic_value = selection_result$ic_value,
    aic = stats::AIC(final_model),
    bic = stats::BIC(final_model),
    multipliers = multipliers,
    asymmetry_tests = asymmetry_tests,
    bounds_test = bounds_test,
    diagnostics = diagnostics,
    dynamic_multipliers = dynamic_mult,
    nobs = nrow(final_model$model),
    df.residual = final_model$df.residual,
    r.squared = summary(final_model)$r.squared,
    adj.r.squared = summary(final_model)$adj.r.squared,
    sigma = summary(final_model)$sigma,
    fstatistic = summary(final_model)$fstatistic,
    call = call,
    formula = formula,
    model = final_model$model,
    decomposed_vars = decompose,
    control_vars = control_vars,
    type = type,
    level = level,
    ic = ic,
    fourier = fourier,
    ssr_by_k = selection_result$ssr_by_k,
    lm_object = final_model
  )

  class(result) <- "fbnardl"

  # ============================================================================
  # Generate plots if requested
  # ============================================================================
  if (graph) {
    plot.fbnardl(result, type = "all")
  }

  return(result)
}


#' @export
print.fbnardl <- function(x, ...) {
  cat("\nFourier Bootstrap Nonlinear ARDL Model\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Type:", toupper(x$type), "\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Decomposed:", paste(x$decomposed_vars, collapse = ", "), "\n")
  if (length(x$control_vars) > 0) {
    cat("Controls:", paste(x$control_vars, collapse = ", "), "\n")
  }
  cat("\nSelected Model:\n")
  cat("  p (dep. var lags):", x$lags$p, "\n")
  for (i in seq_along(x$decomposed_vars)) {
    cat("  q_", x$decomposed_vars[i], " lags: ", x$lags$q[i], "\n", sep = "")
  }
  if (length(x$control_vars) > 0) {
    for (i in seq_along(x$control_vars)) {
      cat("  r_", x$control_vars[i], " lags: ", x$lags$r[i], "\n", sep = "")
    }
  }
  if (x$fourier) {
    cat("  k* (Fourier freq):", sprintf("%.2f", x$kstar), "\n")
  }
  cat("\n", toupper(x$ic), ":", sprintf("%.4f", x$ic_value), "\n")
  cat("Observations:", x$nobs, "\n")
  cat("R-squared:", sprintf("%.4f", x$r.squared), "\n")
  cat("Adj. R-squared:", sprintf("%.4f", x$adj.r.squared), "\n")
  invisible(x)
}


#' @export
summary.fbnardl <- function(object, ...) {
  # Print header
  cat("\n", paste(rep("=", 78), collapse = ""), "\n")
  if (object$type == "fbnardl") {
    cat("  Fourier Bootstrap NARDL (FBNARDL) Estimation Results\n")
  } else {
    cat("  Fourier NARDL (FNARDL) Estimation Results\n")
  }
  cat(paste(rep("=", 78), collapse = ""), "\n\n")

  # Model selection summary
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 1: Model Selection\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Selected lag p (depvar lags)       :", object$lags$p, "\n")
  for (i in seq_along(object$decomposed_vars)) {
    cat("  Selected lag q (", object$decomposed_vars[i], " lags)    : ",
        object$lags$q[i], "\n", sep = "")
  }
  if (length(object$control_vars) > 0) {
    for (i in seq_along(object$control_vars)) {
      cat("  Selected lag r (", object$control_vars[i], " lags)    : ",
          object$lags$r[i], "\n", sep = "")
    }
  }
  if (object$fourier) {
    cat("  Selected Fourier frequency (k*)    :", sprintf("%.2f", object$kstar), "\n")
  }
  cat("  Information criterion (", toupper(object$ic), ")     : ",
      sprintf("%.4f", object$ic_value), "\n", sep = "")
  cat("  AIC                               :", sprintf("%.4f", object$aic), "\n")
  cat("  BIC                               :", sprintf("%.4f", object$bic), "\n")
  cat("  Observations (used)                :", object$nobs, "\n")
  cat("  R-squared                          :", sprintf("%.4f", object$r.squared), "\n")
  cat("  Adjusted R-squared                 :", sprintf("%.4f", object$adj.r.squared), "\n")
  if (!is.null(object$fstatistic)) {
    fstat_p <- stats::pf(object$fstatistic[1], object$fstatistic[2],
                         object$fstatistic[3], lower.tail = FALSE)
    cat("  F-statistic                        : ",
        sprintf("%.4f", object$fstatistic[1]),
        " (p=", sprintf("%.4f", fstat_p), ")\n", sep = "")
  }
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # Coefficient table
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 2: Estimation Results\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")

  coef_summary <- summary(object$lm_object)$coefficients
  print_coef_table(coef_summary)

  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # Multipliers
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 3: Short-Run & Long-Run Multipliers\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  print_multipliers(object$multipliers, object$decomposed_vars, object$control_vars)
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # Asymmetry tests
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 4: Wald Tests for Asymmetry\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  print_asymmetry_tests(object$asymmetry_tests, object$decomposed_vars)
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # Cointegration test
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 5: Cointegration Test\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  print_bounds_test(object$bounds_test, object$type)
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # Diagnostics
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Table 6: Diagnostic Tests\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  print_diagnostics(object$diagnostics)
  cat(paste(rep("-", 78), collapse = ""), "\n\n")

  # References
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  References\n")
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat("  Shin, Yu & Greenwood-Nimmo (2014). Modelling asymmetric cointegration\n")
  cat("    and dynamic multipliers in a nonlinear ARDL framework.\n")
  cat("  Pesaran, Shin & Smith (2001). Bounds testing approaches to the\n")
  cat("    analysis of level relationships. JASA, 16(3), 289-326.\n")
  if (object$type == "fnardl") {
    cat("  Kripfganz & Schneider (2020). Response surface regressions for\n")
    cat("    critical value bounds. Oxford Bull. Econ. Stat., 82(6).\n")
  } else {
    cat("  Bertelli, Vacca & Zoia (2022). Bootstrap cointegration tests\n")
    cat("    in ARDL models. Stat. Methods & Applications, 31, 1231-1268.\n")
    cat("  McNown, Sam & Goh (2018). Bootstrapping the ARDL test for\n")
    cat("    cointegration. Applied Economics, 50(13), 1509-1521.\n")
  }
  if (object$fourier) {
    cat("  Yilanci, Bozoklu & Gorus (2020). Fourier ARDL approach.\n")
  }
  cat(paste(rep("-", 78), collapse = ""), "\n")

  invisible(object)
}


#' Plot method for fbnardl objects
#'
#' @param x An object of class \code{"fbnardl"}
#' @param type Character string specifying the plot type. One of
#'   \code{"multipliers"}, \code{"kstar"}, \code{"residuals"}, or \code{"all"}.
#' @param ... Additional arguments passed to plotting functions
#'
#' @export
plot.fbnardl <- function(x, type = c("all", "multipliers", "kstar", "residuals"), ...) {
  type <- match.arg(type)

  if (type == "all" || type == "kstar") {
    if (!is.null(x$ssr_by_k) && x$fourier) {
      plot_kstar_selection(x$ssr_by_k, x$kstar)
    }
  }

  if (type == "all" || type == "multipliers") {
    plot_dynamic_multipliers(x$dynamic_multipliers, x$decomposed_vars)
  }

  if (type == "all" || type == "residuals") {
    plot_residual_diagnostics(x)
  }

  invisible(x)
}


#' Extract coefficients from fbnardl object
#' @param object An object of class \code{"fbnardl"}
#' @param ... Additional arguments (ignored)
#' @export
coef.fbnardl <- function(object, ...) {
  object$coefficients
}


#' Extract variance-covariance matrix from fbnardl object
#' @param object An object of class \code{"fbnardl"}
#' @param ... Additional arguments (ignored)
#' @export
vcov.fbnardl <- function(object, ...) {
  object$vcov
}


#' Extract residuals from fbnardl object
#' @param object An object of class \code{"fbnardl"}
#' @param ... Additional arguments (ignored)
#' @export
residuals.fbnardl <- function(object, ...) {
  object$residuals
}


#' Extract fitted values from fbnardl object
#' @param object An object of class \code{"fbnardl"}
#' @param ... Additional arguments (ignored)
#' @export
fitted.fbnardl <- function(object, ...) {
  object$fitted.values
}


#' Extract number of observations from fbnardl object
#' @param object An object of class \code{"fbnardl"}
#' @param ... Additional arguments (ignored)
#' @importFrom stats nobs
#' @export
nobs.fbnardl <- function(object, ...) {
  object$nobs
}

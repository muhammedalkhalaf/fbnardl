#' Plot Fourier Frequency Selection (SSR vs k*)
#'
#' @param ssr_by_k Data frame with k and ssr columns
#' @param kstar Optimal k* value
#' @keywords internal
plot_kstar_selection <- function(ssr_by_k, kstar) {
  if (is.null(ssr_by_k) || all(is.na(ssr_by_k$ssr))) {
    return(invisible(NULL))
  }

  # Remove NA values
  plot_data <- ssr_by_k[!is.na(ssr_by_k$ssr), ]

  if (nrow(plot_data) < 2) {
    return(invisible(NULL))
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mar = c(4, 4, 3, 1))

  # Main plot
  graphics::plot(plot_data$k, plot_data$ssr,
                 type = "b", pch = 16, col = "navy",
                 xlab = expression(paste("Fourier Frequency (", k^"*", ")")),
                 ylab = "Sum of Squared Residuals (SSR)",
                 main = "Fourier Frequency Selection (Yilanci et al. 2020)")

  # Mark optimal k*
  opt_idx <- which.min(abs(plot_data$k - kstar))
  if (length(opt_idx) > 0) {
    graphics::points(kstar, plot_data$ssr[opt_idx],
                     pch = 18, col = "red", cex = 2)
    graphics::abline(v = kstar, col = "red", lty = 2)
    graphics::legend("topright",
                     legend = paste0("Optimal k* = ", sprintf("%.2f", kstar)),
                     pch = 18, col = "red", bty = "n")
  }
}


#' Plot Dynamic Multipliers
#'
#' @param dyn_mult List of dynamic multipliers from compute_dynamic_multipliers
#' @param decomposed_vars Names of decomposed variables
#' @keywords internal
plot_dynamic_multipliers <- function(dyn_mult, decomposed_vars) {
  if (is.null(dyn_mult) || length(dyn_mult) == 0) {
    return(invisible(NULL))
  }

  n_vars <- length(decomposed_vars)

  # Set up multi-panel plot
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(n_vars, 2), mar = c(4, 4, 3, 1))

  for (var in decomposed_vars) {
    if (!var %in% names(dyn_mult)) next

    dm <- dyn_mult[[var]]

    # Panel 1: Cumulative multipliers (positive and negative)
    ylim <- range(c(dm$H_pos, dm$H_neg), na.rm = TRUE)
    graphics::plot(dm$horizon, dm$H_pos,
                   type = "l", col = "blue", lwd = 2,
                   xlab = "Horizon", ylab = "Cumulative Multiplier",
                   main = paste0(var, ": Cumulative Dynamic Multipliers"),
                   ylim = ylim)
    graphics::lines(dm$horizon, dm$H_neg, col = "red", lwd = 2)
    graphics::abline(h = 0, col = "gray", lty = 2)
    graphics::legend("topright",
                     legend = c("Positive shock", "Negative shock"),
                     col = c("blue", "red"), lwd = 2, bty = "n", cex = 0.8)

    # Panel 2: Asymmetry (H+ - H-)
    graphics::plot(dm$horizon, dm$H_diff,
                   type = "l", col = "purple", lwd = 2,
                   xlab = "Horizon", ylab = "Asymmetry",
                   main = paste0(var, ": Dynamic Asymmetry (H+ - H-)"))
    graphics::abline(h = 0, col = "gray", lty = 2)
  }
}


#' Plot Residual Diagnostics
#'
#' @param x fbnardl object
#' @keywords internal
plot_residual_diagnostics <- function(x) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  resid <- x$residuals
  fitted <- x$fitted.values

  # Panel 1: Residuals vs Fitted
  graphics::plot(fitted, resid,
                 pch = 16, col = grDevices::rgb(0, 0, 0.5, 0.5),
                 xlab = "Fitted Values", ylab = "Residuals",
                 main = "Residuals vs Fitted")
  graphics::abline(h = 0, col = "red", lty = 2)
  graphics::lines(stats::lowess(fitted, resid), col = "blue", lwd = 2)

  # Panel 2: Histogram of residuals with normal curve
  graphics::hist(resid, breaks = 20, freq = FALSE,
                 col = "lightblue", border = "white",
                 xlab = "Residuals", main = "Distribution of Residuals")
  x_seq <- seq(min(resid), max(resid), length.out = 100)
  graphics::lines(x_seq, stats::dnorm(x_seq, mean(resid), stats::sd(resid)),
                  col = "red", lwd = 2)

  # Panel 3: Q-Q plot
  stats::qqnorm(resid, pch = 16, col = grDevices::rgb(0, 0, 0.5, 0.5),
                main = "Normal Q-Q Plot")
  stats::qqline(resid, col = "red", lwd = 2)

  # Panel 4: ACF of residuals
  stats::acf(resid, main = "ACF of Residuals", na.action = stats::na.pass)
}


#' Print Coefficient Table with Significance Stars
#'
#' @param coef_summary Coefficient summary matrix from summary.lm
#' @keywords internal
print_coef_table <- function(coef_summary) {
  cat(sprintf("  %-25s %10s %10s %8s %8s\n",
              "Variable", "Estimate", "Std.Err.", "t-stat", "Pr(>|t|)"))
  cat(paste(rep("-", 78), collapse = ""), "\n")

  for (i in seq_len(nrow(coef_summary))) {
    var_name <- rownames(coef_summary)[i]
    est <- coef_summary[i, 1]
    se <- coef_summary[i, 2]
    t_val <- coef_summary[i, 3]
    p_val <- coef_summary[i, 4]

    stars <- get_stars(p_val)

    cat(sprintf("  %-25s %10.4f %10.4f %8.3f %8.4f %s\n",
                substr(var_name, 1, 25), est, se, t_val, p_val, stars))
  }
}


#' Get Significance Stars
#'
#' @param p_val P-value
#' @return Character string with stars
#' @keywords internal
get_stars <- function(p_val) {
  if (is.na(p_val)) return("")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  if (p_val < 0.1) return(".")
  return("")
}


#' Print Multipliers Table
#'
#' @param multipliers List of multipliers
#' @param decomposed_vars Names of decomposed variables
#' @param control_vars Names of control variables
#' @keywords internal
print_multipliers <- function(multipliers, decomposed_vars, control_vars) {

  for (var in decomposed_vars) {
    if (!var %in% names(multipliers)) next

    m <- multipliers[[var]]

    cat("\n  Variable:", var, "(decomposed)\n")
    cat("  ", paste(rep("-", 68), collapse = ""), "\n")
    cat(sprintf("  %-18s %10s %10s %8s %8s\n",
                "Component", "Estimate", "Std.Err.", "t-stat", "p-value"))
    cat("  ", paste(rep("-", 68), collapse = ""), "\n")

    cat(sprintf("  %-18s %10.4f\n", "Short-Run (+)", m$sr_pos))
    cat(sprintf("  %-18s %10.4f\n", "Short-Run (-)", m$sr_neg))
    cat("  ", paste(rep("-", 58), collapse = ""), "\n")

    if (!is.na(m$lr_pos)) {
      cat(sprintf("  %-18s %10.4f %10.4f %8.3f %8.4f %s\n",
                  "Long-Run (+)", m$lr_pos, m$lr_pos_se, m$lr_pos_t, m$lr_pos_p,
                  get_stars(m$lr_pos_p)))
    }
    if (!is.na(m$lr_neg)) {
      cat(sprintf("  %-18s %10.4f %10.4f %8.3f %8.4f %s\n",
                  "Long-Run (-)", m$lr_neg, m$lr_neg_se, m$lr_neg_t, m$lr_neg_p,
                  get_stars(m$lr_neg_p)))
    }
    cat("  ", paste(rep("-", 68), collapse = ""), "\n")

    # Asymmetry ratios
    if (m$sr_neg != 0) {
      cat(sprintf("  SR Asymmetry |SR(+)/SR(-)| = %.3f\n", abs(m$sr_pos / m$sr_neg)))
    }
    if (!is.na(m$lr_neg) && m$lr_neg != 0) {
      cat(sprintf("  LR Asymmetry |LR(+)/LR(-)| = %.3f\n", abs(m$lr_pos / m$lr_neg)))
    }
  }

  for (var in control_vars) {
    if (!var %in% names(multipliers)) next

    m <- multipliers[[var]]

    cat("\n  Variable:", var, "(control)\n")
    cat("  ", paste(rep("-", 68), collapse = ""), "\n")
    cat(sprintf("  %-18s %10.4f\n", "Short-Run", m$sr))
    if (!is.na(m$lr)) {
      cat(sprintf("  %-18s %10.4f %10.4f %8.3f %8.4f %s\n",
                  "Long-Run", m$lr, m$lr_se, m$lr_t, m$lr_p,
                  get_stars(m$lr_p)))
    }
    cat("  ", paste(rep("-", 68), collapse = ""), "\n")
  }
}


#' Print Asymmetry Tests
#'
#' @param tests List of asymmetry test results
#' @param decomposed_vars Names of decomposed variables
#' @keywords internal
print_asymmetry_tests <- function(tests, decomposed_vars) {

  for (var in decomposed_vars) {
    if (!var %in% names(tests)) next

    t <- tests[[var]]

    cat("\n  Variable:", var, "\n")
    cat("  ", paste(rep("-", 60), collapse = ""), "\n")

    if (!is.na(t$sr_f)) {
      cat(sprintf("  Short-run asymmetry: F = %8.4f  p-value = %6.4f %s\n",
                  t$sr_f, t$sr_p, get_stars(t$sr_p)))
    } else {
      cat("  Short-run asymmetry: not estimable\n")
    }

    if (!is.na(t$lr_chi2)) {
      cat(sprintf("  Long-run asymmetry:  Chi2 = %8.4f  p-value = %6.4f %s\n",
                  t$lr_chi2, t$lr_p, get_stars(t$lr_p)))
    } else {
      cat("  Long-run asymmetry:  not estimable\n")
    }

    cat("  ", paste(rep("-", 60), collapse = ""), "\n")
  }
}


#' Print Bounds Test Results
#'
#' @param test List with bounds test results
#' @param type "fnardl" or "fbnardl"
#' @keywords internal
print_bounds_test <- function(test, type) {

  cat("\n  Method:", test$method, "\n")
  cat("  Critical Values:", test$cv_source, "\n")
  cat("  Number of long-run forcing variables (k) =", test$k, "\n\n")

  if (type == "fnardl") {
    # PSS bounds test
    cat("  ", paste(rep("-", 65), collapse = ""), "\n")
    cat(sprintf("  %-20s %10s %10s %10s\n",
                "Test", "Statistic", "p-value", ""))
    cat("  ", paste(rep("-", 65), collapse = ""), "\n")

    cat(sprintf("  %-20s %10.4f %10.4f\n",
                "F_overall (Fov)", test$Fov, test$Fov_p))
    cat(sprintf("  %-20s %10.4f %10.4f\n",
                "t_dependent", test$t_dep, test$t_dep_p))
    cat(sprintf("  %-20s %10.4f %10.4f\n",
                "F_independent (Find)", test$Find, test$Find_p))
    cat("  ", paste(rep("-", 65), collapse = ""), "\n\n")

    # Critical values table
    cat("  PSS Bounds Test Critical Values (Case III):\n")
    cat("  ", paste(rep("-", 65), collapse = ""), "\n")
    cat(sprintf("  %-10s %12s %12s %15s\n",
                "Signif.", "I(0) Bound", "I(1) Bound", "Decision"))
    cat("  ", paste(rep("-", 65), collapse = ""), "\n")

    for (i in seq_len(nrow(test$cv))) {
      lb <- test$cv$I0_bound[i]
      ub <- test$cv$I1_bound[i]

      if (!is.na(test$Fov)) {
        if (test$Fov > ub) {
          decision <- "Reject H0"
        } else if (test$Fov < lb) {
          decision <- "Fail to Reject"
        } else {
          decision <- "Inconclusive"
        }
      } else {
        decision <- "N/A"
      }

      cat(sprintf("  %-10s %12.3f %12.3f %15s\n",
                  test$cv$level[i], lb, ub, decision))
    }
    cat("  ", paste(rep("-", 65), collapse = ""), "\n")

  } else {
    # Bootstrap test
    cat("  Bootstrap replications:", test$reps, "\n\n")
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")
    cat(sprintf("  %-15s %8s %9s %9s %9s %8s\n",
                "Test", "Stat", "CV(10%)", "CV(5%)", "CV(1%)", "p-val"))
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")

    cat(sprintf("  %-15s %8.3f %9.3f %9.3f %9.3f %8.4f %s\n",
                "F_overall", test$Fov, test$Fov_cv10, test$Fov_cv05,
                test$Fov_cv01, test$Fov_pval, get_stars(test$Fov_pval)))
    cat(sprintf("  %-15s %8.3f %9.3f %9.3f %9.3f %8.4f %s\n",
                "t_dependent", test$t_dep, test$t_cv10, test$t_cv05,
                test$t_cv01, test$t_pval, get_stars(test$t_pval)))
    cat(sprintf("  %-15s %8.3f %9.3f %9.3f %9.3f %8.4f %s\n",
                "F_independent", test$Find, test$Find_cv10, test$Find_cv05,
                test$Find_cv01, test$Find_pval, get_stars(test$Find_pval)))
    cat("  ", paste(rep("-", 70), collapse = ""), "\n\n")

    # Decision
    cat("  Decision at 5% level:\n")
    if (!is.na(test$Fov_pval) && test$Fov_pval < 0.05) {
      cat("    Fov: Reject H0 => Evidence of a long-run relationship\n")
    } else {
      cat("    Fov: Fail to reject H0 => No evidence of long-run relationship\n")
    }
    if (!is.na(test$t_pval) && test$t_pval < 0.05) {
      cat("    t:   Reject H0 => Dependent variable participates in ECM\n")
      cat("         (rules out degenerate case #1)\n")
    } else {
      cat("    t:   Fail to reject H0 => Possible degenerate case #1\n")
    }
    if (!is.na(test$Find_pval) && test$Find_pval < 0.05) {
      cat("    Find: Reject H0 => Independent variables enter long-run\n")
      cat("          (rules out degenerate case #2)\n")
    } else {
      cat("    Find: Fail to reject H0 => Possible degenerate case #2\n")
    }
  }
}


#' Print Diagnostic Tests
#'
#' @param diag List of diagnostic test results
#' @keywords internal
print_diagnostics <- function(diag) {

  cat("\n")
  cat(sprintf("  %-30s %10s %8s %10s\n",
              "Test", "Statistic", "df", "p-value"))
  cat("  ", paste(rep("-", 65), collapse = ""), "\n")

  # Serial correlation
  if (!is.null(diag$bg1)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$bg1$name, diag$bg1$statistic, diag$bg1$df,
                diag$bg1$p_value, get_stars(diag$bg1$p_value)))
  }
  if (!is.null(diag$bg4)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$bg4$name, diag$bg4$statistic, diag$bg4$df,
                diag$bg4$p_value, get_stars(diag$bg4$p_value)))
  }

  # ARCH
  if (!is.null(diag$arch1)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$arch1$name, diag$arch1$statistic, diag$arch1$df,
                diag$arch1$p_value, get_stars(diag$arch1$p_value)))
  }
  if (!is.null(diag$arch4)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$arch4$name, diag$arch4$statistic, diag$arch4$df,
                diag$arch4$p_value, get_stars(diag$arch4$p_value)))
  }

  # Heteroskedasticity
  if (!is.null(diag$bp)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$bp$name, diag$bp$statistic, diag$bp$df,
                diag$bp$p_value, get_stars(diag$bp$p_value)))
  }

  # Normality
  if (!is.null(diag$jb)) {
    cat(sprintf("  %-30s %10.4f %8d %10.4f %s\n",
                diag$jb$name, diag$jb$statistic, diag$jb$df,
                diag$jb$p_value, get_stars(diag$jb$p_value)))
    cat(sprintf("    (Skewness: %.3f, Kurtosis: %.3f)\n",
                diag$jb$skewness, diag$jb$kurtosis))
  }

  # RESET
  if (!is.null(diag$reset)) {
    cat(sprintf("  %-30s %10.4f %3d,%3d %10.4f %s\n",
                diag$reset$name, diag$reset$statistic,
                diag$reset$df1, diag$reset$df2,
                diag$reset$p_value, get_stars(diag$reset$p_value)))
  }

  cat("  ", paste(rep("-", 65), collapse = ""), "\n")
  cat("  Note: Large p-values indicate no evidence against null hypothesis.\n")
}

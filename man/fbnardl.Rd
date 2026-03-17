\name{fbnardl}
\alias{fbnardl}
\title{Fourier Bootstrap Nonlinear ARDL Estimation}
\usage{
fbnardl(y, x_decomp, x_ctrl = NULL,
  type = c("fnardl", "fbnardl"),
  max_lag = 4L, max_k = 3,
  ic = c("aic", "bic"),
  reps = 499L,
  no_fourier = FALSE)
}
\arguments{
\item{y}{Numeric vector. Dependent variable.}

\item{x_decomp}{Numeric matrix or vector. Variable(s) to decompose into
  positive and negative partial sums.}

\item{x_ctrl}{Numeric matrix or vector, or \code{NULL}. Additional
  control regressors. Default is \code{NULL}.}

\item{type}{Character. \code{"fnardl"} (default) for Fourier NARDL with
  PSS/Kripfganz-Schneider critical values, or \code{"fbnardl"} for Fourier
  Bootstrap NARDL.}

\item{max_lag}{Integer. Maximum lag order. Default is \code{4}.}

\item{max_k}{Numeric. Maximum Fourier frequency. Default is \code{3}.}

\item{ic}{Character. Information criterion: \code{"aic"} (default) or
  \code{"bic"}.}

\item{reps}{Integer. Bootstrap replications for \code{type = "fbnardl"}.
  Default is \code{499}.}

\item{no_fourier}{Logical. If \code{TRUE}, suppress Fourier terms.
  Default is \code{FALSE}.}
}
\value{
A list of class \code{"fbnardl"} with components including
\code{coefficients}, \code{se}, \code{tstat}, \code{pval},
\code{best_p}, \code{best_q}, \code{best_kstar}, \code{r2},
\code{Fov}, \code{t_dep}, \code{Find}, \code{lr_pos}, \code{lr_neg},
\code{bootstrap} (if type is \code{"fbnardl"}), and \code{residuals}.
}
\description{
Estimates a nonlinear ARDL (NARDL) model augmented with flexible Fourier
terms and tests for cointegration using either PSS bounds-testing critical
values (FNARDL) or residual-based bootstrap critical values (FBNARDL).
}
\details{
Each variable in \code{x_decomp} is split into positive and negative
partial sums. A two-step procedure selects the optimal Fourier frequency
\eqn{k^*} by minimum SSR (Yilanci et al. 2020), then searches over all lag
combinations using the specified information criterion.
}
\examples{
set.seed(42)
n <- 80
x <- cumsum(rnorm(n))
y <- 0.6 * x + rnorm(n, sd = 0.5)
res <- fbnardl(y, x, type = "fnardl", max_lag = 2, max_k = 2)
print(res)
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework.
\doi{10.1007/978-1-4899-8008-3_9}

Bertelli, S., Vacca, G., & Zoia, M. G. (2022). Bootstrap cointegration
tests in ARDL models. \emph{Statistical Methods and Applications}, 31,
1231--1268. \doi{10.1007/s10260-021-00571-7}

Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions
for critical value bounds. \emph{Oxford Bulletin of Economics and
Statistics}, 82(6), 1456--1481. \doi{10.1111/obes.12376}
}

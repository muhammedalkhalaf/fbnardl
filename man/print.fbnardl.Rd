\name{print.fbnardl}
\alias{print.fbnardl}
\title{Print Method for fbnardl Objects}
\usage{
\method{print}{fbnardl}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"fbnardl"}.}

\item{...}{Further arguments passed to or from other methods (unused).}
}
\value{
Invisibly returns \code{x}.
}
\description{
Prints a formatted summary of Fourier Bootstrap Nonlinear ARDL results.
}
\examples{
set.seed(1)
n <- 60
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n)
res <- fbnardl(y, x, type = "fnardl", max_lag = 1, max_k = 1)
print(res)
}

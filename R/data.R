#' Oil Price and Inflation Example Dataset
#'
#' A simulated quarterly time series dataset for demonstrating the fbnardl
#' function. Contains data on inflation, oil prices, and GDP growth, designed
#' to exhibit asymmetric effects of oil price changes on inflation.
#'
#' @format A data frame with 120 rows and 4 columns:
#' \describe{
#'   \item{date}{Date, quarterly observations from 1990Q1 to 2019Q4}
#'   \item{inflation}{Numeric, consumer price inflation rate in percent}
#'   \item{oil_price}{Numeric, crude oil price index with base equals 100}
#'   \item{gdp_growth}{Numeric, real GDP growth rate in percent}
#' }
#'
#' @details
#' This dataset is simulated to illustrate the use of the FBNARDL model for
#' analyzing asymmetric effects. The data exhibits the following properties:
#' oil price increases have a larger effect on inflation than decreases,
#' a cointegrating relationship exists between inflation and oil prices,
#' and smooth structural change is present (suitable for Fourier approximation).
#'
#' @source Simulated data for package demonstration
#'
#' @examples
#' data(oil_inflation)
#' head(oil_inflation)
#' summary(oil_inflation)
#'
"oil_inflation"

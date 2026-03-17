# fbnardl

**Fourier Bootstrap Nonlinear ARDL Estimation**

[![CRAN status](https://www.r-pkg.org/badges/version/fbnardl)](https://CRAN.R-project.org/package=fbnardl)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`fbnardl` estimates asymmetric (nonlinear) ARDL models augmented with flexible
Fourier terms to account for smooth structural breaks. Two estimators are
provided:

| Estimator | Critical Values | Reference |
|-----------|----------------|-----------|
| **FNARDL** | PSS / Kripfganz-Schneider | Kripfganz & Schneider (2020) |
| **FBNARDL** | Residual bootstrap | Bertelli, Vacca & Zoia (2022) |

The nonlinear decomposition of regressors into positive and negative partial
sums follows Shin, Yu & Greenwood-Nimmo (2014).

## Installation

```r
install.packages("fbnardl")
```

## Quick Start

```r
library(fbnardl)

set.seed(42)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.6 * x + rnorm(n, sd = 0.5)

# FNARDL with PSS bounds
res <- fbnardl(y, x, type = "fnardl", max_lag = 3, max_k = 3)
print(res)

# FBNARDL with bootstrap critical values
res_bs <- fbnardl(y, x, type = "fbnardl", max_lag = 2, max_k = 2, reps = 499)
print(res_bs)
```

## References

Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework. In W. C. Horrace & R. C. Sickles (Eds.), *Festschrift in Honor of Peter Schmidt*. Springer. <https://doi.org/10.1007/978-1-4899-8008-3_9>

Bertelli, S., Vacca, G., & Zoia, M. G. (2022). Bootstrap cointegration tests in ARDL models. *Statistical Methods and Applications*, 31, 1231–1268. <https://doi.org/10.1007/s10260-021-00571-7>

Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions for critical value bounds and approximate p-values in equilibrium correction models. *Oxford Bulletin of Economics and Statistics*, 82(6), 1456–1481. <https://doi.org/10.1111/obes.12376>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>

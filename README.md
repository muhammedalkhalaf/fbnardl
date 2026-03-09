# fbnardl: Fourier Bootstrap Nonlinear ARDL for R

## Overview

The **fbnardl** package implements the Fourier Bootstrap Nonlinear Autoregressive Distributed Lag (FBNARDL) model for analyzing asymmetric cointegration relationships in R. It combines:

1. **Nonlinear ARDL (NARDL)** framework of Shin, Yu & Greenwood-Nimmo (2014) for asymmetric effects
2. **Fourier approximation** for smooth structural breaks (Yilanci, Bozoklu & Gorus, 2020)
3. **Bootstrap inference** for cointegration testing (Bertelli, Vacca & Zoia, 2022)

## Installation

```r
# Install from CRAN (when available)
install.packages("fbnardl")

# Or install development version from GitHub
# devtools::install_github("muhammedalkhalaf/fbnardl")
```

## Usage

```r
library(fbnardl)

# Load example data
data(oil_inflation)

# Estimate FNARDL model (analytical critical values)
model <- fbnardl(
  inflation ~ oil_price + gdp_growth,
  data = oil_inflation,
  decompose = "oil_price",
  type = "fnardl",
  maxlag = 4,
  maxk = 3
)

# View summary
summary(model)

# Plot dynamic multipliers
plot(model, type = "multipliers")
```

## Features

- **Asymmetric decomposition**: Decomposes specified variables into positive and negative partial sums to capture asymmetric effects
- **Fourier approximation**: Adds sine and cosine terms to approximate smooth structural breaks without specifying break dates
- **Two-step model selection**: 
  - Step 1: Select optimal Fourier frequency k* by minimum SSR
  - Step 2: Select optimal lag orders (p, q, r) by AIC/BIC
- **Two inference methods**:
  - `type = "fnardl"`: Analytical critical values (Kripfganz & Schneider, 2020)
  - `type = "fbnardl"`: Bootstrap critical values (Bertelli, Vacca & Zoia, 2022)
- **Comprehensive output**:
  - Short-run and long-run multipliers
  - Wald tests for asymmetry
  - PSS bounds test or bootstrap cointegration test
  - Diagnostic tests (serial correlation, ARCH, heteroskedasticity, normality, RESET)
  - Dynamic multiplier paths and plots

## Model Specification

The FBNARDL model is specified as:

```
Δy_t = α + Σ φ_j Δy_{t-j} + Σ (θ⁺_j Δx⁺_{t-j} + θ⁻_j Δx⁻_{t-j})
       + ρ y_{t-1} + β⁺ x⁺_{t-1} + β⁻ x⁻_{t-1}
       + γ₁ sin(2πk*t/T) + γ₂ cos(2πk*t/T) + ε_t
```

where:
- `x⁺` and `x⁻` are positive and negative partial sums
- `k*` is the optimal Fourier frequency
- `ρ` is the error correction coefficient (speed of adjustment)

## References

- Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework. In *Festschrift in Honor of Peter Schmidt* (pp. 281-314). Springer. [DOI: 10.1007/978-1-4899-8008-3_9](https://doi.org/10.1007/978-1-4899-8008-3_9)

- Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches to the analysis of level relationships. *Journal of Applied Econometrics*, 16(3), 289-326. [DOI: 10.1002/jae.616](https://doi.org/10.1002/jae.616)

- Yilanci, V., Bozoklu, S., & Gorus, M. S. (2020). Are BRICS countries pollution havens? Evidence from a bootstrap ARDL bounds testing approach with a Fourier function. *Sustainable Cities and Society*, 55, 102035. [DOI: 10.1016/j.scs.2020.102035](https://doi.org/10.1016/j.scs.2020.102035)

- Kripfganz, S., & Schneider, D. C. (2020). Response surface regressions for critical value bounds and approximate p-values in equilibrium correction models. *Oxford Bulletin of Economics and Statistics*, 82(6), 1456-1481. [DOI: 10.1111/obes.12377](https://doi.org/10.1111/obes.12377)

- McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the autoregressive distributed lag test for cointegration. *Applied Economics*, 50(13), 1509-1521. [DOI: 10.1080/00036846.2017.1366643](https://doi.org/10.1080/00036846.2017.1366643)

## Author

Independent Researcher

R port by OpenClaw

## License

GPL-3

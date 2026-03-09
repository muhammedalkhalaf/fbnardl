# fbnardl 1.0.0

* Initial CRAN release
* Main function: `fbnardl()` for Fourier Bootstrap Nonlinear ARDL estimation
* Features:
  - Asymmetric effects via positive/negative partial sum decomposition
  - Fourier terms for smooth structural breaks
  - Two-step model selection (Yilanci et al. 2020)
  - PSS bounds testing with Kripfganz & Schneider (2020) critical values
  - Bootstrap cointegration test (Bertelli, Vacca & Zoia, 2022)
  - Short-run and long-run multipliers with delta method inference
  - Wald tests for short-run and long-run asymmetry
  - Comprehensive diagnostic tests
  - Dynamic multiplier plots
* S3 methods: `print`, `summary`, `plot`, `coef`, `vcov`, `residuals`, `fitted`, `nobs`
* Example dataset: `oil_inflation`

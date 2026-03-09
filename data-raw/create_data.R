# Script to create example dataset for fbnardl package
# Run this script to regenerate data/oil_inflation.rda

set.seed(42)

n <- 120  # 30 years of quarterly data

# Generate dates
dates <- seq(as.Date("1990-01-01"), by = "quarter", length.out = n)

# Generate oil price with trend, cycles, and breaks
t <- 1:n
trend <- 50 + 0.3 * t
cycle <- 15 * sin(2 * pi * t / 40) + 10 * sin(2 * pi * 2 * t / n)
noise <- cumsum(rnorm(n, 0, 2))
oil_price <- trend + cycle + noise
oil_price <- pmax(oil_price, 20)  # Floor at 20

# Generate GDP growth (AR(1) process with drift)
gdp_growth <- numeric(n)
gdp_growth[1] <- 0.5
for (i in 2:n) {
  gdp_growth[i] <- 0.2 + 0.6 * gdp_growth[i-1] + rnorm(1, 0, 0.3)
}

# Generate inflation with asymmetric response to oil
# Positive oil changes have larger effect than negative
d_oil <- c(0, diff(oil_price))
d_oil_pos <- pmax(d_oil, 0)
d_oil_neg <- pmin(d_oil, 0)

# Cumulative partial sums
oil_pos <- cumsum(d_oil_pos)
oil_neg <- cumsum(d_oil_neg)

# Inflation follows NARDL-type DGP
inflation <- numeric(n)
inflation[1:2] <- c(2.5, 2.6)

for (i in 3:n) {
  # Short-run effects
  sr_pos <- 0.15 * d_oil_pos[i]  # Larger positive effect
  sr_neg <- 0.05 * abs(d_oil_neg[i])  # Smaller negative effect

  # Long-run equilibrium (asymmetric)
  lr_pos <- 0.03 * oil_pos[i-1]
  lr_neg <- 0.01 * oil_neg[i-1]

  # ECM
  ecm <- -0.15 * (inflation[i-1] - 2 - lr_pos - lr_neg)

  # GDP effect
  gdp_effect <- 0.2 * gdp_growth[i-1]

  # Fourier component (smooth structural change)
  fourier_effect <- 0.5 * sin(2 * pi * 1.5 * i / n)

  # AR component
  ar_effect <- 0.4 * (inflation[i-1] - inflation[i-2])

  # Combine
  d_inflation <- ar_effect + sr_pos + sr_neg + ecm + gdp_effect +
                 fourier_effect + rnorm(1, 0, 0.2)

  inflation[i] <- inflation[i-1] + d_inflation
}

# Create data frame
oil_inflation <- data.frame(
  date = dates,
  inflation = round(inflation, 2),
  oil_price = round(oil_price, 2),
  gdp_growth = round(gdp_growth, 2)
)

# Save to data directory
if (!dir.exists("../data")) {
  dir.create("../data")
}

save(oil_inflation, file = "../data/oil_inflation.rda", compress = "xz")

message("Dataset created: data/oil_inflation.rda")
message("Observations: ", nrow(oil_inflation))
message("Variables: ", paste(names(oil_inflation), collapse = ", "))

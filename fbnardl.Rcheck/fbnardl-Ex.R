pkgname <- "fbnardl"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "fbnardl-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('fbnardl')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("fbnardl")
### * fbnardl

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fbnardl
### Title: Fourier Bootstrap Nonlinear ARDL Estimation
### Aliases: fbnardl

### ** Examples

## No test: 
# Load example data
data(oil_inflation)

# Estimate FNARDL model (analytical critical values)
model_fnardl <- fbnardl(
  inflation ~ oil_price + gdp_growth,
  data = oil_inflation,
  decompose = "oil_price",
  type = "fnardl",
  maxlag = 4,
  maxk = 3
)
summary(model_fnardl)

# Estimate FBNARDL model (bootstrap critical values)
model_fbnardl <- fbnardl(
  inflation ~ oil_price + gdp_growth,
  data = oil_inflation,
  decompose = "oil_price",
  type = "fbnardl",
  maxlag = 4,
  reps = 999
)
summary(model_fbnardl)

# Extract dynamic multipliers
plot(model_fnardl, type = "multipliers")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fbnardl", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("oil_inflation")
### * oil_inflation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: oil_inflation
### Title: Oil Price and Inflation Example Dataset
### Aliases: oil_inflation
### Keywords: datasets

### ** Examples

data(oil_inflation)
head(oil_inflation)
summary(oil_inflation)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("oil_inflation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

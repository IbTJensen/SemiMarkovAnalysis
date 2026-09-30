library(spatstat)
library(data.table)
library(SeMaCoPE)

# Loading data -----------------------------------------------------------------
load("2_analysis_bank_data/data/pointsAndPoly.RData")
Covar_dt <- fread("2_analysis_bank_data/data/Covariates_processed.csv")

# Grid search for sat and R ----------------------------------------------------
Sys.time()
SS <- SeMaCoPE::SemiMarkov(
  X = ppp.banks,
  covariate = Covar_dt,
  edgecorrection = 0.2,
  R_within = seq(0.002, 0.1, 0.002),
  R_between = seq(0.001, 0.02, 0.001),
  sat = NULL,
  sat_within = c(1:10, Inf),
  sat_between = c(1:10, Inf),
  ncores = 12
)
Sys.time()

SS$R_within #0.004
SS$R_between #0.011
SS$maximum_log_likelihood
SS$sat_within #1
SS$sat_between #1
library(spatstat)
library(data.table)
library(SemiMarkov)

# Loading and covariates -------------------------------------------------------
load("Data/covariatesfct_alr.rdata")
load("Data/pointsAndPoly.RData")

Covar_dt <- data.table(xcoord = ppp.type$x,
                       ycoord = ppp.type$y)

for(i in 1:length(covariatesfct)){
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)

plot(covariatesfct[[1]])

# Grid search for sat and R ----------------------------------------------------
SS <- SemiMarkov(X = ppp.type,
                 covariate = Covar_dt[,-3],
                 edgecorrection = NULL,
                 R_within = seq(0.002, 0.01, 0.002),
                 R_between = seq(0.002, 0.01, 0.002),
                 sat = Inf)

SS$R_within # 0.35
SS$R_between # 0.005
SS$sat # 6

SS$CI

# a <- Sys.time()
SS <- SemiMarkov(X = ppp.type,
                 covariate = Covar_dt,
                 edgecorrection = NULL,
                 R_within = 0.35,
                 R_between = 0.005,
                 sat = 6,
                 standardize = T,
                 Poisson = F)
# Sys.time() - a

SS$CI



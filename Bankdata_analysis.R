library(spatstat)
library(data.table)
library(SemiMarkov)

# Loading and covariates -------------------------------------------------------
load("Data/covariatesfct_alr.rdata")
load("Data/pointsAndPoly.RData")
do_grid_search <- F

Covar_dt <- data.table(xcoord = ppp.type$x,
                       ycoord = ppp.type$y)

for(i in 1:length(covariatesfct)){
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)

plot(covariatesfct[[1]])

# Grid search for sat and R ----------------------------------------------------
if(do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt[,-3],
                   edgecorrection = NULL,
                   R_within = seq(0.001, 0.01, 0.001),
                   R_between = seq(0.001, 0.01, 0.001),
                   sat = seq(2, 10, 2))
}

SS$R_within # 0.006
SS$R_between # 0.004
SS$maximum_log_likelihood
SS$sat # 4

# Fit model and interpret results ----------------------------------------------
if(!do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt[,-3],
                   edgecorrection = NULL,
                   R_within = 0.008,
                   R_between = 0.004,
                   sat = 4,
                   standardize = T,
                   Poisson = F)
}

SS$CI



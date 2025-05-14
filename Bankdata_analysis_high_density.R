library(spatstat)
library(data.table)
library(SemiMarkov)

# Loading and covariates -------------------------------------------------------
load("Data/covariatesfct_alr_high_density.rdata")
load("Data/pointsAndPoly.RData")
load("Data/high_density_regions_polygon.rdata")

# Create high density point process
poly <- as.data.table(poly)
cities <- unique(poly[density > 500, city])
high_density <- inside.owin(ppp.type$x, ppp.type$y, tt)
ppp.type <- ppp.type[high_density]

Covar_dt <- data.table(xcoord = ppp.type$x,
                       ycoord = ppp.type$y)

for(i in 1:length(covariatesfct)){
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)

# SS <- SemiMarkov(X = ppp.type,
#                  covariate = Covar_dt,
#                  edgecorrection = NULL,
#                  R_within = c(0.005, 0.01, 0.015, 0.02),
#                  R_between = c(0.0025, 0.005, 0.01, 0.015),
#                  sat = c(2, 5, 10, 20))

SS <- SemiMarkov(X = ppp.type,
                 covariate = Covar_dt,
                 edgecorrection = NULL,
                 R_within = seq(0.004, 0.008, 0.001),
                 R_between = seq(0.002, 0.01, 0.001),
                 sat = Inf)

# SS <- SemiMarkov(X = ppp.type,
#                  covariate = Covar_dt,
#                  edgecorrection = NULL,
#                  R_within = seq(0.33, 0.4, 0.01),
#                  R_between = seq(0.01, 0.006, 0.001),
#                  sat = seq(2,10,2))

SS$R_within # 0.35
SS$R_between # 0.005
SS$sat # 6

SS$R_within # 0.006
SS$R_between # 0.003
SS$sat # Inf

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



library(spatstat)
library(data.table)
library(SemiMarkov)
library(ggplot2)

# Loading and covariates -------------------------------------------------------
load("Data/covariatesfct_alr.rdata")
load("Data/pointsAndPoly.RData")
load("Data/regions.rdata")
do_grid_search <- F

Covar_dt <- data.table(xcoord = ppp.type$x,
                       ycoord = ppp.type$y)

for(i in 3:length(covariatesfct)){
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)[-(1:2)]

# plot(covariatesfct[[3]])

CC <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")
# plot strong correlations
plot(Covar_dt$density, Covar_dt$decile)
plot(Covar_dt[,c("prop0.15", "prop16.24", "prop25.64")])
plot(Covar_dt[,c("activity", "median", "poverty")])
plot(Covar_dt[,c("density", "decile", "proptrade")])
plot(Covar_dt$proptrade, Covar_dt$proppublic)

# Covar_dt[,":="(prop16.24 = NULL, prop25.64 = NULL, median = NULL,
#                poverty = NULL, proptrade = NULL, decile = NULL)]
Covar_dt[,log_density:=log(density)]
Covar_dt[,":="(prop0.15 = NULL, prop25.64 = NULL, density = NULL,
               median = NULL, poverty = NULL)]
CC2 <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")

plot(Covar_dt$log_density, Covar_dt$prop16.24)
plot(Covar_dt$activity, Covar_dt$proppublic)

# Grid search for sat and R ----------------------------------------------------
if(do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = NULL,
                   R_within = seq(0.001, 0.01, 0.001),
                   R_between = seq(0.001, 0.01, 0.001),
                   sat = 2:10)
}

# Fit model and interpret results ----------------------------------------------
if(!do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = 0.016,
                   R_within = 0.008,
                   R_between = 0.004,
                   sat = 4,
                   standardize = T,
                   Poisson = F)
}

SS$R_within # 0.008
SS$R_between # 0.004
SS$maximum_log_likelihood
SS$sat # 4

SS$CI

SS_CI <- data.table(SS$CI)
SS_CI[,Covariate:=gsub(":1", "", Covariate)]
SS_CI[11:13,Covariate:=gsub("1", "Luc", Covariate)]
SS_CI[11:13,Covariate:=gsub("2", "Coop", Covariate)]
SS_CI[,"95% Confidence interval":=paste0(
  "[", round(Lower_CI, 2), ", ", round(Upper_CI, 2), "]")]
SS_CI[,":="(Lower_CI = NULL, Upper_CI = NULL)]
SS_CI$Covariate[1] <- "Intercept"
xt <- xtable::xtable(SS_CI, digits = 2)
print(xt, include.rownames=FALSE)

logit <- function(x) log(x/(1-x))
pred <- SS$pred[type_pred == 1 & !is.na(prob)]

# Plots ------------------------------------------------------------------------
pred[,logit_prob:=logit(prob)]
poly_dt <- data.table(poly)

ggplot(poly_dt, aes(x, y, group = code))+
  geom_polygon(colour = "white", fill = "white")+
  labs(x = "Longitude", y = "Latitude", colour = "Logit-probability")+
  geom_point(data = pred, mapping = aes(x = xcoord, y = ycoord,
                                        group = 1, colour = logit_prob),
             size = 0.5)+
  scale_colour_gradientn(colours = c("yellow", "orange", "red",
                                     "purple", "blue"))+
  theme(legend.position = "bottom")+
  NULL -> Points_in_France

ggsave("Figures/logit_prob_points.pdf", Points_in_France,
       width = 165, height = 140, units = "mm")

# Phi0 kernel estiamte ---------------------------------------------------------
bandwidth <- 0.25

# density_fct <- covariatesfct[[3]]
# bandwidth <- function(x, y){
#   dns <- density_fct(x, y)
#   if(is.na(dns)){
#     bdw <- 0.1
#     return(bdw)
#   }
#   bdw <- NULL
#   if(dns < 150){
#     bdw <- 0.25
#   }
#   if(dns >= 150 & dns < 400){
#     bdw <- 0.1
#   }
#   if(dns >= 400 & dns < 700){
#     bdw <- 0.075
#   }
#   if(dns >= 700){
#     bdw <- 0.05
#   }
#   return(bdw)
# }

x_pts <- seq(-5.1, 9.5, bandwidth)
y_pts <- seq(41.3, 51, bandwidth)
points_grid <- expand.grid(x = x_pts, y = y_pts)
in_France <- inside.owin(points_grid$x, points_grid$y, ppp.type$window)
points_grid <- points_grid[in_France,]

points_grid_ppp <- ppp(points_grid$x, points_grid$y, ppp.type$window)
ccc <- crosspairs(points_grid_ppp, ppp.type, rmax = bandwidth)

phi0 <- rep(0, nrow(points_grid))
N <- rep(0, nrow(points_grid))
for(i in 1:length(phi0)){
  cat(i, "\r")
  w <- which(ccc$i == i)
  p1 <- ppp.type[ccc$j[w]]
  idx <- (1:ncol(SS$w))[-c(1:4, ncol(SS$w))]
  A <- SS$w[xcoord %in% p1$x & ycoord %in% p1$y & type_obs == j, idx, with = F]
  gamma_v <- as.matrix(A) %*% SS$betahat
  if(nrow(A)>0){
    phi0[i] <- 1/2*sum(1/exp(gamma_v), na.rm = T)
  }
  N[i] <- nrow(A)
}

phi0_dt <- data.table(points_grid, phi0 = phi0, N = N)

phi0_mat <- matrix(NA, length(y_pts), length(x_pts))
colnames(phi0_mat) <- as.character(x_pts)
rownames(phi0_mat) <- as.character(y_pts)

for(i in 1:nrow(phi0_dt)){
  cat(i, "\r")
  i1 <- as.character(phi0_dt$x[i])
  i2 <- as.character(phi0_dt$y[i])
  phi0_mat[i2, i1] <- phi0_dt$phi0[i]
}

phi0_im <- im(phi0_mat)

# png(filename = "Figures/Kernel_est_0.25.png", width = 660, height = 435)
# plot(log(phi0_im), main = expression(paste("Log of Kernel estimate of ", phi[0])))
# dev.off()

pdf(file = "Figures/Kernel_est_0.25.pdf", width = 165/25.4, height = 140/25.4)
plot(log(phi0_im),
     main = expression(paste("Log of Kernel estimate of ", phi[0])))
dev.off()

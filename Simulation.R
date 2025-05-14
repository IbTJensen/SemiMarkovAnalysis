library(spatstat)
library(data.table)
library(geoR)
library(SemiMarkov)
# devtools::document("~/Git/SemiMarkov")
# devtools::load_all("~/Git/SemiMarkov")

# Setting up spatial covariates ------------------------------------------------
set.seed(123)
grid_pts <- expand.grid(x = seq(0, 2, 0.02), y = seq(0, 2, 0.02))

X <- grf(1, grid = grid_pts, cov.model = "exponential",
         cov.pars = c(sigmasq = 0.2, phi = 0.1))

Y <- grf(1, grid = grid_pts, cov.model = "exponential",
         cov.pars = c(sigmasq = 30^2, phi = 0.1), mean = 350)

X_mat <- matrix(NA, 101, 101)
Y_mat <- matrix(NA, 101, 101)
rownames(Y_mat) <- seq(0, 2, 0.02)
rownames(X_mat) <- seq(0, 2, 0.02)
colnames(Y_mat) <- seq(0, 2, 0.02)
colnames(X_mat) <- seq(0, 2, 0.02)
for(i in 1:101^2){
  iname <- as.character(grid_pts[i,1])
  jname <- as.character(grid_pts[i,2])
  X_mat[iname,jname] <- X$data[i]
  Y_mat[iname,jname] <- Y$data[i]
}

X_im <- im(X_mat, xcol = seq(0, 2, 0.02), yrow = seq(0, 2, 0.02))
phi0 <- im(Y_mat, xcol = seq(0, 2, 0.02), yrow = seq(0, 2, 0.02))

Spat_covar <- X_im
unit_window <- owin(c(0,1), c(0,1))

# Poisson process --------------------------------------------------------------
lambda_large1 <- phi0*exp(0.5*Spat_covar)
lambda_large2 <- phi0*exp(-0.5*Spat_covar)
lambda_large3 <- phi0

lambda_small1 <- lambda_large1[unit_window]
lambda_small2 <- lambda_large2[unit_window]
lambda_small3 <- lambda_large3[unit_window]

target_param <- c(0, 0, 0.5 ,-0.5, 0, 0, 0, 0, 0, 0)

nsim=1800
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA, n = NA)[-1]
for(i in 1:nsim){
  set.seed(i)
  cat(i, "\r")
  X_small1=rpoispp(lambda_small1)
  X_small2=rpoispp(lambda_small2)
  X_small3=rpoispp(lambda_small3)
  X_small=superimpose(X_small1, X_small2, X_small3)
  marks=rep(1,X_small$n)
  marks[(X_small1$n+1):(X_small1$n+X_small2$n)]=2
  marks[(X_small1$n+X_small2$n+1):X_small$n]=3
  X_small$marks=factor(marks)

  X_large1=rpoispp(lambda_large1)
  X_large2=rpoispp(lambda_large2)
  X_large3=rpoispp(lambda_large3)
  X_large=superimpose(X_large1, X_large2, X_large3)
  marks=rep(1,X_large$n)
  marks[(X_large1$n+1):(X_large1$n+X_large2$n)]=2
  marks[(X_large1$n+X_large2$n+1):X_large$n]=3
  X_large$marks=factor(marks)

  suppressWarnings(
    temp_small <- SemiMarkov(X = X_small,
                             covariate = Spat_covar[unit_window],
                             R_within = 0.02,
                             R_between = 0.04)
  )

  suppressWarnings(
    temp_large <- SemiMarkov(X = X_large,
                             covariate = Spat_covar,
                             R_within = 0.02,
                             R_between = 0.04)
  )

  data.table(Parameter = names(temp_small$betahat),
             Estimate = temp_small$betahat,
             Std.error = temp_small$std_err,
             within_CI = target_param > temp_small$CI$Lower_CI &
               target_param < temp_small$CI$Upper_CI,
             window_size = 1, n = X_small$n) -> res_small

  data.table(Parameter = names(temp_large$betahat),
             Estimate = temp_large$betahat,
             Std.error = temp_large$std_err,
             within_CI = target_param > temp_large$CI$Lower_CI &
               target_param < temp_large$CI$Upper_CI,
             window_size = 4, n = X_large$n) -> res_large

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,.(Estimate = mean(Estimate), Std.error = mean(Std.error),
                      MC_std.error = sd(Estimate), Coverage = mean(within_CI),
                      n = mean(n)),
                   by = list(Parameter, window_size)]

data.table(Parameter = res$Parameter,
           Target_param = target_param,
           Estimate = res$Estimate,
           Lower = res$Estimate - 1.96*res$MC_std.error/sqrt(nsim),
           Upper = res$Estimate + 1.96*res$MC_std.error/sqrt(nsim)) -> dt

fwrite(dt, "Poisson results/Sim_Poisson_parameter_estimate_MC.csv")
fwrite(res, "Poisson results/Poisson_res_summary.csv")
fwrite(Results_sim, "Poisson results/Poisson_full_results.csv")

# Multi-Strauss ----------------------------------------------------------------
set.seed(123)
ns=1000
nr=1e7
nv=1e6

phi0 <- 1.6*phi0

trend1 <- function(x, y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts]*exp(0.5*Spat_covar[pts]))
}

trend2 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts]*exp(-0.5*Spat_covar[pts]))
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts])
}

btw_int <- c(0.9, 0.9, 0.9)
wtn_int <- c(0.8, 0.8, 0.8)
gmma <- matrix(0, 3, 3)
gmma[lower.tri(gmma)] <- btw_int
gmma <- gmma + t(gmma)
diag(gmma) <- wtn_int

r    <- matrix(0.04, 3, 3)
diag(r) <- 0.02
mod_small <- list(cif = "straussm",
                  par = list(beta = rep(1, 3), gamma = gmma, radii = r),
                  w = c(0, 1, 0, 1),
                  trend = list(trend1, trend2, trend3))

mod_large <- list(cif = "straussm",
                  par = list(beta = rep(1, 3), gamma = gmma, radii = r),
                  w = c(0, 2, 0, 2),
                  trend = list(trend1, trend2, trend3))

X_small <- rmh(model = mod_small, start = list(n.start = ns),
               control = list(ptypes = rep(1/3, 3),
                              nrep = nr, nverb = nv, track = T))
X_large <- rmh(model = mod_large, start = list(n.start = ns),
               control = list(ptypes = rep(1/3, 3),
                              nrep = nr, nverb = nv, track = T))

X_small$marks <- factor(X_small$marks, levels = 1:3)
X_large$marks <- factor(X_large$marks, levels = 1:3)
target_param <- c(0, 0, 0.5, -0.5, log(gmma[lower.tri(gmma, diag = T)])/2)

nsim=1800
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA)[-1]
for (i in 1:nsim){
  set.seed(i)
  print(i)
  X_small <- rmh(model = mod_small, start = list(x.start = X_small),
                 control = list(ptypes = rep(1/3, 3),
                                nrep = 1e+5))
  X_large <- rmh(model = mod_large, start = list(x.start = X_large),
                 control = list(ptypes = rep(1/3, 3),
                                nrep = 1e+5))

  suppressWarnings(
    temp_small <- SemiMarkov(X = X_small,
                             covariate = Spat_covar[unit_window],
                             R_within = 0.02,
                             R_between = 0.04,
                             standardize = T,
                             sat = Inf)
  )

  suppressWarnings(
    temp_large <- SemiMarkov(X = X_large,
                             covariate = Spat_covar,
                             R_within = 0.02,
                             R_between = 0.04,
                             standardize = T,
                             sat = Inf)
  )

  data.table(Parameter = names(temp_small$betahat),
             Estimate = temp_small$betahat,
             Std.error = temp_small$std_err,
             within_CI = target_param > temp_small$CI$Lower_CI &
                         target_param < temp_small$CI$Upper_CI,
             window_size = 1) -> res_small

  data.table(Parameter = names(temp_large$betahat),
             Estimate = temp_large$betahat,
             Std.error = temp_large$std_err,
             within_CI = target_param > temp_large$CI$Lower_CI &
               target_param < temp_large$CI$Upper_CI,
             window_size = 4) -> res_large

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,.(Estimate = mean(Estimate), Std.error = mean(Std.error),
                      MC_std.error = sd(Estimate), Coverage = mean(within_CI)),
                   by = list(Parameter, window_size)]

data.table(Parameter = res$Parameter,
           Target_param = target_param,
           Estimate = res$Estimate,
           Lower = res$Estimate - 1.96*res$MC_std.error/sqrt(nsim),
           Upper = res$Estimate + 1.96*res$MC_std.error/sqrt(nsim)) -> dt

fwrite(dt, "Strauss results/Sim_Strauss_parameter_estimate_MC.csv")
fwrite(res, "Strauss results/Strauss_res_summary.csv")
fwrite(Results_sim, "Strauss results/Strauss_full_results.csv")

# Geyer saturation -------------------------------------------------------------
set.seed(123)
ns=1000
nr=1e+7
nv=1e+6

phi0 <- phi0/1.6

trend1 <- function(x, y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts]*exp(-log(3) + 0.5*Spat_covar[pts]))
}

trend2 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts]*exp(-log(5) - 0.5*Spat_covar[pts]))
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts])
}

# Saturation process 1
mod1_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend1)
X1_small <- rmh(model = mod1_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod1_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 5),
                   w = c(0,2,0,2),
                   trend = trend1)
X1_large <- rmh(model = mod1_large, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))
X1_small$marks <- rep(1, X1_small$n)
X1_large$marks <- rep(1, X1_large$n)

# Saturation process 2
mod2_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend2)
X2_small <- rmh(model = mod2_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod2_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 5),
                   w = c(0,2,0,2),
                   trend = trend2)
X2_large <- rmh(model = mod2_large, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))
X2_small$marks <- rep(2, X2_small$n)
X2_large$marks <- rep(2, X2_large$n)

# Saturation process 3
mod3_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend3)
X3_small <- rmh(model = mod3_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod3_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 5),
                   w = c(0,2,0,2),
                   trend = trend3)
X3_large <- rmh(model = mod3_large, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))
X3_small$marks <- rep(3, X3_small$n)
X3_large$marks <- rep(3, X3_large$n)

X_small <- superimpose(X1_small, X2_small, X3_small)
X_small$marks <- factor(X_small$marks)
X_large <- superimpose(X1_large, X2_large, X3_large)
X_large$marks <- factor(X_large$marks)

target_param <- c( c(-log(3), -log(5), 0.5, -0.5), log(c(1.1, 1, 1, 1.2, 1, 0.8)) )

nsim=1800
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA)[-1]
for (i in 1:nsim){
  set.seed(i)
  print(i)

  X1_small <- rmh(model = mod1_small, start = list(x.start = X1_small),
                 control = list(nrep = 1e+5))
  X1_large <- rmh(model = mod1_large, start = list(x.start = X1_large),
                 control = list(nrep = 1e+5))
  X1_small$marks <- rep(1, X1_small$n)
  X1_large$marks <- rep(1, X1_large$n)

  X2_small <- rmh(model = mod2_small, start = list(x.start = X2_small),
                  control = list(nrep = 1e+5))
  X2_large <- rmh(model = mod2_large, start = list(x.start = X2_large),
                  control = list(nrep = 1e+5))
  X2_small$marks <- rep(2, X2_small$n)
  X2_large$marks <- rep(2, X2_large$n)

  X3_small <- rmh(model = mod3_small, start = list(x.start = X3_small),
                  control = list(nrep = 1e+5))
  X3_large <- rmh(model = mod3_large, start = list(x.start = X3_large),
                  control = list(nrep = 1e+5))
  X3_small$marks <- rep(3, X3_small$n)
  X3_large$marks <- rep(3, X3_large$n)

  X_small <- superimpose(X1_small, X2_small, X3_small)
  X_small$marks <- factor(X_small$marks)
  X_large <- superimpose(X1_large, X2_large, X3_large)
  X_large$marks <- factor(X_large$marks)

  suppressWarnings(
    temp_small <- SemiMarkov(X = X_small,
                             covariate = Spat_covar[unit_window],
                             R_within = 0.02,
                             R_between = 0.04,
                             standardize = T,
                             sat = 5)
  )

  suppressWarnings(
    temp_large <- SemiMarkov(X = X_large,
                             covariate = Spat_covar,
                             R_within = 0.02,
                             R_between = 0.04,
                             standardize = T,
                             sat = 5)
  )

  data.table(Parameter = names(temp_small$betahat),
             Estimate = temp_small$betahat,
             Std.error = temp_small$std_err,
             within_CI = target_param > temp_small$CI$Lower_CI &
               target_param < temp_small$CI$Upper_CI,
             window_size = 1) -> res_small

  data.table(Parameter = names(temp_large$betahat),
             Estimate = temp_large$betahat,
             Std.error = temp_large$std_err,
             within_CI = target_param > temp_large$CI$Lower_CI &
               target_param < temp_large$CI$Upper_CI,
             window_size = 4) -> res_large

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,.(Estimate = mean(Estimate, na.rm = T),
                      Std.error = mean(Std.error, na.rm = T),
                      MC_std.error = sd(Estimate),
                      Coverage = mean(within_CI, na.rm = T)),
                   by = list(Parameter, window_size)]

data.table(Parameter = res$Parameter,
           Target_param = target_param,
           Estimate = res$Estimate,
           Lower = res$Estimate - 1.96*res$MC_std.error/sqrt(nsim),
           Upper = res$Estimate + 1.96*res$MC_std.error/sqrt(nsim)) -> dt

fwrite(dt, "Geyer results/Sim_Geyer_parameter_estimate_MC.csv")
fwrite(res, "Geyer results/Geyer_res_summary.csv")
fwrite(Results_sim, "Geyer results/Geyer_full_results.csv")


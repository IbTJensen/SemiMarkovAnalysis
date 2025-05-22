library(spatstat)
library(data.table)
library(geoR)
library(SemiMarkov)
# devtools::document("~/Git/SemiMarkov")
# devtools::load_all("~/Git/SemiMarkov")

# Setting up spatial covariates ------------------------------------------------
set.seed(123)
grid_pts <- expand.grid(x = seq(0, 4, 0.04), y = seq(0, 4, 0.04))

X <- grf(1, grid = grid_pts, cov.model = "exponential",
         cov.pars = c(sigmasq = 0.2, phi = 0.1))

Y <- grf(1, grid = grid_pts, cov.model = "exponential",
         cov.pars = c(sigmasq = 30^2, phi = 0.1), mean = 350)

X_mat <- matrix(NA, 101, 101)
Y_mat <- matrix(NA, 101, 101)
rownames(Y_mat) <- seq(0, 4, 0.04)
rownames(X_mat) <- seq(0, 4, 0.04)
colnames(Y_mat) <- seq(0, 4, 0.04)
colnames(X_mat) <- seq(0, 4, 0.04)
for(i in 1:101^2){
  iname <- as.character(grid_pts[i,1])
  jname <- as.character(grid_pts[i,2])
  X_mat[iname,jname] <- X$data[i]
  Y_mat[iname,jname] <- Y$data[i]
}

X_im <- im(X_mat, xcol = seq(0, 4, 0.04), yrow = seq(0, 4, 0.04))
phi0 <- im(Y_mat, xcol = seq(0, 4, 0.04), yrow = seq(0, 4, 0.04))

Spat_covar <- X_im

# Geyer saturation -------------------------------------------------------------
set.seed(123)
ns=1000
nr=1e+7
nv=1e+6

phi0 <- phi0

trend1 <- function(x, y){
  pts <- ppp(x = x, y = y, owin(c(0,4), c(0,4)))
  return(phi0[pts]*exp(-log(3) + 0.5*Spat_covar[pts]))
}

trend2 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,4), c(0,4)))
  return(phi0[pts]*exp(-log(5) - 0.5*Spat_covar[pts]))
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,4), c(0,4)))
  return(phi0[pts])
}

# Saturation process 1
mod1 <- list(cif = "geyer",
             par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 5),
             w = c(0,4,0,4),
             trend = trend1)
X1 <- rmh(model = mod1, start = list(n.start = ns),
          control = list(nrep = nr, nverb = nv, track = T))
X1$marks <- rep(1, X1$n)

# Saturation process 2
mod2 <- list(cif = "geyer",
             par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 5),
             w = c(0,4,0,4),
             trend = trend2)
X2 <- rmh(model = mod2, start = list(n.start = ns),
          control = list(nrep = nr, nverb = nv, track = T))
X2$marks <- rep(2, X2$n)

# Saturation process 3
mod3 <- list(cif = "geyer",
             par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 5),
             w = c(0,4,0,4),
             trend = trend3)
X3 <- rmh(model = mod3, start = list(n.start = ns),
          control = list(nrep = nr, nverb = nv, track = T))
X3$marks <- rep(3, X3$n)

X <- superimpose(X1, X2, X3)
X$marks <- factor(X$marks)

target_param <- c( c(-log(3), -log(5), 0.5, -0.5),
                   log(c(1.1, 1, 1, 1.2, 1, 0.8)) )

nsim=1800
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA)[-1]
for (i in 1:nsim){
  set.seed(i)
  print(i)

  X1 <- rmh(model = mod1, start = list(x.start = X1),
            control = list(nrep = 1e+6))
  X1$marks <- rep(1, X1$n)

  X2 <- rmh(model = mod2, start = list(x.start = X2),
            control = list(nrep = 1e+6))
  X2$marks <- rep(2, X2$n)

  X3 <- rmh(model = mod3, start = list(x.start = X3),
            control = list(nrep = 1e+6))
  X3$marks <- rep(3, X3$n)

  X <- superimpose(X1, X2, X3)
  X$marks <- factor(X$marks)

  suppressWarnings(
    temp <- SemiMarkov(X = X,
                       covariate = Spat_covar,
                       R_within = 0.02,
                       R_between = 0.04,
                       standardize = T,
                       sat = 5)
  )

  data.table(Parameter = names(temp$betahat),
             Estimate = temp$betahat,
             Std.error = temp$std_err,
             within_CI = target_param > temp$CI$Lower_CI &
               target_param < temp$CI$Upper_CI,
             window_size = 16) -> res

  Results_sim <- rbind(Results_sim, res)
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

fwrite(dt, "Geyer results/Sim_Geyer_parameter_estimate_MC_large.csv")
fwrite(res, "Geyer results/Geyer_res_summary_large.csv")
fwrite(Results_sim, "Geyer results/Geyer_full_results_large.csv")

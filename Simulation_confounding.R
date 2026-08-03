library(spatstat)
library(data.table)
library(geoR)
library(SemiMarkov)

# Setting up spatial covariates ------------------------------------------------
set.seed(123)
grid_pts <- expand.grid(x = seq(0, 2, 0.02), y = seq(0, 2, 0.02))

X <- grf(
  1,
  grid = grid_pts,
  cov.model = "exponential",
  cov.pars = c(sigmasq = 0.2, phi = 0.1)
)

Y <- grf(
  1,
  grid = grid_pts,
  cov.model = "exponential",
  cov.pars = c(sigmasq = 30^2, phi = 0.1),
  mean = 350
)

Z <- grf(
  1,
  grid = grid_pts,
  cov.model = "exponential",
  cov.pars = c(sigmasq = 0.1, phi = 0.2)
)


X_mat <- matrix(NA, 101, 101)
Y_mat <- matrix(NA, 101, 101)
Z_mat <- matrix(NA, 101, 101)
rownames(Y_mat) <- seq(0, 2, 0.02)
rownames(X_mat) <- seq(0, 2, 0.02)
rownames(Z_mat) <- seq(0, 2, 0.02)
colnames(Y_mat) <- seq(0, 2, 0.02)
colnames(X_mat) <- seq(0, 2, 0.02)
colnames(Z_mat) <- seq(0, 2, 0.02)
for (i in 1:101^2) {
  iname <- as.character(grid_pts[i, 1])
  jname <- as.character(grid_pts[i, 2])
  X_mat[iname, jname] <- X$data[i]
  Y_mat[iname, jname] <- Y$data[i]
  Z_mat[iname, jname] <- Z$data[i]
}

X_im <- im(X_mat, xcol = seq(0, 2, 0.02), yrow = seq(0, 2, 0.02))
phi0_init <- im(Y_mat, xcol = seq(0, 2, 0.02), yrow = seq(0, 2, 0.02))
Z_im <- im(Z_mat, xcol = seq(0, 2, 0.02), yrow = seq(0, 2, 0.02))

Spat_covar <- X_im
spat_covar_con_init <- Z_im
unit_window <- owin(c(0,1), c(0,1))
eps <- 1.5

# Poisson process --------------------------------------------------------------
phi0 <- phi0_init
phi0_log <- log(phi0)
phi0_log_std <- (phi0_log - mean(phi0_log))/sd(phi0_log)
spat_covar_con <- spat_covar_con_init + eps * phi0_log_std

intercepts <- c(0, 0)
param_spat <- c(0.25, -0.25)
param_spat_con <- c(-0.5, 0.5)

lambda_large1 <- phi0 *
  exp(
    intercepts[1] +
    param_spat[1] * Spat_covar +
    param_spat_con[1] * spat_covar_con
  )

lambda_large2 <- phi0 *
  exp(
    intercepts[2] +
    param_spat[2] * Spat_covar +
    param_spat_con[2] * spat_covar_con
  )

lambda_large3 <- phi0

lambda_small1 <- lambda_large1[unit_window]
lambda_small2 <- lambda_large2[unit_window]
lambda_small3 <- lambda_large3[unit_window]

target_param <- c(intercepts, param_spat, param_spat_con, 0, 0, 0, 0, 0, 0)

nsim <- 1800
Results_sim <- data.table(
  Parameter = NA,
  Estimate = NA,
  True_value = NA,
  Std.error = NA,
  within_CI = NA,
  window_size = NA,
  n1 = NA,
  n2 = NA,
  n3 = NA
)[-1]

for (i in 1:nsim) {
  set.seed(i)
  cat(i, "\r")
  X_small1 <- rpoispp(lambda_small1)
  X_small2 <- rpoispp(lambda_small2)
  X_small3 <- rpoispp(lambda_small3)
  X_small <- superimpose(X_small1, X_small2, X_small3)
  marks <- rep(1, X_small$n)
  marks[(X_small1$n + 1):(X_small1$n + X_small2$n)] <- 2
  marks[(X_small1$n + X_small2$n + 1):X_small$n] <- 3
  X_small$marks <- factor(marks)

  X_large1 <- rpoispp(lambda_large1)
  X_large2 <- rpoispp(lambda_large2)
  X_large3 <- rpoispp(lambda_large3)
  X_large <- superimpose(X_large1, X_large2, X_large3)
  marks <- rep(1, X_large$n)
  marks[(X_large1$n + 1):(X_large1$n + X_large2$n)] <- 2
  marks[(X_large1$n + X_large2$n + 1):X_large$n] <- 3
  X_large$marks <- factor(marks)

  covar_mat_small <- data.table(
    xcoord = X_small$x,
    ycoord = X_small$y,
    Covariate = Spat_covar[X_small],
    Covariate_con = spat_covar_con[X_small]
  )

  covar_mat_large <- data.table(
    xcoord = X_large$x,
    ycoord = X_large$y,
    Covariate = Spat_covar[X_large],
    Covariate_con = spat_covar_con[X_large]
  )

  suppressWarnings(
    temp_small <- SemiMarkov(
      X = X_small,
      covariate = covar_mat_small,
      R_within = 0.02,
      R_between = 0.04
    )
  )

  suppressWarnings(
    temp_large <- SemiMarkov(
      X = X_large,
      covariate = covar_mat_large,
      R_within = 0.02,
      R_between = 0.04
    )
  )

  n_small <- table(X_small$marks)
  n_large <- table(X_large$marks)

  res_small <- data.table(
    Parameter = names(temp_small$betahat),
    Estimate = temp_small$betahat,
    True_value = target_param,
    Std.error = temp_small$std_err,
    within_CI = target_param > temp_small$CI$Lower_CI &
      target_param < temp_small$CI$Upper_CI,
    window_size = 1,
    n1 = n_small[1],
    n2 = n_small[2],
    n3 = n_small[3]
  )

  res_large <- data.table(
    Parameter = names(temp_large$betahat),
    Estimate = temp_large$betahat,
    True_value = target_param,
    Std.error = temp_large$std_err,
    within_CI = target_param > temp_large$CI$Lower_CI &
      target_param < temp_large$CI$Upper_CI,
    window_size = 4,
    n1 = n_large[1],
    n2 = n_large[2],
    n3 = n_large[3]
  )

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,
  .(
    Estimate = mean(Estimate),
    Std.error = mean(Std.error),
    MC_std.error = sd(Estimate),
    Coverage = mean(within_CI),
    n1 = mean(n1),
    n2 = mean(n2),
    n3 = mean(n3)
  ),
  by = list(Parameter, window_size, True_value)
]

dt <- data.table(
  Parameter = res$Parameter,
  Target_param = target_param,
  Estimate = res$Estimate,
  Lower = res$Estimate - 1.96 * res$MC_std.error / sqrt(nsim),
  Upper = res$Estimate + 1.96 * res$MC_std.error / sqrt(nsim)
)

nm1 <- paste(
  "Sim", "Poisson", "parameter", "estimate", "MC", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm2 <- paste(
  "Sim", "Poisson", "res", "summary", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm3 <- paste(
  "Sim", "Poisson", "full", "results", "conf",
  as.character(Sys.Date()),
  sep = "_"
)

fwrite(dt, paste0("Poisson results/", nm1, ".csv") )
fwrite(res, paste0("Poisson results/", nm2, ".csv"))
fwrite(Results_sim, paste0("Poisson results/", nm3, ".csv"))

# Multi-Strauss ----------------------------------------------------------------
set.seed(123)
ns <- 1000
nr <- 1e7
nv <- 1e6

# Non-parametric baseline
phi0 <- 1.6*phi0_init

# Confounded covariate
phi0_log <- log(phi0)
phi0_log_std <- (phi0_log - mean(phi0_log))/sd(phi0_log)
spat_covar_con <- spat_covar_con_init + eps * phi0_log_std

# First-order parameters
intercepts <- c(0, 0)
param_spat <- c(0.25, -0.25)
param_spat_con <- c(-0.5, 0.5)

# Spatial trends
trend1 <- function(x, y) {
  pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
  trend <- phi0[pts] * exp(
    intercepts[1] +
    param_spat[1] * Spat_covar[pts] +
    param_spat_con[1] * spat_covar_con[pts]
  )
  return(trend)
}

trend2 <- function(x, y) {
  pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
  trend <- phi0[pts] * exp(
    intercepts[2] +
    param_spat[2] * Spat_covar[pts] +
    param_spat_con[2] * spat_covar_con[pts]
  )
  return(trend)
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts])
}

# Interaction parameters
btw_int <- c(0.9, 0.9, 0.9)
wtn_int <- c(0.8, 0.8, 0.8)
gmma <- matrix(0, 3, 3)
gmma[lower.tri(gmma)] <- btw_int
gmma <- gmma + t(gmma)
diag(gmma) <- wtn_int

# Initiazing Markov-chain
r <- matrix(0.04, 3, 3)
diag(r) <- 0.02
mod_small <- list(
  cif = "straussm",
  par = list(beta = rep(1, 3), gamma = gmma, radii = r),
  w = c(0, 1, 0, 1),
  trend = list(trend1, trend2, trend3)
)

mod_large <- list(
  cif = "straussm",
  par = list(beta = rep(1, 3), gamma = gmma, radii = r),
  w = c(0, 2, 0, 2),
  trend = list(trend1, trend2, trend3)
)

X_small <- rmh(
  model = mod_small,
  start = list(n.start = ns),
  control = list(ptypes = rep(1 / 3, 3), nrep = nr, nverb = nv, track = T)
)
X_large <- rmh(
  model = mod_large,
  start = list(n.start = ns),
  control = list(ptypes = rep(1 / 3, 3), nrep = nr, nverb = nv, track = T)
)

X_small$marks <- factor(X_small$marks, levels = 1:3)
X_large$marks <- factor(X_large$marks, levels = 1:3)
target_param <- c(
  intercepts, param_spat, param_spat_con, 
  log(gmma[lower.tri(gmma, diag = T)])/2
)

nsim <- 1800
Results_sim <- data.table(
  Parameter = NA,
  Estimate = NA,
  True_value = NA,
  Std.error = NA,
  within_CI = NA,
  window_size = NA,
  n1 = NA,
  n2 = NA,
  n3 = NA
)[-1]

for (i in 1:nsim) {
  set.seed(i)
  print(i)
  X_small <- rmh(
    model = mod_small,
    start = list(x.start = X_small),
    control = list(ptypes = rep(1 / 3, 3), nrep = 1e+5)
  )
  X_large <- rmh(
    model = mod_large,
    start = list(x.start = X_large),
    control = list(ptypes = rep(1 / 3, 3), nrep = 1e+5)
  )

  covar_mat_small <- data.table(
    xcoord = X_small$x,
    ycoord = X_small$y,
    Covariate = Spat_covar[X_small],
    Covariate_con = spat_covar_con[X_small]
  )

  covar_mat_large <- data.table(
    xcoord = X_large$x,
    ycoord = X_large$y,
    Covariate = Spat_covar[X_large],
    Covariate_con = spat_covar_con[X_large]
  )

  suppressWarnings(
    temp_small <- SemiMarkov(
      X = X_small,
      covariate = covar_mat_small,
      R_within = 0.02,
      R_between = 0.04,
      standardize = T,
      sat = Inf
    )
  )

  suppressWarnings(
    temp_large <- SemiMarkov(
      X = X_large,
      covariate = covar_mat_large,
      R_within = 0.02,
      R_between = 0.04,
      standardize = T,
      sat = Inf
    )
  )

  n_small <- table(X_small$marks)
  n_large <- table(X_large$marks)

  res_small <- data.table(
    Parameter = names(temp_small$betahat),
    Estimate = temp_small$betahat,
    True_value = target_param,
    Std.error = temp_small$std_err,
    within_CI = target_param > temp_small$CI$Lower_CI &
      target_param < temp_small$CI$Upper_CI,
    window_size = 1,
    n1 = n_small[1],
    n2 = n_small[2],
    n3 = n_small[3]
  )

  res_large <- data.table(
    Parameter = names(temp_large$betahat),
    Estimate = temp_large$betahat,
    True_value = target_param,
    Std.error = temp_large$std_err,
    within_CI = target_param > temp_large$CI$Lower_CI &
      target_param < temp_large$CI$Upper_CI,
    window_size = 4,
    n1 = n_large[1],
    n2 = n_large[2],
    n3 = n_large[3]
  )

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,
  .(
    Estimate = mean(Estimate),
    Std.error = mean(Std.error),
    MC_std.error = sd(Estimate),
    Coverage = mean(within_CI),
    n1 = mean(n1),
    n2 = mean(n2),
    n3 = mean(n3)
  ),
  by = list(Parameter, window_size, True_value)
]

dt <- data.table(
  Parameter = res$Parameter,
  Target_param = target_param,
  Estimate = res$Estimate,
  Lower = res$Estimate - 1.96 * res$MC_std.error / sqrt(nsim),
  Upper = res$Estimate + 1.96 * res$MC_std.error / sqrt(nsim)
)

nm1 <- paste(
  "Sim", "Strauss", "parameter", "estimate", "MC", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm2 <- paste(
  "Sim", "Strauss", "res", "summary", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm3 <- paste(
  "Sim", "Strauss", "full", "results", "conf",
  as.character(Sys.Date()),
  sep = "_"
)

fwrite(dt, paste0("Strauss results/", nm1, ".csv") )
fwrite(res, paste0("Strauss results/", nm2, ".csv"))
fwrite(Results_sim, paste0("Strauss results/", nm3, ".csv"))

# Geyer saturation -------------------------------------------------------------
set.seed(123)
ns <- 1000
nr <- 1e+7
nv <- 1e+6

# Non-parametric baseline
phi0 <- phi0_init*1.3

# Confounded covariate
phi0_log <- log(phi0)
phi0_log_std <- (phi0_log - mean(phi0_log))/sd(phi0_log)
spat_covar_con <- spat_covar_con_init + eps * phi0_log_std

# First-order parameters
intercepts <- c(-log(1.4), -log(3.2))
param_spat <- c(0.25, -0.25)
param_spat_con <- c(-0.5, 0.5)

# Spatial trends
trend1 <- function(x, y) {
  pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
  trend <- phi0[pts] * exp(
    intercepts[1] +
    param_spat[1] * Spat_covar[pts] +
    param_spat_con[1] * spat_covar_con[pts]
  )
  return(trend)
}

trend2 <- function(x, y) {
  pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
  trend <- phi0[pts] * exp(
    intercepts[2] +
    param_spat[2] * Spat_covar[pts] +
    param_spat_con[2] * spat_covar_con[pts]
  )
  return(trend)
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(phi0[pts])
}

# Saturation process 1
mod1_small <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 10),
  w = c(0, 1, 0, 1),
  trend = trend1
)
X1_small <- rmh(
  model = mod1_small,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)

mod1_large <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 10),
  w = c(0, 2, 0, 2),
  trend = trend1
)
X1_large <- rmh(
  model = mod1_large,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)
X1_small$marks <- rep(1, X1_small$n)
X1_large$marks <- rep(1, X1_large$n)

# Saturation process 2
mod2_small <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 10),
  w = c(0, 1, 0, 1),
  trend = trend2
)
X2_small <- rmh(
  model = mod2_small,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)

mod2_large <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 10),
  w = c(0, 2, 0, 2),
  trend = trend2
)
X2_large <- rmh(
  model = mod2_large,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)
X2_small$marks <- rep(2, X2_small$n)
X2_large$marks <- rep(2, X2_large$n)

# Saturation process 3
mod3_small <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 10),
  w = c(0, 1, 0, 1),
  trend = trend3
)
X3_small <- rmh(
  model = mod3_small,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)

mod3_large <- list(
  cif = "geyer",
  par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 10),
  w = c(0, 2, 0, 2),
  trend = trend3
)
X3_large <- rmh(
  model = mod3_large,
  start = list(n.start = ns),
  control = list(nrep = nr, nverb = nv, track = T)
)
X3_small$marks <- rep(3, X3_small$n)
X3_large$marks <- rep(3, X3_large$n)

X_small <- superimpose(X1_small, X2_small, X3_small)
X_small$marks <- factor(X_small$marks)
X_large <- superimpose(X1_large, X2_large, X3_large)
X_large$marks <- factor(X_large$marks)

target_param <- c( 
  c(intercepts, param_spat, param_spat_con),
  log( c(1.1, 1, 1, 1.2, 1, 0.8) ) )

nsim <- 1800
Results_sim <- data.table(
  Parameter = NA,
  Estimate = NA,
  True_value = NA,
  Std.error = NA,
  within_CI = NA,
  window_size = NA,
  n1 = NA,
  n2 = NA,
  n3 = NA
)[-1]

for (i in 1:nsim) {
  set.seed(i)
  print(i)

  X1_small <- rmh(
    model = mod1_small,
    start = list(x.start = X1_small),
    control = list(nrep = 1e+5)
  )
  X1_large <- rmh(
    model = mod1_large,
    start = list(x.start = X1_large),
    control = list(nrep = 1e+5)
  )
  X1_small$marks <- rep(1, X1_small$n)
  X1_large$marks <- rep(1, X1_large$n)

  X2_small <- rmh(
    model = mod2_small,
    start = list(x.start = X2_small),
    control = list(nrep = 1e+5)
  )
  X2_large <- rmh(
    model = mod2_large,
    start = list(x.start = X2_large),
    control = list(nrep = 1e+5)
  )
  X2_small$marks <- rep(2, X2_small$n)
  X2_large$marks <- rep(2, X2_large$n)

  X3_small <- rmh(
    model = mod3_small,
    start = list(x.start = X3_small),
    control = list(nrep = 1e+5)
  )
  X3_large <- rmh(
    model = mod3_large,
    start = list(x.start = X3_large),
    control = list(nrep = 1e+5)
  )
  X3_small$marks <- rep(3, X3_small$n)
  X3_large$marks <- rep(3, X3_large$n)

  X_small <- superimpose(X1_small, X2_small, X3_small)
  X_small$marks <- factor(X_small$marks)
  X_large <- superimpose(X1_large, X2_large, X3_large)
  X_large$marks <- factor(X_large$marks)

  covar_mat_small <- data.table(
    xcoord = X_small$x,
    ycoord = X_small$y,
    Covariate = Spat_covar[X_small],
    Covariate_con = spat_covar_con[X_small]
  )

  covar_mat_large <- data.table(
    xcoord = X_large$x,
    ycoord = X_large$y,
    Covariate = Spat_covar[X_large],
    Covariate_con = spat_covar_con[X_large]
  )

  suppressWarnings(
    temp_small <- SemiMarkov(
      X = X_small,
      covariate = covar_mat_small,
      R_within = 0.02,
      R_between = 0.04,
      standardize = T,
      sat = 10
    )
  )

  suppressWarnings(
    temp_large <- SemiMarkov(
      X = X_large,
      covariate = covar_mat_large,
      R_within = 0.02,
      R_between = 0.04,
      standardize = T,
      sat = 10
    )
  )

  n_small <- table(X_small$marks)
  n_large <- table(X_large$marks)

  res_small <- data.table(
    Parameter = names(temp_small$betahat),
    Estimate = temp_small$betahat,
    True_value = target_param,
    Std.error = temp_small$std_err,
    within_CI = target_param > temp_small$CI$Lower_CI &
      target_param < temp_small$CI$Upper_CI,
    window_size = 1,
    n1 = n_small[1],
    n2 = n_small[2],
    n3 = n_small[3]
  )

  res_large <- data.table(
    Parameter = names(temp_large$betahat),
    Estimate = temp_large$betahat,
    True_value = target_param,
    Std.error = temp_large$std_err,
    within_CI = target_param > temp_large$CI$Lower_CI &
      target_param < temp_large$CI$Upper_CI,
    window_size = 4,
    n1 = n_large[1],
    n2 = n_large[2],
    n3 = n_large[3]
  )

  Results_sim <- rbind(Results_sim, res_small, res_large)
}
Sys.time()

res <- Results_sim[,
  .(
    Estimate = mean(Estimate, na.rm = T),
    Std.error = mean(Std.error, na.rm = T),
    MC_std.error = sd(Estimate),
    Coverage = mean(within_CI, na.rm = T),
    n1 = mean(n1),
    n2 = mean(n2),
    n3 = mean(n3)
  ),
  by = list(Parameter, True_value, window_size)
]

dt <- data.table(
  Parameter = res$Parameter,
  Target_param = target_param,
  Estimate = res$Estimate,
  Lower = res$Estimate - 1.96 * res$MC_std.error / sqrt(nsim),
  Upper = res$Estimate + 1.96 * res$MC_std.error / sqrt(nsim)
)

nm1 <- paste(
  "Sim", "Geyer", "parameter", "estimate", "MC", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm2 <- paste(
  "Sim", "Geyer", "res", "summary", "conf",
  as.character(Sys.Date()),
  sep = "_"
)
nm3 <- paste(
  "Sim", "Geyer", "full", "results", "conf",
  as.character(Sys.Date()),
  sep = "_"
)

fwrite(dt, paste0("Geyer results/", nm1, ".csv"))
fwrite(res, paste0("Geyer results/", nm2, ".csv"))
fwrite(Results_sim, paste0("Geyer results/", nm3, ".csv"))

# Int_22_est <- Results_sim[Parameter == "2-2" & window_size == 1, Estimate]
# plot(density(Int_22_est))


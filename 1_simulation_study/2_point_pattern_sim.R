library(spatstat)
library(data.table)
library(geoR)

# Parameters -------------------------------------------------------------------
intercepts <- c(0, 0)
intercepts_geyer <- c(-log(1.4), -log(1.6))
param_spat <- c(0.5, -0.5)
param_spat_con <- c(-0.5, 0.5)

# Interaction parameters Poisson
int_pois <- rep(0, 6)

# Interaction parameters Strauss
btw_int <- c(0.9, 0.9, 0.9)
wtn_int <- c(0.8, 0.8, 0.8)
gmma <- matrix(0, 3, 3)
gmma[lower.tri(gmma)] <- btw_int
gmma <- gmma + t(gmma)
diag(gmma) <- wtn_int
int_strauss <- log(gmma[lower.tri(gmma, diag = T)])/2

# Interaction parameters Geyer
int_geyer <- c(log(1.1), 0, 0, log(1.2), 0, log(0.8))

# Spatial covariates -----------------------------------------------------------
# Loading random fields
spat_covar <- readRDS("1_simulation_study/simulated_data/X_im.rds")
spat_covar_con_init <- readRDS("1_simulation_study/simulated_data/Z_im.rds")
phi0_init <- readRDS("1_simulation_study/simulated_data/phi0_init.rds")

# Simulation settings
eps <- c(0.5, 4, 0)
int_model <- c("Poisson", "Strauss", "Geyer")
win_size <- c(1, 4)
settings <- expand.grid(interaction = int_model, win_size = win_size, eps = eps)
nsim <- 1800

# Names
settings_names <- paste(
  settings$interaction,
  "win_size",
  settings$win_size,
  "eps",
  settings$eps,
  sep = "_"
)

# Settings up lists of simulation settings
settings_list <- list()
for (i in 1:nrow(settings)) {
  common_fac <- fcase(
    settings$interaction[i] == "Poisson" , 1   ,
    settings$interaction[i] == "Strauss" , 1.6 ,
    settings$interaction[i] == "Geyer"   , 1.3
  )

  phi0 <- common_fac * phi0_init
  phi0_log <- log(phi0)
  phi0_log_std <- (phi0_log - mean(phi0_log)) / max(phi0_log - mean(phi0_log))
  # phi0_log_std <- (phi0_log - mean(phi0_log)) / sd(phi0_log)
  if(settings$eps[i] == -1){
    spat_covar_con <- NULL
  } else{
    spat_covar_con <- spat_covar_con_init + settings$eps[i] * phi0_log_std
  }

  settings_list[[i]] <- list(
    phi0 = phi0,
    spat_covar = spat_covar,
    spat_covar_con = spat_covar_con,
    eps = settings$eps[i],
    win_size = settings$win_size[i],
    intercepts = ifelse(
      rep(settings$interaction[i], 2) == "Geyer",
      intercepts_geyer,
      intercepts
    ),
    param_spat = param_spat,
    param_spat_con = if(settings$eps[i] == -1){NULL} else{param_spat_con},
    interaction = settings$interaction[i],
    setting_ids = settings_names[i],
    nsim = nsim
  )
}
names(settings_list) <- settings_names
saveRDS(
  settings_list,
  paste0(
    "1_simulation_study/simulated_data/sim_covar_and_settings_",
    Sys.Date(),
    ".rds"
  )
)

# Poisson process --------------------------------------------------------------
poisson_sim <- function(lst, n = NULL) {
  # Extracting variables from input list
  phi0 <- lst$phi0
  spat_covar <- lst$spat_covar
  spat_covar_con <- lst$spat_covar_con
  eps <- lst$eps
  win_size <- lst$win_size
  intercepts <- lst$intercepts
  param_spat <- lst$param_spat
  param_spat_con <- lst$param_spat_con
  nsim <- lst$nsim
  setting_ids <- lst$setting_ids

  window <- owin(c(0, sqrt(win_size)), c(0, sqrt(win_size)))
  target_param <- c(intercepts, param_spat, param_spat_con, int_pois)

  # Intensities per type
  lambda1 <- phi0[window] *
    exp(
      intercepts[1] +
        param_spat[1] * spat_covar[window] +
        if(eps == -1){0} else{param_spat_con[1] * spat_covar_con[window]}
    )

  lambda2 <- lst$phi0[window] *
    exp(
      intercepts[2] +
        param_spat[2] * spat_covar[window] +
        if(eps == -1){0} else{param_spat_con[2] * spat_covar_con[window]}
    )

  lambda3 <- lst$phi0[window]

  # Simulate point patterns
  if(!is.null(n)) cat(paste("Run", n), "\r")
  out <- list()
  for (i in 1:nsim) {
    if(!is.null(n)) cat(paste("Run", n, " - Iteration", i, "     "), "\r")
    X1 <- rpoispp(lambda1)
    X2 <- rpoispp(lambda2)
    X3 <- rpoispp(lambda3)
    X <- superimpose(X1, X2, X3)
    marks <- rep(1, X$n)
    marks[(X1$n + 1):(X1$n + X2$n)] <- 2
    marks[(X1$n + X2$n + 1):X$n] <- 3
    X$marks <- factor(marks)
    
    out[[i]] <- list(
      point_pattern = X,
      win_size = win_size,
      true_val = target_param,
      process = "Poisson",
      eps = eps,
      setting_ids = setting_ids
    )
  }

  return(out)
}

# Multi-Strauss ----------------------------------------------------------------
strauss_sim <- function(lst, n = NULL) {
  # Extracting variables from input list
  phi0 <- lst$phi0
  spat_covar <- lst$spat_covar
  spat_covar_con <- lst$spat_covar_con
  eps <- lst$eps
  win_size <- lst$win_size
  intercepts <- lst$intercepts
  param_spat <- lst$param_spat
  param_spat_con <- lst$param_spat_con
  nsim <- lst$nsim
  setting_ids <- lst$setting_ids

  window <- owin(c(0, sqrt(win_size)), c(0, sqrt(win_size)))
  target_param <- c(intercepts, param_spat, param_spat_con, int_strauss)

  # Spatial trend for type 1
  trend1 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    trend <- phi0[pts] *
      exp(
        intercepts[1] +
          param_spat[1] * spat_covar[pts] +
          if(eps == -1){0} else{param_spat_con[1] * spat_covar_con[pts]}
      )
    return(trend)
  }

  # Spatial trend for type 2
  trend2 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    trend <- phi0[pts] *
      exp(
        intercepts[2] +
          param_spat[2] * spat_covar[pts] +
          if(eps == -1){0} else{param_spat_con[2] * spat_covar_con[pts]}
      )
    return(trend)
  }

  # Spatial trend for type 3 (reference type)
  trend3 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    return(phi0[pts])
  }

  # Settings up model to be sampled from
  r <- matrix(0.04, 3, 3)
  diag(r) <- 0.02
  mod <- list(
    cif = "straussm",
    par = list(beta = rep(1, 3), gamma = gmma, radii = r),
    w = c(0, sqrt(win_size), 0, sqrt(win_size)),
    trend = list(trend1, trend2, trend3)
  )

  # MCMC Burn-in
  if(!is.null(n)) cat(paste("Run", n, "- Burn-in            "), "\r")
  X <- rmh(
    model = mod,
    start = list(n.start = 1000),
    control = list(ptypes = rep(1 / 3, 3), nrep = 1e+7),
    verbose = FALSE
  )
  
  X <- ppp(x = X$x, y = X$y, window = X$window, marks = X$marks)
  gc()

  # Simulating point patterns
  out <- list()
  for (i in 1:nsim) {
    if(!is.null(n)) cat(paste("Run", n, "- Iteration", i, "     "), "\r")
    X <- rmh(
      model = mod,
      start = list(x.start = X),
      control = list(ptypes = rep(1 / 3, 3), nrep = 1e+5),
      verbose = FALSE
    )

    X <- ppp(x = X$x, y = X$y, window = X$window, marks = X$marks)

    out[[i]] <- list(
      point_pattern = X,
      win_size = win_size,
      true_val = target_param,
      process = "Strauss",
      eps = eps,
      setting_ids = setting_ids
    )
  }

  return(out)
}

# Geyer saturation -------------------------------------------------------------
geyer_sim <- function(lst, n = NULL) {
  # Extracting variables from input list
  phi0 <- lst$phi0
  spat_covar <- lst$spat_covar
  spat_covar_con <- lst$spat_covar_con
  eps <- lst$eps
  win_size <- lst$win_size
  intercepts <- lst$intercepts
  param_spat <- lst$param_spat
  param_spat_con <- lst$param_spat_con
  nsim <- lst$nsim
  setting_ids <- lst$setting_ids

  window <- owin(c(0, sqrt(win_size)), c(0, sqrt(win_size)))
  target_param <- c(intercepts, param_spat, param_spat_con, int_geyer)

  # Spatial trend for type 1
  trend1 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    trend <- phi0[pts] *
      exp(
        intercepts[1] +
          param_spat[1] * spat_covar[pts] +
          if(eps == -1){0} else{param_spat_con[1] * spat_covar_con[pts]}
      )
    return(trend)
  }

  # Spatial trend for type 2
  trend2 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    trend <- phi0[pts] *
      exp(
        intercepts[2] +
          param_spat[2] * spat_covar[pts] +
          if(eps == -1){0} else{param_spat_con[2] * spat_covar_con[pts]}
      )
    return(trend)
  }

  # Spatial trend for type 3 (reference type)
  trend3 <- function(x, y) {
    pts <- ppp(x = x, y = y, owin(c(0, 2), c(0, 2)))
    return(phi0[pts])
  }

  # Settings up model to be sampled from (type 1)
  mod1 <- list(
    cif = "geyer",
    par = list(beta = 1, gamma = 1.1, r = 0.02, sat = 10),
    w = c(0, sqrt(win_size), 0, sqrt(win_size)),
    trend = trend1
  )
  
  # MCMC burn-in type 1
  if(!is.null(n)) cat(paste("Run", n, "- Burn-in 1            "), "\r")
  X1 <- rmh(
    model = mod1,
    start = list(n.start = 1000),
    control = list(nrep = 1e+7),
    verbose = FALSE
  )
  
  X1 <- ppp(x = X1$x, y = X1$y, window = X1$window, marks = X1$marks)
  gc()

  X1$marks <- rep(1, X1$n)

  # Settings up model to be sampled from (type 2)
  mod2 <- list(
    cif = "geyer",
    par = list(beta = 1, gamma = 1.2, r = 0.02, sat = 10),
    w = c(0, sqrt(win_size), 0, sqrt(win_size)),
    trend = trend2
  )
  
  # MCMC burn-in type 2
  if(!is.null(n)) cat(paste("Run", n, "- Burn-in 2            "), "\r")
  X2 <- rmh(
    model = mod2,
    start = list(n.start = 1000),
    control = list(nrep = 1e+7),
    verbose = FALSE
  )
  
  X2 <- ppp(x = X2$x, y = X2$y, window = X2$window, marks = X2$marks)
  gc()

  X2$marks <- rep(2, X2$n)

  # Settings up model to be sampled from (type 3)
  mod3 <- list(
    cif = "geyer",
    par = list(beta = 1, gamma = 0.8, r = 0.02, sat = 10),
    w = c(0, sqrt(win_size), 0, sqrt(win_size)),
    trend = trend3
  )

  # MCMC burn-in type 3
  if(!is.null(n)) cat(paste("Run", n, "- Burn-in 3            "), "\r")
  X3 <- rmh(
    model = mod3,
    start = list(n.start = 1000),
    control = list(nrep = 1e+7),
    verbose = FALSE
  )
  
  X3 <- ppp(x = X3$x, y = X3$y, window = X3$window, marks = X3$marks)
  gc()

  X3$marks <- rep(3, X3$n)

  # Simulate point patterns
  out <- list()
  for (i in 1:nsim) {
    if(!is.null(n)) cat(paste("Run", n, "- Iteration", i, "     "), "\r")
    X1 <- rmh(
      model = mod1,
      start = list(x.start = X1),
      control = list(nrep = 1e+5),
      verbose = FALSE
    )
    X1$marks <- rep(1, X1$n)

    X2 <- rmh(
      model = mod2,
      start = list(x.start = X2),
      control = list(nrep = 1e+5),
      verbose = FALSE
    )
    X2$marks <- rep(2, X2$n)

    X3 <- rmh(
      model = mod3,
      start = list(x.start = X3),
      control = list(nrep = 1e+5),
      verbose = FALSE
    )
    X3$marks <- rep(3, X3$n)

    X <- superimpose(X1, X2, X3)
    X$marks <- factor(X$marks)

    X <- ppp(
      x = X$x,
      y = X$y,
      window = X$window,
      marks = X$marks
    )

    out[[i]] <- list(
      point_pattern = X,
      win_size = win_size,
      true_val = target_param,
      process = "Geyer",
      eps = eps,
      setting_ids = setting_ids
    )
  }

  return(out)
}

# Simulation -------------------------------------------------------------------
point_pattern_sim <- function(lst, n){
  if(lst$interaction == "Poisson"){
    out <- poisson_sim(lst, n)
  }
  if(lst$interaction == "Strauss"){
    out <- strauss_sim(lst, n)
  }
  if(lst$interaction == "Geyer"){
    out <- geyer_sim(lst, n)
  }
  return(out)
}

start_time <- Sys.time(); start_time
all_res_list <- list()
for(i in 1:length(settings_list)){
  set.seed(i)
  all_res_list[[i]] <- point_pattern_sim(settings_list[[i]], i)
}
end_time <- Sys.time(); end_time
end_time - start_time

res_list <- Reduce(c, all_res_list)
print("Saving simulated point patterns")
saveRDS(
  res_list,
  paste0(
    "1_simulation_study/simulated_data/simulated_point_patterns_",
    Sys.Date(),
    ".rds"
  )
)

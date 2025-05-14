library(spatstat)
library(data.table)
library(SemiMarkov)
library(geoR)
# devtools::document("~/Git/SemiMarkov")
# devtools::load_all("~/Git/SemiMarkov")

# Setting up spatial covariates ------------------------------------------------
bei.extra$elev <- bei.extra$elev[owin(c(0,999), c(0,499))]
bei.extra$elev <- shift.im(bei.extra$elev, vec = c(2.5, 2.5))
bei.extra$elev <- rescale.im(bei.extra$elev, 250, "units")
bei.extra$elev <- bei.extra$elev[owin(c(0,2), c(0,2))]

bei.extra$grad <- bei.extra$grad[owin(c(0,999), c(0,499))]
bei.extra$grad <- shift.im(bei.extra$grad, vec = c(2.5, 2.5))
bei.extra$grad <- rescale.im(bei.extra$grad, 250, "units")
bei.extra$grad <- bei.extra$grad[owin(c(2,4), c(0,2))]
bei.extra$grad <- shift.im(bei.extra$grad, vec = c(-2, 0))


# Poisson process --------------------------------------------------------------
unit_elev <- bei.extra$elev[owin(c(0,1), c(0,1))]
unit_grad <- bei.extra$grad[owin(c(0,1), c(0,1))]

ll0=(unit_elev - mean(unit_elev))
ll0=ll0/max(ll0)

exp_points <- 0.02^2*sum(exp(ll0))
ll0 <- ll0 - log(exp_points)
ll0 <- ll0 + log(350)

lambda_small1=exp(ll0)
lambda_small2=exp(ll0 - 1 + unit_grad)
lambda_small3=exp(ll0 + 1 - unit_grad)

lambda_large1=exp(ll0 + log(4))
lambda_large2=exp(ll0 - 1 + log(4) + unit_grad)
lambda_large3=exp(ll0 + 1 + log(4) - unit_grad)

target_param <- c(0, 0, 1 ,2, 0, 0, 0, 0, 0, 0)

nsim=1800
# betahat=matrix(0,nsim,10)
# coverage=matrix(0,nsim,10)
# std_errors <- matrix(0,nsim,10)
# std_errors_S <- matrix(0,nsim,10)
# score_fct <- matrix(0,nsim,10)
# second_order_list <- list()
# sensitivity_list <- list()
# score_mat <- matrix(NA, nsim, 10)
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA, n = NA)[-1]
for (i in 1:nsim){
  set.seed(i)
  # a <- Sys.time()
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
    temp_small <- SemiMarkov_fixed_R(X = X_small,
                                     covariate = unit_grad,
                                     R_within = 0.05,
                                     R_between = 0.02)
  )

  suppressWarnings(
    temp_large <- SemiMarkov_fixed_R(X = X_large,
                                     covariate = unit_grad,
                                     R_within = 0.05,
                                     R_between = 0.02)
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

  # coverage[i,] <- target_param > temp$CI$Lower_CI & target_param < temp$CI$Upper_CI
  # std_errors[i,] <- temp$std_err
  # std_errors_S[i,] <- diag(solve(temp$Sensitivity))
  #
  # w <- temp$w[,1:14]
  # log_lambda <- as.matrix(temp$w[,5:14]) %*% target_param
  # w[,lambda:=exp(log_lambda)]
  # w_Lambda <- w[,.(Lambda = sum(lambda)),list(xcoord, ycoord, type_obs)]
  # w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  # w[,prob:=lambda/Lambda]
  #
  # h <- w[,lapply(.SD, function(x) x - sum(x*prob)),
  #        by = list(xcoord, ycoord, type_obs),
  #        .SDcols = 5:14]
  # h <- cbind(w[,1:4], h[,-(1:3)])
  # score_fct[i,] <- colSums(h[type_obs == j, -(1:4)]) # Score function in the true value
  #
  # Check <- SemiMarkov:::Covariance_true_val(X = X,
  #                                           covariate = covariate,
  #                                           edgecorrection = 10,
  #                                           R_within = 5,
  #                                           R_between = 2,
  #                                           symmetric = T,
  #                                           standardize = F,
  #                                           sat = Inf,
  #                                           Poisson = F,
  #                                           verbose = F,
  #                                           true.param = target_param)
  #
  # betahat[i,]=temp$betahat
  # second_order_list[[i]] <- Check$Sigma_term # Second order term in the true value
  # sensitivity_list[[i]] <- Check$S # Sensitivity in the true value
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

# data.table(Parameter = names(temp$betahat),
#            Estimate = colMeans(betahat),
#            "True value" = target_param,
#            Coverage = colMeans(coverage),
#            "Estimated std. error" = colMeans(std_errors),
#            "Monte-Carlo std. error" = apply(betahat, 2, sd) ) -> res
#
# abs(res$Estimate - target_param)/(res$`Estimated std. error`/sqrt(nsim))
#
# score_fct_var <- array(NA, dim = c(10, 10, nsim))
# for(i in 1:nsim){
#   score_fct_var[,,i] <- score_fct[i,] %*% t(score_fct[i,])
# }
# apply(score_fct_var, 1:2, mean)
# apply(score_fct_var, 1:2, sd)/sqrt(1000)
#
# Sens_array <- array(NA, dim = c(10, 10, nsim))
# for(i in 1:nsim){
#   Sens_array[,,i] <- sensitivity_list[[i]]
# }
#
# cov_mat_array <- array(NA, dim = c(10, 10, nsim))
# for(i in 1:nsim){
#   Sigma <- sensitivity_list[[i]] + second_order_list[[i]]
#   S <- sensitivity_list[[i]]
#   S_inv <- solve(S)
#   cov_mat_array[,,i] <- S_inv %*% Sigma %*% S_inv
# }
# apply(cov_mat_array, 1:2, mean)
# apply(cov_mat_array, 1:2, sd)
# cov(betahat)
# # cov(score_fct)
#
# A <- array(NA, dim = c(10, 10, nsim))
# for(i in 1:nsim){
#   A[,,i] <- second_order_list[[i]]
# }
#
# A_mean <- apply(A, 1:2, mean)
# A_se <- apply(A, 1:2, sd)/sqrt(nsim)
# abs(A_mean) > 1.96*A_se
# A_mean/A_se
#
# fwrite(res, "Poisson results/Poisson_res_summary.csv")
#
# fwrite(betahat, "Poisson results/betahat_poisson.csv")
# fwrite(coverage, "Poisson results/coverage_poisson.csv")
# fwrite(coverage, "Poisson results/std_errors_poisson.csv")
# fwrite(coverage, "Poisson results/std_errors_S_poisson.csv")
# fwrite(coverage, "Poisson results/score_fct_true_value_poisson.csv")
# saveRDS(second_order_list, "Poisson results/Second_order_term_Poisson.rds")
# saveRDS(sensitivity_list, "Poisson results/Sensitivity_Poisson.rds")
#
# data.table(Mean = colMeans(betahat),
#            Lower = colMeans(betahat) - 1.96*apply(betahat, 2, sd)/sqrt(nsim),
#            Upper = colMeans(betahat) + 1.96*apply(betahat, 2, sd)/sqrt(nsim)) -> dt
#
# fwrite(dt, "Poisson results/Sim_poisson_parameter_estimate_MC.csv")

# Multi-Strauss ----------------------------------------------------------------
set.seed(123)
ns=10000
nr=1e7
nv=1e6
beta <- c(375, 450, 300)*1.5
# gmma <- matrix(c(0.7,0.80,0.85,0.80,0.75,0.9,0.85,0.9,0.8),3,3)
btw_int <- c(0.8, 0.85, 0.9)
wtn_int <- c(0.5, 0.6, 0.7)
gmma <- matrix(0, 3, 3)
gmma[lower.tri(gmma)] <- btw_int
gmma <- gmma + t(gmma)
diag(gmma) <- wtn_int
# gmma <- matrix(c(0.6,0.80,0.9,0.8,0.5,0.85,0.9,0.85,0.7),3,3)
r    <- matrix(0.04,3,3)
diag(r) <- 0.02
mod_small <- list(cif = "straussm",
                  par = list(beta = beta, gamma = gmma, radii = r),
                  w = c(0,1,0,1))
mod_large <- list(cif = "straussm",
                  par = list(beta = beta, gamma = gmma, radii = r),
                  w = c(0,2,0,2))
X_small <- rmh(model = mod_small, start = list(n.start = ns),
               control = list(ptypes = beta/sum(beta),
                              nrep = nr, nverb = nv, track = T))
X_large <- rmh(model = mod_large, start = list(n.start = ns),
               control = list(ptypes = beta/sum(beta),
                              nrep = nr, nverb = nv, track = T))
X_small$marks <- factor(X_small$marks, levels = 1:3)
X_large$marks <- factor(X_large$marks, levels = 1:3)
target_param <- c(log(beta[-3]/beta[3]), log(gmma[lower.tri(gmma, diag = T)])/2)
# target_param <- log(c(beta[1]/beta[2], 0.6, 0.8, 0.5))/c(1,2,2,2)

nsim=1800
# betahat=matrix(NA,nsim,8)
# coverage=matrix(NA,nsim,8)
# std_errors <- matrix(NA,nsim,8)
# std_errors_S <- matrix(NA,nsim,8)
# score_fct <- matrix(NA,nsim,8)
# second_order_list <- list()
# sensitivity_list <- list()
# score_mat <- matrix(NA, nsim, 8)
Results_sim <- data.table(Parameter = NA, Estimate = NA, Std.error = NA,
                          within_CI = NA, window_size = NA)[-1]
for (i in 1:nsim){
  set.seed(i)
  print(i)
  X_small <- rmh(model = mod_small, start = list(x.start = X_small),
                 control = list(ptypes = beta/sum(beta),
                                nrep = 1e+5))
  X_large <- rmh(model = mod_large, start = list(x.start = X_large),
                 control = list(ptypes = beta/sum(beta),
                                nrep = 1e+5))

  suppressWarnings(
    temp_small <- SemiMarkov_fixed_R(X = X_small,
                                     covariate = NULL,
                                     R_within = 0.02,
                                     R_between = 0.04,
                                     standardize = T,
                                     sat = Inf)
  )

  suppressWarnings(
    temp_large <- SemiMarkov_fixed_R(X = X_large,
                                     covariate = NULL,
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

  # coverage[i,] <- target_param > temp$CI$Lower_CI & target_param < temp$CI$Upper_CI
  # std_errors[i,] <- temp$std_err
  # std_errors_S[i,] <- diag(solve(temp$Sensitivity))
  # betahat[i,]=temp$betahat

  # w <- temp$w[,1:12]
  # log_lambda <- as.matrix(temp$w[,5:12]) %*% target_param
  # w[,lambda:=exp(log_lambda)]
  # w_Lambda <- w[,.(Lambda = sum(lambda)),list(xcoord, ycoord, type_obs)]
  # w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  # w[,prob:=lambda/Lambda]
  #
  # h <- w[,lapply(.SD, function(x) x - sum(x*prob)),
  #        by = list(xcoord, ycoord, type_obs),
  #        .SDcols = 5:12]
  # h <- cbind(w[,1:4], h[,-(1:3)])
  # score_fct[i,] <- colSums(h[type_obs == j, -(1:4)]) # Score function in the true value
  #
  # Check <- SemiMarkov:::Covariance_true_val(X = X,
  #                                           covariate = NULL,
  #                                           R_within = 0.02,
  #                                           R_between = 0.04,
  #                                           sat = Inf,
  #                                           true.param = target_param)
  #
  # second_order_list[[i]] <- Check$Sigma_term # Second order term in the true value
  # sensitivity_list[[i]] <- Check$S # Sensitivity in the true value
}
Sys.time()

res <- Results_sim[,.(Estimate = mean(Estimate), Std.error = mean(Std.error),
                      MC_std.error = sd(Estimate), Coverage = mean(within_CI)),
                   by = list(Parameter, window_size)]

# abs(res$Estimate - rep(target_param,2))/(res$MC_std.error/sqrt(nsim))

data.table(Parameter = res$Parameter,
           Target_param = target_param,
           Estimate = res$Estimate,
           Lower = res$Estimate - 1.96*res$MC_std.error/sqrt(nsim),
           Upper = res$Estimate + 1.96*res$MC_std.error/sqrt(nsim)) -> dt

fwrite(dt, "Strauss results/Sim_Strauss_parameter_estimate_MC.csv")
fwrite(res, "Strauss results/Strauss_res_summary.csv")
fwrite(Results_sim, "Strauss results/Strauss_full_results.csv")

# data.table(Parameter = names(temp$betahat),
#            Estimate = colMeans(betahat),
#            "True value" = target_param,
#            Coverage = colMeans(coverage),
#            "Estimated std. error" = colMeans(std_errors),
#            "Monte-Carlo std. error" = apply(betahat, 2, sd) ) -> res

# score_fct_var <- array(NA, dim = c(10, 10, nsim))
# for(i in 1:nsim){
#   score_fct_var[,,i] <- score_fct[i,] %*% t(score_fct[i,])
# }
# apply(score_fct_var, 1:2, mean)
# apply(score_fct_var, 1:2, sd)/sqrt(1000)
#
# Sens_array <- array(NA, dim = c(4, 4, nsim))
# for(i in 1:nsim){
#   Sens_array[,,i] <- sensitivity_list[[i]] + second_order_list[[i]]
# }
#
# cov_mat_array <- array(NA, dim = c(4, 4, nsim))
# for(i in 1:nsim){
#   Sigma <- sensitivity_list[[i]] + second_order_list[[i]]
#   S <- sensitivity_list[[i]]
#   S_inv <- solve(S)
#   cov_mat_array[,,i] <- S_inv %*% Sigma %*% S_inv
# }
# apply(cov_mat_array, 1:2, mean)
# apply(cov_mat_array, 1:2, sd)
# cov(betahat)
# cov(score_fct)

# fwrite(res, "Strauss results/Strauss_res_summary.csv")
#
# fwrite(betahat, "Strauss results/betahat_Strauss.csv")
# fwrite(coverage, "Strauss results/coverage_Strauss.csv")
# fwrite(coverage, "Strauss results/std_errors_Strauss.csv")
# fwrite(coverage, "Strauss results/std_errors_S_Strauss.csv")
# fwrite(coverage, "Strauss results/score_fct_true_value_Strauss.csv")
# saveRDS(second_order_list, "Strauss results/Second_order_term_Strauss.rds")
# saveRDS(sensitivity_list, "Strauss results/Sensitivity_Strauss.rds")
#
# data.table(Mean = colMeans(betahat),
#            Lower = colMeans(betahat) - 1.96*apply(betahat, 2, sd)/sqrt(nsim),
#            Upper = colMeans(betahat) + 1.96*apply(betahat, 2, sd)/sqrt(nsim)) -> dt
#
# fwrite(dt, "Strauss results/Sim_Strauss_parameter_estimate_MC.csv")


# Geyer saturation -------------------------------------------------------------
set.seed(123)
ns=1000
nr=1e+6
nv=1e+5

ll1=(bei.extra$elev - mean(bei.extra$elev))
ll1=ll1/max(ll1)

exp_points <- 0.02^2*sum(exp(ll1))
ll1 <- ll1 - log(exp_points)
ll1 <- ll1 + log(2000)

trend1 <- function(x, y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(exp(ll1[pts]))
}

trend2 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(exp(ll1[pts] - 1 - bei.extra$grad[pts]))
}

trend3 <- function(x,y){
  pts <- ppp(x = x, y = y, owin(c(0,2), c(0,2)))
  return(exp(ll1[pts] + 2 + bei.extra$grad[pts]))
}

# Saturation process 1
mod1_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.2, r = 0.05, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend1)
X1_small <- rmh(model = mod1_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod1_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.2, r = 0.05, sat = 5),
                   w = c(0,2,0,2),
                   trend = trend1)
X1_large <- rmh(model = mod1_large, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))
X1_small$marks <- rep(1, X1_small$n)
X1_large$marks <- rep(1, X1_large$n)

# Saturation process 2
mod2_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.4, r = 0.05, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend2)
X2_small <- rmh(model = mod2_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod2_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 1.4, r = 0.05, sat = 5),
                   w = c(0,2,0,2),
                   trend = trend2)
X2_large <- rmh(model = mod2_large, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))
X2_small$marks <- rep(2, X2_small$n)
X2_large$marks <- rep(2, X2_large$n)

# Saturation process 3
mod3_small <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 0.8, r = 0.05, sat = 5),
                   w = c(0,1,0,1),
                   trend = trend3)
X3_small <- rmh(model = mod3_small, start = list(n.start = ns),
                control = list(nrep = nr, nverb = nv, track = T))

mod3_large <- list(cif = "geyer",
                   par = list(beta = 1, gamma = 0.8, r = 0.05, sat = 5),
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

# X_small$marks <- factor(X_small$marks, levels = 1:3)
# X_large$marks <- factor(X_large$marks, levels = 1:3)
target_param <- c( c(0, 0, 1, 2), log(c(1.2, 1, 1, 1.4, 1, 0.8)) )

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
    temp_small <- SemiMarkov_fixed_R(X = X_small,
                                     covariate = bei.extra$grad,
                                     R_within = 0.05,
                                     R_between = 0.02,
                                     standardize = T,
                                     sat = 5)
  )

  suppressWarnings(
    temp_large <- SemiMarkov_fixed_R(X = X_large,
                                     covariate = bei.extra$grad,
                                     R_within = 0.05,
                                     R_between = 0.02,
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

# abs(res$Estimate - rep(target_param,2))/(res$MC_std.error/sqrt(nsim))

data.table(Parameter = res$Parameter,
           Target_param = target_param,
           Estimate = res$Estimate,
           Lower = res$Estimate - 1.96*res$MC_std.error/sqrt(nsim),
           Upper = res$Estimate + 1.96*res$MC_std.error/sqrt(nsim)) -> dt

fwrite(dt, "Geyer results/Sim_Geyer_parameter_estimate_MC.csv")
fwrite(res, "Geyer results/Geyer_res_summary.csv")
fwrite(Results_sim, "Geyer results/Geyer_full_results.csv")


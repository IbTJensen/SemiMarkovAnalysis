library(spatstat)
library(data.table)
# library(VGAM)
devtools::document("~/Git/SemiMarkov")
devtools::load_all("~/Git/SemiMarkov")
set.seed(1234)
#Test on tri-variate inhomogeneous Poisson
#use bei.extra covariates
ll0=(bei.extra$elev - 144.2534)
# ll0=ll0/max(ll0)-log(1000)
ll0=ll0/max(ll0)-log(150)
lambda1=exp(ll0)
lambda2=exp(ll0+bei.extra$grad)
lambda3=exp(ll0-bei.extra$grad)
covariate=bei.extra$grad

X1=rpoispp(lambda1)
X2=rpoispp(lambda2)
X3=rpoispp(lambda3)
X=superimpose(X1,X2,X3)
marks=rep(1,X$n)
marks[(X1$n+1):(X1$n+X2$n)]=2
marks[(X1$n+X2$n+1):X$n]=3
X$marks=factor(marks)

# R_within = 10
# R_within_grid = 50
# R_between = 5
# R_between_grid = 20
# sat = Inf
# sat_grid = Inf
# symmetric = T
# edgecorrection = 50
# standardize = T
# n.cores <- 1
# Poisson <- F

# Parameter.value <- c(0, 0, 1, 2, 0, 0, 0, 0, 0, 0)

# S <- SemiMarkov(X = X, covariate = covariate, edgecorrection = 50,
#                 # R_within_grid = c(5, 10, 20, 50, 100),
#                 R_within_grid = 50,
#                 # R_between_grid = c(2, 5, 10, 20, 50),
#                 R_between_grid = 20,
#                 sat_grid = Inf)
# S$R_within; S$R_between; S$sat
# S$CI

# a <- Sys.time()
# S <- SemiMarkov_fixed_Rij(X = X,
#                           covariate = covariate,
#                           edgecorrection = 10,
#                           R_within = 5,
#                           R_between = 2.5,
#                           Poisson = T,
#                           symmetric = T,
#                           standardize = T)
# Sys.time() - a


# R_within <- 100
# R_between <- 50
# sat <- 20

target_param <- c(0, 0, 1 ,2, 0, 0, 0, 0, 0, 0)
# target_param <- c(0, 0, 1 ,2)
# names(target_param) <- S$CI$Covariate

nsim=1000
betahat=matrix(0,nsim,10)
coverage=matrix(0,nsim,10)
std_errors <- matrix(0,nsim,10)
std_errors_S <- matrix(0,nsim,10)
score_fct <- matrix(0,nsim,10)
second_order_list <- list()
sensitivity_list <- list()
score_mat <- matrix(NA, nsim, 10)
# betahat2=matrix(0,nsim,4)
for (i in 1:nsim){
  set.seed(i)
  # a <- Sys.time()
  cat(i, "\r")
  X1=rpoispp(lambda1)
  X2=rpoispp(lambda2)
  X3=rpoispp(lambda3)
  X=superimpose(X1,X2,X3)
  marks=rep(1,X$n)
  marks[(X1$n+1):(X1$n+X2$n)]=2
  marks[(X1$n+X2$n+1):X$n]=3
  X$marks=factor(marks)

  a <- Sys.time()
  suppressWarnings(
  temp <- SemiMarkov_fixed_Rij_new(X = X,
                               covariate = covariate,
                               edgecorrection = 10,
                               R_within = 5,
                               R_between = 2,
                               symmetric = T,
                               standardize = T,
                               sat = Inf,
                               Poisson = F,
                               verbose = F,
                               true.param = NULL)
  )
  Sys.time() - a

  coverage[i,] <- target_param > temp$CI$Lower_CI & target_param < temp$CI$Upper_CI
  std_errors[i,] <- temp$std_err
  std_errors_S[i,] <- diag(solve(temp$Sensitivity))

  w <- temp$w[,1:14]
  log_lambda <- as.matrix(temp$w[,5:14]) %*% target_param
  w[,lambda:=exp(log_lambda)]
  w_Lambda <- w[,.(Lambda = sum(lambda)),list(xcoord, ycoord, type_obs)]
  w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  w[,prob:=lambda/Lambda]

  h <- w[,lapply(.SD, function(x) x - sum(x*prob)),
         by = list(xcoord, ycoord, type_obs),
         .SDcols = 5:14]
  h <- cbind(w[,1:4], h[,-(1:3)])
  score_fct[i,] <- colSums(h[type_obs == j, -(1:4)]) # Score function in the true value

  Check <- Check_second_order_term(X = X,
                                   covariate = covariate,
                                   edgecorrection = 10,
                                   R_within = 5,
                                   R_between = 2,
                                   symmetric = T,
                                   standardize = F,
                                   sat = Inf,
                                   Poisson = F,
                                   verbose = F,
                                   true.param = target_param)

  # suppressWarnings(
  # temp2 <- SemiMarkov_fixed_Rij(X = X,
  #                               covariate = covariate,
  #                               edgecorrection = 50,
  #                               R_within = 50,
  #                               R_between = 20,
  #                               symmetric = T,
  #                               standardize = F,
  #                               sat = Inf)
  # )
  # temp2 <- SemiMarkov(X = X, covariate = NULL, edgecorrection = 0.1,
  #                     # R_within_grid = c(5, 10, 20, 50, 100),
  #                     R_within_grid = 0.05,
  #                     # R_between_grid = c(2, 5, 10, 20, 50),
  #                     R_between_grid = 0.08,
  #                     symmetric = T,
  #                     sat_grid = Inf)

  # conv[i]=temp$convergence
  betahat[i,]=temp$betahat
  # second_order_list[[i]] <- temp$Second_order_term
  second_order_list[[i]] <- Check$Sigma_term # Second order term in the true value
  sensitivity_list[[i]] <- Check$S # Sensitivity in the true value
  # score_mat[i,] <- temp$score
  # betahat2[i,]=temp2$betahat
  # betahat2[i,]=temp2
  # Sys.time() - a
}
beepr::beep()
Sys.time()

data.table(Parameter = names(temp$betahat),
           Estimate = colMeans(betahat),
           "True value" = target_param,
           Coverage = colMeans(coverage),
           "Estimated std. error" = colMeans(std_errors),
           "Monte-Carlo std. error" = apply(betahat, 2, sd) ) -> res

abs(res$Estimate - target_param)/(res$`Estimated std. error`/sqrt(nsim))

score_fct_var <- array(NA, dim = c(10, 10, nsim))
for(i in 1:nsim){
  score_fct_var[,,i] <- score_fct[i,] %*% t(score_fct[i,])
}
apply(score_fct_var, 1:2, mean)
apply(score_fct_var, 1:2, sd)/sqrt(1000)

Sens_array <- array(NA, dim = c(10, 10, nsim))
for(i in 1:nsim){
  Sens_array[,,i] <- sensitivity_list[[i]]
}

cov_mat_array <- array(NA, dim = c(10, 10, nsim))
for(i in 1:nsim){
  Sigma <- sensitivity_list[[i]] + second_order_list[[i]]
  S <- sensitivity_list[[i]]
  S_inv <- solve(S)
  cov_mat_array[,,i] <- S_inv %*% Sigma %*% S_inv
}
apply(cov_mat_array, 1:2, mean)
apply(cov_mat_array, 1:2, sd)
cov(betahat)
# cov(score_fct)

A <- array(NA, dim = c(10, 10, nsim))
for(i in 1:nsim){
  A[,,i] <- second_order_list[[i]]
}

A_mean <- apply(A, 1:2, mean)
A_se <- apply(A, 1:2, sd)/sqrt(nsim)
abs(A_mean) > 1.96*A_se
A_mean/A_se

fwrite(res, "Poisson results/Poisson_res_summary.csv")

fwrite(betahat, "Poisson results/betahat_poisson.csv")
fwrite(coverage, "Poisson results/coverage_poisson.csv")
fwrite(coverage, "Poisson results/std_errors_poisson.csv")
fwrite(coverage, "Poisson results/std_errors_S_poisson.csv")
fwrite(coverage, "Poisson results/score_fct_true_value_poisson.csv")
saveRDS(second_order_list, "Poisson results/Second_order_term_Poisson.rds")
saveRDS(sensitivity_list, "Poisson results/Sensitivity_Poisson.rds")

data.table(Mean = colMeans(betahat),
           Lower = colMeans(betahat) - 1.96*apply(betahat, 2, sd)/sqrt(nsim),
           Upper = colMeans(betahat) + 1.96*apply(betahat, 2, sd)/sqrt(nsim)) -> dt

fwrite(dt, "Poisson results/Sim_poisson_parameter_estimate_MC.csv")

# Multi-Strauss ----------------------------------------------------------------
set.seed(123)
ns=5000
nr=1e7
nv=1e6
beta <- c(125,150)
gmma <- matrix(c(0.6,0.80,0.80,0.5),2,2)
r    <- matrix(c(0.05,0.08,0.08,0.05),2,2)
mod08 <- list(cif = "straussm",
              par = list(beta = beta, gamma = gmma, radii = r),
              w = c(0,5,0,5))
X <- rmh(model = mod08, start = list(n.start = ns),
         control = list(ptypes = c(0.75, 0.25),
                        nrep = nr, nverb = nv, track = T))
X$marks <- factor(X$marks, levels = 1:2)
target_param <- log(c(beta[1]/beta[2], 0.6, 0.8, 0.5))/c(1,2,2,2)

nsim=1000
betahat=matrix(0,nsim,4)
coverage=matrix(0,nsim,4)
std_errors <- matrix(0,nsim,4)
std_errors_S <- matrix(0,nsim,4)
score_fct <- matrix(0,nsim,4)
second_order_list <- list()
sensitivity_list <- list()
score_mat <- matrix(NA, nsim, 4)
# betahat2=matrix(0,nsim,4)
for (i in 1:nsim){
  set.seed(i)
  print(i)
  X <- rmh(model = mod08, start = list(x.start = X),
           control = list(ptypes = c(0.75, 0.25),
                          nrep = 1e+5))

  suppressWarnings(
    temp <- SemiMarkov_fixed_Rij_new(X = X,
                                     covariate = NULL,
                                     R_within = 0.05,
                                     R_between = 0.08,
                                     symmetric = T,
                                     standardize = F,
                                     sat = Inf,
                                     true.param = NULL)
  )

  coverage[i,] <- target_param > temp$CI$Lower_CI & target_param < temp$CI$Upper_CI
  std_errors[i,] <- temp$std_err
  std_errors_S[i,] <- diag(solve(temp$Sensitivity))

  w <- temp$w[,1:8]
  log_lambda <- as.matrix(temp$w[,5:8]) %*% target_param
  w[,lambda:=exp(log_lambda)]
  w_Lambda <- w[,.(Lambda = sum(lambda)),list(xcoord, ycoord, type_obs)]
  w <- merge(w, w_Lambda, by = c("xcoord", "ycoord", "type_obs"))
  w[,prob:=lambda/Lambda]

  h <- w[,lapply(.SD, function(x) x - sum(x*prob)),
         by = list(xcoord, ycoord, type_obs),
         .SDcols = 5:8]
  h <- cbind(w[,1:4], h[,-(1:3)])
  score_fct[i,] <- colSums(h[type_obs == j, -(1:4)]) # Score function in the true value

  Check <- Check_second_order_term(X = X,
                                   covariate = NULL,
                                   R_within = 0.05,
                                   R_between = 0.08,
                                   symmetric = T,
                                   standardize = F,
                                   sat = Inf,
                                   true.param = target_param)

  betahat[i,]=temp$betahat
  second_order_list[[i]] <- Check$Sigma_term # Second order term in the true value
  sensitivity_list[[i]] <- Check$S # Sensitivity in the true value
}
beepr::beep()
Sys.time()

data.table(Parameter = names(temp$betahat),
           Estimate = colMeans(betahat),
           "True value" = target_param,
           Coverage = colMeans(coverage),
           "Estimated std. error" = colMeans(std_errors),
           "Monte-Carlo std. error" = apply(betahat, 2, sd) ) -> res

abs(res$Estimate - target_param)/(res$`Estimated std. error`/sqrt(nsim))

score_fct_var <- array(NA, dim = c(10, 10, nsim))
for(i in 1:nsim){
  score_fct_var[,,i] <- score_fct[i,] %*% t(score_fct[i,])
}
apply(score_fct_var, 1:2, mean)
apply(score_fct_var, 1:2, sd)/sqrt(1000)

Sens_array <- array(NA, dim = c(4, 4, nsim))
for(i in 1:nsim){
  Sens_array[,,i] <- sensitivity_list[[i]] + second_order_list[[i]]
}

cov_mat_array <- array(NA, dim = c(4, 4, nsim))
for(i in 1:nsim){
  Sigma <- sensitivity_list[[i]] + second_order_list[[i]]
  S <- sensitivity_list[[i]]
  S_inv <- solve(S)
  cov_mat_array[,,i] <- S_inv %*% Sigma %*% S_inv
}
apply(cov_mat_array, 1:2, mean)
apply(cov_mat_array, 1:2, sd)
cov(betahat)
# cov(score_fct)

fwrite(res, "Strauss results/Strauss_res_summary.csv")

fwrite(betahat, "Strauss results/betahat_Strauss.csv")
fwrite(coverage, "Strauss results/coverage_Strauss.csv")
fwrite(coverage, "Strauss results/std_errors_Strauss.csv")
fwrite(coverage, "Strauss results/std_errors_S_Strauss.csv")
fwrite(coverage, "Strauss results/score_fct_true_value_Strauss.csv")
saveRDS(second_order_list, "Strauss results/Second_order_term_Strauss.rds")
saveRDS(sensitivity_list, "Strauss results/Sensitivity_Strauss.rds")

data.table(Mean = colMeans(betahat),
           Lower = colMeans(betahat) - 1.96*apply(betahat, 2, sd)/sqrt(nsim),
           Upper = colMeans(betahat) + 1.96*apply(betahat, 2, sd)/sqrt(nsim)) -> dt

fwrite(dt, "Strauss results/Sim_Strauss_parameter_estimate_MC.csv")

# Old stuff ----

#Let us try multitype Strauss
set.seed(123)
ns=5000
nr=1e7
nv=1e6
beta <- c(125,150)
gmma <- matrix(c(0.6,0.80,0.80,0.5),2,2)
r    <- matrix(c(0.05,0.08,0.08,0.05),2,2)
mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),w=c(0,5,0,5))
X <- rmh(model=mod08,start=list(n.start=ns),control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv,track=T))
X$marks <- factor(X$marks, levels = 1:2)

# R_within = 0.05
# R_within_grid = 0.05
# R_between = 0.08
# R_between_grid = 0.08
# sat = Inf
# sat_grid = Inf
# covariate = NULL
# symmetric = T
# edgecorrection = 0.2
# standardize = T
# Poisson <- F

nsim=1000
betahat=matrix(0,nsim,4)
coverage=matrix(0,nsim,4)
coverage_S=matrix(0,nsim,4)
std_errors <- matrix(0,nsim,4)
std_errors_S <- matrix(0,nsim,4)
score_mat <- matrix(NA, nsim, 4)
S <- list()
Second_order <- list()
target_param <- log(c(beta[1]/beta[2], 0.6, 0.8, 0.5))/c(1,2,2,2)
# betahat2=matrix(0,nsim,4)
for(i in 1:nsim){
  set.seed(i)
  print(i)
  # X <- rmh(model=mod08,start=list(n.start=ns),control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
  X <- rmh(model=mod08,start=list(x.start=X),control=list(ptypes=c(0.75,0.25),nrep=1e+5))
  temp <- SemiMarkov_fixed_Rij_new(X = X,
                               covariate = NULL,
                               R_within = 0.05,
                               R_between = 0.08,
                               symmetric = T,
                               standardize = F,
                               sat = Inf,
                               true.param = NULL)
  betahat[i,]=temp$betahat

  std_errors[i,] <- temp$std_err
  # std_errors_S[i,] <- temp$std_err_S
  coverage[i,] <- target_param > temp$CI$Lower_CI & target_param < temp$CI$Upper_CI
  # coverage_S[i,] <- target_param > temp$betahat - 1.96*temp$std_err_S & target_param < temp$betahat + 1.96*temp$std_err_S
  S[[i]] <- temp$Sensitivity
  Second_order[[i]] <- temp$Second_order_term
  # score_mat[i,] <- temp$score
}

data.table(Parameter = c("Intercept", "1-1", "1-2", "2-2"),
           Estimate = colMeans(betahat),
           "True value" = target_param,
           Coverage = colMeans(coverage),
           "Estimated std. error" = colMeans(std_errors),
           "Monte-Carlo std. error" = apply(betahat, 2, sd) ) -> res

fwrite(res, "Strauss_res.csv")

data.table("Target paramter" = target_param,
           Mean = colMeans(betahat),
           Lower = colMeans(betahat) - 1.96*apply(betahat, 2, sd)/sqrt(nsim),
           Upper = colMeans(betahat) + 1.96*apply(betahat, 2, sd)/sqrt(nsim)) -> dt

fwrite(dt, "Sim_Strauss_500.csv")

#with 1e5
#true
log(c((300/250),0.5,0.8,0.6))
(est <- apply(betahat,2,mean)*c(1,2,2,2))
# (est2 <- apply(betahat2,2,mean)*c(1,2,2,2))
data.table(True = log(c((250/300),0.6,0.8,0.5)),
           # Estimated_new = est,
           Estimated_new = c(1,2,2,2)*est,
           Estimated_old = c(1,2,2,2)*est2)

# est <- c(1,2,2,2)*est
sd <- c(1,2,2,2)*apply(betahat,2,sd)/sqrt(nsim)
data.table(True = log(c((300/250),0.5,0.8,0.6)),
           Estimated = est,
           lower_95 = est - 1.96*sd,
           upper_95 = est + 1.96*sd)

[1] -0.1823216 -0.5108256 -0.6931472 -0.2231436
> [1] -0.1683157 -0.5457130 -0.7145818 -0.2188887
> [1] 0.02106053 0.02760113 0.02574501 0.01447793
#looks good. Differences within Monte Carlo error

##Model with unknown common trend
nr=nv=1e6
ns=200

beta  <- c(150,200)
#gmma <- matrix(c(0.6,0.9,0.9,0.7),2,2)
gmma=matrix(c(0.7,0.9,0.9,0.7),2,2)
ri    <- matrix(c(0.04,0.08,0.08,0.04),2,2)
rhc  <- matrix(c(0.01,0.01,0.01,0.01),2,2)
#tr3   <- function(x,y){x <- x; y <- y; exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6) }
tr3   <- function(x,y){x <- x/2; y <- y/2; exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6) }
# log quadratic trend
#tr4   <- function(x,y){x <- x/250; y <- y/250;
#                     exp(-0.6*x+0.5*y)}
#                   # log linear trend
win=c(0,2,0,2)
#win=c(0,1,0,1)
mod10 <- rmhmodel(cif="straushm",par=list(beta=beta,gamma=gmma,
                                          iradii=ri,hradii=rhc),w=win,
                  trend=list(tr3,tr3))#common unknown trend
nsim=1000
betahat=matrix(0,nsim,4)
conv=rep(0,nsim)
interaction=list(StraussHard(0.04,0.01),StraussHard(0.08,0.01),StraussHard(0.08,0.01),StraussHard(0.04,0.01))
for (i in 1:nsim){
  print(i)
  X <- rmh(model=mod10,start=list(n.start=ns),control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
  temp=FirstOrderCCL(X,covariate=NULL,interaction,symmetric=T)
  conv[i]=temp$converg
  betahat[i,]=temp$betahat
}

apply(betahat,2,mean)
#true 100 sim.
log(c(150/200,0.6,0.7,1.1))
apply(betahat,2,mean)
apply(betahat,2,sd)/sqrt(nsim)#deviation for last parameter.
[1] -0.28768207 -0.51082562 -0.35667494  0.09531018
> [1] -0.2987835 -0.5615678 -0.3761232  0.1177104
> [1] 0.019782440 0.025597786 0.022567651 0.009538402

#without between interaction. 100 sim.
log(c(150/200,0.6,0.7,1))
apply(betahat,2,mean)
apply(betahat,2,sd)/sqrt(nsim)#that is fine.

#for 1000 simulations
simout(betahat,log(c(150/200,0.6,0.7,1.1))
       [1] -0.28768207 -0.51082562 -0.35667494  0.09531018
       [1] -0.3039229 -0.5432143 -0.3762612  0.1175274
       [1] 0.006880843 0.007505677 0.006225720 0.003069169
       [1] 1.204236 2.201648 1.605118 3.693284#bias for all parameters

       simout(betahat,log(c(150/200,0.6,0.7,0.9))
              [1] -0.2876821 -0.5108256 -0.3566749 -0.1053605
              [1] -0.30327644 -0.56637792 -0.38682845 -0.09297216
              [1] 0.006706344 0.012301340 0.009346002 0.004088689
              [1] 1.186385 2.304059 1.646099 1.545872#bias for all parameters.

              #try larger window. 2 x 2.
              simout(betahat,log(c(150/200,0.6,0.7,0.9))
                     1] -0.2876821 -0.5108256 -0.3566749 -0.1053605
[1] -0.2904768 -0.5352045 -0.3684985 -0.1011396
[1] 0.003051094 0.005792538 0.004352685 0.001815987
[1] 0.4673344 2.1472760 1.3859108 1.1858774#Bias for three last parameters.
hist(betahat[,2])
summary(betahat[,2])
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
-1.1559 -0.6486 -0.5288 -0.5352 -0.4064 -0.0593
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
-0.85631 -0.46273 -0.36212 -0.36850 -0.27887  0.03103
hist(betahat[,3])
qqnorm(betahat[,3])
qqline(betahat[,3])#looks very normal
qqnorm(betahat[,2])
qqline(betahat[,2])#looks pretty normal
hist(betahat[,2])#looks a bit skew. But qq looks quite fine.

#with weaker repulsion for type 1
simout(betahat,log(c(150/200,0.7,0.7,0.9))
       [1] -0.2876821 -0.3566749 -0.3566749 -0.1053605
       [1] -0.30356512 -0.39139590 -0.37956266 -0.09718099
       [1] 0.006018390 0.009558878 0.007841555 0.003231743
       [1] 1.346472 1.853227 1.489170 1.291325
       [1] 1.346472 1.853227 1.489170 1.291325#bias for all but seems moderate.

       #with 1e5.
       simout(betahat,log(c(150/200,0.7,0.7,0.9))
              [1] -0.2876821 -0.3566749 -0.3566749 -0.1053605
              [1] -0.28730169 -0.44150670 -0.39468573 -0.09976983
              [1] 0.006656535 0.012240551 0.009877366 0.004017643
              [1] 0.05714478 6.93038737 3.84827108 1.39153370
              [1]  TRUE FALSE FALSE  TRUE#Strong bias due to lack of MCMC convergence ?


              #with 1e6
              simout(betahat,log(c(150/200,0.7,0.7,0.9)))
              [1] -0.2876821 -0.3566749 -0.3566749 -0.1053605
              [1] -0.30304122 -0.39707476 -0.38257569 -0.09494279
              [1] 0.006441417 0.011206715 0.009419882 0.003792208
              [1] 2.384436 3.604965 2.749583 2.747138
              [1] FALSE FALSE FALSE FALSE

              #with 1e6 2x2
              simout(betahat,log(c(150/200,0.7,0.7,0.9)))
              [1] -0.2876821 -0.3566749 -0.3566749 -0.1053605
              [1] -0.2948093 -0.3698203 -0.3640818 -0.1019107
              [1] 0.003173496 0.005177526 0.004435933 0.001830099
              [1] 2.245858 2.538928 1.669742 1.885067#still bias, but not much.
              [1] FALSE FALSE  TRUE  TRUE



              #Look at tr3
              x=y=seq(0,1,len=100)

              xy=expand.grid(x,y)

              z=tr3(xy[,1],xy[,2])
              image(t(matrix(z,ncol=100,byrow=T)))
              summary(z)
              Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
              0.1353  0.7451  1.0687  0.9863  1.2505  1.4549


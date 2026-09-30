library(spatstat)
library(data.table)
library(geoR)

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

saveRDS(X_im, "1_simulation_study/simulated_data/X_im.rds")
saveRDS(Z_im, "1_simulation_study/simulated_data/Z_im.rds")
saveRDS(phi0_init, "1_simulation_study/simulated_data/phi0_init.rds")

library(data.table)
library(spatstat)
library(ggplot2)
library(ggforce)
library(VGAM)

# X1 <- ppp(x = runif(3), y = runif(3), marks = rep(1, 3))
# X2 <- ppp(x = runif(5), y = runif(5), marks = rep(2, 5))
# X3 <- ppp(x = runif(4), y = runif(4), marks = rep(3, 4))
#
# X <- superimpose(X1, X2, X3)
# X$marks <- factor(X$marks)
# plot(X)

X1 <- ppp(x = c(0.79, 0.4175, 0.426), y = c(0.50495, 0.83862, 0.00816), marks = rep(1, 3))
X2 <- ppp(x = c(0.1565, 0.6737, 0.7016, 0.954, 0.5002), y = c(0.2549, 0.9651, 0.01358, 0.73809, 0.94365), marks = rep(2,5))
X3 <- ppp(x = c(0.51, 0.8234, 0.8361, 0.0339), y = c(0.80402, 0.85902, 0.46108, 0.61231), marks = rep(3, 4))

X <- superimpose(X1, X2, X3)
X$marks <- factor(X$marks)

dt <- data.table(x = X$x, y = X$y, mark = X$marks)
plot(X)

ggplot(data = dt, mapping = aes(x = x, y = y, col = mark))+
  geom_point()+
  geom_circle(data = dt[mark == 3], aes(x0 = x, y0 = y, r = 0.5))

# Strauss
Neighboors <- data.table(dt,
                         s_11 = c(0, 0, 0, 1, 2, 2, 1, 1, 2, 2, 1, 1),
                         s_12 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         s_13 = c(1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0),
                         s_21 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         s_22 = c(3, 2, 2, 0, 2, 0, 2, 2, 3, 3, 2, 1),
                         s_23 = c(1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0),
                         s_31 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         s_32 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         s_33 = c(3, 3, 0, 1, 2, 1, 3, 2, 2, 2, 2, 0),
                         DeltaS_11 = c(0, 0, 0, 1, 2, 2, 1, 1, 2, 2, 1, 1),
                         DeltaS_12 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         DeltaS_13 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         DeltaS_21 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         DeltaS_22 = c(3, 2, 2, 0, 2, 0, 2, 2, 3, 3, 2, 1),
                         DeltaS_23 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         DeltaS_31 = c(1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0),
                         DeltaS_32 = c(1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 0, 0),
                         DeltaS_33 = c(3, 3, 0, 1, 2, 1, 3, 2, 2, 2, 2, 0)
                         )


# Saturation
sat_i <- 1.8*c(3,5,4)/4
Neighboor2 <- data.table(dt,
                         s_11 = c(0, 0, 0, 1, 1.35, 1.35, 1, 1, 1.35, 1.35, 1, 1),
                         s_12 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         s_13 = c(1, 1, 0, 0, 1.80, 0, 1, 1, 0, 0, 0, 0),
                         s_21 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         s_22 = c(2.25, 2, 2, 0, 2, 0, 2, 2, 2.25, 2.25, 2, 1),
                         s_23 = c(1, 1, 0, 0, 1.80, 0, 1, 1, 0, 0, 0, 0),
                         s_31 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         s_32 = c(0, 1, 0, 0, 1, 0, 0, 1, 2, 2, 0, 0),
                         s_33 = c(1.80, 1.80, 0, 1, 1.80, 1, 1.80, 1.80, 1.80, 1.80, 1.80, 0),
                         DeltaS_11 = c(0, 0, 0, 1, 2, 2, 1, 1, 2, 2, 1, 1),
                         DeltaS_12 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                         DeltaS_13 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         DeltaS_21 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0),
                         DeltaS_22 = c(1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1),
                         DeltaS_23 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         DeltaS_31 = c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0),
                         DeltaS_32 = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         DeltaS_33 = c(0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
                         )

covariate <- NULL
symmetric <- T
edgecorrection <- NULL
R_within <- 0.5
R_between <- 0.25
standardize <- T
sat <- 1.8

# prelim_dt[match(Neighboor2$x,prelim_dt$xcoord)] == Neighboor2

# Chekcs ----
plot(X)

crosspairs(X1, X1, rmax = 0.5)
crosspairs(X1, X2, rmax = 0.25)
crosspairs(X1, X3, rmax = 0.25)
crosspairs(X2, X2, rmax = 0.5)
crosspairs(X2, X3, rmax = 0.25)
crosspairs(X3, X3, rmax = 0.5)

Neigbour_test <- data.table(xcoord = X$x, ycoord = X$y, type_obs = X$marks,
                            Neigh_type1 = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0),
                            Neigh_type2 = c(0, 1, 0, 0, 2, 0, 2, 2, 2, 2, 0, 0),
                            Neigh_type3 = c(1, 1, 0, 0, 2, 0, 1, 1, 2, 2, 2, 0) )

SemiMarkov_fixed_Rij(X, covariate, symmetric = T, edgecorrection = NULL,
                     R_within, R_between, standardize = T, sat = Inf)







library(spatstat)
library(data.table)
library(SemiMarkov)
library(ggplot2)
library(ggpubr)
library(sf)

# Loading data -----------------------------------------------------------------
load("Data/covariatesfct_alr.rdata")
load("Data/pointsAndPoly.RData")
load("Data/regions.rdata")
do_grid_search <- F

# Plotting point pattern -------------------------------------------------------
ppp.dt <- data.table(x = ppp.type$x, y = ppp.type$y,
                     Type = ifelse(ppp.type$marks == "cooperative",
                                   "Cooperative banks", "Lucrative banks")
                     )
http_comp <- c("https://raw.githubusercontent.com", "gregoiredavid",
               "france-geojson", "master",
               "regions-version-simplifiee.geojson")
region_url <- paste(http_comp, collapse = "/")
France <- st_read(region_url)
france_outer <- st_union(France)

ggplot(data = france_outer) +
  geom_sf(fill = "white", color = "black") +
  geom_point(data = ppp.dt, mapping = aes(x = x, y = y, color = Type),
             size = 0.2)+
  facet_wrap(~Type)+
  guides(color = "none")+
  scale_color_manual(values = c("#bc272d", "#0000a2"))+
  theme_minimal() +
  labs(x = NULL, y = NULL)+
  scale_x_continuous(breaks = seq(-4,12,4)) +
  NULL -> ppp.plot

ggsave("Figures/Banks_plot_ppp.pdf", ppp.plot,
       width = 165, height = 120, units = "mm")

# Loading an managing covariates -----------------------------------------------
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

colnames(Covar_dt) <- gsub("prop16.24", "prop.age", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("evolution", "pop.growth", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("decile", "income.decile", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("prop19", "prop.firms", colnames(Covar_dt))

# Grid search for sat and R ----------------------------------------------------
if(do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = 0.02,
                   R_within = seq(0.001, 0.01, 0.001),
                   R_between = seq(0.001, 0.01, 0.001),
                   sat = 2:10)
}

# Fit model and interpret results ----------------------------------------------
if(!do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = 0.02,
                   R_within = 0.006,
                   R_between = 0.004,
                   sat = 3,
                   standardize = T,
                   Poisson = F)
}

SS$R_within # 0.006
SS$R_between # 0.004
SS$maximum_log_likelihood
SS$sat # 3

SS$CI

SS_CI <- data.table(SS$CI)
SS_CI[,Covariate:=gsub(":1", "", Covariate)]
SS_CI[10,Covariate:=gsub("_", "-", Covariate)]
SS_CI[11:13,Covariate:=gsub("1", "Luc", Covariate)]
SS_CI[11:13,Covariate:=gsub("2", "Coop", Covariate)]

SS_plot <- data.table(SS_CI)
SS_plot[1,Covariate:="        Intercept"]
SS_plot[,Covariate:=factor(Covariate, levels = rev(Covariate))]
ggplot(data = SS_plot[-1],
       mapping = aes(x = Estimate,
                     y = Covariate))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI))+
  theme_bw()+
  labs(y = "Covariate")+
  NULL -> CI_plot

# ggplot(data = SS_plot[1],
#        mapping = aes(x = Estimate,
#                      y = Covariate))+
#   geom_point()+
#   geom_vline(xintercept = 0, linetype = "dashed")+
#   geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI))+
#   theme_bw()+
#   labs(y = "", x = NULL)+
#   NULL -> CI_plot2
#
# CI_plots <- ggarrange(CI_plot2, CI_plot, ncol = 1, heights = c(0.175, 0.825))
# ggsave(filename = "Figures/CI_plot.pdf",
#        width = 170, height = 120, units = "mm")
ggsave(
  filename = "Figures/CI_plot.pdf",
  plot = CI_plot, width = 170, height = 66, units = "mm"
)

# Cooperative probabilities - supplementary ------------------------------------
logit <- function(x) log(x/(1-x))
pred <- SS$pred[type_pred == 1 & !is.na(prob)]
pred[,logit_prob:=logit(prob)]
poly_dt <- data.table(poly)

ewbrks <- seq(-4, 12, 4)
nsbrks <- seq(42, 50, 2)
ewlbls <- sapply(ewbrks, function(x) fcase(x < 0, paste(x, "°E"),
                                           x > 0, paste(x, "°W"),
                                           x == 0, paste(0, "°")))
nslbls <- sapply(nsbrks, function(x) fcase(x < 0, paste(x, "°S"),
                                           x > 0, paste(x, "°N"),
                                           x == 0, paste(0, "°")))
ggplot(data = france_outer) +
  geom_sf(fill = "white", color = "black") +
  geom_point(data = pred, mapping = aes(x = xcoord, y = ycoord,
                                        group = 1, colour = logit_prob),
             size = 0.5)+
  scale_colour_gradientn(colours = c("yellow", "orange", "red",
                                     "purple", "blue"))+
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(x = NULL, y = NULL, colour = "Logit-probability")+
  scale_x_continuous(breaks = seq(-4,12,4)) +
  NULL -> Points_in_France

ggsave("Figures/logit_prob_points.pdf", Points_in_France,
       width = 165, height = 140, units = "mm")

# Phi0 kernel estiamte ---------------------------------------------------------
bandwidth <- 0.25

x_pts <- seq(-5.1, 9.5, 0.05)
y_pts <- seq(41.3, 51, 0.05)
points_grid <- expand.grid(x = x_pts, y = y_pts)
in_France <- inside.owin(points_grid$x, points_grid$y, ppp.type$window)
points_grid <- points_grid[in_France,]

points_grid_ppp <- ppp(points_grid$x, points_grid$y, ppp.type$window)
ccc <- crosspairs(points_grid_ppp, ppp.type, rmax = bandwidth)

# Epanechnikov kernel
K <- function(u){
  2/pi*(1-(u/bandwidth)^2)*(u <= bandwidth)
}
# Uniform kernel
# K <- function(u){
#   return(1/(pi*bandwidth^2))
# }

phi0 <- rep(0, nrow(points_grid))
N <- rep(0, nrow(points_grid))
Admin_unit <- rep(0, nrow(points_grid))
gm_save <- rep(0, nrow(points_grid))
Log_dens <- rep(0, nrow(points_grid))
for(i in 1:length(phi0)){
  cat(i, "\r")
  w <- which(ccc$i == i)
  p1 <- ppp.type[ccc$j[w]]
  idx <- (1:ncol(SS$w))[-c(1:4, ncol(SS$w))]
  A <- SS$w[xcoord %in% p1$x & ycoord %in% p1$y & type_obs == j]
  gamma_v <- as.numeric(as.matrix(A[,idx, with = F]) %*% SS$betahat)
  if(nrow(A)>0){
    vec_diffs <- t(t(as.matrix(A[,1:2])) - c(points_grid_ppp$x[i],
                                             points_grid_ppp$y[i]))
    distances <- apply(vec_diffs, 1, function(x) sqrt(sum(x^2)))
    phi0[i] <- 1/2*sum(K(distances)/exp(gamma_v), na.rm = T)
    # gm_save[i] <- gamma_v
  }
  N[i] <- nrow(A)
  Admin_unit[i] <- covariatesfct[[2]](points_grid_ppp$x[i],
                                      points_grid_ppp$y[i])
  Log_dens[i] <- covariatesfct[[3]](points_grid_ppp$x[i],
                                    points_grid_ppp$y[i])
}

phi0_dt <- data.table(points_grid, phi0 = phi0, N = N, Admin_unit, Log_dens)

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

splot(log(phi0_im), main = NULL)

ewbrks <- seq(-4, 12, 4)
nsbrks <- seq(42, 50, 2)
ewlbls <- sapply(ewbrks, function(x) fcase(x < 0, paste(x, "°E"),
                                           x > 0, paste(x, "°W"),
                                           x == 0, paste(0, "°")))
nslbls <- sapply(nsbrks, function(x) fcase(x < 0, paste(x, "°S"),
                                           x > 0, paste(x, "°N"),
                                           x == 0, paste(0, "°")))
by_unit <- phi0_dt[,.(phi0 = mean(phi0, na.rm = T),
                      Log_dens = mean(Log_dens, na.rm = T)),
                   Admin_unit]
colnames(by_unit)[1] <- "city"
poly_dt2 <- merge(poly_dt, by_unit, by = "city")
poly_dt2[,":="(phi0 = log(phi0), Log_dens = log(Log_dens))]

ggplot(poly_dt2, aes(x, y, group = code, fill = phi0))+
  geom_polygon(colour = "black", linewidth = 0.1)+
  labs(x = NULL, y = NULL,
       fill = expression(paste("log ", hat(phi)[0])))+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours = rev(c("yellow", "orange", "red",
                                   "purple", "blue", "blue", "blue")))+
  coord_quickmap()+
  scale_x_continuous(breaks = ewbrks, labels = ewlbls) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls) +
  NULL -> phi0_plot

ggplot(poly_dt2, aes(x, y, group = code, fill = Log_dens))+
  geom_polygon(colour = "black", linewidth = 0.1)+
  labs(x = NULL, y = NULL, fill = "log-density")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_fill_gradientn(colours = rev(c("yellow", "orange", "red",
                                   "purple", "blue")))+
  coord_quickmap()+
  scale_x_continuous(breaks = ewbrks, labels = ewlbls) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls) +
  NULL -> Dens_plot

phi0_dens_plot <- ggarrange(phi0_plot, Dens_plot)
ggsave("Figures/phi0_dens_plot.pdf", phi0_dens_plot,
       width = 165, height = 140, units = "mm")

cor(log(by_unit$phi0[by_unit$phi0 != 0]),
    log(by_unit$Log_dens[by_unit$phi0 != 0]),
    use = "pairwise.complete.obs")

plot(log(by_unit$phi0[by_unit$phi0 != 0]),
     log(by_unit$Log_dens[by_unit$phi0 != 0]))

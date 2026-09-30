library(spatstat)
library(data.table)
library(SeMaCoPE)
library(ggplot2)
library(ggpubr)
library(sf)

# Loading data -----------------------------------------------------------------
load("2_analysis_bank_data/data/covariatesfct_alr.rdata")
load("2_analysis_bank_data/data/pointsAndPoly.RData")
load("2_analysis_bank_data/data/regions.rdata")
Covar_dt <- fread("2_analysis_bank_data/data/Covariates_processed.csv")
ppp.type$marks <- factor(ppp.type$marks, levels = c("lucrative", "cooperative"))

# Plotting point pattern -------------------------------------------------------
ppp.dt <- data.table(
  x = ppp.type$x,
  y = ppp.type$y,
  Type = ifelse(
    ppp.type$marks == "cooperative",
    "Cooperative banks",
    "Lucrative banks"
  )
)
http_comp <- c(
  "https://raw.githubusercontent.com",
  "gregoiredavid",
  "france-geojson",
  "master",
  "regions-version-simplifiee.geojson"
)
region_url <- paste(http_comp, collapse = "/")
France <- st_read(region_url)
france_outer <- st_union(France)

ggplot(data = france_outer) +
  geom_sf(fill = "white", color = "black") +
  geom_point(
    data = ppp.dt,
    mapping = aes(x = x, y = y, color = Type),
    size = 0.2
  ) +
  facet_wrap(~Type) +
  guides(color = "none") +
  scale_color_manual(values = c("#bc272d", "#0000a2")) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = seq(-4, 12, 4)) +
  NULL -> ppp.plot

ggsave(
  "2_analysis_bank_data/figures/Banks_plot_ppp.pdf",
  ppp.plot,
  width = 165,
  height = 120,
  units = "mm"
)

# Fit model and interpret results ----------------------------------------------
SS <- SeMaCoPE:::SemiMarkov(
  X = ppp.type,
  covariate = Covar_dt,
  edgecorrection = 0.064,
  R_within = 0.032,
  R_between = 0.006,
  sat = NULL,
  sat_within = 1,
  sat_between = 2,
  standardize = TRUE,
  Poisson = FALSE
)

SS_CI <- data.table(SS$CI)
SS_CI[,
  Covariate := fcase(
    Covariate == "1-1" , "Coop-Coop" ,
    Covariate == "1-2" , "Coop-Luc"  ,
    Covariate == "2-2" , "Luc-Luc"   ,
    Covariate == "log_density:1" , "log-density"   ,
    default = gsub(":1", "", Covariate)
  )
]

SS_plot <- data.table(SS_CI)
SS_plot <- SS_plot[Covariate != "(Intercept)"]
SS_plot[,Covariate:=factor(Covariate, levels = rev(Covariate))]
ggplot(data = SS_plot, mapping = aes(x = Estimate, y = Covariate)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI)) +
  theme_bw() +
  labs(y = "Covariate") +
  NULL -> CI_plot

ggsave(
  filename = paste0(
    "2_analysis_bank_data/figures/CI_plot_",
    Sys.Date(),
    ".pdf"
  ),
  plot = CI_plot,
  width = 170,
  height = 66,
  units = "mm"
)

# Cooperative probabilities - supplementary ------------------------------------
logit <- function(x) log(x/(1-x))
pred <- SS$pred[type_pred == 1 & !is.na(prob)]
pred[,logit_prob:=logit(prob)]
poly_dt <- data.table(poly)

ewbrks <- seq(-4, 12, 4)
nsbrks <- seq(42, 50, 2)
ewlbls <- sapply(ewbrks, function(x) {
  fcase(
    x < 0  , paste(x, "°E") ,
    x > 0  , paste(x, "°W") ,
    x == 0 , paste(0, "°")
  )
})
nslbls <- sapply(nsbrks, function(x) {
  fcase(
    x < 0  , paste(x, "°S") ,
    x > 0  , paste(x, "°N") ,
    x == 0 , paste(0, "°")
  )
})
ggplot(data = france_outer) +
  geom_sf(fill = "white", color = "black") +
  geom_point(
    data = pred,
    mapping = aes(x = xcoord, y = ycoord, group = 1, colour = logit_prob),
    size = 0.5
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "orange", "red", "purple", "blue")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = NULL, colour = "Logit-probability") +
  scale_x_continuous(breaks = seq(-4, 12, 4)) +
  NULL -> Points_in_France

ggsave(
  paste0("2_analysis_bank_data/figures/logit_prob_points_", Sys.Date(), ".pdf"),
  Points_in_France,
  width = 165,
  height = 140,
  units = "mm"
)

# Hosmer-Lemeshow type plot
type_pred_dt <- SS$pred[type_pred == 1]
type_pred_dt <- type_pred_dt[!is.na(prob)]
type_pred_dt <- type_pred_dt[order(prob)]

x <- seq(0.05, 1, 0.1)
y <- rep(NA, length(x))
for (i in 1:length(x)) {
  btw_dt <- type_pred_dt[between(prob, x[i] - 0.05, x[i] + 0.05)]
  y[i] <- mean(btw_dt$type_obs == 1)
}
dt <- data.table(prob = x, pred_prob = y)

ggplot(data = dt, aes(x = prob, y = pred_prob)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  labs(x = "Cooperative probability-bin midpoint", y = "Observed proportion of cooperative")+
  NULL -> prob_bins

ggsave(
  paste0("2_analysis_bank_data/figures/prob_bins_", Sys.Date(), ".pdf"),
  prob_bins,
  height = 90,
  width = 160,
  units = "mm"
)

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

phi0 <- rep(0, nrow(points_grid))
N <- rep(0, nrow(points_grid))
Admin_unit <- rep(0, nrow(points_grid))
gm_save <- rep(0, nrow(points_grid))
Log_dens <- rep(0, nrow(points_grid))
for (i in 1:length(phi0)) {
  cat(i, "\r")
  w <- which(ccc$i == i)
  p1 <- ppp.type[ccc$j[w]]
  idx <- (1:ncol(SS$w))[-c(1:4, ncol(SS$w))]
  A <- SS$w[xcoord %in% p1$x & ycoord %in% p1$y & type_obs == j]
  gamma_v <- as.numeric(as.matrix(A[, idx, with = F]) %*% SS$betahat)
  if (nrow(A) > 0) {
    vec_diffs <- t(
      t(as.matrix(A[, 1:2])) - c(points_grid_ppp$x[i], points_grid_ppp$y[i])
    )
    distances <- apply(vec_diffs, 1, function(x) sqrt(sum(x^2)))
    phi0[i] <- 1 / 2 * sum(K(distances) / exp(gamma_v), na.rm = T)
  }
  N[i] <- nrow(A)
  Admin_unit[i] <- covariatesfct[[2]](
    points_grid_ppp$x[i],
    points_grid_ppp$y[i]
  )
  Log_dens[i] <- covariatesfct[[3]](points_grid_ppp$x[i], points_grid_ppp$y[i])
}

phi0_dt <- data.table(points_grid, phi0 = phi0, N = N, Admin_unit, Log_dens)

dens_fine <- fread("Data/France_pop_density_finescale.csv")
dens_fine[,admin_unit:=covariatesfct[[2]](X, Y)]

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

plot(log(phi0_im), main = NULL)

ewbrks <- seq(-4, 12, 4)
nsbrks <- seq(42, 50, 2)
ewlbls <- sapply(ewbrks, function(x) {
  fcase(
    x < 0  , paste(x, "°E") ,
    x > 0  , paste(x, "°W") ,
    x == 0 , paste(0, "°")
  )
})
nslbls <- sapply(nsbrks, function(x) {
  fcase(
    x < 0  , paste(x, "°S") ,
    x > 0  , paste(x, "°N") ,
    x == 0 , paste(0, "°")
  )
})
by_unit <- phi0_dt[,
  .(phi0 = mean(phi0, na.rm = T), Log_dens = mean(Log_dens, na.rm = T)),
  Admin_unit
]
colnames(by_unit)[1] <- "city"
poly_dt2 <- merge(poly_dt, by_unit, by = "city")
poly_dt2[,":="(phi0 = log(phi0), Log_dens = log(Log_dens))]

ggplot(poly_dt2, aes(x, y, group = code, fill = phi0)) +
  geom_polygon(colour = "black", linewidth = 0.1) +
  labs(x = NULL, y = NULL, fill = expression(paste("log ", hat(phi)[0]))) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_gradientn(
    colours = rev(c(
      "yellow",
      "orange",
      "red",
      "purple",
      "blue",
      "blue",
      "blue"
    ))
  ) +
  coord_quickmap() +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls) +
  NULL -> phi0_plot

ggplot(poly_dt2, aes(x, y, group = code, fill = Log_dens)) +
  geom_polygon(colour = "black", linewidth = 0.1) +
  labs(x = NULL, y = NULL, fill = "log-density") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_gradientn(
    colours = rev(c("yellow", "orange", "red", "purple", "blue"))
  ) +
  coord_quickmap() +
  scale_x_continuous(breaks = ewbrks, labels = ewlbls) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls) +
  NULL -> Dens_plot

phi0_dens_plot <- ggarrange(phi0_plot, Dens_plot)
ggsave(
  paste0("2_analysis_bank_data/figures/phi0_dens_plot_", Sys.Date(), ".pdf"),
  phi0_dens_plot,
  width = 165,
  height = 140,
  units = "mm"
)

cor(
  log(by_unit$phi0[by_unit$phi0 != 0]),
  log(by_unit$Log_dens[by_unit$phi0 != 0]),
  use = "pairwise.complete.obs"
)

plot(
  log(by_unit$phi0[by_unit$phi0 != 0]),
  log(by_unit$Log_dens[by_unit$phi0 != 0])
)

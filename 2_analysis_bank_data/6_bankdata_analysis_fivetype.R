library(spatstat)
library(data.table)
library(SeMaCoPE)
library(ggplot2)
library(ggpubr)
library(sf)
library(xtable)

# Loading data -----------------------------------------------------------------
load("2_analysis_bank_data/data/covariatesfct_alr.rdata")
load("2_analysis_bank_data/data/pointsAndPoly.RData")
load("2_analysis_bank_data/data/regions.rdata")
Covar_dt <- fread("2_analysis_bank_data/data/Covariates_processed.csv")
ppp.banks$marks <- factor(
  ppp.banks$marks,
  levels = c("sg", "bp", "bnp", "ca", "cm")
)

# Plotting point pattern -------------------------------------------------------
ppp.dt <- data.table(
  x = ppp.banks$x,
  y = ppp.banks$y,
  Type = as.character(ppp.banks$marks)
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
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = seq(-4, 12, 4)) +
  NULL -> ppp.plot

ggsave(
  "2_analysis_bank_data/figures/Banks_plot_ppp_fivetype.pdf",
  ppp.plot,
  width = 165,
  height = 120,
  units = "mm"
)

# Fit model and interpret results ----------------------------------------------
SS <- SemiMarkov(
  X = ppp.banks,
  covariate = Covar_dt,
  edgecorrection = 0.022,
  R_within = 0.004,
  R_between = 0.011,
  sat = 1,
  standardize = TRUE,
  Poisson = FALSE
)

SS_CI <- data.table(SS$CI)
mark.pp <- ppp.banks$marks |> unique() |> sort() |> as.character()

SS_CI[, Covariate := gsub("1", mark.pp[1], Covariate)]
SS_CI[, Covariate := gsub("2", mark.pp[2], Covariate)]
SS_CI[, Covariate := gsub("3", mark.pp[3], Covariate)]
SS_CI[, Covariate := gsub("4", mark.pp[4], Covariate)]
SS_CI[, Covariate := gsub("5", mark.pp[5], Covariate)]

SS_CI_spat <- SS_CI[1:40]
SS_CI_int <- SS_CI[41:55]

xtable(SS_CI_spat) |> print(include.rownames = FALSE)
xtable(SS_CI_int) |> print(include.rownames = FALSE)

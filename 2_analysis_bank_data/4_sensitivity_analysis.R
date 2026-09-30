library(spatstat)
library(data.table)
library(SeMaCoPE)
library(ggplot2)
library(ggpubr)

# Loading data -----------------------------------------------------------------
load("2_analysis_bank_data/data/covariatesfct_alr.rdata")
load("2_analysis_bank_data/data/pointsAndPoly.RData")
load("2_analysis_bank_data/data/regions.rdata")
Covar_dt <- fread("2_analysis_bank_data/data/Covariates_processed.csv")

# Sensivity analysis R_within --------------------------------------------------
r_within <- seq(0.028, 0.036, 0.002)
sens_within <- lapply(1:5, function(i) {
  SS <- SemiMarkov(
    X = ppp.type,
    covariate = Covar_dt,
    edgecorrection = 0.064,
    R_within = r_within[i],
    R_between = 0.006,
    sat = NULL,
    sat_within = 1,
    sat_between = 2,
    standardize = TRUE,
    Poisson = FALSE
  )
})

CIs <- lapply(sens_within, function(x) x$CI)
CIs <- lapply(1:5, function(i) data.table(CIs[[i]], R_within = r_within[i]))
CIs <- rbindlist(CIs)

SS_CI <- data.table(CIs)
SS_CI[,Covariate:=gsub(":1", "", Covariate)]
SS_CI[,
  Covariate := fcase(
    Covariate == "1-1" , "Coop-Coop" ,
    Covariate == "1-2" , "Coop-Luc"  ,
    Covariate == "2-2" , "Luc-Luc"   ,
    default = Covariate
  )
]

SS_plot <- data.table(SS_CI)
SS_plot <- SS_plot[Covariate != "(Intercept)"]
SS_plot[, ":="(
  Covariate = factor(Covariate, levels = Covariate[1:12]),
  R_within = factor(R_within, levels = r_within)
)]

ggplot(data = SS_plot, mapping = aes(x = Estimate, y = R_within)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI)) +
  theme_bw() +
  labs(y = "R_within") +
  facet_grid(Covariate ~ .) +
  theme(
    strip.text.y.right = element_text(angle = 0, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  ) +
  NULL -> CI_plot

ggsave(
  filename = paste0(
    "2_analysis_bank_data/figures/sens_analysis_rwithin_",
    Sys.Date(),
    ".pdf"
  ),
  plot = CI_plot,
  width = 170,
  height = 260,
  units = "mm"
)

# Sensivity analysis R_between -------------------------------------------------
r_between <- seq(0.004, 0.008, 0.001)
sens_between <- lapply(1:5, function(i) {
  SS <- SemiMarkov(
    X = ppp.type,
    covariate = Covar_dt,
    edgecorrection = 0.064,
    R_within = 0.032,
    R_between = r_between[i],
    sat = NULL,
    sat_within = 1,
    sat_between = 2,
    standardize = TRUE,
    Poisson = FALSE
  )
})

CIs <- lapply(sens_between, function(x) x$CI)
CIs <- lapply(1:5, function(i) data.table(CIs[[i]], R_between = r_between[i]))
CIs <- rbindlist(CIs)

SS_CI <- data.table(CIs)
SS_CI[, Covariate := gsub(":1", "", Covariate)]
SS_CI[,
  Covariate := fcase(
    Covariate == "1-1" , "Coop-Coop" ,
    Covariate == "1-2" , "Coop-Luc"  ,
    Covariate == "2-2" , "Luc-Luc"   ,
    default = Covariate
  )
]

SS_plot <- data.table(SS_CI)
SS_plot <- SS_plot[Covariate != "(Intercept)"]
SS_plot[, ":="(
  Covariate = factor(Covariate, levels = Covariate[1:12]),
  R_between = factor(R_between, levels = r_between)
)]

ggplot(data = SS_plot, mapping = aes(x = Estimate, y = R_between)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI)) +
  theme_bw() +
  labs(y = "R_between") +
  facet_grid(Covariate ~ .) +
  theme(
    strip.text.y.right = element_text(angle = 0, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  ) +
  NULL -> CI_plot

ggsave(
  filename = paste0(
    "2_analysis_bank_data/figures/sens_analysis_rbetween_",
    Sys.Date(),
    ".pdf"
  ),
  plot = CI_plot,
  width = 170,
  height = 260,
  units = "mm"
)

# Sensivity analysis saturation within -----------------------------------------
sat_within <- 1:5
sens_between <- lapply(1:5, function(i) {
  SS <- SemiMarkov(
    X = ppp.type,
    covariate = Covar_dt,
    edgecorrection = 0.064,
    R_within = 0.032,
    R_between = 0.06,
    sat = NULL,
    sat_within = sat_within[i],
    sat_between = 2,
    standardize = TRUE,
    Poisson = FALSE
  )
})

CIs <- lapply(sens_between, function(x) x$CI)
CIs <- lapply(1:5, function(i) data.table(CIs[[i]], sat_within = sat_within[i]))
CIs <- rbindlist(CIs)

SS_CI <- data.table(CIs)
SS_CI[, Covariate := gsub(":1", "", Covariate)]
SS_CI[,
  Covariate := fcase(
    Covariate == "1-1" , "Coop-Coop" ,
    Covariate == "1-2" , "Coop-Luc"  ,
    Covariate == "2-2" , "Luc-Luc"   ,
    default = Covariate
  )
]

SS_plot <- data.table(SS_CI)
SS_plot <- SS_plot[Covariate != "(Intercept)"]
SS_plot[, ":="(
  Covariate = factor(Covariate, levels = Covariate[1:12]),
  sat_within = factor(sat_within, levels = 5:1)
)]

ggplot(data = SS_plot, mapping = aes(x = Estimate, y = sat_within)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI)) +
  theme_bw() +
  labs(y = "Saturation parameter") +
  facet_grid(Covariate ~ .) +
  theme(
    strip.text.y.right = element_text(angle = 0, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  ) +
  NULL -> CI_plot

ggsave(
  filename = paste0(
    "2_analysis_bank_data/figures/sens_analysis_sat_within_,",
    Sys.Date(),
    ".pdf"
  ),
  plot = CI_plot,
  width = 170,
  height = 260,
  units = "mm"
)

# Sensivity analysis saturation between ----------------------------------------
sat_between <- 1:5
sens_between <- lapply(1:5, function(i) {
  SS <- SemiMarkov(
    X = ppp.type,
    covariate = Covar_dt,
    edgecorrection = 0.064,
    R_within = 0.032,
    R_between = 0.06,
    sat = NULL,
    sat_within = 1,
    sat_between = sat_between[i],
    standardize = TRUE,
    Poisson = FALSE
  )
})

CIs <- lapply(sens_between, function(x) x$CI)
CIs <- lapply(1:5, function(i) data.table(CIs[[i]], sat_between = sat_between[i]))
CIs <- rbindlist(CIs)

SS_CI <- data.table(CIs)
SS_CI[, Covariate := gsub(":1", "", Covariate)]
SS_CI[,
  Covariate := fcase(
    Covariate == "1-1" , "Coop-Coop" ,
    Covariate == "1-2" , "Coop-Luc"  ,
    Covariate == "2-2" , "Luc-Luc"   ,
    default = Covariate
  )
]

SS_plot <- data.table(SS_CI)
SS_plot <- SS_plot[Covariate != "(Intercept)"]
SS_plot[, ":="(
  Covariate = factor(Covariate, levels = Covariate[1:12]),
  sat_between = factor(sat_between, levels = 5:1)
)]

ggplot(data = SS_plot, mapping = aes(x = Estimate, y = sat_between)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(xmin = Lower_CI, xmax = Upper_CI)) +
  theme_bw() +
  labs(y = "Saturation parameter") +
  facet_grid(Covariate ~ .) +
  theme(
    strip.text.y.right = element_text(angle = 0, colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  ) +
  NULL -> CI_plot

ggsave(
  filename = paste0(
    "2_analysis_bank_data/figures/sens_analysis_sat_between_,",
    Sys.Date(),
    ".pdf"
  ),
  plot = CI_plot,
  width = 170,
  height = 260,
  units = "mm"
)
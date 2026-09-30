library(spatstat)
library(data.table)
library(SeMaCoPE)

# Point pattern
settings_list <- readRDS(
  "1_simulation_study/simulated_data/sim_covar_and_settings_2026-09-19.rds"
)
sim_list <- readRDS(
  "1_simulation_study/simulated_data/simulated_point_patterns_2026-09-19.rds"
)

run_semacope <- function(X) {
  settings <- settings_list[[X$setting_ids]]

  spat_covar <- settings$spat_covar
  spat_covar_con <- settings$spat_covar_con

  covar_mat <- data.table(
    xcoord = X$point_pattern$x,
    ycoord = X$point_pattern$y,
    Covariate = spat_covar[X$point_pattern]
  )

  if (X$eps != -1) {
    covar_mat[, Covariate_con := spat_covar_con[X$point_pattern]]
  }

  temp <- SemiMarkov(
    X = X$point_pattern,
    covariate = covar_mat,
    R_within = 0.02,
    R_between = 0.04,
    sat = ifelse(X$process == "Geyer", 10, Inf),
    quiet = TRUE
  )

  n <- table(X$point_pattern$marks)

  res <- data.table(
    Parameter = names(temp$betahat),
    Estimate = temp$betahat,
    True_value = X$true_val,
    Std.error = temp$std_err,
    within_CI = between(X$true_val, temp$CI$Lower_CI, temp$CI$Upper_CI),
    window_size = X$win_size,
    n1 = n[1],
    n2 = n[2],
    n3 = n[3],
    Interaction = X$process,
    eps = X$eps
  )

  return(res)
}

rs <- list()
for(i in 1:length(sim_list)){
  cat(i, "\r")
  rs[[i]] <- run_semacope(sim_list[[i]])
}

all_res <- rbindlist(rs)
fwrite(
  all_res,
  paste0("1_simulation_study/results/sim_res_all_", Sys.Date(), ".csv")
)

res <- all_res[,
  .(
    Estimate = mean(Estimate, na.rm = T),
    Std.error = mean(Std.error, na.rm = T),
    MC_std.error = sd(Estimate),
    Coverage = mean(within_CI, na.rm = T),
    n1 = mean(n1),
    n2 = mean(n2),
    n3 = mean(n3)
  ),
  by = list(eps, Interaction, Parameter, True_value, window_size)
]


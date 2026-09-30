library(data.table)
library(ggplot2)
library(ggh4x)

res <- fread("1_simulation_study/results/sim_res_all_2026-09-20.csv")

# Plots of coverage rate -------------------------------------------------------
res_sum <- res[
  eps == 0,
  .(Coverage = mean(within_CI, na.rm = T)),
  list(eps, Interaction, window_size, Parameter)
]
res_sum[, group := paste(Interaction, window_size, sep = "_")]
res_sum[, window_size := factor(window_size)]
min_MC <- 0.95 - 1.96 * sqrt(0.95 * 0.05 / 1800)
max_MC <- 0.95 + 1.96 * sqrt(0.95 * 0.05 / 1800)

res_sum[, window_size := ifelse(window_size == 1, "W[1]", "W[2]")]
res_sum[,
  Parameter := factor(
    Parameter,
    levels = c(
      "(Intercept):1",
      "(Intercept):2",
      "Covariate:1",
      "Covariate:2",
      "Covariate_con:1",
      "Covariate_con:2",
      "1-1",
      "1-2",
      "1-3",
      "2-2",
      "2-3",
      "3-3"
    )
  )
]
res_sum[,
  Interaction := factor(Interaction, levels = c("Poisson", "Strauss", "Geyer"))
]

ggplot(
  data = res_sum,
  aes(x = Parameter, y = Coverage, group = group, color = Interaction)
) +
  geom_point(shape = 21) +
  geom_abline(intercept = min_MC, slope = 0, linetype = "dashed") +
  geom_abline(intercept = 0.95, slope = 0) +
  geom_abline(intercept = max_MC, slope = 0, linetype = "dashed") +
  # geom_line(aes(linetype = window_size)) +
  geom_line() +
  facet_wrap(
    ~window_size,
    nrow = 2,
    strip.position = "right",
    labeller = label_parsed,
    scales = "free_y"
  ) +
  scale_color_manual(values = c("#bc272d", "#0000a2", "#50ad9f")) +
  labs(y = "Coverage rate") +
  scale_x_discrete(
    labels = c(
      expression("\u03b2"["10,1"]),
      expression("\u03b2"["20,1"]),
      expression("\u03b2"["10,2"]),
      expression("\u03b2"["20,2"]),
      expression("\u03b2"["10,3"]),
      expression("\u03b2"["20,3"]),
      expression("\u03b3"["11"]),
      expression("\u03b3"["12"]),
      expression("\u03b3"["13"]),
      expression("\u03b3"["22"]),
      expression("\u03b3"["23"]),
      expression("\u03b3"["33"])
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text.y.right = element_text(angle = 0)
  ) +
  ylim(0.92, 0.98) +
  NULL -> Coverage_plot

filename <- paste0("1_simulation_study/figures/Coverage_plot_main_", Sys.Date(), ".pdf")
ggsave(
  filename = filename,
  plot = Coverage_plot,
  width = 160,
  height = 100,
  units = "mm",
  device = cairo_pdf
)

# Standard error plot ----------------------------------------------------------
res[,
  Parameter := factor(
    Parameter,
    levels = c(
      "(Intercept):1",
      "(Intercept):2",
      "Covariate:1",
      "Covariate:2",
      "Covariate_con:1",
      "Covariate_con:2",
      "1-1",
      "1-2",
      "1-3",
      "2-2",
      "2-3",
      "3-3"
    )
  )
]
res[,
  Interaction := factor(Interaction, levels = c("Poisson", "Strauss", "Geyer"))
]
A1 <- res[
  eps == 0,
  .(`Standard error` = mean(Std.error)),
  .(Parameter, eps, window_size, Interaction)
]
A2 <- res[
  eps == 0,
  .(`Standard error` = sd(Estimate)),
  .(Parameter, eps, window_size, Interaction)
]

A1[,Estimate:="Model esimtate"]
A2[,Estimate:="Monte-Carlo"]
A <- rbind(A1, A2)
A[, group := paste(window_size, Interaction, Estimate, sep = "_")]
A[, window_size := factor(window_size)]

B <- merge(A1, A2, by = c("Parameter", "eps", "window_size", "Interaction"))
B[, div := (`Standard error.x` - `Standard error.y`)/`Standard error.y`]
B[,":="(`Standard error.x` = NULL, `Standard error.y` = NULL)]
B[, window_size := ifelse(window_size == 1, "W[1]", "W[2]")]

ggplot(
  data = B,
  mapping = aes(
    x = Parameter,
    y = div,
    group = Interaction,
    color = Interaction
  )
) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#bc272d", "#0000a2", "#50ad9f")) +
  # scale_shape_manual(values = c(24, 21)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  facet_wrap(
    ~window_size,
    nrow = 2,
    strip.position = "right",
    labeller = label_parsed,
    scales = "free_y"
  ) +
  # facet_wrap(~window_size, labeller = label_parsed, nrow = 2) +
  scale_x_discrete(
    labels = c(
      expression("\u03b2"["10,1"]),
      expression("\u03b2"["20,1"]),
      expression("\u03b2"["10,2"]),
      expression("\u03b2"["20,2"]),
      expression("\u03b2"["10,3"]),
      expression("\u03b2"["20,3"]),
      expression("\u03b3"["11"]),
      expression("\u03b3"["12"]),
      expression("\u03b3"["13"]),
      expression("\u03b3"["22"]),
      expression("\u03b3"["23"]),
      expression("\u03b3"["33"])
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text.y.right = element_text(angle = 0)
  ) +
  labs(y = "Relative difference est vs. MC std. error") +
  # ylim(c(-0.06, 0.06))+
  # expand_limits(y = 0) +
  NULL -> Std_plot

filename <- paste0("1_simulation_study/figures/Std_plot_main_", Sys.Date(), ".pdf")
ggsave(
  filename,
  Std_plot,
  width = 160,
  height = 120,
  units = "mm",
  device = cairo_pdf
)

# Kernel density plots ---------------------------------------------------------
St_dt <- res[eps == 0 & Interaction == "Strauss"]

display_covar <- c("(Intercept):1", "Covariate:1", "1-1", "1-2")
res2 <- St_dt[Parameter %in% display_covar & window_size == 1]
res3 <- res2[,.(Mean = mean(Estimate)), Parameter]
res2 <- merge(res2, res3, all.x = T)

res2[Parameter == "(Intercept):1", Parameter := "beta['10,1']"]
res2[Parameter == "Covariate:1", Parameter := "beta['10,2']"]
res2[Parameter == "1-1", Parameter := "gamma[11]"]
res2[Parameter == "1-2", Parameter := "gamma[12]"]

ggplot(data = res2, aes(x = Estimate)) +
  geom_density() +
  facet_wrap2(
    ~Parameter,
    strip = strip_vanilla(),
    labeller = label_parsed,
    scales = "free"
  ) +
  geom_vline(aes(xintercept = Mean), color = "#bc272d") +
  geom_vline(aes(xintercept = True_value), color = "#0000a2") +
  NULL -> Kernel_density_plot

file_name <- paste0(
  "1_simulation_study/figures/Kernel_density_plot_main_",
  Sys.Date(),
  ".pdf"
)
ggsave(
  filename = file_name,
  plot = Kernel_density_plot,
  width = 160,
  height = 120,
  units = "mm",
  device = cairo_pdf
)

suppl_table <- res[
  eps == 0,
  .(
    True_value = mean(True_value),
    Estimate = mean(Estimate),
    Model_std_error = mean(Std.error),
    MC_std_error = sd(Estimate),
    Coverage = mean(within_CI)
  ),
  .(Parameter, window_size, Interaction)
]

suppl_table[,Parameter:=gsub("Covariate_con", "z2", Parameter)]
suppl_table[,Parameter:=gsub("Covariate", "z1", Parameter)]
setnames(
  suppl_table,
  c("True_value", "Model_std_error", "MC_std_error"),
  c("True value", "Model std. error", "MC std. error")
)
suppl_table_w1 <- suppl_table[window_size == 1]
suppl_table_w2 <- suppl_table[window_size == 4]
suppl_table_w1[,window_size:=NULL]
suppl_table_w2[,window_size:=NULL]

xtable::xtable(suppl_table_w1, digits = 3) |> print(include.rownames = FALSE)
xtable::xtable(suppl_table_w2, digits = 3) |> print(include.rownames = FALSE)

sds_w1 <- suppl_table[window_size == 1, Std.error]
sds_w2 <- suppl_table[window_size == 4, Std.error]
quantile(sds_w2/sds_w1)

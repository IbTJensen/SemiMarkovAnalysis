library(data.table)
library(ggplot2)
library(ggh4x)

Ps_dt <- fread("Poisson results/Sim_Poisson_full_results_2025-07-25.csv")
Pt_dt2 <- fread("Poisson results/Sim_Poisson_res_summary_2025-07-25.csv")
St_dt <- fread("Strauss results/Sim_Strauss_full_results_2025-07-25.csv")
St_dt2 <- fread("Strauss results/Sim_Strauss_res_summary_2025-07-25.csv")
Gs_dt <- fread("Geyer results/Sim_Geyer_full_results_2025-07-26.csv")
Gs_dt2 <- fread("Geyer results/Sim_Geyer_res_summary_2025-07-26.csv")

# Plots of coverage rate -------------------------------------------------------
Ps_dt[,Model:="Poisson"]
St_dt[,Model:="Strauss"]
Gs_dt[,Model:="Geyer"]

res <- rbind(Ps_dt, St_dt, Gs_dt)
res_sum <- res[,.(Coverage = mean(within_CI, na.rm = T)),
               list(Model, window_size, Parameter)]
res_sum[,group:=paste(Model, window_size, sep = "_")]
res_sum[,window_size:=factor(window_size)]
min_MC <- 0.95 - 1.96*sqrt(0.95*0.05/1800)
max_MC <- 0.95 + 1.96*sqrt(0.95*0.05/1800)

res_sum[,window_size:=ifelse(window_size == 1, "W[1]", "W[2]")]
res_sum[, Parameter:=factor(Parameter, levels = c("(Intercept):1",
                                                  "(Intercept):2",
                                                  "Covariate:1", "Covariate:2",
                                                  "1-1", "1-2", "1-3", "2-2",
                                                  "2-3", "3-3"))]
res_sum[,Model:=factor(Model, levels = c("Poisson", "Strauss", "Geyer"))]

ggplot(data = res_sum, aes(x = Parameter, y = Coverage,
                           group = group, color = Model))+
  geom_point(shape = 21) +
  geom_abline(intercept = min_MC, slope = 0, linetype = "dashed") +
  geom_abline(intercept = 0.95, slope = 0) +
  geom_abline(intercept = max_MC, slope = 0, linetype = "dashed") +
  # geom_line(aes(linetype = window_size)) +
  geom_line() +
  facet_wrap(~window_size, nrow = 2, strip.position = "right",
             labeller = label_parsed, scales = "free_y")+
  scale_color_manual(values = c("#bc272d", "#0000a2", "#50ad9f"))+
  labs(y = "Coverage rate")+
  scale_x_discrete(labels = c(expression("\u03b2"["01,1"]),
                              expression("\u03b2"["02,1"]),
                              expression("\u03b2"["01,2"]),
                              expression("\u03b2"["02,2"]),
                              expression("\u03b3"["11"]), expression("\u03b3"["12"]),
                              expression("\u03b3"["13"]), expression("\u03b3"["22"]),
                              expression("\u03b3"["23"]), expression("\u03b3"["33"])))+
  theme_bw()+
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0))+
  ylim(0.925, 0.975)+
  NULL -> Coverage_plot

filename <- paste0("Figures/Coverage_plot_sep_", Sys.Date(), ".pdf")
ggsave(filename = filename, plot = Coverage_plot,
       width = 160, height = 100, units = "mm", device = cairo_pdf)

# Standard error plot ----------------------------------------------------------
St_dt2[,Model:="Strauss"]
Pt_dt2[,Model:="Poisson"]
Gs_dt2[,Model:="Geyer"]
res_sum <- rbind(Gs_dt2, St_dt2, Pt_dt2)

res_sum[, Parameter:=factor(Parameter, levels = c("(Intercept):1",
                                                  "(Intercept):2",
                                                  "Covariate:1", "Covariate:2",
                                                  "1-1", "1-2", "1-3", "2-2",
                                                  "2-3", "3-3"))]
res_sum[,Model:=factor(Model, levels = c("Poisson", "Strauss", "Geyer"))]
A1 <- res_sum[,c("Parameter", "window_size", "Std.error", "Model")]
A2 <- res_sum[,c("Parameter", "window_size", "MC_std.error", "Model")]
colnames(A2)[3] <- "Std.error"
A1[,Estimate:="Model esimtate"]
A2[,Estimate:="Monte-Carlo"]
A <- rbind(A1, A2)
A[,group:=paste(window_size, Model, Estimate, sep = "_")]
A[,window_size:=factor(window_size)]
colnames(A)[5] <- "Standard error"

ggplot(data = A, mapping = aes(x = Parameter, y = Std.error, group = group,
                               color = Model, shape = `Standard error`,
                               linetype = `Standard error`))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("#bc272d", "#0000a2", "#50ad9f"))+
  scale_shape_manual(values = c(24, 21))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_x_discrete(labels = c(expression("\u03b2"["01,1"]),
                              expression("\u03b2"["02,1"]),
                              expression("\u03b2"["01,2"]),
                              expression("\u03b2"["02,2"]),
                              expression("\u03b3"["11"]),
                              expression("\u03b3"["12"]),
                              expression("\u03b3"["13"]),
                              expression("\u03b3"["22"]),
                              expression("\u03b3"["23"]),
                              expression("\u03b3"["33"])))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y = "Standard error")+
  expand_limits(y = 0)+
  NULL -> Std_plot

filename <- paste0("Figures/Std_plot_", Sys.Date(), ".pdf")
ggsave(filename, Std_plot, width = 160, height = 70,
       units = "mm", device = cairo_pdf)

# Kernel density plots ---------------------------------------------------------
St_dt <- merge(St_dt, St_dt2[,c("Parameter", "window_size")],
               by = c("Parameter", "window_size"))

display_covar <- c("(Intercept):1", "Covariate:1", "1-1", "1-2")
res2 <- St_dt[Parameter %in% display_covar & window_size == 1]
res3 <- res2[,.(Mean = mean(Estimate)), Parameter]
res2 <- merge(res2, res3, all.x = T)

res2[Parameter == "(Intercept):1", Parameter := "beta['01,1']"]
res2[Parameter == "Covariate:1", Parameter := "beta['01,2']"]
res2[Parameter == "1-1", Parameter := "gamma[11]"]
res2[Parameter == "1-2", Parameter := "gamma[12]"]
ggplot(data = res2, aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free")+
  geom_vline(aes(xintercept = Mean), color = "#bc272d")+
  geom_vline(aes(xintercept = True_value), color = "#0000a2")+
  NULL -> Kernel_density_plot

ggsave(filename = "Figures/Kernel_density_plot.pdf", plot = Kernel_density_plot,
       width = 160, height = 120, units = "mm", device = cairo_pdf)

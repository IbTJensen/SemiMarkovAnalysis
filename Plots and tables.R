library(data.table)
library(ggplot2)
library(ggh4x)

Ps_dt <- fread("Poisson results/Sim_Poisson_full_results_2025-07-25.csv")
St_dt <- fread("Strauss results/Sim_Strauss_full_results_2025-07-25.csv")
St_dt2 <- fread("Strauss results/Sim_Strauss_res_summary_2025-07-25.csv")
Gs_dt <- fread("Geyer results/Sim_Geyer_full_results_2025-07-26.csv")

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

# Standard error plot Poisson/Strauss ------------------------------------------
St_dt2 <- fread("Strauss results/Sim_Strauss_res_summary_2025-07-25.csv")
St_dt2[,Model:="Strauss"]
Pt_dt2 <- fread("Poisson results/Sim_Poisson_res_summary_2025-07-25.csv")
Pt_dt2[,Model:="Poisson"]
GS_dt2 <- fread("Geyer results/Sim_Geyer_res_summary_2025-07-26.csv")
GS_dt2[,Model:="Geyer"]
res_sum <- rbind(GS_dt2, St_dt2, Pt_dt2)

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
  NULL -> Std_plot_poi_strauss

filename <- paste0("Figures/Std_plot_", Sys.Date(), ".pdf")
ggsave(filename, Std_plot_poi_strauss, width = 160, height = 70,
       units = "mm", device = cairo_pdf)

# Standard error plot Geyer ----------------------------------------------------
Gs_dt2 <- fread("Geyer results/Sim_Geyer_res_summary_2025-07-26.csv")

Gs_dt2[, Parameter:=factor(Parameter, levels = c("(Intercept):1",
                                                  "(Intercept):2",
                                                  "Covariate:1", "Covariate:2",
                                                  "1-1", "1-2", "1-3", "2-2",
                                                  "2-3", "3-3"))]
A1 <- Gs_dt2[,c("Parameter", "window_size", "Std.error")]
A2 <- Gs_dt2[,c("Parameter", "window_size", "MC_std.error")]
colnames(A2)[3] <- "Std.error"
A1[,Estimate:="Model esimtate"]
A2[,Estimate:="Monte-Carlo"]
A <- rbind(A1, A2)
A[,group:=paste(window_size, Estimate, sep = "_")]
A[,window_size:=factor(window_size)]
colnames(A)[4] <- "Standard error"
colnames(A)[2] <- "Window"
A[,Window:=ifelse(Window == 1, "W[1]", "W[2]")]

ggplot(data = A, mapping = aes(x = Parameter, y = Std.error, group = group,
                               color = Window, shape = `Standard error`,
                               linetype = `Standard error`))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c("#bc272d", "#0000a2"),
                     labels = c(expression(W[1], W[2])))+
  scale_shape_manual(values = c(24, 21))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_x_discrete(labels = c(expression("\u03b2"["01,1"]),
                              expression("\u03b2"["02,1"]),
                              expression("\u03b2"["01,2"]),
                              expression("\u03b2"["02,2"]),
                              expression("\u03b3"["11"]), expression("\u03b3"["12"]),
                              expression("\u03b3"["13"]), expression("\u03b3"["22"]),
                              expression("\u03b3"["23"]), expression("\u03b3"["33"])))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y = "Standard error")+
  # coord_cartesian(ylim = c(0,2), clip = "off")+
  NULL -> Std_plot_Geyer

filename <- paste0("Figures/Std_plot_Geyer_", Sys.Date(), ".pdf")
ggsave(filename, Std_plot_Geyer,
       width = 160, height = 70, units = "mm", device = cairo_pdf)

# True vs estimated Geyer ------------------------------------------------------
A1 <- Gs_dt2[,c("Parameter", "window_size", "Estimate")]
A2 <- Gs_dt2[,c("Parameter", "window_size", "True_value")]
colnames(A2)[3] <- "Estimate"
A1[,":="(Value = "Estimate", Model = "Geyer")]
A2[,":="(Value = "True value", Model = "Geyer")]

B1 <- Pt_dt2[,c("Parameter", "window_size", "Estimate")]
B2 <- Pt_dt2[,c("Parameter", "window_size", "True_value")]
colnames(B2)[3] <- "Estimate"
B1[,":="(Value = "Estimate", Model = "Poisson")]
B2[,":="(Value = "True value", Model = "Poisson")]

C1 <- St_dt2[,c("Parameter", "window_size", "Estimate")]
C2 <- St_dt2[,c("Parameter", "window_size", "True_value")]
colnames(C2)[3] <- "Estimate"
C1[,":="(Value = "Estimate", Model = "Strauss")]
C2[,":="(Value = "True value", Model = "Strauss")]

A <- rbind(A1, A2, B1, B2, C1, C2)
A[,group:=paste(window_size, Value, sep = "_")]
A[,window_size:=factor(window_size)]
colnames(A)[2] <- "Window"
A <- A[group != "4_True value"]
A[,Model:=factor(Model, levels = c("Poisson", "Strauss", "Geyer"))]

ggplot(data = A, mapping = aes(x = Parameter, y = Estimate, group = group,
                               color = group, shape = group,
                               linetype = Window))+
  facet_wrap(~Model, nrow = 3, scales = "free_y", strip.position = "right")+
  geom_abline(intercept = 0, slope = 0)+
  geom_line(data = A[Value == "Estimate" & Window == 1], linetype = "solid")+
  geom_line(data = A[Value == "Estimate" & Window == 4], linetype = "dashed")+
  # geom_line()+
  geom_point()+
  # scale_color_manual(values = c("#bc272d", "#0000a2"))+
  scale_color_manual(values = c("#f55f74", "#4a2377", "#0d7d87"))+
  scale_shape_manual(values = c(21, 22, 24))+
  scale_linetype_manual(values = c("dashed", "solid"),
                        labels = c(expression(W[1], W[2])))+
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
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))+
  labs(y = "Parameter value")+
  NULL -> Est_plot_all

filename <- paste0("Figures/Est_plot_", Sys.Date(), ".pdf")
ggsave(filename, Est_plot_all, width = 160,
       height = 120, units = "mm", device = cairo_pdf)

# Kernel density plots ---------------------------------------------------------
St_dt <- merge(St_dt, St_dt2[,c("Parameter", "window_size", "True_value")],
               by = c("Parameter", "window_size"))

plot(density(Ps_dt[Parameter == "1-1" & window_size == 4, Estimate]),
     main = "Estimate of interaction 1-1 on Poisson process")
abline(v = mean(Ps_dt[Parameter == "1-1" & window_size == 4, Estimate]),
       lty = 2)

plot(density(St_dt[Parameter == "1-1" & window_size == 4, Estimate]),
     main = "Estimate of interaction 1-1 on Strauss process")
abline(v = mean(St_dt[Parameter == "1-1" & window_size == 4, Estimate]),
       lty = 2)

plot(density(Gs_dt[Parameter == "1-1" & window_size == 4, Estimate]),
     main = "Estimate of interaction 1-1 on Geyer saturation process")
abline(v = mean(Gs_dt[Parameter == "1-1" & window_size == 4, Estimate]),
       lty = 2)

# res2 <- res[Parameter %in% c("Covariate:1", "1-1", "1-2") & window_size == 1]
res2 <- St_dt[Parameter %in% c("(Intercept):1", "Covariate:1", "1-1", "1-2") & window_size == 1]
res3 <- res2[,.(Mean = mean(Estimate)), Parameter]
res2 <- merge(res2, res3, all.x = T)

res2[Parameter == "(Intercept):1", Parameter := "beta['01,1']"]
res2[Parameter == "Covariate:1", Parameter := "beta['01,2']"]
res2[Parameter == "1-1", Parameter := "gamma[11]"]
res2[Parameter == "1-2", Parameter := "gamma[12]"]
# res2[Parameter == "Covariate:1", Parameter := paste0("(beta[","0","1","])","[2]")]
# res2[Parameter == "Covariate:1", Parameter := expression(paste("(\u03b2"["01"],")"["2"]))]
# res2[Parameter == "1-1", Parameter:="Interaction 1-1"]
ggplot(data = res2, aes(x = Estimate))+
  geom_density()+
  # facet_grid(Parameter~Model, labeller = label_bquote(alpha[.(label)]))+
  # facet_grid2(Model~Parameter, strip = strip_vanilla(), labeller = label_parsed,
  #             scales = "free_y")+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free")+
  geom_vline(aes(xintercept = Mean), color = "#bc272d")+
  geom_vline(aes(xintercept = True_value), color = "#0000a2")+
  NULL -> Kernel_density_plot

ggsave(filename = "Figures/Kernel_density_plot.pdf", plot = Kernel_density_plot,
       width = 160, height = 120, units = "mm", device = cairo_pdf)

res_geyer <- res[Model == "Geyer"]
res_geyer <- res_geyer[Parameter %in% c("1-2", "2-2")]

ggplot(data = res_geyer, aes(x = Estimate))+
  geom_density()+
  # facet_grid(Parameter~Model, labeller = label_bquote(alpha[.(label)]))+
  facet_grid2(window_size~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free_y")+
  geom_vline(xintercept = log(1.2))+
  # geom_vline(aes(xintercept = Mean), linetype = "dashed")+
  NULL -> Kernel_density_plot

# Tables -----------------------------------------------------------------------
res[window_size == 4, .(Std.error = mean(Std.error)), list(Model, Parameter)]


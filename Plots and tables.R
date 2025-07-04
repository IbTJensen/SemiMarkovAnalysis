library(data.table)
library(ggplot2)
library(ggh4x)

Ps_dt <- fread("Poisson results/Poisson_full_results.csv")
St_dt <- fread("Strauss results/Strauss_full_results.csv")
Gs_dt <- fread("Geyer results/Geyer_full_results.csv")

# Plots of coverage rate -------------------------------------------------------
Ps_dt[,":="(Model = "Poisson", n = NULL)]
St_dt[,Model:="Strauss"]
Gs_dt[,Model:="Geyer"]

res <- rbind(Ps_dt, St_dt, Gs_dt)
res_sum <- res[,.(Coverage = mean(within_CI, na.rm = T)),
               list(Model, window_size, Parameter)]
res_sum[,group:=paste(Model, window_size, sep = "_")]
res_sum[,window_size:=factor(window_size)]
min_MC <- 0.95 - 1.96*sqrt(0.95*0.05/1800)
max_MC <- 0.95 + 1.96*sqrt(0.95*0.05/1800)

res_sum[, Parameter:=factor(Parameter, levels = c("(Intercept):1",
                                                  "(Intercept):2",
                                                  "Covariate:1", "Covariate:2",
                                                  "1-1", "1-2", "1-3", "2-2",
                                                  "2-3", "3-3"))]
res_sum[,Model:=factor(Model, levels = c("Poisson", "Strauss", "Geyer"))]

ggplot(data = res_sum, aes(x = Parameter, y = Coverage,
                           group = group, color = Model))+
  geom_point() +
  geom_abline(intercept = min_MC, slope = 0, linetype = "dashed") +
  geom_abline(intercept = 0.95, slope = 0) +
  geom_abline(intercept = max_MC, slope = 0, linetype = "dashed") +
  geom_line(aes(linetype = window_size)) +
  scale_color_manual(values = c("#bc272d", "#50ad9f", "#0000a2"))+
  labs(y = "Coverage rate")+
  scale_x_discrete(labels = c(expression(paste("(\u03b2"["01"],")"["1"])),
                              expression(paste("(\u03b2"["02"],")"["1"])),
                              expression(paste("(\u03b2"["01"],")"["2"])),
                              expression(paste("(\u03b2"["02"],")"["2"])),
                              expression("\u03b3"["11"]), expression("\u03b3"["12"]),
                              expression("\u03b3"["13"]), expression("\u03b3"["22"]),
                              expression("\u03b3"["23"]), expression("\u03b3"["33"])))+
  NULL -> Coverage_plot

ggsave(filename = "Figures/Coverage_plot.pdf", plot = Coverage_plot,
       width = 160, height = 80, units = "mm", device = cairo_pdf)

# Kernel density plots ---------------------------------------------------------
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

res2 <- res[Parameter %in% c("Covariate:1", "1-1", "1-2") & window_size == 4]
res3 <- res2[,.(Mean = mean(Estimate)), list(Model, Parameter)]
res2 <- merge(res2, res3, all.x = T)

res2[Parameter == "Covariate:1", Parameter := "(beta['01'])[2]"]
res2[Parameter == "1-1", Parameter := "gamma[11]"]
res2[Parameter == "1-2", Parameter := "gamma[12]"]
# res2[Parameter == "Covariate:1", Parameter := paste0("(beta[","0","1","])","[2]")]
# res2[Parameter == "Covariate:1", Parameter := expression(paste("(\u03b2"["01"],")"["2"]))]
# res2[Parameter == "1-1", Parameter:="Interaction 1-1"]
ggplot(data = res2, aes(x = Estimate))+
  geom_density()+
  # facet_grid(Parameter~Model, labeller = label_bquote(alpha[.(label)]))+
  facet_grid2(Parameter~Model, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free")+
  geom_vline(aes(xintercept = Mean), linetype = "dashed")+
  NULL -> Kernel_density_plot

ggsave(filename = "Figures/Kernel_density_plot.pdf", plot = Kernel_density_plot,
       width = 160, height = 160, units = "mm", device = cairo_pdf)

# Tables -----------------------------------------------------------------------
res[window_size == 4, .(Std.error = mean(Std.error)), list(Model, Parameter)]


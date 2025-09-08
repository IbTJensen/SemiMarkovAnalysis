library(data.table)
library(ggplot2)
library(ggh4x)

Ps_dt <- fread("Poisson results/Sim_Poisson_full_results_2025-07-25.csv")
St_dt <- fread("Strauss results/Sim_Strauss_full_results_2025-07-25.csv")
Gs_dt <- fread("Geyer results/Sim_Geyer_full_results_2025-07-26.csv")

# Poisson densities ------------------------------------------------------------
Ps_dt[Parameter == "(Intercept):1", Parameter := "beta['01,1']"]
Ps_dt[Parameter == "(Intercept):2", Parameter := "beta['02,1']"]
Ps_dt[Parameter == "Covariate:1", Parameter := "beta['01,2']"]
Ps_dt[Parameter == "Covariate:2", Parameter := "beta['02,2']"]
Ps_dt[Parameter == "1-1", Parameter := "gamma[11]"]
Ps_dt[Parameter == "1-2", Parameter := "gamma[12]"]
Ps_dt[Parameter == "2-2", Parameter := "gamma[22]"]
Ps_dt[Parameter == "2-3", Parameter := "gamma[23]"]
Ps_dt[Parameter == "1-3", Parameter := "gamma[13]"]
Ps_dt[Parameter == "3-3", Parameter := "gamma[33]"]

ggplot(data = Ps_dt[window_size == 1], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
             scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Poisson simulations on W"[1]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Poisson_est_W1

ggplot(data = Ps_dt[window_size == 4], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Poisson simulations on W"[2]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Poisson_est_W2

ggsave(plot = Poisson_est_W1, filename = "Figures/Poisson_all_densities_W1.pdf",
       width = 165, height = 230, units = "mm")

ggsave(plot = Poisson_est_W2, filename = "Figures/Poisson_all_densities_W2.pdf",
       width = 165, height = 230, units = "mm")

# Strauss densities ------------------------------------------------------------
St_dt[Parameter == "(Intercept):1", Parameter := "beta['01,1']"]
St_dt[Parameter == "(Intercept):2", Parameter := "beta['02,1']"]
St_dt[Parameter == "Covariate:1", Parameter := "beta['01,2']"]
St_dt[Parameter == "Covariate:2", Parameter := "beta['02,2']"]
St_dt[Parameter == "1-1", Parameter := "gamma[11]"]
St_dt[Parameter == "1-2", Parameter := "gamma[12]"]
St_dt[Parameter == "2-2", Parameter := "gamma[22]"]
St_dt[Parameter == "2-3", Parameter := "gamma[23]"]
St_dt[Parameter == "1-3", Parameter := "gamma[13]"]
St_dt[Parameter == "3-3", Parameter := "gamma[33]"]

ggplot(data = St_dt[window_size == 1], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Strauss simulations on W"[1]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Strauss_est_W1

ggplot(data = St_dt[window_size == 4], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Strauss simulations on W"[2]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Strauss_est_W2

ggsave(plot = Strauss_est_W1, filename = "Figures/Strauss_all_densities_W1.pdf",
       width = 165, height = 230, units = "mm")

ggsave(plot = Strauss_est_W2, filename = "Figures/Strauss_all_densities_W2.pdf",
       width = 165, height = 230, units = "mm")

# Geyer densities ------------------------------------------------------------
Gs_dt[Parameter == "(Intercept):1", Parameter := "beta['01,1']"]
Gs_dt[Parameter == "(Intercept):2", Parameter := "beta['02,1']"]
Gs_dt[Parameter == "Covariate:1", Parameter := "beta['01,2']"]
Gs_dt[Parameter == "Covariate:2", Parameter := "beta['02,2']"]
Gs_dt[Parameter == "1-1", Parameter := "gamma[11]"]
Gs_dt[Parameter == "1-2", Parameter := "gamma[12]"]
Gs_dt[Parameter == "2-2", Parameter := "gamma[22]"]
Gs_dt[Parameter == "2-3", Parameter := "gamma[23]"]
Gs_dt[Parameter == "1-3", Parameter := "gamma[13]"]
Gs_dt[Parameter == "3-3", Parameter := "gamma[33]"]

ggplot(data = Gs_dt[window_size == 1], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Geyer simulations on W"[1]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Geyer_est_W1

ggplot(data = Gs_dt[window_size == 4], mapping = aes(x = Estimate))+
  geom_density()+
  facet_wrap2(~Parameter, strip = strip_vanilla(), labeller = label_parsed,
              scales = "free", nrow = 5)+
  geom_vline(aes(xintercept = True_value), linetype = "dashed")+
  labs(title = expression("Densities of parameter estimates for Geyer simulations on W"[2]))+
  theme(strip.text = element_text(size = 12))+
  NULL -> Geyer_est_W2

ggsave(plot = Geyer_est_W1, filename = "Figures/Geyer_all_densities_W1.pdf",
       width = 165, height = 230, units = "mm")

ggsave(plot = Geyer_est_W2, filename = "Figures/Geyer_all_densities_W2.pdf",
       width = 165, height = 230, units = "mm")

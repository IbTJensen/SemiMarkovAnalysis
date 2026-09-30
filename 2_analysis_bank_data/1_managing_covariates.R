library(spatstat)
library(data.table)

load("2_analysis_bank_data/data/covariatesfct_alr.rdata")
load("2_analysis_bank_data/data/pointsAndPoly.RData")
Covar_dt <- data.table(xcoord = ppp.type$x, ycoord = ppp.type$y)

for (i in 3:length(covariatesfct)) {
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)[-(1:2)]

CC <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")

Covar_dt[,log_density:=log(density)]
Covar_dt[, ":="(
  prop0.15 = NULL,
  prop25.64 = NULL,
  density = NULL,
  median = NULL,
  poverty = NULL
)]
CC2 <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")

colnames(Covar_dt) <- gsub("prop16.24", "prop.age", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("evolution", "pop.growth", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("decile", "income.decile", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("prop19", "prop.firms", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("proppublic", "prop.public", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("propindustry", "prop.industry", colnames(Covar_dt))
colnames(Covar_dt) <- gsub("proptrade", "prop.trade", colnames(Covar_dt))

fwrite(Covar_dt, file = "2_analysis_bank_data/data/Covariates_processed.csv")

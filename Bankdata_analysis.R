library(spatstat)
library(data.table)
library(SemiMarkov)
library(ggplot2)

# Loading and covariates -------------------------------------------------------
load("Data/covariatesfct_alr.rdata")
load("Data/pointsAndPoly.RData")
load("Data/regions.rdata")
do_grid_search <- F

Covar_dt <- data.table(xcoord = ppp.type$x,
                       ycoord = ppp.type$y)

for(i in 3:length(covariatesfct)){
  cat(i, "\r")
  Covar_dt <- cbind(Covar_dt, covariatesfct[[i]](ppp.type$x, ppp.type$y))
}

colnames(Covar_dt)[-(1:2)] <- names(covariatesfct)[-(1:2)]

# plot(covariatesfct[[3]])

CC <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")
# plot strong correlations
plot(Covar_dt$density, Covar_dt$decile)
plot(Covar_dt[,c("prop0.15", "prop16.24", "prop25.64")])
plot(Covar_dt[,c("activity", "median", "poverty")])
plot(Covar_dt[,c("density", "decile", "proptrade")])
plot(Covar_dt$proptrade, Covar_dt$proppublic)

Covar_dt[,":="(prop16.24 = NULL, prop25.64 = NULL, median = NULL,
               poverty = NULL, proptrade = NULL, decile = NULL)]
CC2 <- cor(Covar_dt[,-(1:2)], use = "pairwise.complete.obs")

plot(Covar_dt$density, Covar_dt$prop0.15)
plot(Covar_dt$activity, Covar_dt$proppublic)

# Grid search for sat and R ----------------------------------------------------
if(do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = NULL,
                   R_within = seq(0.001, 0.01, 0.001),
                   R_between = seq(0.001, 0.01, 0.001),
                   sat = 2:10)
}

# Fit model and interpret results ----------------------------------------------
if(!do_grid_search){
  SS <- SemiMarkov(X = ppp.type,
                   covariate = Covar_dt,
                   edgecorrection = NULL,
                   R_within = 0.009,
                   R_between = 0.004,
                   sat = 4,
                   standardize = T,
                   Poisson = F)
}

SS$R_within # 0.009
SS$R_between # 0.004
SS$maximum_log_likelihood
SS$sat # 4

SS$CI

xt <- xtable::xtable(SS$CI, digits = 3)
print(xt, include.rownames=FALSE)

logit <- function(x) log(x/(1-x))
pred <- SS$pred_no_int[type_pred == 1 & !is.na(prob)]
pred[,code:=covariatesfct[[1]](xcoord, ycoord)]
# pred <- pred[,.(logit_prob = mean(logit(prob)), prob = mean(prob)),code]

# Split in urban and rural -----------------------------------------------------
# Make tessalation of French regions
data <- data[data$code %in% pred$code,]
codes=unique(data$code)
length(codes)
cities <- unique(data$city)

polygs=list()
for (i in 1:length(codes)){
  cat(i, "\r")
  admunit=poly[poly$code==codes[i],]

  polygs[[i]]=owin(poly=list(x=admunit$x,y=admunit$y))
}

tt <- tess(tiles = polygs)

# Urban
urban_codes <- unique(poly[poly$density>500,"code"])
tt_urban <- tess(tiles = polygs[codes %in% urban_codes])

high_density <- inside.owin(ppp.type$x, ppp.type$y, tt_urban)
ppp.urban <- ppp.type[high_density]

Covar_dt[,bank_type:=ppp.type$marks]
Covar_urban <- Covar_dt[density>500]

all(Covar_urban$xcoord == ppp.urban$x)
all(Covar_urban$ycoord == ppp.urban$y)

# Rural
rural_codes <- unique(poly[poly$density<500,"code"])
tt_rural <- tess(tiles = polygs[codes %in% rural_codes])

low_density <- inside.owin(ppp.type$x, ppp.type$y, tt_rural)
ppp.rural <- ppp.type[low_density]

Covar_rural <- Covar_dt[density<500]
all(Covar_rural$xcoord == ppp.rural$x)
all(Covar_rural$ycoord == ppp.rural$y)

# Grid search for high density -------------------------------------------------
# SS_urban <- SemiMarkov(X = ppp.urban,
#                        covariate = Covar_urban,
#                        edgecorrection = NULL,
#                        R_within = seq(0.005, 0.015, 0.001),
#                        R_between = seq(0.003, 0.007, 0.001),
#                        sat = c(2,4,6))

# SS_urban <- SemiMarkov(X = ppp.urban,
#                        covariate = Covar_urban[,-ncol(Covar_urban),with=F],
#                        edgecorrection = NULL,
#                        R_within = seq(0.005, 0.014, 0.001),
#                        R_between = seq(0.003, 0.012, 0.001),
#                        sat = 2:6)

SS_urban <- SemiMarkov(X = ppp.urban,
                       covariate = Covar_urban[,-ncol(Covar_urban),with=F],
                       edgecorrection = NULL,
                       R_within = 0.006,
                       R_between = 0.007,
                       sat = 3)

SS_urban$sat # 3
SS_urban$R_within # 0.006
SS_urban$R_between # 0.007
SS_urban$CI
pred_urban <- SS_urban$pred_no_int[type_pred == 1]
pred_urban2 <- SS_urban$pred[type_pred == 1]

# SS_rural <- SemiMarkov(X = ppp.rural,
#                        covariate = Covar_rural[,-ncol(Covar_rural),with=F],
#                        edgecorrection = NULL,
#                        R_within = seq(0.005, 0.014, 0.001),
#                        R_between = seq(0.003, 0.012, 0.001),
#                        sat = 2:6)

SS_rural <- SemiMarkov(X = ppp.rural,
                       covariate = Covar_rural[,-ncol(Covar_rural),with=F],
                       edgecorrection = NULL,
                       R_within = 0.009,
                       R_between = 0.004,
                       sat = 3)

SS_rural$sat # 3
SS_rural$R_within # 0.009
SS_rural$R_between # 0.004
SS_rural$CI
pred_rural <- SS_rural$pred_no_int[type_pred == 1]
pred_rural2 <- SS_rural$pred[type_pred == 1]

pred_full <- rbind(pred_rural, pred_urban)
pred_full <- pred_full[!is.na(prob)]
pred_full[,code:=covariatesfct[[1]](xcoord, ycoord)]
pred_full[,logit_prob:=logit(prob)]

pred_full2 <- rbind(pred_rural2, pred_urban2)
pred_full2 <- pred_full2[!is.na(prob)]
pred_full2[,code:=covariatesfct[[1]](xcoord, ycoord)]
pred_full2[,logit_prob:=logit(prob)]
# pred_full <- pred_full[match(codes, pred_full$code)]

# Separate rural and urban regions ---------------------------------------------
# Covar_urban <- as.data.frame(Covar_dt[,-(1:2)]) * as.numeric(Covar_dt$density>500)
# Covar_rural <- as.data.frame(Covar_dt[,-(1:2)]) * as.numeric(Covar_dt$density<500)
# colnames(Covar_urban) <- paste(colnames(Covar_urban), "urban", sep = "_")
# colnames(Covar_rural) <- paste(colnames(Covar_rural), "rural", sep = "_")
# Covar_full <- cbind(Covar_dt[,1:2], Covar_urban, Covar_rural)


# Plot -------------------------------------------------------------------------
data <- data[data$code %in% pred$code,]
codes=unique(data$code)
length(codes)
cities <- unique(data$city)
# pred <- pred[match(codes, pred$code)]

polygs=list()
for (i in 1:length(codes)){
  cat(i, "\r")
  admunit=poly[poly$code==codes[i],]

  polygs[[i]]=owin(poly=list(x=admunit$x,y=admunit$y))
}
tt <- tess(tiles = polygs)

# pred_fct <- as.function(tt, values = pred$logit_prob)

# pdf("Figures/pred_logit.pdf", width = 180/25.4, height = 180/25.4)
# plot(pred_fct, main = "Predicted logit-proportion of cooperative banks")
# dev.off()

# Exploratory ------------------------------------------------------------------
# summary_ratios <- poly_dt[,.(log_ratio = mean(log(cooperative/lucrative)),
#                              lucrative = mean(lucrative),
#                              cooperative = mean(cooperative)), city]
#
# summary_ratios <- summary_ratios[abs(log_ratio) != Inf & !is.na(log_ratio)]
# summary_ratios[order(log_ratio)]

# model checks -----------------------------------------------------------------
pred[,logit_prob:=logit(prob)]
pred <- merge(pred, SS$h[,-4])
# pred <- cbind(pred, SS$h)
poly_dt <- data.table(poly)
pred_to_merge <- pred[,-(1:5)]
pred_to_merge <- pred_to_merge[,lapply(.SD, mean), code]
poly_dt <- merge(poly_dt, pred_to_merge, by = "code", all = T)

model_check <- merge(Covar_dt, pred[,c(1:2,8)])
# model_check[, lapply(.SD, mean, na.rm = T), `Intercept:1`<0,
#             .SDcols = c("density", "prop0.15", "prop16.24", "prop25.64",
#                         "evolution", "poverty", "activity", "median",
#                         "decile", "prop19", "proppublic", "propindustry",
#                         "proptrade")]

model_check[, mean(`Intercept:1`), `Intercept:1`<0]

ggplot(poly_dt, aes(x, y, group = code, fill = `density:1`))+
  geom_polygon(colour = "black")+
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))

ggplot(poly_dt, aes(x, y, group = code))+
  geom_polygon(colour = "black", fill = "grey80", linewidth = 0.1)+
  geom_point(data = pred, aes(x = xcoord, y = ycoord, colour = `density:1`),
             size = 1)+
  scale_colour_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))+
  NULL

# Brier scores
SS_pred <- SS$pred[type_pred == 1]
SS_obs <- data.table(xcoord = ppp.type$x, ycoord = ppp.type$y,
                     coop = ifelse(ppp.type$marks == "cooperative", 1, 0))
SS_dt <- merge(SS_obs, SS_pred[,c(1,2,5)], by = c("xcoord", "ycoord"))
SS_dt <- SS_dt[!is.na(prob)]
mean(abs(SS_dt$coop - SS_dt$prob)^2) # Overall Brier score
mean((1 - SS_dt[coop == 1, prob])^2) # Brier score for cooperative
mean((0 - SS_dt[coop == 0, prob])^2) # Brier score for lucrative
# Aberage Brier score of the two types
0.5*(mean((1 - SS_dt[coop == 1, prob])^2) + mean((0 - SS_dt[coop == 0, prob])^2))

SS_comb_pred <- pred_full2[type_pred == 1]
SS_comb_obs <- data.table(xcoord = ppp.type$x, ycoord = ppp.type$y,
                          coop = ifelse(ppp.type$marks == "cooperative", 1, 0))
SS_comb_dt <- merge(SS_comb_obs, SS_comb_pred[,c(1,2,5)],
                    by = c("xcoord", "ycoord"))
SS_comb_dt <- SS_comb_dt[!is.na(prob)]
mean(abs(SS_comb_dt$coop - SS_comb_dt$prob)^2) # Overall Brier score
mean((1 - SS_comb_dt[coop == 1, prob])^2) # Brier score for cooperative
mean((0 - SS_comb_dt[coop == 0, prob])^2) # Brier score for lucrative
# Aberage Brier score of the two types
0.5*(mean((1 - SS_comb_dt[coop == 1, prob])^2) + mean((0 - SS_comb_dt[coop == 0, prob])^2))

# Plot with interactions -------------------------------------------------------
# france_map <- map_data("france")

# pred2 <- SS$pred[type_pred == 1 & !is.na(prob)]
# pred2[,code:=covariatesfct[[1]](xcoord, ycoord)]
# pred2[,logit_prob:=logit(prob)]

# One model for whole France
pred <- pred[,1:7]
poly_dt <- data.table(poly)
pred_to_merge <- pred[,c("code", "logit_prob")]
pred_to_merge <- pred_to_merge[,.(logit_prob = mean(logit_prob)), code]
poly_dt <- merge(poly_dt, pred_to_merge, by = "code", all = T)

ggplot(poly_dt, aes(x, y, group = code, fill = logit_prob))+
  geom_polygon(colour = "black", linewidth = 0.1)+
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"),
                       name = "Logit-probability")+
  labs(x = "Longitude", y = "Latitude")+
  NULL -> g
ggsave("Figures/Pred_plot.pdf", width = 150, height = 110, units = "mm")

# Urban-rural model
pred_full[,logit_prob:=logit(prob)]
poly_dt <- data.table(poly)
pred_to_merge <- pred_full[,c("code", "logit_prob")]
pred_to_merge <- pred_to_merge[,.(logit_prob = mean(logit_prob)), code]
poly_dt <- merge(poly_dt, pred_to_merge, by = "code", all = T)

ggplot(poly_dt, aes(x, y, group = code, fill = logit_prob))+
  geom_polygon(colour = "black")+
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))

ggplot(poly_dt[code %in% urban_codes], aes(x, y, group = code, fill = logit_prob))+
  geom_polygon(colour = "black")+
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))

ggplot(poly_dt[code %in% rural_codes], aes(x, y, group = code, fill = logit_prob))+
  geom_polygon(colour = "black")+
  scale_fill_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))

# Old
data_obs <- data.table(x = ppp.type$x, y = ppp.type$y, type = ppp.type$marks)

ggplot(poly_dt, aes(x, y, group = code))+
  geom_polygon(colour = "black", fill = "grey90")+
  geom_point(data = data_obs, aes(x = x, y = y, group = 1, colour = type),
             size = 0.1, alpha = 0.5)

ggplot(poly_dt, aes(x, y, group = code))+
  geom_polygon(colour = "black", fill = "grey90")+
  geom_point(data = pred_full2, aes(x = xcoord, y = ycoord,
                                group = 1, colour = logit_prob), size = 0.1)+
  # scale_colour_gradientn(colours = rainbow(8))+
  scale_colour_gradientn(colours = c("yellow", "orange", "red", "purple", "blue"))+
  NULL

ggplot(poly_dt, aes(x, y, group = code))+
  geom_polygon(colour = "black", fill = "grey90")+
  geom_point(data = pred2, aes(x = xcoord, y = ycoord,
                               group = 1, colour = prob_logit), size = 0.1)+
  scale_colour_gradientn(colours = rainbow(8))

ggplot(france_map, aes(long, lat, group = group)) +
  geom_polygon(fill = "grey30", colour = NA) +
  geom_point(data = pred2, aes(x = xcoord, y = ycoord,
                               group = 1, colour = prob_logit),
             size = 0.1)+
  scale_color_gradientn(colours = rainbow(8))+
  coord_quickmap()



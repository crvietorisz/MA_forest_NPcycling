# MA_forest_NPcycling 
# Stepwise model selection for generalized linear model explaining net nitrification
# 11/21/24 Corinne Vietorisz

library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(ggfortify)

setwd("MA_forest_NPcycling/")

# read in full data
sub_dist <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)

#make list of nitrification predictors

Nitr_pred_list <- c("Nitrification", "pH", "litter_depth",
                    "total_floor_litter",
                    "Ferns", "hw_floor_litter", "AM_basal",
                    "num_stems",
                    "AM_abundance", "med_dist_expl", "Nitrifier_abundance",
                    "Copiotroph_abundance",
                    "denitrification_ECbac_abund",
                    "bac_nitr_pos_module_abundance",
                    "bac_shannon", "fun_shannon", "ITS_copy_number", "EMF_Shannon"
)


#create subsetted dataset with just Nitrification and predictors
Nitr_preds <- sub_dist[,c(Nitr_pred_list)]

#visualize covariance of the predictors in PCA space
pca.nit <- princomp(na.omit(Nitr_preds), cor = TRUE)

barplot(pca.nit$loadings[,1]) 
barplot(pca.nit$loadings[,2])

p <- autoplot(pca.nit, data = na.omit(Nitr_preds),
              loadings = TRUE, loadings.colour = 'blue',
              loadings.label = FALSE, loadings.label.size = 3, loadings.label.colour = 'red3',
)
p + coord_cartesian(xlim = c(-0.3, 0.35))

#When variables that correlate strongly in PC space, choose only the strongest predictor from the cluster


########################### Make full model ####################################

mod.full <- glm((Nitrification+1) ~ litter_depth + denitrification_ECbac_abund + med_dist_expl + bac_shannon + Ferns + Nitrifier_abundance + pH + Copiotroph_abundance + num_stems + bac_nitr_pos_module_abundance + ITS_copy_number, 
                    data = sub_dist,
                    family = Gamma(link = "log"))
summary(mod.full)
MuMIn::r.squaredGLMM(mod.full)
AIC(mod.full) 


########################### Model selection rounds ####################################

######## ROUND 1 ###########
vars_sel_1 <- c(
  "litter_depth",
  "denitrification_ECbac_abund",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Nitrifier_abundance", 
  "pH", 
  "Copiotroph_abundance",
  "num_stems",
  "bac_nitr_pos_module_abundance",
  "ITS_copy_number"
)

# create a dataframe to hold the AIC info
aic_df_1 <- data.frame(
  Model = character(length(vars_sel_1)),
  AIC = numeric(length(vars_sel_1)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_1)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_1[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_1$Model[i] <- paste("Model without", vars_sel_1[i])
  aic_df_1$AIC[i] <- AIC(model)
}

# REMOVING NUM STEMS


######## ROUND 2 ###########

vars_sel_2 <- c(
  "litter_depth",
  "denitrification_ECbac_abund",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Nitrifier_abundance", 
  "pH", 
  "Copiotroph_abundance",
  "bac_nitr_pos_module_abundance",
  "ITS_copy_number"
)


# create a dataframe to hold the AIC info
aic_df_2 <- data.frame(
  Model = character(length(vars_sel_2)),
  AIC = numeric(length(vars_sel_2)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_2)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_2[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_2$Model[i] <- paste("Model without", vars_sel_2[i])
  aic_df_2$AIC[i] <- AIC(model)
}

# REMOVING BAC MODULE

######## ROUND 3 ###########

vars_sel_3 <- c(
  "litter_depth",
  "denitrification_ECbac_abund",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Nitrifier_abundance", 
  "pH", 
  "Copiotroph_abundance",
  "ITS_copy_number"
)


# create a dataframe to hold the AIC info
aic_df_3 <- data.frame(
  Model = character(length(vars_sel_3)),
  AIC = numeric(length(vars_sel_3)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_3)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_3[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_3$Model[i] <- paste("Model without", vars_sel_3[i])
  aic_df_3$AIC[i] <- AIC(model)
}

# REMOVING DENITRIFICATION EC

######## ROUND 4 ###########

vars_sel_4 <- c(
  "litter_depth",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Nitrifier_abundance", 
  "pH", 
  "Copiotroph_abundance",
  "ITS_copy_number"
)

# create a dataframe to hold the AIC info
aic_df_4 <- data.frame(
  Model = character(length(vars_sel_4)),
  AIC = numeric(length(vars_sel_4)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_4)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_4[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_4$Model[i] <- paste("Model without", vars_sel_4[i])
  aic_df_4$AIC[i] <- AIC(model)
}

# REMOVING NITRIFIER ABUNDANCE

######## ROUND 5 ###########

vars_sel_5 <- c(
  "litter_depth",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "pH", 
  "Copiotroph_abundance",
  "ITS_copy_number"
)

# create a dataframe to hold the AIC info
aic_df_5 <- data.frame(
  Model = character(length(vars_sel_5)),
  AIC = numeric(length(vars_sel_5)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_5)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_5[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_5$Model[i] <- paste("Model without", vars_sel_5[i])
  aic_df_5$AIC[i] <- AIC(model)
}

# REMOVING PH

######## ROUND 6 ###########

vars_sel_6 <- c(
  "litter_depth",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Copiotroph_abundance",
  "ITS_copy_number"
)

# create a dataframe to hold the AIC info
aic_df_6 <- data.frame(
  Model = character(length(vars_sel_6)),
  AIC = numeric(length(vars_sel_6)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_6)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_6[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_6$Model[i] <- paste("Model without", vars_sel_6[i])
  aic_df_6$AIC[i] <- AIC(model)
}

# REMOVING ITS COPY NUMBER


######## ROUND 7 ###########

vars_sel_7 <- c(
  "litter_depth",
  "med_dist_expl",
  "bac_shannon",
  "Ferns",
  "Copiotroph_abundance"
)


# create a dataframe to hold the AIC info
aic_df_7 <- data.frame(
  Model = character(length(vars_sel_7)),
  AIC = numeric(length(vars_sel_7)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_7)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_7[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_7$Model[i] <- paste("Model without", vars_sel_7[i])
  aic_df_7$AIC[i] <- AIC(model)
}

# REMOVING BAC SHANNON

######## ROUND 8 ###########

vars_sel_8 <- c(
  "litter_depth",
  "med_dist_expl",
  "Ferns",
  "Copiotroph_abundance"
)


# create a dataframe to hold the AIC info
aic_df_8 <- data.frame(
  Model = character(length(vars_sel_8)),
  AIC = numeric(length(vars_sel_8)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_8)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_8[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_8$Model[i] <- paste("Model without", vars_sel_8[i])
  aic_df_8$AIC[i] <- AIC(model)
}

# REMOVING FERNS

######## ROUND 9 ###########

vars_sel_9 <- c(
  "litter_depth",
  "med_dist_expl",
  "Copiotroph_abundance"
)

# create a dataframe to hold the AIC info
aic_df_9 <- data.frame(
  Model = character(length(vars_sel_9)),
  AIC = numeric(length(vars_sel_9)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_9)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_9[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
                   data = sub_dist,
                   family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_9$Model[i] <- paste("Model without", vars_sel_9[i])
  aic_df_9$AIC[i] <- AIC(model)
}

# REMOVING LITTER DEPTH

######## ROUND 10 ###########

vars_sel_10 <- c(
  "med_dist_expl",
  "Copiotroph_abundance"
)

# create a dataframe to hold the AIC info
aic_df_10 <- data.frame(
  Model = character(length(vars_sel_10)),
  AIC = numeric(length(vars_sel_10)),
  stringsAsFactors = FALSE
)

# set response variable
response <- "(Nitrification+1)"

# run for loop that drops one variable at a time and records the AIC for each model
for (i in 1:length(vars_sel_10)) {
  # Select predictors excluding the i-th predictor
  current_predictors <- vars_sel_10[-i]
  
  # Create the formula for the model
  formula_str <- paste(response, "~", paste(current_predictors, collapse = " + "))
  
  # Fit the linear model
  model <- glm(as.formula(formula_str),
               data = sub_dist,
               family = Gamma(link = "log"))
  
  # Put AIC of each model into df
  aic_df_10$Model[i] <- paste("Model without", vars_sel_10[i])
  aic_df_10$AIC[i] <- AIC(model)
}



# YAY WE ARE DONE!! so fun


############### BUILD FINAL MODEL ####################

# Pick the model from all rounds with an AIC within 3 of the lowest AIC (251) with the fewest parameters
mod.final <- glm((Nitrification+1) ~  litter_depth + Ferns + med_dist_expl + Copiotroph_abundance + bac_shannon + ITS_copy_number,
                     data = sub_dist,
                     family = Gamma(link = "log"))
summary(mod.final)
AIC(mod.final) # 253.6
MuMIn::r.squaredGLMM(mod.final) # R2 = 0.47

# ITS_copy_number is not significant, so removing that parameter
mod.final2 <- glm((Nitrification+1) ~  litter_depth + Ferns + med_dist_expl + Copiotroph_abundance + bac_shannon,
                 data = sub_dist,
                 family = Gamma(link = "log"))
summary(mod.final2)
AIC(mod.final2) # 258.6
anova(mod.final2)
MuMIn::r.squaredGLMM(mod.final2) # R2 = 0.47
confint(mod.final2, 'litter_depth', level = 0.95)
confint(mod.final2, 'Ferns', level = 0.95)
confint(mod.final2, 'med_dist_expl', level = 0.95)
confint(mod.final2, 'Copiotroph_abundance', level = 0.95)
confint(mod.final2, 'bac_shannon', level = 0.95)


# use the DHARMa package to test for over-dispersion
testDispersion(mod.final2) # no over-dispersion

### compare to model without microbial predictors

mod.nomic <- glm((Nitrification+1) ~  litter_depth + Ferns,
                     data = sub_dist,
                     family = Gamma(link = "log"))
summary(mod.nomic) 
AIC(mod.nomic) # 391.2
MuMIn::r.squaredGLMM(mod.nomic) # R2 = 0.09


# MA_forest_NPcycling 
# Multivariate linear model selection for transect-aggregated models explaining net ammonification and net phosphate change
# 9/13/24 C. Vietorisz

library(lme4)
library(MuMIn)
library(leaps)

setwd("MA_forest_NPcycling/")

##############################################################################################################################
# NET AMMONIFICATION
##############################################################################################################################

# read in data
amm_data <- read.csv("Data/Transect aggregated/ammon_data_transect_12Sep24.csv", row.names=1)
amm_data$Site <- c(rep("HF1",4), rep("HF2",4), rep("HF3",4), rep("US1",4), rep("US2",4), rep("US3",4)) #add in site info

hist(amm_data$Ammonification)

# select out input predictors for model selection
Ammon_sig_preds <- amm_data[c("Ammonification", "soil_temp", "soil_moisture", "SOM", 'med_dist_expl', 'NH4uptake_ECfun_abund', 'all_fungal_Ndecomp_ECfun', 'all_EMF_Ndecomp_ECfun', 'all_bac_Ndecomp_ECbac', 'all_copio_Ndecomp_ECbac')]

# format for input into leaps
Ammon_sig_preds_leaps <- na.omit(Ammon_sig_preds)
Ammon_sig_leaps_names <- colnames(Ammon_sig_preds_leaps[2:ncol(Ammon_sig_preds_leaps)]) #change this range to the size of your input dataframe, minus your response variable

# Use the leaps package to test all subsets of our data using Mallow's Cp as our selection criteria
Ammon_sig_leaps <- leaps(x=Ammon_sig_preds_leaps[c(2:ncol(Ammon_sig_preds_leaps))], y=Ammon_sig_preds_leaps$Ammonification,method="Cp")

#plot all possible models by the number of parameters vs. the Cp score
plot(Ammon_sig_leaps$size,Ammon_sig_leaps$Cp,
     xlab="Number of Parameters",
     ylab="Cp",pch=16,col="blue")
abline(0,1)

# create a dataframe that we can use to interpret results
cp.dat <- cbind(Ammon_sig_leaps$size,
                abs(Ammon_sig_leaps$Cp-(Ammon_sig_leaps$size)),
                Ammon_sig_leaps$which)                           
rownames(cp.dat) <- 1:length(Ammon_sig_leaps$size)                 # add rownames
colnames(cp.dat)[1:2] <- c("N.p","Cp-p")                      # add colnames for p, Cp-p-1
colnames(cp.dat)[3:ncol(cp.dat)] <- Ammon_sig_leaps_names

# pick the model(s) where the difference between p and Cp is the smallest
cp3 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==3,])
cp4 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==4,])

#### CP3 ####

# create lists of variables in each model combination
result <- lapply(1:nrow(cp3), function(i) {
  row_data <- cp3[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 2)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:ncol(new_df))
  new_df
})

# create empty matrix to hold AIC values
aic_matrix <- matrix(nrow = 10, ncol =1)
colnames(aic_matrix) <- "CP3 AIC"

# run an lmer model with each variable combination and input the AICs into the matrix
for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("Ammonification", "~", paste(row, collapse = " + "), "+ (1|Site)")) # include Site as a random effect
  model <- lmer(formula, data = amm_data)
  aic_matrix[x,1] <- AIC(model)
}

#set x equal to the lowest AIC row then run through the code within the for loop to determine which variables are included in that combination
x <- 1
summary(model)

#### CP4 ####

# create lists of variables in each model combination
result <- lapply(1:nrow(cp4), function(i) {
  row_data <- cp4[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 3)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:3) # change to ncol
  new_df
})

# create empty matrix to hold AIC values
aic_matrix <- matrix(nrow = 10, ncol =1)
colnames(aic_matrix) <- "CP4 AIC"

# run an lmer model with each variable combination and input the AICs into the matrix
for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("Ammonification", "~", paste(row, collapse = " + "), "+ (1|Site)"))
  model <- lmer(formula, data = amm_data)
  aic_matrix[x,1] <- AIC(model)
}

#set x equal to the lowest AIC row then run through the code within the for loop to determine which variables are included in that combination
x <- 1
summary(model)

##### Build the lowest AIC model that has a Cp score under 2 ########

Ammon.lmer.1 <- lmer(Ammonification ~ soil_temp + soil_moisture + med_dist_expl + (1|Site), data = amm_data)
summary(Ammon.lmer.1)
pt(q=1.395, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #0.08
pt(q=3.481, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #0.0003
pt(q=2.713, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #0.003
MuMIn::r.squaredGLMM(Ammon.lmer.1) #R2m = 0.51 - this is the amount of variance explained by predictor variables by themselves (excluding the random effect)
AIC(Ammon.lmer.1) # 130

######## assess linear independence of predictors
plot(amm_data$soil_moisture, amm_data$soil_temp) #not highly correlated
plot(amm_data$soil_moisture, amm_data$med_dist_expl) #not highly correlated
plot(amm_data$med_dist_expl, amm_data$soil_temp) #not highly correlated


##############################################################################################################################
# NET NITRIFICATION
##############################################################################################################################

# read in data
nitr_data <- read.csv("Data/Transect aggregated/nitrif_data_transect_12Sep24.csv", row.names=1)
nitr_data$Site <- c(rep("HF1",4), rep("HF2",4), rep("HF3",4), rep("US1",4), rep("US2",4), rep("US3",4))

hist(log(nitr_data$Nitrification+0.01))

# select out input predictors for model selection
Nitr_sig_preds <- nitr_data[c("Nitrification", 'pH', 'litter_depth', 'SOM', 'AM_basal', 'total_basal', 'no_pine_saplings', 'ECM_abundance', 'Nitrifier_abundance', 'Copiotroph_abundance', 'fun_nitr_neg_module_abundance', 'bac_richness', 'fun_shannon', 'denitrification_ECbac_abund', 'Copios_denitrification_ECbac', 'S_litterfall', 'C_litterfall')] 

# format for input into leaps
Nitr_sig_preds_leaps <- na.omit(Nitr_sig_preds)
Nitr_sig_leaps_names <- colnames(Nitr_sig_preds_leaps[2:ncol(Nitr_sig_preds_leaps)]) #change this range to the size of your input dataframe, minus your response variable

# Use the leaps package to test all subsets of our data using Mallow's Cp as our selection criteria
Nitr_sig_leaps <- leaps(x=Nitr_sig_preds_leaps[c(2:ncol(Nitr_sig_preds_leaps))], y=log(nitr_data$Nitrification+0.01), method="Cp")

#plot all possible models by the number of parameters vs. the Cp score
plot(Nitr_sig_leaps$size,Nitr_sig_leaps$Cp,
     xlab="Number of Parameters",
     ylab="Cp",pch=16,col="blue")
abline(0,1)

# create a dataframe that we can use to interpret results
cp.dat <- cbind(Nitr_sig_leaps$size,
                abs(Nitr_sig_leaps$Cp-(Nitr_sig_leaps$size)),
                Nitr_sig_leaps$which)                           
rownames(cp.dat) <- 1:length(Nitr_sig_leaps$size)                 # add rownames
colnames(cp.dat)[1:2] <- c("N.p","Cp-p")                      # add colnames for p, Cp-p-1
colnames(cp.dat)[3:ncol(cp.dat)] <- Nitr_sig_leaps_names

# pick the model(s) where the difference between p and Cp is the smallest
cp3 <- cp.dat[cp.dat[,"N.p"]==3,]

#### CP3 ####

# create lists of variables in each model combination
result <- lapply(1:nrow(cp3), function(i) {
  row_data <- cp3[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 2)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:ncol(new_df))
  new_df
})

# create empty matrix to hold AIC values
nitr_cp3_aic_matrix <- matrix(nrow = 10, ncol =1)
colnames(nitr_cp3_aic_matrix) <- "CP3 AIC"

# run an lmer model with each variable combination and input the AICs into the matrix
for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("log(Nitrification +0.01)", "~", paste(row, collapse = " + "), "+ (1|Site)")) # include Site as a random effect
  model <- lmer(formula, data = nitr_data)
  nitr_cp3_aic_matrix[x,1] <- AIC(model)
}

#set x equal to the lowest AIC row then run through the code within the for loop to determine which variables are included in that combination
x <- 4
summary(model)

##### Build the lowest AIC model where all predictors are significant ########
# NOTE!! All of the Cp = 3 models that had a Cp score under 2 had one predictor that was not significant. So, we chose the next lowest AIC model with all predictors significant from all Cp=3 options 

Nitr.lmer.1 <- lmer(log(Nitrification +0.01) ~ Copiotroph_abundance + fun_shannon + (1|Site), data = nitr_data)
summary(Nitr.lmer.1)
pt(q=3.377, df=length(na.omit(nitr_data$Nitrification)-2), lower.tail=FALSE) #0.001 - Copiotroph_abundance
pt(q=1.861, df=length(na.omit(nitr_data$Nitrification)-2), lower.tail=FALSE) #0.04 - fun_shannon
MuMIn::r.squaredGLMM(Nitr.lmer.1) #R2m = 0.55 - this is the amount of variance explained by predictor variables by themselves (excluding the random effect)
AIC(Nitr.lmer.1) # 69

######## assess linear independence of predictors
plot(nitr_data$Copiotroph_abundance, nitr_data$fun_shannon) #not highly correlated


##############################################################################################################################
# NET PO4 CHANGE
##############################################################################################################################

# read in data
pmin_data <- read.csv("Data/Transect aggregated/pmin_data_transect_12Sep24.csv", row.names=1)
pmin_data$Site <- c(rep("HF1",4), rep("HF2",4), rep("HF3",4), rep("US1",4), rep("US2",4), rep("US3",4))

# select out input predictors for model selection
pmin_sig_preds <- pmin_data[c('PO4_release', 'hw_floor_litter', 'pine_basal', 'bac_pmin_neg_module_abundance', 'bac_copy_number', 'soil_temp2022', 'B_litterfall', 'Fe_litterfall', 'Mg_litterfall', 'Mo_litterfall')]

# format for input into leaps
pmin_sig_preds_leaps <- na.omit(pmin_sig_preds)
pmin_sig_leaps_names <- colnames(pmin_sig_preds_leaps[2:ncol(pmin_sig_preds_leaps)]) #change this range to the size of your input dataframe, minus your response variable

# Use the leaps package to test all subsets of our data using Mallow's Cp as our selection criteria
pmin_sig_leaps <- leaps(x=pmin_sig_preds_leaps[c(2:ncol(pmin_sig_preds_leaps))], y=pmin_data$PO4_release, method="Cp")

#plot all possible models by the number of parameters vs. the Cp score
plot(pmin_sig_leaps$size,pmin_sig_leaps$Cp,
     xlab="Number of Parameters",
     ylab="Cp",pch=16,col="blue")
abline(0,1)

# create a dataframe that we can use to interpret results
cp.dat <- cbind(pmin_sig_leaps$size,
                abs(pmin_sig_leaps$Cp-(pmin_sig_leaps$size)),
                pmin_sig_leaps$which)                           
rownames(cp.dat) <- 1:length(pmin_sig_leaps$size)                 # add rownames
colnames(cp.dat)[1:2] <- c("N.p","Cp-p")                      # add colnames for p, Cp-p-1
colnames(cp.dat)[3:ncol(cp.dat)] <- pmin_sig_leaps_names

# pick the model(s) where the difference between p and Cp is the smallest
cp4 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==4,])
cp5 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==5,])

#### CP4 ####

# create lists of variables in each model combination
result <- lapply(1:nrow(cp4), function(i) {
  row_data <- cp4[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 3)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:3)
  new_df
})

# create empty matrix to hold AIC values
pmin_cp4_aic_matrix <- matrix(nrow = 10, ncol =1)
colnames(pmin_cp4_aic_matrix) <- "CP4 AIC"

# run an lm model with each variable combination and input the AICs into the matrix
# note: not including site as a random effect here because it does not explain any variation when included
for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("PO4_release", "~", paste(row, collapse = " + ")))
  model <- lm(formula, data = pmin_data)
  pmin_cp4_aic_matrix[x,1] <- AIC(model)
}

#set x equal to the lowest AIC row then run through the code within the for loop to determine which variables are included in that combination
x <- 1
summary(model)

#### CP5 ####
result <- lapply(1:nrow(cp5), function(i) {
  row_data <- cp5[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 4)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:4)
  new_df
})

pmin_cp5_aic_matrix <- matrix(nrow = 10, ncol =1)
colnames(pmin_cp5_aic_matrix) <- "CP5 AIC"

for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("PO4_release", "~", paste(row, collapse = " + ")))
  model <- lm(formula, data = pmin_data)
  pmin_cp5_aic_matrix[x,1] <- AIC(model)
}

#set x equal to the lowest AIC row then run through the code within the for loop to determine which variables are included in that combination
x <- 7
summary(model)

##### Build the lowest AIC model that has a Cp score under 2 ########

Pmin.lm1 <- lm(PO4_release ~ bac_pmin_neg_module_abundance + Mg_litterfall + Mo_litterfall, data = pmin_data)
summary(Pmin.lm1)
AIC(Pmin.lm1) #-110

######## assess linear independence of predictors
plot(pmin_data$Mg_litterfall, pmin_data$bac_pmin_neg_module_abundance) #not highly correlated
plot(pmin_data$Mg_litterfall, pmin_data$Mo_litterfall) #not highly correlated
plot(pmin_data$Mo_litterfall, pmin_data$bac_pmin_neg_module_abundance) #not highly correlated


##############################################################################################################################
# TEST LINEAR MODEL ASSUMPTIONS
##############################################################################################################################

#RESIDUALS PLOT

mod <- Pmin.lm1 # change to your model

#get unstandardized predicted and residual values
unstandardizedPredicted <- predict(mod)
unstandardizedResiduals <- resid(mod)

#get standardized values
standardizedPredicted <- (unstandardizedPredicted - mean(unstandardizedPredicted)) / sd(unstandardizedPredicted)
standardizedResiduals <- (unstandardizedResiduals - mean(unstandardizedResiduals)) / sd(unstandardizedResiduals)

#create standardized residuals plot
plot(standardizedPredicted, standardizedResiduals, main = "Standardized Residuals Plot", xlab = "Standardized Predicted Values", ylab = "Standardized Residuals")

#add horizontal line
abline(0,0)

#RESIDUALS HISTOGRAM
#create residuals histogram
hist(standardizedResiduals, freq = FALSE)

#add normal curve
curve(dnorm, add = TRUE)

#PP PLOT
#get probability distribution for residuals
probDist <- pnorm(standardizedResiduals)

#create PP plot
plot(ppoints(length(standardizedResiduals)), sort(probDist), main = "PP Plot", xlab = "Observed Probability", ylab = "Expected Probability")

#add diagonal line
abline(0,1)



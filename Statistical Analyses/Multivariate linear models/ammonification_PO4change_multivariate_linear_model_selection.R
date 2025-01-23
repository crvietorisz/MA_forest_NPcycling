# MA_forest_NPcycling 
# Multivariate linear model selection for models explaining net ammonification and net phosphate change
# 9/13/24 C. Vietorisz

library(lme4)
library(MuMIn)
library(leaps)

setwd("MA_forest_NPcycling/")

# read in full data
sub_dist <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)

##############################################################################################################################
# NET AMMONIFICATION
##############################################################################################################################


Ammon_sig_pred_list <- c('Ammonification', 
                          'edaphic_PC1',
                          'Shrubs',
                          'num_stems',
                          'root_density',
                          'ECM_abundance',
                          'Oligotroph_abundance',
                          'chitinase_ECbac_abund',
                          'bac_amm_pos_module_abundance',
                          'bac_shannon',
                          'fun_richness',
                          'bac_copy_number',
                          'fun_amm_pos_module_abundance',
                          'EMF_Shannon',
                          'glycosidaseOS_ECbac_abund',
                          'all_oligo_Ndecomp_ECbac',
                          'glycosidaseOS_ECfun_abund',
                          'all_ECM_Ndecomp_ECfun')

# format for input into leaps
Ammon_sig_preds <- sub_dist[,c(Ammon_sig_pred_list)]
Ammon_sig_preds <- Ammon_sig_preds [!is.na(Ammon_sig_preds$Ammonification),] 
Ammon_sig_preds_leaps <- na.omit(Ammon_sig_preds)
Ammon_sig_leaps_names <- colnames(Ammon_sig_preds_leaps[2:ncol(Ammon_sig_preds_leaps)]) #change this range to the size of your input dataframe, minus your response variable

# Use the leaps package to test all subsets of our data using Mallow's Cp as our selection criteria
library(leaps)
Ammon_sig_leaps <- leaps(x=Ammon_sig_preds_leaps[c(2:ncol(Ammon_sig_preds_leaps))], y=sqrt(Ammon_sig_preds_leaps$Ammonification),method="Cp")

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

# from plot, Cp=3 is the lowest number of parameters that is close to 1:1 line
cp3 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==3,])

#### CP3 ####
result <- lapply(1:nrow(cp3), function(i) {
  row_data <- cp3[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 4)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:3)
  new_df
})

# create empty matrix to hold AIC values
amm_cp3_aic_matrix <- matrix(nrow = 10, ncol =2) #change nrow to length of cp matrix
colnames(amm_cp3_aic_matrix) <- c("cp3 AIC", "Cp-p")
# merge with Cp data and predictor combos
amm_cp3_aic_matrix[,2] <- cp3$`Cp-p`
amm_cp3_aic_df <- as.data.frame(amm_cp3_aic_matrix)

for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("sqrt(Ammonification)", "~", paste(row, collapse = " + "), "+ (1|Site)"))
  model <- lmer(formula, data = sub_dist)
  amm_cp3_aic_df[x,1] <- AIC(model)
}
amm_cp3_aic_df <- merge(amm_cp3_aic_df, cp3, by = "Cp-p")


##### Build the lowest AIC model that has a Cp score under 2 ########

Ammon.lmer.1 <- lmer(sqrt(Ammonification) ~ edaphic_PC1 + ECM_abundance + (1|Site), data = sub_dist)
summary(Ammon.lmer.1)
MuMIn::r.squaredGLMM(Ammon.lmer.1) #R2m = 0.38, R2c = 0.42
AIC(Ammon.lmer.1) # 328.8
anova(Ammon.lmer.1)
confint(Ammon.lmer.1, 'edaphic_PC1', level = 0.95)
confint(Ammon.lmer.1, 'ECM_abundance', level = 0.95)

Ammon.lmer.2 <- lmer(sqrt(Ammonification) ~ edaphic_PC1 + fun_amm_pos_module_abundance + (1|Site), data = sub_dist)
summary(Ammon.lmer.2)
MuMIn::r.squaredGLMM(Ammon.lmer.2) #R2m = 0.38, R2c = 0.42
AIC(Ammon.lmer.2) # 328.3
anova(Ammon.lmer.2)
confint(Ammon.lmer.2, 'edaphic_PC1', level = 0.95)
confint(Ammon.lmer.2, 'fun_amm_pos_module_abundance', level = 0.95)


######## assess linear independence of predictors
plot(sub_dist$edaphic_PC1, sub_dist$ECM_abundance) #not highly correlated
plot(sub_dist$edaphic_PC1, sub_dist$fun_amm_pos_module_abundance) #not highly correlated


##############################################################################################################################
# NET PHOSPHATE CHANGE
##############################################################################################################################

# select out input predictors for model selection
Pmin_pred_list <- c("PO4_release", "no_pine_saplings", "ITS_copy_number", "soil_temp2022", "CPhydrolase_ECbac_abund",
                    "Oligotroph_Phydrolase_ECbac", "Poxidoreductase_ECfun_abund")


# format for input into leaps
Pmin_preds <- sub_dist[,c(Pmin_pred_list)]
Pmin_preds <- Pmin_preds [!is.na(Pmin_preds$PO4_release),] 
Pmin_preds_leaps <- na.omit(Pmin_preds)
Pmin_leaps_names <- colnames(Pmin_preds_leaps[2:ncol(Pmin_preds_leaps)])

# Use the leaps package to test all subsets of our data using Mallow's Cp as our selection criteria
library(leaps)
Pmin_leaps <- leaps(x=Pmin_preds_leaps[2:ncol(Pmin_preds_leaps)], y=Pmin_preds_leaps$PO4_release,method="Cp")

#plot all possible models by the number of parameters vs. the Cp score
plot(Pmin_leaps$size,Pmin_leaps$Cp,
     xlab="Number of Parameters",
     ylab="Cp",pch=16,col="blue")
# make it a bit easier to visualize by limiting y axis scale
plot(Pmin_leaps$size,Pmin_leaps$Cp,
     xlab="Number of Parameters",
     ylab="Cp",
     pch=16,col="blue",
     ylim=c(0,50))
abline(0,1)

# create a dataframe that we can use to interpret results
cp.dat <- cbind(Pmin_leaps$size,
                abs(Pmin_leaps$Cp-(Pmin_leaps$size)),
                Pmin_leaps$which)                           
rownames(cp.dat) <- 1:length(Pmin_leaps$size)                 # add rownames
colnames(cp.dat)[1:2] <- c("N.p","Cp-p")                      # add colnames for p, Cp-p-1
colnames(cp.dat)[3:ncol(cp.dat)] <- Pmin_leaps_names

# from plot, Cp=5 or 6 is the lowest number of parameters that is close to 1:1 line
cp5 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==5,])
cp6 <- as.data.frame(cp.dat[cp.dat[,"N.p"]==6,]) 

#### CP5 ####
result <- lapply(1:nrow(cp5), function(i) {
  row_data <- cp5[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 4)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:4)
  new_df
})

# create empty matrix to hold AIC values
pmin_cp5_aic_matrix <- matrix(nrow = 10, ncol =2) #change nrow to length of cp matrix
colnames(pmin_cp5_aic_matrix) <- c("cp5 AIC", "Cp-p")
# merge with Cp data and predictor combos
pmin_cp5_aic_matrix[,2] <- cp5$`Cp-p`
pmin_cp5_aic_df <- as.data.frame(pmin_cp5_aic_matrix)

for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("PO4_release", "~", paste(row, collapse = " + ")))
  model <- lm(formula, data = sub_dist)
  pmin_cp5_aic_df[x,1] <- AIC(model)
}
pmin_cp5_aic_df <- merge(pmin_cp5_aic_df, cp5, by = "Cp-p")

#### CP6 ####
result <- lapply(1:nrow(cp6), function(i) {
  row_data <- cp6[i, ]
  col_names <- names(row_data)[row_data == 1]
  new_df <- data.frame(matrix(col_names, nrow = 1, ncol = 5)) #change ncol to the number of predictors included in each model
  names(new_df) <- paste0("Column", 1:5) # change to ncol
  new_df
})

# create empty matrix to hold AIC values
pmin_cp6_aic_matrix <- matrix(nrow = 6, ncol =2)
colnames(pmin_cp6_aic_matrix) <- c("CP6 AIC", "Cp-p")
# merge with Cp data and predictor combos
pmin_cp6_aic_matrix[,2] <- cp6$`Cp-p`
pmin_cp6_aic_df <- as.data.frame(pmin_cp6_aic_matrix)

for (x in 1:length(result)) {
  result_df <- as.data.frame(result[x])
  row <- as.list(result_df[1,])
  formula <- as.formula(paste("PO4_release", "~", paste(row, collapse = " + ")))
  model <- lm(formula, data = sub_dist)
  pmin_cp6_aic_df[x,1] <- AIC(model)
}
pmin_cp6_aic_df <- merge(pmin_cp6_aic_df, cp6, by = "Cp-p")

##### Build the lowest AIC model that has a Cp score under 2 ########

Pmin.lm1 <- lm(PO4_release ~ no_pine_saplings + soil_temp2022 + CPhydrolase_ECbac_abund + Poxidoreductase_ECfun_abund, data = sub_dist)
summary(Pmin.lm1)
AIC(Pmin.lm1) #-330
anova(Pmin.lm1)
confint(Pmin.lm1, 'no_pine_saplings', level = 0.95)
confint(Pmin.lm1, 'soil_temp2022', level = 0.95)
confint(Pmin.lm1, 'CPhydrolase_ECbac_abund', level = 0.95)
confint(Pmin.lm1, 'Poxidoreductase_ECfun_abund', level = 0.95)


######## assess linear independence of predictors
plot(sub_dist$soil_temp2022, sub_dist$no_pine_saplings) #not highly correlated
plot(sub_dist$soil_temp2022, sub_dist$CPhydrolase_ECbac_abund) #not highly correlated
plot(sub_dist$soil_temp2022, sub_dist$Poxidoreductase_ECfun_abund) #slightly correlated? but not highly
plot(sub_dist$no_pine_saplings, sub_dist$CPhydrolase_ECbac_abund) #not highly correlated
plot(sub_dist$no_pine_saplings, sub_dist$Poxidoreductase_ECfun_abund) #not highly correlated
plot(sub_dist$CPhydrolase_ECbac_abund, sub_dist$Poxidoreductase_ECfun_abund) #not highly correlated


##############################################################################################################################
# TEST LINEAR MODEL ASSUMPTIONS
##############################################################################################################################

#RESIDUALS PLOT

mod <- Ammon.lmer.1 # change to your model

#get unstandardized predicted and residual values
unstandardizedPredicted <- predict(mod)
unstandardizedResiduals <- resid(mod)

#get standardized values
standardizedPredicted <- (unstandardizedPredicted - mean(unstandardizedPredicted)) / sd(unstandardizedPredicted)
standardizedResiduals <- (unstandardizedResiduals - mean(unstandardizedResiduals)) / sd(unstandardizedResiduals)

#create standardized residuals plot
plot(standardizedPredicted, standardizedResiduals, main = "Standardized Residuals Plot", xlab = "Standardized Predicted Values", ylab = "Standardized Residuals"); text(standardizedPredicted, standardizedResiduals, SampleID, cex = 0.5)

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





# White pine ch1
# Plots and stats for microbial metrics data
# C Vietorisz 8/18/24

library(lme4)
library(MuMIn)
library(ggplot2)
library(tidyverse)
library(DHARMa)
library(lmerTest)

setwd("MA_forest_NPcycling/")

# read in full data
sub_dist <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)

########################## AMMONIFICATION #########################################################
###################################################################################################

############ EMF relative abundance ##########

# Figure
EMFab.plot <- ggplot(sub_dist, aes(x = ECM_abundance, y = sqrt(Ammonification))) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, color = "black", lwd = 1) +
  labs(x = "EMF relative abundance", y = "Net ammonification")+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=18))
EMFab.plot

sum(!is.na(sub_dist$Ammonification) & !is.na(sub_dist$ECM_abundance))

# Stats
EMF.amm <- lmer(sqrt(Ammonification) ~ ECM_abundance + (1|Site), data = sub_dist)
summary(EMF.amm)
anova(EMF.amm, ddf = "Satterthwaite") # p = 0.0004, F = 13.3, numDF = 1, denDF = 109
MuMIn::r.squaredGLMM(EMF.amm) # R2m = 0.11

############ fungal high ammonification ICM relative abundance ##########

# Figure
fun.high.amm.plot <- ggplot(sub_dist, aes(x = fun_amm_pos_module_abundance, y = sqrt(Ammonification))) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, color = "black", lwd = 1) +
  labs(x = "Relative abundance of fungal high ammonification ICM", y = "Net ammonification")+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=18))
fun.high.amm.plot 

# Stats
funICM.amm <- lmer(sqrt(Ammonification) ~ fun_amm_pos_module_abundance + (1|Site), data = sub_dist)
summary(funICM.amm)
anova(funICM.amm, ddf = "Satterthwaite") # p = 0.0006, F = 12.5, numDF = 1, denDF = 109
MuMIn::r.squaredGLMM(funICM.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

############ Bacterial N-decomposition genes relative abundance ##########

### Stats

# All bacterial chitinases
chitinaseEC.amm <- lmer(sqrt(Ammonification) ~ sqrt(chitinase_ECbac_abund) + (1|Site), data = sub_dist)
summary(chitinaseEC.amm)
anova(chitinaseEC.amm, ddf = "Satterthwaite") # p = 7e-07, F = 33.9, numDF = 1, denDF = 43
MuMIn::r.squaredGLMM(chitinaseEC.amm) #R2m = 0.24

# All bacterial glycosidases (hydrolysing O- and S- glycosyl compounds)
glycosOS.EC.amm <- lmer(sqrt(Ammonification) ~ sqrt(glycosidaseOS_ECbac_abund) + (1|Site), data = sub_dist)
summary(glycosOS.EC.amm)
anova(glycosOS.EC.amm, ddf = "Satterthwaite") # p = 3e-06, F = 31.4, numDF = 1, denDF = 35
MuMIn::r.squaredGLMM(glycosOS.EC.amm) #R2m = 0.22


# Oligotrophic bacterial chitinases
olig.chitinaseEC.amm <- lmer(sqrt(Ammonification) ~ sqrt(Oligotroph_chitinase_ECbac) + (1|Site), data = sub_dist)
summary(olig.chitinaseEC.amm)
anova(olig.chitinaseEC.amm, ddf = "Satterthwaite") # p = 4e-05, F = 17.7, numDF = 1, denDF = 117.5
MuMIn::r.squaredGLMM(olig.chitinaseEC.amm) #R2m = 0.13

# Oligotrophic bacterial glycosidases (hydrolysing O- and S- glycosyl compounds)
olig.glycosOS.EC.amm <- lmer(sqrt(Ammonification) ~ sqrt(Oligotroph_glycosidaseOS_ECbac) + (1|Site), data = sub_dist)
summary(olig.glycosOS.EC.amm)
anova(olig.glycosOS.EC.amm, ddf = "Satterthwaite") # p = 0.0002, F = 14.7, numDF = 1, denDF = 115
MuMIn::r.squaredGLMM(olig.glycosOS.EC.amm) #R2m = 0.11

### Figure
glycosOS_chitinase.ECbac.amm.olig.plot <- ggplot(sub_dist) +
  geom_jitter(aes(Ndecomp_EC_abund, sqrt(Ammonification), color = Oligotroph_Ndecomp_EC)) +
  geom_smooth(aes(Ndecomp_EC_abund,sqrt(Ammonification)), method=lm, se=FALSE, color = "black", lwd = 0.5, linetype = "longdash") +
  geom_jitter(aes(sqrt(glycosidaseOS_EC_abund), sqrt(Ammonification), color = Oligotroph_glycosidaseOS_EC), shape = 17) +
  scale_color_gradient(low = "lightblue1", high = "deepskyblue4") +
  geom_smooth(aes(sqrt(glycosidaseOS_EC_abund), sqrt(Ammonification)), method=lm, se=FALSE, color = "black", lwd = 0.5) +
  labs(x = "Relative abundance of 
  bacterial N-decomposition genes", y = "sqrt(Net ammonification)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
glycosOS_chitinase.EC.amm.olig.plot   

########################## NITRIFICATION #########################################################
###################################################################################################

############ Copiotrophic bacteria relative abundance ##########

# Stats
copios.nitrif <- glm((Nitrification+1) ~ Copiotroph_abundance, # when running a glmer with site as an error term, the error term did not explain any variance
                         data = sub_dist,
                         family = Gamma(link = "log"))
summary(copios.nitrif) 
anova(copios.nitrif) # p = 9e-11, F = 50.3, numDF - 1, denDF = 124
MuMIn::r.squaredGLMM(copios.nitrif) # R2m (lognormal) = 0.32

# use the DHARMa package to test for over-dispersion
testDispersion(copios.nitrif) # over-dispersion present

# Figure
copios.nitrif.plot <- ggplot(sub_dist, aes(x = Copiotroph_abundance, y = log(Nitrification+1))) +
  geom_point() +
  geom_abline(slope=7.688, intercept=-1.74) +
  labs(x = "Copiotroph abundance", y = "log(Nitrification+1)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
copios.nitrif.plot               

############ Nitrifying bacteria relative abundance ##########

# Stats
nitrifiers.nitrif <- glm((Nitrification+1) ~ Nitrifier_abundance,
                             data = sub_dist,
                             family = Gamma(link = "log"))
summary(nitrifiers.nitrif) 
anova(nitrifiers.nitrif) # p = 0.002, F = 9.8, numDF = 1, denDF = 124
MuMIn::r.squaredGLMM(nitrifiers.nitrif) # R2m (lognormal) = 0.09

# use the DHARMa package to test for over-dispersion
testDispersion(nitrifiers.nitrif) # over-dispersion present

# Figure
nitrifiers.nitrif.plot <- ggplot(sub_dist, aes(x = Nitrifier_abundance, y = log(Nitrification+1))) +
  geom_point() +
  geom_abline(slope=42.07165, intercept=0.52033) +
  labs(x = "Nitrifier abundance", y = "log(Nitrification+1)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
nitrifiers.nitrif.plot       


############ Bacterial high nitrification ICM relative abundance ##########

# Stats
bacmodule.nitrif <- glm((Nitrification+1) ~ bac_nitr_pos_module_abundance,
                            data = sub_dist,
                        family = Gamma(link = "log"))
summary(bacmodule.nitrif)
anova(bacmodule.nitrif) # p = 2e-06, F = 25.1, DFnum = 1, DFden = 124
MuMIn::r.squaredGLMM(bacmodule.nitrif) # R2m (lognormal) = 0.15


# use the DHARMa package to test for over-dispersion
testDispersion(bacmodule.nitrif) # over-dispersion is present

# Figure
bacmodule.nitrif.plot <- ggplot(sub_dist, aes(x = bac_nitr_pos_module_abundance, y = log(Nitrification+1))) +
  geom_point() +
  geom_abline(slope=4.5707, intercept=-2.2608) +
  labs(x = "Bacterial nitrification module", y = "log(Nitrification+0.01)")+
  xlim(0.43, 0.71)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
bacmodule.nitrif.plot               

########################## PO4 RELEASE #########################################################
###################################################################################################

# Stats
Poxidored.pmin <- lmer(PO4_release ~ Poxidoreductase_ECfun_abund + (1|Site), data = sub_dist)
summary(Poxidored.pmin)
anova(Poxidored.pmin, ddf = "Satterthwaite") # p = 0.002, F = 10.3, numDF = 1, denDF = 96
MuMIn::r.squaredGLMM(Poxidored.pmin) #R2m = 0.08

# Figure
Poxidored.plot <- ggplot(sub_dist, aes(x = Poxidoreductase_ECfun_abund, y = PO4_release)) +
  geom_point() +
  geom_abline(slope=0.006978, intercept=-0.039836) +
  labs(x = "Relative abundance of fungal 
  P-cycling oxidoreductases", y = "PO4 change")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=18))
Poxidored.plot               

##############################################################################################################################
# TEST LINEAR MODEL ASSUMPTIONS
##############################################################################################################################

# use for models explaining net ammonification and net phosphate change (non-zero inflated models)

#RESIDUALS PLOT

mod <- Poxidored.pmin # change to the model you want

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


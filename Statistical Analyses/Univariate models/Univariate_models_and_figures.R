# White pine ch1
# Plots and stats for microbial metrics data
# C Vietorisz 8/18/24

library(lme4)
library(MuMIn)
library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(DHARMa)

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

# Stats
EMF.amm <- lmer(sqrt(Ammonification) ~ ECM_abundance + (1|Site), data = sub_dist)
summary(EMF.amm)
pt(q=3.642, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(EMF.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

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
pt(q=3.533, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(funICM.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

############ Bacterial N-decomposition genes relative abundance ##########

### Stats

# All bacterial chitinases
chitinaseEC.amm <- lmer(sqrt(Ammonification) ~ sqrt(chitinase_ECbac_abund) + (1|Site), data = sub_dist)
summary(chitinaseEC.amm)
pt(q=5.823, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(chitinaseEC.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

# All bacterial glycosidases (hydrolysing O- and S- glycosyl compounds)
glycosOS.EC.amm <- lmer(sqrt(Ammonification) ~ sqrt(glycosidaseOS_EC_abund) + (1|Site), data = sub_dist)
summary(glycosOS.EC.amm)
pt(q=5.60, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(glycosOS.EC.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model
# p = 6e-08, R2m = 0.23

# Oligotrophic bacterial chitinases
olig.chitinaseEC.amm <- lmer(sqrt(Ammonification) ~ sqrt(Oligotroph_chitinase_ECbac) + (1|Site), data = sub_dist)
summary(olig.chitinaseEC.amm)
pt(q=4.212, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(olig.chitinaseEC.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

# Oligotrophic bacterial glycosidases (hydrolysing O- and S- glycosyl compounds)
olig.glycosOS.EC.amm <- lmer(sqrt(Ammonification) ~ sqrt(Oligotroph_glycosidaseOS_ECbac) + (1|Site), data = sub_dist)
summary(olig.glycosOS.EC.amm)
pt(q=3.829, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(olig.glycosOS.EC.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

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
ggsave("GlycosidaseOS_chitinase_ECs_oligotroph_ammonification.png", path = "/Users/moniquegagnon/Desktop/BU/PhD/White pine/Figures/Nutr cycling vs. microbes/Gene abunds", width = 5.1, height = 3.5, dpi=300)

########################## NITRIFICATION #########################################################
###################################################################################################

############ Copiotrophic bacteria relative abundance ##########

# Stats
copios.nitrif <- glmmTMB((Nitrification+0.01) ~ Copiotroph_abundance + (1|Site),
                         data = sub_dist,
                         ziformula = ~ 1,
                         family = Gamma(link = "log"))
summary(copios.nitrif)
performance::r2(copios.nitrif) # p = 2e-12, R2m = 0.34

# use the DHARMa package to test if the data is zero inflated (to confirm that a zero-inflated model is appropriate)
mod.res<- simulateResiduals(copios.nitrif)
testZeroInflation(mod.res)

# Figure
copios.nitrif.plot <- ggplot(sub_dist, aes(x = Copiotroph_abundance, y = log(Nitrification+0.01))) +
  geom_point() +
  geom_abline(slope=20.6556, intercept=-6.7707) +
  labs(x = "Copiotroph abundance", y = "log(Nitrification+0.01)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
copios.nitrif.plot               

############ Nitrifying bacteria relative abundance ##########

# Stats
nitrifiers.nitrif <- glmmTMB((Nitrification+1) ~ Nitrifier_abundance + (1|Site),
                             data = sub_dist,
                             ziformula = ~ 1,
                             family = Gamma(link = "log"))
summary(nitrifiers.nitrif)
performance::r2(nitrifiers.nitrif) # p = 3e-4, R2m = 0.11

# use the DHARMa package to test if the data is zero inflated (to confirm that a zero-inflated model is appropriate)
mod.res<- simulateResiduals(nitrifiers.nitrif)
testZeroInflation(nitrifiers.nitrif)

# Figure
nitrifiers.nitrif.plot <- ggplot(sub_dist, aes(x = Nitrifier_abundance, y = log(Nitrification+1))) +
  geom_point() +
  geom_abline(slope=38.06525, intercept=0.52442) +
  scale_x_continuous(breaks=c(0, 0.01, 0.02), limits = c(0,0.021)) +
  scale_y_continuous(breaks=c(0, 0.5, 1, 2)) +
  labs(x = "Nitrifier abundance", y = "log(Nitrification+1)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
nitrifiers.nitrif.plot               


############ Bacterial high nitrification ICM relative abundance ##########

# Stats
bacmodule.nitrif <- glmmTMB((Nitrification+0.01) ~ bac_nitr_pos_module_abundance + (1|Site),
                            data = sub_dist,
                            ziformula = ~ 1,
                            family = Gamma(link = "log"))
summary(bacmodule.nitrif)
performance::r2(bacmodule.nitrif) # p = 1e-05, R2m = 0.08

# use the DHARMa package to test if the data is zero inflated (to confirm that a zero-inflated model is appropriate)
mod.res<- simulateResiduals(bacmodule.nitrif)
testZeroInflation(bacmodule.nitrif)

# Figure
bacmodule.nitrif.plot <- ggplot(sub_dist, aes(x = bac_nitr_pos_module_abundance, y = log(Nitrification+0.01))) +
  geom_point() +
  geom_abline(slope=7.685, intercept=-5.118) +
  labs(x = "Bacterial nitrification module", y = "log(Nitrification+0.01)")+
  xlim(0.43, 0.71)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=13),axis.title.y = element_text(size=13), axis.text=element_text(size=12))
bacmodule.nitrif.plot               

########################## PO4 RELEASE #########################################################
###################################################################################################

# Stats
Poxidored.amm <- lmer(PO4_release ~ Poxidoreductase_ECfun_abund + (1|Site), data = sub_dist)
summary(Poxidored.amm)
pt(q=3.207, df=length(na.omit(sub_dist$Ammonification)-2), lower.tail=FALSE) #get p value - MUST INPUT T STAT!
MuMIn::r.squaredGLMM(Poxidored.amm) #get Rsquared value - R2c is the conditional R2 and is variance explained by entire model

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

# use for models explaining net ammonification and net phosphate change

#RESIDUALS PLOT

mod <- Poxidored.amm # change to the model you want

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


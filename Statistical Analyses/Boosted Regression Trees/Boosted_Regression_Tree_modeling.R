# MA_forest_NPcycling
# Running Boosted Regression Tree models to explain nutrient cycling rates
# 7/9/24 C. Vietorisz

library(dismo)

setwd("MA_forest_NPcycling/")

# read in full data
sub_dist <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)

######################### AMMONIFICATION ###############################################################

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

amm_data <- sub_dist[Ammon_sig_pred_list]
amm_data <- amm_data[! is.na(amm_data$Ammonification),]
#turn ammonification into a positive integer, as required by the gbm.step function
amm_data$Ammonification <- round(amm_data$Ammonification*10, digits = 0)
hist(amm_data$Ammonification)

### Run model
amm.1 <- dismo::gbm.step(data=amm_data,
                                  gbm.x = c(2:18),
                                  gbm.y = 1,
                                  family = "poisson", # to view all family options, run: fix(gbm.step)
                                  tree.complexity = 2,
                                  learning.rate = 0.001,
                                  bag.fraction = 0.75,
                                  step.size = 50) 
amm.1$contributions

########## Run a for loop to run the model 1000 times

contributions <- amm.1$contributions
contributions_ord <- contributions[order(contributions$var),]

# create empty matrix to hold relative contribution and n.trees data
amm.boosted.contribs <- matrix(nrow=17, ncol=1000)
rownames(amm.boosted.contribs) <- rownames(contributions_ord)
amm.boosted.trees <- matrix(nrow=1000, ncol=1)
colnames(amm.boosted.trees) <- "Number of trees"

for (x in 1:1000) {
amm.boost.mod <- dismo::gbm.step(data=amm_sig3_data,
                             gbm.x = c(2:18),
                             gbm.y = 1,
                             family = "poisson", # to view all family options, run: fix(gbm.step)
                             tree.complexity = 2,
                             learning.rate = 0.001,
                             bag.fraction = 0.75,
                             step.size = 50) 
amm.contribs <- amm.boost.mod$contributions
amm.contribs <- amm.contribs[order(amm.contribs$var),]
amm.boosted.amm.contribs[,x] <- amm.contribs$rel.inf
amm.boosted.trees[x,] <- boost.mod$n.trees
}

#look at mean number of trees
mean(amm.boosted.trees[1])


############ Make figure #############

# read in variable categories
amm.var.types <- read.csv("WP21_ammon_variable_categories.csv")

# format dataframe for ggplot
library(tidyverse)

amm.boosted.contribs.t <- as.data.frame(t(amm.contribs))
amm.contribs.gg <- amm.boosted.contribs.t %>% 
  pivot_longer(col= everything(), 
               values_to = "contributions", 
               names_to = "variable")
amm.contribs.gg <- merge(amm.contribs.gg, amm.var.types, by = "variable")

### with axis labels
contribs.plot.amm = ggplot(amm.contribs.gg, aes(x=contributions, y=reorder(variable, contributions, FUN = median), fill = var_type, color = var_type))+
  geom_boxplot(outlier.size = 0.75, outlier.alpha = 0.5) +
  labs(x="Relative influence", y="Variable", fill = "Variable type")+
  scale_fill_manual(values=c("deepskyblue2", "sienna3", "green3")) +
  scale_color_manual(values=c("blue3", "sienna4", "darkgreen")) +
  guides(
    fill = guide_legend(override.aes = list(
      fill = c("deepskyblue2", "sienna3", "green3"),
      color = c("blue3", "sienna4", "darkgreen")
    )),
    color = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),legend.background = element_blank(), legend.key = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15), axis.text=element_text(size=12), legend.key.size = unit(2, "lines"), legend.title = element_text(size=15), legend.text = element_text(size=12))
contribs.plot.amm


######################### NITRIFICATION ###############################################################

hist(sub_dist$Nitrification, breaks = 50)

# select out variables that are p<0.05 significant with nitrification
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

nitr_data <- sub_dist[Nitr_pred_list]
nitr_data <- nitr_data[! is.na(nitr_data$Nitrification),]
#turn nitrification into a positive integer, as required by the gbm.step function
nitr_data$Nitrification <- round(nitr_data$Nitrification*10, digits = 0)
hist(nitr_data$Nitrification, breaks = 50)

## run model
nitr.1 <- dismo::gbm.step(data=nitr_data,
                         gbm.x = c(2:18),
                         gbm.y = 1,
                         family = "poisson", # to view all family options, run: ?gbm.step
                         tree.complexity = 2,
                         learning.rate = 0.001,
                         bag.fraction = 0.75,
                         step.size = 50) 
nitr.1$contributions

########## Run a for loop to run the model many times

contributions <- nitr.3$contributions
contributions_ord <- contributions[order(contributions$var),]

# create empty matrix to hold relative contribution and n.trees data
nitr.boosted.contribs <- matrix(nrow=17, ncol=1000)
rownames(nitr.boosted.contribs) <- rownames(contributions_ord)
nitr.boosted.trees <- matrix(nrow=1000, ncol=1)
colnames(nitr.boosted.trees) <- "Number of trees"

for (x in 1:1000) {
  boost.mod <- dismo::gbm.step(data=nitr_data,
                               gbm.x = c(2:18),
                               gbm.y = 1,
                               family = "poisson", # to view all family options, run: ?gbm.step
                               tree.complexity = 2,
                               learning.rate = 0.001,
                               bag.fraction = 0.75,
                               step.size = 50) 
  nitr.contribs <- boost.mod$contributions
  nitr.contribs <- nitr.contribs[order(nitr.contribs$var),]
  nitr.boosted.contribs[,x] <- nitr.contribs$rel.inf
  nitr.boosted.trees[x,] <- boost.mod$n.trees
}

#look at mean number of trees
mean(nitr.boosted.trees[1])


############### Make figures #############

# read in variable categories
nitr.var.types <- read.csv("WP21_nitrif_variable_categories.csv")

# format dataframe for ggplot
library(tidyverse)

nitr.boosted.contribs.t <- as.data.frame(t(nitr.contribs))
nitr.contribs.gg <- nitr.boosted.contribs.t %>% 
  pivot_longer(col= everything(), 
               values_to = "contributions", 
               names_to = "variable")
nitr.contribs.gg <- merge(nitr.contribs.gg, nitr.var.types, by = "variable")

contribs.plot.nitr = ggplot(nitr.contribs.gg, aes(x=contributions, y=reorder(variable, contributions, FUN = median), fill = var_type, color = var_type))+
  geom_boxplot(outlier.size = 0.75, outlier.alpha = 0.25) +
  labs(x="Relative influence", y="Variable", fill = "Variable type")+
  scale_fill_manual(values=c("deepskyblue2", "sienna3", "green3")) +
  scale_color_manual(values=c("blue3", "sienna4", "darkgreen")) +
  guides(
    fill = guide_legend(override.aes = list(
      fill = c("deepskyblue2", "sienna3", "green3"),
      color = c("blue3", "sienna4", "darkgreen")
    )),
    color = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),legend.background = element_blank(), legend.key = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15), axis.text=element_text(size=12), legend.key.size = unit(2, "lines"), legend.title = element_text(size=15), legend.text = element_text(size=12))
contribs.plot.nitr
ggsave("nitrif_boosted_reg_plot_1000reps_v2.png", path = "/Users/moniquegagnon/Desktop/BU/PhD/White pine/Figures/Boosted regressions", width=12, height=7, dpi=300)


######################### PO4 CHANGE ###############################################################

hist(sub_dist$PO4_release, breaks = 50)
#test if this is a laplace distribution - it is not!
library(goft)
laplace_test(na.omit(sub_dist$PO4_release), method = "transf", N = 10^5)


# select out variables that are univariate significant with PO4 change
Pmin_pred_list <- c("PO4_release", "no_pine_saplings", "ITS_copy_number", "soil_temp2022", "CPhydrolase_ECbac_abund",
                    "Oligotroph_Phydrolase_ECbac", "Poxidoreductase_ECfun_abund")


pmin_data <- sub_dist[Pmin_pred_list]
pmin_data <- pmin_data[! is.na(pmin_data$PO4_release),]

#turn PO4 release into a positive integer, as required by the gbm.step function
pmin_data$PO4_release <- pmin_data$PO4_release + 1
pmin_data$PO4_release <- round((pmin_data$PO4_release*100), digits = 0)
hist(pmin_data$PO4_release, breaks = 50)

# manually run the gbm.step function - set trees to minimum recommended of 1000 because otherwise gbm.step will run with fewer than 1000 trees
pmin.gbm.fix <- gbm.fixed(data=pmin_data, 
                              gbm.x = c(2:7), 
                              gbm.y = 1,
                              family = "gaussian",
                              learning.rate = 0.001, 
                              tree.complexity = 2,
                              bag.fraction = 0.75,
                              n.trees = 1000) 
pmin.gbm.fix$contributions


########## Run a for loop to run the model many times

contributions <- pmin.gbm.fix$contributions
contributions_ord <- contributions[order(contributions$var),]

# create empty matrix to hold relative contribution and n.trees data
pmin.boosted.contribs <- matrix(nrow=6, ncol=1000)
rownames(pmin.boosted.contribs) <- rownames(contributions_ord)
pmin.boosted.trees <- matrix(nrow=1000, ncol=1)
colnames(pmin.boosted.trees) <- "Number of trees"

for (x in 1:1000) {
  boost.mod <- dismo::gbm.fixed(data=pmin_data, 
                                gbm.x = c(2:7), 
                                gbm.y = 1,
                                family = "gaussian",
                                learning.rate = 0.001, 
                                tree.complexity = 2,
                                bag.fraction = 0.75,
                                n.trees = 1000) 
  pmin.contribs <- boost.mod$contributions
  pmin.contribs <- pmin.contribs[order(pmin.contribs$var),]
  pmin.boosted.contribs[,x] <- pmin.contribs$rel.inf
  pmin.boosted.trees[x,] <- boost.mod$n.trees
}

write.csv(pmin.boosted.contribs, "/Users/moniquegagnon/Desktop/BU/PhD/White pine/Amplicon chapter/Analyses/Boosted regression/WP21_Pmin_BoostedReg_1000reps_25Jul24.csv")

############### Make figures #############

pmin.boosted.contribs <- read.csv("/Users/moniquegagnon/Desktop/BU/PhD/White pine/Amplicon chapter/Analyses/Boosted regression/WP21_Pmin_BoostedReg_1000reps_25Jul24.csv", row.names=1)

# read in variable categories
pmin.var.types <- read.csv("WP21_pmin_var_categories.csv")

# format dataframe for ggplot
library(tidyverse)

pmin.boosted.contribs.t <- as.data.frame(t(pmin.boosted.contribs))
pmin.contribs.gg <- pmin.boosted.contribs.t %>% 
  pivot_longer(col= everything(), 
               values_to = "contributions", 
               names_to = "variable")
pmin.contribs.gg <- merge(pmin.contribs.gg, pmin.var.types, by = "variable")

contribs.plot.pmin = ggplot(pmin.contribs.gg, aes(x=contributions, y=reorder(variable, contributions, FUN = median), fill = var_type, color = var_type))+
  geom_boxplot(outlier.size = 0.75, outlier.alpha = 0.25) +
  labs(x="Relative influence", y="Variable", fill = "Variable type")+
  scale_fill_manual(values=c("deepskyblue2", "sienna3", "green3")) +
  scale_color_manual(values=c("blue3", "sienna4", "darkgreen")) +
  guides(
    fill = guide_legend(override.aes = list(
      fill = c("deepskyblue2", "sienna3", "green3"),
      color = c("blue3", "sienna4", "darkgreen")
    )),
    color = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),legend.background = element_blank(), legend.key = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15), axis.text=element_text(size=12), legend.key.size = unit(2, "lines"), legend.title = element_text(size=15), legend.text = element_text(size=12))
contribs.plot.pmin

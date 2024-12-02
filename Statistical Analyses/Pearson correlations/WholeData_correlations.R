# MA_forest_NPcycling 
# Run Pearson correlations between all variables in our dataset
# 11/25/24 C. Vietorisz

setwd("MA_forest_NPcycling/")

# read in full data
sub_dist <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)
# subset out just the continuous variables
contin <- sub_dist[5:186]

# run Pearson correlations
library(Hmisc)
cor_value <- rcorr(as.matrix(contin))$r # the correlation value
p_value <- signif(rcorr(as.matrix(contin))$P,2) # the P value for each correlation

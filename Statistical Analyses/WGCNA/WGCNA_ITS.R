# MA_forest_NPcycling 
# Running WGCNA to identify networks of ASVs associated with high or low nutrient cycling rates
# ITS data
# 6/8/23 C. Vietorisz

# Information on WGCNA package and workflow here:
# Langfelder and Horvath 2008: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559
# Package info: https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf

BiocManager::install("WGCNA")

library(tidyverse)
library(WGCNA)
library(flashClust)

setwd("MA_forest_NPcycling/")

################### Read in data and format ################################

#read in nutrient cycling data
fulldat <- read.csv("Data/MA_forest_NPcycling_fullData.csv", row.names=1)
nutr <- fulldat[c(6,7,40)]

#format data
rownames(nutr) <- nutr$ID

#read in ASV table
its_asv <- read.table("WP21_ITS_ASV_table_nodups_ratio_normalized.txt")

#format sample names to match nutrient data
rownames(its_asv) <- gsub("_", " ", rownames(its_asv))
rownames(its_asv) <- gsub(" ITS", "", rownames(its_asv))

#remove any rows in the nutrient cycling data that are not present in the ASV table
nutr <- merge(nutr, its_asv, by = 0)
nutr <- nutr[1:4]
rownames(nutr) <- nutr[,1]
nutr <- nutr[,c(2:4)]

#### Filter out rare ASVs ####

#create a plot showing what proportion of samples each ASV appears in
asv_pres <- as.data.frame(t(its_asv))
asv_pres[asv_pres>0] <- 1 #convert each cell to presence vs. absence
asv_pres$pres_sum <- rowSums(asv_pres[1:176]) #total number of samples the ASV is present in
hist(asv_pres$pres_sum, breaks = 100)

#Select out all ASVs present in more than 1 sample
asv_more1 <- asv_pres[asv_pres$pres_sum > 1 ,] # 2688 ASVs retained
over1.names <- rownames(asv_more1)
its_asv_1 <- t(its_asv)
its_asv_1 <- as.data.frame(its_asv_1[rownames(its_asv_1) %in% over1.names ,])
its_asv_1 <- as.data.frame(t(its_asv_1))

############################ Run WGCNA ############################

options(stringsAsFactors=FALSE)
allowWGCNAThreads()

datExpr0 = its_asv_1
datTraits = nutr
table(rownames(datTraits)==rownames(datExpr0)) #check that rownames are the same in both 

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0)) #type back to default, no direction in ITS counts
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.
# this calculates the whole network connectivity we choose signed (this is the default) because we are dealing with OTU counts, not gene expression
k=as.numeric(apply(A,2,sum))-1 #Summing columns of adjacency matrix (-1 to account for self correlation) 
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average") 
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
#there are no outliers! if there were, the "outlierC" band would highlight the outliers.

save(datExpr0,datTraits, file="WP21_Samples_Traits_ALL2_filter_less1_sample.RData")

######################### NETWORK CONSTRUCTION AND MODULE DETECTION ######################
################Moving on!  Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset
library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads use this in base R
allowWGCNAThreads() 
#lnames = load(file="WP21_Samples_Traits_ALL2_filter5perc.RData")

###################### SOFT POWER THRESHOLD #####################
#Figure out proper SFT
#**notes from tutorial: The function adjacency calculates the adjacency matrix from expression data.
#*Adjacency functions for both weighted and unweighted networks require the user to choose threshold parameters, for example by applying the approximate scale-free topology criterion 
#*general framework for soft thresholding that weighs each connection, suggested to look at scale-free topology index in WGCNA faq to decide soft threshold based on number of samples
#*dont want to include all possible data (not be too high) bc then its overfitting wgcna model
#*emphasizing network on stronger association/larger correlation coefficient and reduce the noise

# Choose a set of soft-thresholding powers 
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

#soft threshold powers are the power to which co-expression similarity is raised to calculate adjacency 
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)
#Printing table to help decide estimate for soft power threshold 
#Soft power threshold chosen based on SFT.R.sq being over .8 and mean k being below the hundreds 
#see info here!! https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#look at actual values for each power
sft

#a power of 7 gets us the closest to 0.9 without going over, so we choose 7
softPower=7
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.

############################### INITIAL DENDROGRAM ##############################
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average") 

sizeGrWindow(10,6)
pdf(file="WP21_dendrogram_thresh12_signed_filt_less1samp_26Jun2023.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()
#each leaf corresponds to an ASV, branches grouping together densely are interconnected, highly co-occurring ASVs

#*each module gets a number and a size
#*assign color to each module, color represents clusters of coexpressed genes

minModuleSize=4 #arbitrarily set the minimum number of ASVs in a module. 4 feels reasonable? play around with this...
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
#DeepSplit: For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced.

table(dynamicMods) #lists the number of ASVs in each module 

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#save Rdata up to this point
save.image(file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/premerge_WGCNA_filt_less1samp_mod4_26Jun2023.RData")

################ MERGING ##############################
#Merge modules whose expression profiles are very similar,
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 7) #*will need to change this to the soft threshold decided earlier
MEs= MEList$eigengenes

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs) 
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_filt_less1samp_mod4_WP21_26Jun2023.RData")

#lnames = load(file = "Network_filt_less1samp_mod4_WP21_26Jun2023")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.9 #*we can change the threshold of our height, merges them different based on value, dont want to overmerge things that arent that similar
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =1)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork_0.9.pdf", width=20, height=20) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save.image(file= "Rdata_signed_merge0.9_filt_less1samp_18Jul2023.RData")

#############################################################
#Relating modules to traits and finding important taxa
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
#lnames = load(file = "Network_invasion_nomerge.RData");
# Load network data saved in the second part.

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)
#*can pick which colors you want to further explore

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

####################### MODULE HEATMAP #####################################
#represent module trait correlations as a heatmap

# Will display correlations and their p-values
pdf(file="WP21_Module_heatmap_filt_less1samp_mod4_merge0.9_18Jul2023.pdf", width=20, height=60) 
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
#*when see groups of modules that are all hot or all cold, they should be merged. this will have it more likely detect enrichment
#* can see tight correlations, what percentage of variation is explained by certain relationships and look for modules that are doing the same thing

######################### RELATING MODULES TO TRAITS #######################

# Run the following code for only one trait at a time (eg, ammonification), then go directly to the module correlation plots! 

############## Ammonification #########################

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Ammonification); #change to your trait name
names(weight) = "Ammonification"

#remove rows where Ammonification is NA
ammon.na <- which(is.na(weight))
datExpr0.amm <- datExpr0[-c(ammon.na),]
MEs.amm <- MEs[-c(ammon.na),]

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0.amm, MEs.amm, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0.amm, na.omit(weight), use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

############## Nitrification #########################

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Nitrification); #change to your trait name
names(weight) = "Nitrification"

#remove rows where nitrification is NA
nitr.na <- which(is.na(weight))
datExpr0.nitr <- datExpr0[-c(nitr.na),]
MEs.nitr <- MEs[-c(nitr.na),]

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0.nitr, MEs.nitr, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0.nitr, na.omit(weight), use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

############## PO4 release #########################

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$PO4_release); #change to your trait name
names(weight) = "Phosphate release"

#remove rows where PO4 release is NA
Pmin.na <- which(is.na(weight))
datExpr0.Pmin <- datExpr0[-c(Pmin.na),]
MEs.Pmin <- MEs[-c(Pmin.na),]

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0.Pmin, MEs.Pmin, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0.Pmin, na.omit(weight), use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")


##################### MODULE CORRELATION PLOT ####################
#*how well things belong to module
#Gene-trait significance correlation plots
# a strong correlation here indicates that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait.

# par(mfrow=c(2,3))
module = "darkviolet" #*change this to the module we're going to look at
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem in", module, "module"),
                     ylab = "Gene Sig for Ammonification",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

  
##################### MODULE SPECIFIC FILES AND PLOTS #############################
####################################################################################

############### palevioletred3 ###########################################################

######## VSD FILES BY MODULE 
  #Making VSD files by module - so we can tell what genera are in each module
  
setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")
  
  vs=t(datExpr0)
  cands=names(datExpr0[moduleColors=="palevioletred3"]) #*change this also to the color of module we're looking at
  
  #*subsetting the genes in this module
  c.vsd=vs[rownames(vs) %in% cands,]
  nrow(c.vsd) #should correspond to module size
  table(moduleColors) #check module sizes here
  write.csv(c.vsd,"rlog_MMpalevioletred3_27jun2023.csv",quote=F)
  
########## KNOW WHICH TAXA ARE IN WHICH MODULE 

  ##########fisher of module vs whole dataset
  #*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
  #*kME is how well that gene belongs to the module 
  
  #Know which genera are in a given module
  library(WGCNA)
  options(stringsAsFactors=FALSE)
  
  vsd <- t(its_asv_1)
  allkME =as.data.frame(signedKME(its_asv_1, MEs))
  
  whichModule="palevioletred3" #*name your color and execute to the end
    
  length(moduleColors)
  inModule=data.frame("module"=rep(0,nrow(vsd)))
  row.names(inModule)=row.names(vsd)
  genes=row.names(vsd)[moduleColors == whichModule]
  inModule[genes,1]=1
  sum(inModule[,1])
  head(inModule)
  write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
  #*sum is sanity check, should be the same number that was in the module
  
  #know which ASVs are in each module
  palevioletred3 <- subset(inModule, module=="1")
  palevioletred3 <- row.names(palevioletred3)
  palevioletred3
  
# assign ASVs to taxa
asv_taxa <- read.table("/projectnb/talbot-lab-data/cviet/White_pine/soilDNA/WP21_ITS/DADA2_output/WP21_ITS_ASV_table.txt")
asv_taxa <- asv_taxa[193:199]
#replace the "x__" in front of each taxonomy
asv_taxa$Kingdom <- gsub("k__","",asv_taxa$Kingdom)
asv_taxa$Phylum <- gsub("p__","",asv_taxa$Phylum)
asv_taxa$Class <- gsub("c__","",asv_taxa$Class)
asv_taxa$Order <- gsub("o__","",asv_taxa$Order)
asv_taxa$Family <- gsub("f__","",asv_taxa$Family)
asv_taxa$Genus <- gsub("g__","",asv_taxa$Genus)
asv_taxa$Species <- gsub("s__","",asv_taxa$Species)
write.csv(asv_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/soilDNA/WP21_ITS/DADA2_output/WP21_ITS_ASV_taxonomy_list.csv")

#subset for only the ASVs in the module
palevioletred3_taxa <- asv_taxa[c(palevioletred3), ]
palevioletred3_taxa
write.csv(palevioletred3_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/palevioletred3_taxa.csv")


###################### kMEs 

  #*this gives kME and input for 
  #*series of how well gene belongs in module
  modColName=paste("kME",whichModule,sep="")
  modkME=as.data.frame(allkME[,modColName])
  row.names(modkME)=row.names(allkME)
  names(modkME)=modColName
  write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)
  
  #make a subset of the palevioletred3 kMEs for only the taxa in the palevioletred3 module
  palevioletred3.kmeinput<- paste(palevioletred3, sep=",")
  palevioletred3.kme<- subset(modkME, rownames(modkME) %in% palevioletred3.kmeinput)
  ASV <- rownames(palevioletred3.kme)
  rownames(palevioletred3.kme) <- NULL
  palevioletred3.kme <- cbind(ASV,palevioletred3.kme)
  
  #plot the kMEs to show which genera fit best into the module
palevioletred3_kme_plot <- ggplot(data= palevioletred3.kme, aes(x=reorder(ASV, -kMEpalevioletred3), y=kMEpalevioletred3)) + 
    geom_bar(color= "palevioletred3", fill="palevioletred3", stat="identity") +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x="ASV", y="kME")
palevioletred3_kme_plot
ggsave("palevioletred3_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### darkviolet ###########################################################

######## VSD FILES BY MODULE 
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="darkviolet"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMdarkviolet_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="darkviolet" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
darkviolet <- subset(inModule, module=="1")
darkviolet <- row.names(darkviolet)
darkviolet

#subset for only the ASVs in the module
darkviolet_taxa <- asv_taxa[c(darkviolet), ]
darkviolet_taxa
write.csv(darkviolet_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/darkviolet_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the darkviolet kMEs for only the taxa in the darkviolet module
darkviolet.kmeinput<- paste(darkviolet, sep=",")
darkviolet.kme<- subset(modkME, rownames(modkME) %in% darkviolet.kmeinput)
ASV <- rownames(darkviolet.kme)
rownames(darkviolet.kme) <- NULL
darkviolet.kme <- cbind(ASV,darkviolet.kme)

#plot the kMEs to show which genera fit best into the module
darkviolet_kme_plot <- ggplot(data= darkviolet.kme, aes(x=reorder(ASV, -kMEdarkviolet), y=kMEdarkviolet)) + 
  geom_bar(color= "darkviolet", fill="darkviolet", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
darkviolet_kme_plot
ggsave("darkviolet_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



############### plum3 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="plum3"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMplum3_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="plum3" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
plum3 <- subset(inModule, module=="1")
plum3 <- row.names(plum3)
plum3

#subset for only the ASVs in the module
plum3_taxa <- asv_taxa[c(plum3), ]
plum3_taxa
write.csv(plum3_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/plum3_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the plum3 kMEs for only the taxa in the plum3 module
plum3.kmeinput<- paste(plum3, sep=",")
plum3.kme<- subset(modkME, rownames(modkME) %in% plum3.kmeinput)
ASV <- rownames(plum3.kme)
rownames(plum3.kme) <- NULL
plum3.kme <- cbind(ASV,plum3.kme)

#plot the kMEs to show which genera fit best into the module
plum3_kme_plot <- ggplot(data= plum3.kme, aes(x=reorder(ASV, -kMEplum3), y=kMEplum3)) + 
  geom_bar(color= "plum3", fill="plum3", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
plum3_kme_plot
ggsave("plum3_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



############### darkseagreen3 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="darkseagreen3"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMdarkseagreen3_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="darkseagreen3" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
darkseagreen3 <- subset(inModule, module=="1")
darkseagreen3 <- row.names(darkseagreen3)
darkseagreen3

#subset for only the ASVs in the module
darkseagreen3_taxa <- asv_taxa[c(darkseagreen3), ]
darkseagreen3_taxa
write.csv(darkseagreen3_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/darkseagreen3_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the darkseagreen3 kMEs for only the taxa in the darkseagreen3 module
darkseagreen3.kmeinput<- paste(darkseagreen3, sep=",")
darkseagreen3.kme<- subset(modkME, rownames(modkME) %in% darkseagreen3.kmeinput)
ASV <- rownames(darkseagreen3.kme)
rownames(darkseagreen3.kme) <- NULL
darkseagreen3.kme <- cbind(ASV,darkseagreen3.kme)

#plot the kMEs to show which genera fit best into the module
darkseagreen3_kme_plot <- ggplot(data= darkseagreen3.kme, aes(x=reorder(ASV, -kMEdarkseagreen3), y=kMEdarkseagreen3)) + 
  geom_bar(color= "darkseagreen3", fill="darkseagreen3", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
darkseagreen3_kme_plot
ggsave("darkseagreen3_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)
# why do some ASVs have negative kME scores???



############### magenta4 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="magenta4"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMmagenta4_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="magenta4" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
magenta4 <- subset(inModule, module=="1")
magenta4 <- row.names(magenta4)
magenta4

#subset for only the ASVs in the module
magenta4_taxa <- asv_taxa[c(magenta4), ]
magenta4_taxa
write.csv(magenta4_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/magenta4_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the magenta4 kMEs for only the taxa in the magenta4 module
magenta4.kmeinput<- paste(magenta4, sep=",")
magenta4.kme<- subset(modkME, rownames(modkME) %in% magenta4.kmeinput)
ASV <- rownames(magenta4.kme)
rownames(magenta4.kme) <- NULL
magenta4.kme <- cbind(ASV,magenta4.kme)

#plot the kMEs to show which genera fit best into the module
magenta4_kme_plot <- ggplot(data= magenta4.kme, aes(x=reorder(ASV, -kMEmagenta4), y=kMEmagenta4)) + 
  geom_bar(color= "magenta4", fill="magenta4", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
magenta4_kme_plot
ggsave("magenta4_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### floralwhite ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="floralwhite"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMfloralwhite_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="floralwhite" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
floralwhite <- subset(inModule, module=="1")
floralwhite <- row.names(floralwhite)
floralwhite

#subset for only the ASVs in the module
floralwhite_taxa <- asv_taxa[c(floralwhite), ]
floralwhite_taxa
write.csv(floralwhite_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/floralwhite_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the floralwhite kMEs for only the taxa in the floralwhite module
floralwhite.kmeinput<- paste(floralwhite, sep=",")
floralwhite.kme<- subset(modkME, rownames(modkME) %in% floralwhite.kmeinput)
ASV <- rownames(floralwhite.kme)
rownames(floralwhite.kme) <- NULL
floralwhite.kme <- cbind(ASV,floralwhite.kme)

#plot the kMEs to show which genera fit best into the module
floralwhite_kme_plot <- ggplot(data= floralwhite.kme, aes(x=reorder(ASV, -kMEfloralwhite), y=kMEfloralwhite)) + 
  geom_bar(color= "floralwhite", fill="floralwhite", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
floralwhite_kme_plot
ggsave("floralwhite_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### greenyellow ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="greenyellow"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMgreenyellow_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="greenyellow" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
greenyellow <- subset(inModule, module=="1")
greenyellow <- row.names(greenyellow)
greenyellow

#subset for only the ASVs in the module
greenyellow_taxa <- asv_taxa[c(greenyellow), ]
greenyellow_taxa
write.csv(greenyellow_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/greenyellow_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the greenyellow kMEs for only the taxa in the greenyellow module
greenyellow.kmeinput<- paste(greenyellow, sep=",")
greenyellow.kme<- subset(modkME, rownames(modkME) %in% greenyellow.kmeinput)
ASV <- rownames(greenyellow.kme)
rownames(greenyellow.kme) <- NULL
greenyellow.kme <- cbind(ASV,greenyellow.kme)

#plot the kMEs to show which genera fit best into the module
greenyellow_kme_plot <- ggplot(data= greenyellow.kme, aes(x=reorder(ASV, -kMEgreenyellow), y=kMEgreenyellow)) + 
  geom_bar(color= "greenyellow", fill="greenyellow", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
greenyellow_kme_plot
ggsave("greenyellow_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### darkorange2 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="darkorange2"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMdarkorange2_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="darkorange2" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
darkorange2 <- subset(inModule, module=="1")
darkorange2 <- row.names(darkorange2)
darkorange2

#subset for only the ASVs in the module
darkorange2_taxa <- asv_taxa[c(darkorange2), ]
darkorange2_taxa
write.csv(darkorange2_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/darkorange2_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the darkorange2 kMEs for only the taxa in the darkorange2 module
darkorange2.kmeinput<- paste(darkorange2, sep=",")
darkorange2.kme<- subset(modkME, rownames(modkME) %in% darkorange2.kmeinput)
ASV <- rownames(darkorange2.kme)
rownames(darkorange2.kme) <- NULL
darkorange2.kme <- cbind(ASV,darkorange2.kme)

#plot the kMEs to show which genera fit best into the module
darkorange2_kme_plot <- ggplot(data= darkorange2.kme, aes(x=reorder(ASV, -kMEdarkorange2), y=kMEdarkorange2)) + 
  geom_bar(color= "darkorange2", fill="darkorange2", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
darkorange2_kme_plot
ggsave("darkorange2_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### grey60 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="grey60"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMgrey60_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="grey60" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
grey60 <- subset(inModule, module=="1")
grey60 <- row.names(grey60)
grey60

#subset for only the ASVs in the module
grey60_taxa <- asv_taxa[c(grey60), ]
grey60_taxa
write.csv(grey60_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/grey60_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the grey60 kMEs for only the taxa in the grey60 module
grey60.kmeinput<- paste(grey60, sep=",")
grey60.kme<- subset(modkME, rownames(modkME) %in% grey60.kmeinput)
ASV <- rownames(grey60.kme)
rownames(grey60.kme) <- NULL
grey60.kme <- cbind(ASV,grey60.kme)

#plot the kMEs to show which genera fit best into the module
grey60_kme_plot <- ggplot(data= grey60.kme, aes(x=reorder(ASV, -kMEgrey60), y=kMEgrey60)) + 
  geom_bar(color= "grey60", fill="grey60", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
grey60_kme_plot
ggsave("grey60_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### bisque4 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="bisque4"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMbisque4_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="bisque4" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
bisque4 <- subset(inModule, module=="1")
bisque4 <- row.names(bisque4)
bisque4

#subset for only the ASVs in the module
bisque4_taxa <- asv_taxa[c(bisque4), ]
bisque4_taxa
write.csv(bisque4_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/bisque4_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the bisque4 kMEs for only the taxa in the bisque4 module
bisque4.kmeinput<- paste(bisque4, sep=",")
bisque4.kme<- subset(modkME, rownames(modkME) %in% bisque4.kmeinput)
ASV <- rownames(bisque4.kme)
rownames(bisque4.kme) <- NULL
bisque4.kme <- cbind(ASV,bisque4.kme)

#plot the kMEs to show which genera fit best into the module
bisque4_kme_plot <- ggplot(data= bisque4.kme, aes(x=reorder(ASV, -kMEbisque4), y=kMEbisque4)) + 
  geom_bar(color= "bisque4", fill="bisque4", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
bisque4_kme_plot
ggsave("bisque4_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### yellowgreen ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="yellowgreen"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMyellowgreen_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="yellowgreen" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
yellowgreen <- subset(inModule, module=="1")
yellowgreen <- row.names(yellowgreen)
yellowgreen

#subset for only the ASVs in the module
yellowgreen_taxa <- asv_taxa[c(yellowgreen), ]
yellowgreen_taxa
write.csv(yellowgreen_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/yellowgreen_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the yellowgreen kMEs for only the taxa in the yellowgreen module
yellowgreen.kmeinput<- paste(yellowgreen, sep=",")
yellowgreen.kme<- subset(modkME, rownames(modkME) %in% yellowgreen.kmeinput)
ASV <- rownames(yellowgreen.kme)
rownames(yellowgreen.kme) <- NULL
yellowgreen.kme <- cbind(ASV,yellowgreen.kme)

#plot the kMEs to show which genera fit best into the module
yellowgreen_kme_plot <- ggplot(data= yellowgreen.kme, aes(x=reorder(ASV, -kMEyellowgreen), y=kMEyellowgreen)) + 
  geom_bar(color= "yellowgreen", fill="yellowgreen", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
yellowgreen_kme_plot
ggsave("yellowgreen_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### paleturquoise ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="paleturquoise"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMpaleturquoise_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="paleturquoise" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
paleturquoise <- subset(inModule, module=="1")
paleturquoise <- row.names(paleturquoise)
paleturquoise

#subset for only the ASVs in the module
paleturquoise_taxa <- asv_taxa[c(paleturquoise), ]
paleturquoise_taxa
write.csv(paleturquoise_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/paleturquoise_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the paleturquoise kMEs for only the taxa in the paleturquoise module
paleturquoise.kmeinput<- paste(paleturquoise, sep=",")
paleturquoise.kme<- subset(modkME, rownames(modkME) %in% paleturquoise.kmeinput)
ASV <- rownames(paleturquoise.kme)
rownames(paleturquoise.kme) <- NULL
paleturquoise.kme <- cbind(ASV,paleturquoise.kme)

#plot the kMEs to show which genera fit best into the module
paleturquoise_kme_plot <- ggplot(data= paleturquoise.kme, aes(x=reorder(ASV, -kMEpaleturquoise), y=kMEpaleturquoise)) + 
  geom_bar(color= "paleturquoise", fill="paleturquoise", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
paleturquoise_kme_plot
ggsave("paleturquoise_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### antiquewhite1 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="antiquewhite1"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMantiquewhite1_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="antiquewhite1" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
antiquewhite1 <- subset(inModule, module=="1")
antiquewhite1 <- row.names(antiquewhite1)
antiquewhite1

#subset for only the ASVs in the module
antiquewhite1_taxa <- asv_taxa[c(antiquewhite1), ]
antiquewhite1_taxa
write.csv(antiquewhite1_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/antiquewhite1_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the antiquewhite1 kMEs for only the taxa in the antiquewhite1 module
antiquewhite1.kmeinput<- paste(antiquewhite1, sep=",")
antiquewhite1.kme<- subset(modkME, rownames(modkME) %in% antiquewhite1.kmeinput)
ASV <- rownames(antiquewhite1.kme)
rownames(antiquewhite1.kme) <- NULL
antiquewhite1.kme <- cbind(ASV,antiquewhite1.kme)

#plot the kMEs to show which genera fit best into the module
antiquewhite1_kme_plot <- ggplot(data= antiquewhite1.kme, aes(x=reorder(ASV, -kMEantiquewhite1), y=kMEantiquewhite1)) + 
  geom_bar(color= "antiquewhite1", fill="antiquewhite1", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
antiquewhite1_kme_plot
ggsave("antiquewhite1_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



############### darkmagenta ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="darkmagenta"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMdarkmagenta_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="darkmagenta" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
darkmagenta <- subset(inModule, module=="1")
darkmagenta <- row.names(darkmagenta)
darkmagenta

#subset for only the ASVs in the module
darkmagenta_taxa <- asv_taxa[c(darkmagenta), ]
darkmagenta_taxa
write.csv(darkmagenta_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/darkmagenta_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the darkmagenta kMEs for only the taxa in the darkmagenta module
darkmagenta.kmeinput<- paste(darkmagenta, sep=",")
darkmagenta.kme<- subset(modkME, rownames(modkME) %in% darkmagenta.kmeinput)
ASV <- rownames(darkmagenta.kme)
rownames(darkmagenta.kme) <- NULL
darkmagenta.kme <- cbind(ASV,darkmagenta.kme)

#plot the kMEs to show which genera fit best into the module
darkmagenta_kme_plot <- ggplot(data= darkmagenta.kme, aes(x=reorder(ASV, -kMEdarkmagenta), y=kMEdarkmagenta)) + 
  geom_bar(color= "darkmagenta", fill="darkmagenta", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
darkmagenta_kme_plot
ggsave("darkmagenta_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### orangered4 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="orangered4"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMorangered4_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="orangered4" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
orangered4 <- subset(inModule, module=="1")
orangered4 <- row.names(orangered4)
orangered4

#subset for only the ASVs in the module
orangered4_taxa <- asv_taxa[c(orangered4), ]
orangered4_taxa
write.csv(orangered4_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/orangered4_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the orangered4 kMEs for only the taxa in the orangered4 module
orangered4.kmeinput<- paste(orangered4, sep=",")
orangered4.kme<- subset(modkME, rownames(modkME) %in% orangered4.kmeinput)
ASV <- rownames(orangered4.kme)
rownames(orangered4.kme) <- NULL
orangered4.kme <- cbind(ASV,orangered4.kme)

#plot the kMEs to show which genera fit best into the module
orangered4_kme_plot <- ggplot(data= orangered4.kme, aes(x=reorder(ASV, -kMEorangered4), y=kMEorangered4)) + 
  geom_bar(color= "orangered4", fill="orangered4", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
orangered4_kme_plot
ggsave("orangered4_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)


############### darkseagreen4 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="darkseagreen4"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMdarkseagreen4_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="darkseagreen4" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
darkseagreen4 <- subset(inModule, module=="1")
darkseagreen4 <- row.names(darkseagreen4)
darkseagreen4

#subset for only the ASVs in the module
darkseagreen4_taxa <- asv_taxa[c(darkseagreen4), ]
darkseagreen4_taxa
write.csv(darkseagreen4_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/darkseagreen4_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the darkseagreen4 kMEs for only the taxa in the darkseagreen4 module
darkseagreen4.kmeinput<- paste(darkseagreen4, sep=",")
darkseagreen4.kme<- subset(modkME, rownames(modkME) %in% darkseagreen4.kmeinput)
ASV <- rownames(darkseagreen4.kme)
rownames(darkseagreen4.kme) <- NULL
darkseagreen4.kme <- cbind(ASV,darkseagreen4.kme)

#plot the kMEs to show which genera fit best into the module
darkseagreen4_kme_plot <- ggplot(data= darkseagreen4.kme, aes(x=reorder(ASV, -kMEdarkseagreen4), y=kMEdarkseagreen4)) + 
  geom_bar(color= "darkseagreen4", fill="darkseagreen4", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
darkseagreen4_kme_plot
ggsave("darkseagreen4_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



############### red3 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="red3"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMred3_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="red3" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
red3 <- subset(inModule, module=="1")
red3 <- row.names(red3)
red3

#subset for only the ASVs in the module
red3_taxa <- asv_taxa[c(red3), ]
red3_taxa
write.csv(red3_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/red3_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the red3 kMEs for only the taxa in the red3 module
red3.kmeinput<- paste(red3, sep=",")
red3.kme<- subset(modkME, rownames(modkME) %in% red3.kmeinput)
ASV <- rownames(red3.kme)
rownames(red3.kme) <- NULL
red3.kme <- cbind(ASV,red3.kme)

#plot the kMEs to show which genera fit best into the module
red3_kme_plot <- ggplot(data= red3.kme, aes(x=reorder(ASV, -kMEred3), y=kMEred3)) + 
  geom_bar(color= "red3", fill="red3", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
red3_kme_plot
ggsave("red3_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



############### coral3 ###########################################################

######## VSD FILES BY MODULE
#Making VSD files by module - so we can tell what genera are in each module

setwd("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/")

vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="coral3"]) #*change this also to the color of module we're looking at

#*subsetting the genes in this module
c.vsd=vs[rownames(vs) %in% cands,]
nrow(c.vsd) #should correspond to module size
table(moduleColors) #check module sizes here
write.csv(c.vsd,"rlog_MMcoral3_27jun2023.csv",quote=F)

########## KNOW WHICH TAXA ARE IN WHICH MODULE 

##########fisher of module vs whole dataset
#*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
#*kME is how well that gene belongs to the module 

#Know which genera are in a given module
library(WGCNA)
options(stringsAsFactors=FALSE)

vsd <- t(its_asv_1)
allkME =as.data.frame(signedKME(its_asv_1, MEs))

whichModule="coral3" #*name your color and execute to the end
  
length(moduleColors)
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
genes=row.names(vsd)[moduleColors == whichModule]
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_merge0.9_fisher.csv",sep=""),quote=F)
#*sum is sanity check, should be the same number that was in the module

#know which ASVs are in each module
coral3 <- subset(inModule, module=="1")
coral3 <- row.names(coral3)
coral3

#subset for only the ASVs in the module
coral3_taxa <- asv_taxa[c(coral3), ]
coral3_taxa
write.csv(coral3_taxa, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_taxa/coral3_taxa.csv")


###################### kMEs 

#*this gives kME and input for 
#*series of how well gene belongs in module
modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

#make a subset of the coral3 kMEs for only the taxa in the coral3 module
coral3.kmeinput<- paste(coral3, sep=",")
coral3.kme<- subset(modkME, rownames(modkME) %in% coral3.kmeinput)
ASV <- rownames(coral3.kme)
rownames(coral3.kme) <- NULL
coral3.kme <- cbind(ASV,coral3.kme)

#plot the kMEs to show which genera fit best into the module
coral3_kme_plot <- ggplot(data= coral3.kme, aes(x=reorder(ASV, -kMEcoral3), y=kMEcoral3)) + 
  geom_bar(color= "coral3", fill="coral3", stat="identity") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="ASV", y="kME")
coral3_kme_plot
ggsave("coral3_kme_plot_merge0.9.png", path = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_plots/", width=6, height=4, dpi=300)



################## checking if most ASVs in one module are from the same sites #######################

its_asv_1_palevioletred3 <- its_asv_1[,c(palevioletred3)]
write.csv(its_asv_1_palevioletred3, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_palevioletred3.csv")  

its_asv_1_plum3 <- its_asv_1[,c(plum3)]
write.csv(its_asv_1_plum3, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_plum3.csv")  

its_asv_1_darkseagreen3 <- its_asv_1[,c(darkseagreen3)]
write.csv(its_asv_1_darkseagreen3, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_darkseagreen3.csv")  

its_asv_1_magenta4 <- its_asv_1[,c(magenta4)]
write.csv(its_asv_1_magenta4, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_magenta4.csv")  

its_asv_1_floralwhite <- its_asv_1[,c(floralwhite)]
write.csv(its_asv_1_floralwhite, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_floralwhite.csv")  

its_asv_1_greenyellow <- its_asv_1[,c(greenyellow)]
write.csv(its_asv_1_greenyellow, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_greenyellow.csv")  

its_asv_1_darkorange2 <- its_asv_1[,c(darkorange2)]
write.csv(its_asv_1_darkorange2, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_darkorange2.csv")  

its_asv_1_grey60 <- its_asv_1[,c(grey60)]
write.csv(its_asv_1_grey60, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_grey60.csv")  

its_asv_1_bisque4 <- its_asv_1[,c(bisque4)]
write.csv(its_asv_1_bisque4, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_bisque4.csv")  

its_asv_1_yellowgreen <- its_asv_1[,c(yellowgreen)]
write.csv(its_asv_1_yellowgreen, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_yellowgreen.csv")  

its_asv_1_paleturquoise <- its_asv_1[,c(paleturquoise)]
write.csv(its_asv_1_paleturquoise, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_paleturquoise.csv")  

its_asv_1_antiquewhite1 <- its_asv_1[,c(antiquewhite1)]
write.csv(its_asv_1_antiquewhite1, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_antiquewhite1.csv")  

its_asv_1_darkmagenta <- its_asv_1[,c(darkmagenta)]
write.csv(its_asv_1_darkmagenta, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_darkmagenta.csv")  

its_asv_1_coral3 <- its_asv_1[,c(coral3)]
write.csv(its_asv_1_coral3, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/its_asv_1_coral3.csv")  
#########################################################

mod_list <- c("palevioletred3", "darkviolet", "plum3", "darkseagreen3", "magenta4", "floralwhite", "greenyellow", "darkorange2", "grey60", "bisque4", "yellowgreen", "paleturquoise", "antiquewhite1", "darkmagenta", "yellowgreen", "orangered4", "darkseagreen4", "red3", "coral3")

# read in all kME files
for (x in 1:19) {
  filename <- paste(mod_list[x],".kME", sep="") 
  wd <- paste("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/", mod_list[x], "_kME.csv", sep="")
  assign(filename, read.csv(wd))
}

# initialize empty lists to hold new dataframes
kME_full_list <- list()
kME_subset_list <- list()

# Create a dataframe with all the ASVs by module and select for ASVs with kMEs >0.66
for (x in 1:19) {
wd <- paste("/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/module_files/", mod_list[x], "_kME.csv", sep="")
df <- read.csv(wd, col.names = c("ASV", "kME"))
df$module_name <- mod_list[x] #create a column with the module name
kME_full_list[[x]] <- df #put each new dataframe in a list
kME_subset_list[[x]] <- df[df[,2] > 0.66,] #subset each dataframe for ASVs with kMEs >0.66
}
kME_subset_df <-  as.data.frame(data.table::rbindlist(kME_subset_list, fill = T) ) #put all the subsetted dfs into one df
kME_full_df <-  as.data.frame(data.table::rbindlist(kME_full_list, fill = T) ) #put all the original dfs into one df

# Merge with taxonomy
asv_taxa2 <- asv_taxa
asv_taxa2$ASV <- rownames(asv_taxa2)
kME_taxa_66 <- merge(kME_subset_df, asv_taxa2, by = "ASV", all.x = TRUE)

write.csv(kME_taxa_66, file = "/projectnb/talbot-lab-data/cviet/White_pine/Analyses/WGCNA/Filter_less1_sample/key_taxa_WGCNA_ITS_merge0.9_kME0.66.csv")


kME_full_df[kME_full_df[,2] > .2,]

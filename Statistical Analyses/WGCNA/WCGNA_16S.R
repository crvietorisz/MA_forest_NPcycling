# MA_forest_NPcycling 
# Running WGCNA to identify networks of ASVs associated with high or low nutrient cycling rates
# 16S data
# 7/6/23 C. Vietorisz

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
bac_asv <- read.csv("WP21_16S_ASV_ratio_normalized.csv")
rownames(bac_asv) <- bac_asv$X
bac_asv <- bac_asv[2:192]
bac_asv <- as.data.frame(t(bac_asv))

#format sample names to match nutrient data
rownames(bac_asv) <- gsub("_", " ", rownames(bac_asv))
rownames(bac_asv) <- gsub(" 16S", "", rownames(bac_asv))

#remove any rows in the nutrient cycling data that are not present in the ASV table
nutr <- merge(nutr, bac_asv, by = 0)
nutr <- nutr[1:4]
rownames(nutr) <- nutr[,1]
nutr <- nutr[,c(2:4)]

#### Filter out rare ASVs ####

#create a plot showing what proportion of samples each ASV appears in
asv_pres <- as.data.frame(t(bac_asv))
asv_pres[asv_pres>0] <- 1 #convert each cell to presence vs. absence
asv_pres$pres_sum <- rowSums(asv_pres[1:191]) #total number of samples the ASV is present in
hist(asv_pres$pres_sum, breaks = 100)

#Select out all ASVs present in more than 1 sample
asv_more1 <- asv_pres[asv_pres$pres_sum > 1 ,] # 6459 ASVs retained

over1.names <- rownames(asv_more1)
bac_asv_1 <- t(bac_asv)
bac_asv_1 <- as.data.frame(bac_asv_1[rownames(bac_asv_1) %in% over1.names ,])
bac_asv_1 <- as.data.frame(t(bac_asv_1))

############################ Run WGCNA ############################

options(stringsAsFactors=FALSE)
allowWGCNAThreads()

datExpr0 = bac_asv_1
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

#a power of 7 gets us the closest to 0.9 without going over, so we choose 9
softPower=5
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#Calculates (correlation or distance) network adjacency from given expression data or from a similarity.

############################### INITIAL DENDROGRAM ##############################
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average") 

sizeGrWindow(10,6)
pdf(file="WP21_dendrogram_thresh12_signed_filt_less1samp_06jul2023.pdf", width=20, height=20)
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

################ MERGING ##############################
#Merge modules whose expression profiles are very similar, 
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 5) #*will need to change this to the soft threshold decided earlier
MEs= MEList$eigengenes

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs) 
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_filt_less1samp_mod4_WP21_06Jul2023.RData")

#lnames = load(file = "Network_filt_less1samp_mod4_WP21_06Jul2023")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.8 #*we can change the threshold of our height, merges them different based on value, dont want to overmerge things that arent that similar
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =1)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork_merge0.85.pdf", width=20, height=20) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save.image(file= "Rdata_signed_merge0.85_filt_less1samp_06Jul2023.RData")

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
pdf(file="WP21_Module_heatmap_filt_less1samp_mod4_merge0.85_06Jul2023.pdf", width=20, height=40) 
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
               cex.text = 1.5,
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
module = "green" #*change this to the module we're going to look at
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem in", module, "module"),
                     ylab = "Gene Sig for Trait",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
  
  ##################### MODULE SPECIFIC FILES AND PLOTS #############################
  ####################################################################################

  mod = "green"
  
  ######## VSD FILES BY MODULE 
  #Making VSD files by module - so we can tell what genera are in each module
  
  vs=t(datExpr0)
  cands=names(datExpr0[moduleColors==mod]) #*change this also to the color of module we're looking at
  
  #*subsetting the genes in this module
  c.vsd=vs[rownames(vs) %in% cands,]
  nrow(c.vsd) #should correspond to module size
  table(moduleColors) #check module sizes here
  write.csv(c.vsd, file=paste(mod,"_cvsd_rlog_06Jul2023.csv",sep=""),quote=F)
  
  ########## KNOW WHICH TAXA ARE IN WHICH MODULE 
  
  ##########fisher of module vs whole dataset
  #*fisher is binary value of 1  or 0, (1 if in module or 0 if its not)
  #*kME is how well that gene belongs to the module 
  
  #Know which genera are in a given module
  library(WGCNA)
  options(stringsAsFactors=FALSE)
  
  vsd <- t(bac_asv_1)
  allkME =as.data.frame(signedKME(bac_asv_1, MEs))
  
  whichModule=mod #*name your color and execute to the end
    
  length(moduleColors)
  inModule=data.frame("module"=rep(0,nrow(vsd)))
  row.names(inModule)=row.names(vsd)
  genes=row.names(vsd)[moduleColors == whichModule]
  inModule[genes,1]=1
  sum(inModule[,1])
  head(inModule)
  write.csv(inModule,file=paste(whichModule,"_merge0.8_fisher.csv",sep=""),quote=F)
  #*sum is sanity check, should be the same number that was in the module
  
  #know which ASVs are in each module
  mod_ASVs <- subset(inModule, module=="1")
  mod_ASVs <- row.names(mod_ASVs)
  mod_ASVs
  
  #subset for only the ASVs in the module
  mod_taxa <- asv_taxa[c(mod_ASVs), ]
  mod_taxa
  write.csv(mod_taxa, file = paste("module_taxa/",mod,"_taxa.csv", sep=""))
  
  
  ###################### kMEs 
  
  #*this gives kME and input for 
  #*series of how well gene belongs in module
  modColName=paste("kME",whichModule,sep="")
  modkME=as.data.frame(allkME[,modColName])
  row.names(modkME)=row.names(allkME)
  names(modkME)=modColName
  write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)
  
  #make a subset of this module's kMEs for only the taxa in this module
  mod.kmeinput<- paste(mod_ASVs, sep=",")
  mod.kme<- subset(modkME, rownames(modkME) %in% mod.kmeinput)
  ASV <- rownames(mod.kme)
  rownames(mod.kme) <- NULL
  mod.kme <- cbind(ASV,mod.kme)
  
  #plot the kMEs to show which genera fit best into the module
  mod_kme_plot <- ggplot(data= mod.kme, aes(x=reorder(ASV, -kMEgreen), y=kMEgreen)) + ### CHANGE ME TO MODULE NAME
    geom_bar(color= mod, fill=mod, stat="identity") +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x="ASV", y="kME")
  mod_kme_plot
  ggsave(paste(mod, "_kme_plot_merge0.8.png", sep=""), width=6, height=4, dpi=300)

#############################################################################################
###read in all module taxa and kMEs###
  
white1.kME <- read.csv("white1_kME.csv")
lightsteelblue.kME <- read.csv("lightsteelblue_kME.csv")
saddlebrown.kME <- read.csv("saddlebrown_kME.csv")
antiquewhite4.kME <- read.csv("antiquewhite4_kME.csv")
greenyellow.kME <- read.csv("greenyellow_kME.csv")
cyan.kME <- read.csv("cyan_kME.csv")
red4.kME <- read.csv("red4_kME.csv")
plum2.kME <- read.csv("plum2_kME.csv")
violet.kME <- read.csv("violet_kME.csv")
darkseagreen4.kME <- read.csv("darkseagreen4_kME.csv")
lightcyan1.kME <- read.csv("lightcyan1_kME.csv")
blue2.kME <- read.csv("blue2_kME.csv")
blue.kME <- read.csv("blue_kME.csv")
darkgreen.kME <- read.csv("darkgreen_kME.csv")
floralwhite.kME <- read.csv("blue_kME.csv")
bisque4.kME <- read.csv("floralwhite_kME.csv")
brown.kME <- read.csv("bisque4_kME.csv")
palevioletred2.kME <- read.csv("palevioletred2_kME.csv")
brown4.kME <- read.csv("brown4_kME.csv")
green.kME <- read.csv("green_kME.csv")

white1.taxa <- read.csv("module_taxa/white1_taxa.csv")
lightsteelblue.taxa <- read.csv("module_taxa/lightsteelblue_taxa.csv")
saddlebrown.taxa <- read.csv("module_taxa/saddlebrown_taxa.csv")
antiquewhite4.taxa <- read.csv("module_taxa/antiquewhite4_taxa.csv")
greenyellow.taxa <- read.csv("module_taxa/greenyellow_taxa.csv")
cyan.taxa <- read.csv("module_taxa/cyan_taxa.csv")
red4.taxa <- read.csv("module_taxa/red4_taxa.csv")
plum2.taxa <- read.csv("module_taxa/plum2_taxa.csv")
violet.taxa <- read.csv("module_taxa/violet_taxa.csv")
darkseagreen4.taxa <- read.csv("module_taxa/darkseagreen4_taxa.csv")
lightcyan1.taxa <- read.csv("module_taxa/lightcyan1_taxa.csv")
blue2.taxa <- read.csv("module_taxa/blue2_taxa.csv")
blue.taxa <- read.csv("module_taxa/blue_taxa.csv")
darkgreen.taxa <- read.csv("module_taxa/darkgreen_taxa.csv")
floralwhite.taxa <- read.csv("module_taxa/blue_taxa.csv")
bisque4.taxa <- read.csv("module_taxa/floralwhite_taxa.csv")
brown.taxa <- read.csv("module_taxa/bisque4_taxa.csv")
palevioletred2.taxa <- read.csv("module_taxa/palevioletred2_taxa.csv")
brown4.taxa <- read.csv("module_taxa/brown4_taxa.csv")
green.taxa <- read.csv("module_taxa/green_taxa.csv")

kME.list <- c(white1.kME, lightsteelblue.kME, saddlebrown.kME, antiquewhite4.kME, greenyellow.kME, cyan.kME, red4.kME, plum2.kME, violet.kME, darkseagreen4.kME, lightcyan1.kME, blue2.kME, darkgreen.kME, floralwhite.kME, bisque4.kME, brown.kME, palevioletred2.kME, brown4.kME, green.kME)

#subset for ASVs with a kME of > 0.66
white1.kME.66 <- white1.kME[white1.kME[2] > 0.66,]
lightsteelblue.kME.66 <- lightsteelblue.kME[lightsteelblue.kME[2] > 0.66 ,]
saddlebrown.kME.66 <- saddlebrown.kME[saddlebrown.kME[2] > 0.66 ,]
antiquewhite4.kME.66 <- antiquewhite4.kME[antiquewhite4.kME[2] > 0.66,]
greenyellow.kME.66 <- greenyellow.kME[greenyellow.kME[2] > 0.66,]
cyan.kME.66 <- cyan.kME[cyan.kME[2] > 0.66,]
red4.kME.66 <- red4.kME[red4.kME[2] > 0.66,]
plum2.kME.66 <- plum2.kME[plum2.kME[2] > 0.66,]
violet.kME.66 <- violet.kME[violet.kME[2] > 0.66,]
darkseagreen4.kME.66 <- darkseagreen4.kME[darkseagreen4.kME[2] > 0.66,]
lightcyan1.kME.66 <- lightcyan1.kME[lightcyan1.kME[2] > 0.66,]
blue2.kME.66 <- blue2.kME[blue2.kME[2] > 0.66,]
blue.kME.66 <- blue.kME[blue.kME[2] > 0.66,]
darkgreen.kME.66 <- darkgreen.kME[darkgreen.kME[2] > 0.66,]
floralwhite.kME.66 <- floralwhite.kME[floralwhite.kME[2] > 0.66,]
bisque4.kME.66 <- bisque4.kME[bisque4.kME[2] > 0.66,]
brown.kME.66 <- brown.kME[brown.kME[2] > 0.66,]
palevioletred2.kME.66 <- palevioletred2.kME[palevioletred2.kME[2] > 0.66,]
brown4.kME.66 <- brown4.kME[brown4.kME[2] > 0.66,]
green.kME.66 <- green.kME[green.kME[2] > 0.66,]

#get taxa names for ASVs with kME > 0.66
white1.taxa.66 <- white1.taxa[white1.taxa$X %in% white1.kME.66$X,]
lightsteelblue.taxa.66 <- lightsteelblue.taxa[lightsteelblue.taxa$X %in% lightsteelblue.kME.66$X,]
saddlebrown.taxa.66 <- saddlebrown.taxa[saddlebrown.taxa$X %in% saddlebrown.kME.66$X,]
antiquewhite4.taxa.66 <- antiquewhite4.taxa[antiquewhite4.taxa$X %in% antiquewhite4.kME.66$X,]
greenyellow.taxa.66 <- greenyellow.taxa[greenyellow.taxa$X %in% greenyellow.kME.66$X,]
cyan.taxa.66 <- cyan.taxa[cyan.taxa$X %in% cyan.kME.66$X,]
red4.taxa.66 <- red4.taxa[red4.taxa$X %in% red4.kME.66$X,]
plum2.taxa.66 <- plum2.taxa[plum2.taxa$X %in% plum2.kME.66$X,]
violet.taxa.66 <- violet.taxa[violet.taxa$X %in% violet.kME.66$X,]
darkseagreen4.taxa.66 <- darkseagreen4.taxa[darkseagreen4.taxa$X %in% darkseagreen4.kME.66$X,]
lightcyan1.taxa.66 <- lightcyan1.taxa[lightcyan1.taxa$X %in% lightcyan1.kME.66$X,]
blue2.taxa.66 <- blue2.taxa[blue2.taxa$X %in% blue2.kME.66$X,]
blue.taxa.66 <- blue.taxa[blue.taxa$X %in% blue.kME.66$X,]
darkgreen.taxa.66 <- darkgreen.taxa[darkgreen.taxa$X %in% darkgreen.kME.66$X,]
floralwhite.taxa.66 <- floralwhite.taxa[floralwhite.taxa$X %in% floralwhite.kME.66$X,]
bisque4.taxa.66 <- bisque4.taxa[bisque4.taxa$X %in% bisque4.kME.66$X,]
brown.taxa.66 <- brown.taxa[brown.taxa$X %in% brown.kME.66$X,]
palevioletred2.taxa.66 <- palevioletred2.taxa[palevioletred2.taxa$X %in% palevioletred2.kME.66$X,]
brown4.taxa.66 <- brown4.taxa[brown4.taxa$X %in% brown4.kME.66$X,]
green.taxa.66 <- green.taxa[green.taxa$X %in% green.kME.66$X,]

# add a column to each df indicating module
white1.taxa.66$module <- "white1"
lightsteelblue.taxa.66$module <- "lightsteelblue"
saddlebrown.taxa.66$module <- "saddlebrown"
antiquewhite4.taxa.66$module <- "antiquewhite4"
greenyellow.taxa.66$module <- "greenyellow"
cyan.taxa.66$module <- "cyan"
red4.taxa.66$module <- "ired4"
plum2.taxa.66$module <- "plum2"
violet.taxa.66$module <- "violet"
darkseagreen4.taxa.66$module <- "darkseagreen4"
lightcyan1.taxa.66$module <- "lightcyan1"
blue2.taxa.66$module <- "blue2"
blue.taxa.66$module <- "blue"
darkgreen.taxa.66$module <- "darkgreen"
floralwhite.taxa.66$module <- "floralwhite"
bisque4.taxa.66$module <- "bisque4"
brown.taxa.66$module <- "brown"
palevioletred2.taxa.66$module <- "palevioletred2"
brown4.taxa.66$module <- "brown4"
green.taxa.66$module <- "green"

#combine all into one dataframe
key_taxa <- rbind(white1.taxa.66, lightsteelblue.taxa.66, saddlebrown.taxa.66, antiquewhite4.taxa.66, blue.taxa.66, greenyellow.taxa.66, cyan.taxa.66, red4.taxa.66, plum2.taxa.66, violet.taxa.66, darkseagreen4.taxa.66, lightcyan1.taxa.66, blue2.taxa.66, darkgreen.taxa.66, floralwhite.taxa.66, bisque4.taxa.66, brown.taxa.66, palevioletred2.taxa.66, brown4.taxa.66, green.taxa.66)
colnames(key_taxa)[1] <- "ASV_no"

write.csv(key_taxa, file = "WP21_WGCNA_16S_key_taxa_merge0.8_06Jul23.csv")



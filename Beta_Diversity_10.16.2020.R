#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("phyloseq"))
library(phyloseq)
#install.packages("ggplot2")
library (ggplot2)
library (vegan)
library(reshape2)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("scales")
library(scales)
library(grid)
library(ape)
#install.packages("dplyr")
library(dplyr)
#install.packages("viridis")
library(viridis)
#install.packages("readxl")
library(readxl)
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(dplyr)
library(magrittr)
library(MASS)
library(ade4)
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
# Install
#install.packages("wesanderson")
# Load
library(wesanderson)
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
library(shades)
#install.packages("tidyverse")
#install.packages("tidygraph")
#install.packages("ggraph")
#install.packages("RAM")
#install.packages("FactoMineR")
library(FactoMineR)
#library(tidyverse)
#library(tidygraph)
#library(ggraph)
#library(RAM)

setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Denise_Data")
getwd()

##### Beta-diversity #####

## so there are A LOT of metrics of beta-diversity,  

## betadisper - made for comparing multivariate dispersion among group 
metadata=metadata[rownames(t.otu_counts),] ## reorder metadata to have same rows as relative abundance table * must have same # of rows
head(t.otu_counts)
head(phyla_RelAb)

otu_RA<-data.frame(decostand(t.otu_counts, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples

otu.RA.d<-vegdist(otu_RA, 'bray') ## create Bray-Curtis dissimilarity matrix with Relative Abundance of original OTU table

head(metadata)

ori<-as.factor(metadata$Origin)

ster<-as.factor(metadata$Sterility)

## ** what betadisper does -- multivariate dispersion! 
### aka multivariate homogeneity of group dispersions (variances)
# calculates average distance of group members to group centroid or spatial median in multivariate space

beta1<-betadisper(otu.RA.d, ori, type="centroid") # beta diversity across sites (origin as group)

boxplot(beta1, xlab="Origin")
## so what is our beta-diversity? (aka Average distance to centroid)
# lab = 0.4882, wild = 0.5422 
summary(beta1)

## beta diversity across sites (sterility as group)
beta2<-betadisper(otu.RA.d, ster, type="centroid")
beta2
## so what is our beta-diversity? (aka Average distance to centroid)
# natural = 0.5245, sterile - 0.5139 
boxplot(beta2, xlab="Sterility")

## we can quickly test if there are significant differences
# To test if one or more groups is more variable than the others, ANOVA of the distances to group centroids can be performed 
# and parametric theory used to interpret the significance of F

anova(beta1) # based on origin -- P val =  0.0831, F = 3.23
anova(beta2) # based on sterility -- P val = 0.7195, F = 0.1316
anova(beta1)

#You can also use Tukey's Honest Significant Differences when you have more than 2 groups
## to test what groups are different

mod.HSD = TukeyHSD(beta1)
mod.HSD

mod.HSD2 = TukeyHSD(beta2)
mod.HSD2

#### Checking out some color pallettes before we visualize ####
wes1<-wes_palette("Chevalier1")
wes2<-wes_palette("Moonrise3")
wes3<-wes_palette("IsleofDogs1")
wes4<-wes_palette("GrandBudapest1")
wes5<-wes_palette("GrandBudapest2")
#scale_fill_manual(values = wes_palette("IsleofDogs1"))

SM_pal <- park_palette("SmokyMountains") # create a palette and specify # of colors youw ant
Arc_pal <- park_palette("Arches") # create a palette and specify # of colors youw ant
CL_pal <- park_palette("CraterLake") # create a palette and specify # of colors youw ant
Sag_pal <- park_palette("Saguaro") # create a palette and specify # of colors youw ant
Aca_pal <- park_palette("Acadia") # create a palette and specify # of colors youw ant
DV_pal <- park_palette("DeathValley") # create a palette and specify # of colors youw ant
CI_pal <- park_palette("ChannelIslands") # create a palette and specify # of colors youw ant
Bad_pal <- park_palette("Badlands") # create a palette and specify # of colors youw ant
MR_pal <- park_palette("MtRainier") # create a palette and specify # of colors youw ant
HI_pal <- park_palette("Hawaii") # create a palette and specify # of colors youw ant


#### Beta Diversity Visualization ####

## beta diversity across sites (origin as group)
df1 <- data.frame(Distance_to_centroid=beta1$distances,Group=beta1$group)
# ^ from beta1<-betadisper(otu.RA.d, ori, type="centroid") - beta diversity w/ origin as group

bd.ori<-ggplot(data=df1,aes(x=Group,y=Distance_to_centroid, fill=Group)) + geom_boxplot(alpha=0.5, color="black") +
  labs(title = "Beta Diversity by Treatment Source", x="Treatment Source", y="Beta Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  scale_fill_manual(values = saturation(CL_pal, 1), name ="Treatment Source", labels=c("lab"="Lab", "wild"="Garden")) 

ggsave(bd.ori,filename = "16S_betadiv_trtmntSource_10.15.20.pdf", width=15, height=10, dpi=600)


## beta diversity across sites (sterility as group)

df2 <- data.frame(Distance_to_centroid=beta2$distances,Group=beta2$group)
## ^ from beta2<-betadisper(otu.RA.d, ster, type="centroid") - beta diversity w/ sterility as group


bd.ster<-ggplot(data=df2,aes(x=Group,y=Distance_to_centroid, fill=Group)) + geom_boxplot(alpha=0.5, color="black") +
  labs(title = "Beta Diversity by Treatment Condition", x="Treatment Condition", y="Beta Diversity") + theme_classic() + scale_x_discrete(name ="Treatment Condition", labels=c("N"="Natural", "S"="Sterile")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15)) + 
  scale_fill_manual(values = saturation(Bad_pal, 2), name ="Treatment Condition", labels=c("N"="Natural", "S"="Sterile")) 

ggsave(bd.ster,filename = "16S_betadiv_trtmntCondition_10.15.20.pdf", width=15, height=10, dpi=600)


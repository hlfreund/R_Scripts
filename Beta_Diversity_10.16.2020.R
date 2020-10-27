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

setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/SierraDust/New_Results_September2020")
getwd()


#### Import and reformat the data ####
otus_counts_taxa<- as.data.frame(read.csv("filtered_table_16S_191117_HF.2.csv"))

dim(otus_counts_taxa)
tail(otus_counts_taxa) ## make sure that OTUs are rows and samples are columns!!

otus_counts_taxa$Phylum <- sub("^$", "Unknown", otus_counts_taxa$Phylum) # changing all empty cells to "Unknown"
otus_counts_taxa$Class <- sub("^$", "Unknown", otus_counts_taxa$Class)
otus_counts_taxa$Order <- sub("^$", "Unknown", otus_counts_taxa$Order)
otus_counts_taxa$Family <- sub("^$", "Unknown", otus_counts_taxa$Family)
otus_counts_taxa$Genus <- sub("^$", "Unknown", otus_counts_taxa$Genus)
otus_counts_taxa$Species <- sub("^$", "Unknown", otus_counts_taxa$Species)
head(otus_counts_taxa)

otus_counts_taxa<-subset(otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaea!!!!!!!!!
otus_counts_taxa_all<-otus_counts_taxa # use this to compare Chloroplasts to all else -- total relative abundance, including Chloroplasts

otus_counts_taxa<-subset(otus_counts_taxa, Class!="Chloroplast") ## keep only bacteria!!!!!!!!!
#otus_counts_taxa[is.na(otus_counts_taxa)]<- "Unknown" # labeling NAs as unknowns instead 
archaea<-subset(otus_counts_taxa, Kingdom=="Archaea") ## keep only bacteria and archaea!!!!!!!!!


'Chloroplast' %in% otus_counts_taxa # check if Chloroplast counts are still in df, should be false because they've been removed
'Unknown' %in% otus_counts_taxa$Kingdom # check if Unknown kingdom counts are still in df, should be false because they've been removed
'Archaea' %in% otus_counts_taxa$Kingdom # check if Archaea kingdom counts are still in df, should be true 

# export OTU + taxa table that excludes chloroplasts counts -- for comparison
otus_counts_bac_melted<-melt(otus_counts_taxa) # using this to export out and compare other OTU tables
colnames(otus_counts_bac_melted)[which(names(otus_counts_bac_melted) == "variable")] <- "SampleID"
colnames(otus_counts_bac_melted)[which(names(otus_counts_bac_melted) == "value")] <- "Count"

otu_tax_ID<-subset(otus_counts_taxa, select=c(OTUID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(otu_tax_ID)
otu_counts<-subset(otus_counts_taxa, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(otu_counts) ## these are RAW otu counts!
rownames(otu_counts)<-otu_counts$OTUID
otu_counts$OTUID<-NULL

t.otu_counts<-data.frame(t(otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
head(t.otu_counts)
#t.otu_counts$SampleID<-rownames(t.otu_counts) 
head(t.otu_counts) # should see long list of OTU IDs and then "SampleID" as last column

### Import metadata ####
metadata<-as.data.frame(read_excel("SierraMapEle2SH.xls"))

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, SampType, TubeID, Description)) # subset only part of the metadata we need

rownames(metadata)<-metadata$SampleID
tail(metadata)

metadata=metadata[rownames(t.otu_counts),] ## reorder metadata to have same rows as relative abundance table * must have same # of rows
head(t.otu_counts)

#### Beta Diversity - Step by Step

#### Relative Abundance of Count Data ####
otu_RA<-data.frame(decostand(t.otu_counts, method="total", MARGIN=1, na.rm=TRUE))  
# relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples

otu.RA.d<-vegdist(otu_RA, 'bray') ## create Bray-Curtis dissimilarity matrix with Relative Abundance of original OTU table

head(metadata)

site<-as.factor(metadata$Site)

ele<-as.factor(metadata$Elevation)

yr<-as.factor(metadata$Year)

## ** what betadisper does -- multivariate dispersion! 
### made for comparing multivariate dispersion among group
### aka multivariate homogeneity of group dispersions (variances)
# calculates average distance of group members to group centroid or spatial median in multivariate space

beta1<-betadisper(otu.RA.d, ele, type="centroid") # beta diversity across sites (origin as group)

boxplot(beta1, xlab="Elevation")
## so what is our beta-diversity? (aka Average distance to centroid)
# lab = 0.4882, wild = 0.5422 
summary(beta1)

## beta diversity across sites (sterility as group)
beta2<-betadisper(otu.RA.d, site, type="centroid")
beta2
## so what is our beta-diversity? (aka Average distance to centroid)
# natural = 0.5245, sterile - 0.5139 
boxplot(beta2, xlab="Site")

## we can quickly test if there are significant differences
# To test if one or more groups is more variable than the others, ANOVA of the distances to group centroids can be performed 
# and parametric theory used to interpret the significance of F

anova(beta1) #
anova(beta2) # 
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


#### Notes on Beta Diversity with betadisper() ####
# reducing the distances produced using any dissimilarity coefficient to principal coordinates, which embeds them within a Euclidean space. 
# **** The analysis then proceeds by calculating the Euclidean distances between group members and the group centroid on the basis of the principal coordinate axes rather than the original distances. 
# Non-metric dissimilarity coefficients can produce principal coordinate axes that have negative Eigenvalues. 
# These correspond to the imaginary, non-metric part of the distance between objects. If negative Eigenvalues are produced, we must correct for these imaginary distances.

# To test if one or more groups is more variable than the others, ANOVA of the distances to group centroids can be performed and parametric theory used to interpret the significance of F

### Playing with NMDS ####

#  head(dist_f)

NMDS_ra = metaMDS(otu.RA.d) ##### non-metric multidimensional scaling of distance matrices to understand distribution of taxa across sites
## ^ distance matrix based on Bray-Curtis dissimilarities of relativized abundance values

class(NMDS_ra)

MDS1 = NMDS_ra$points[,1]
MDS1
MDS2 = NMDS_ra$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2)
NMDS$SampleID<-rownames(NMDS)

head(metadata)

chem_meta<-subset(metadata, select=c(SampleID, Cu, Fe, Mg, Mn, Ni, P, S, Zn)) # subset only part of the metadata we need
rownames(chem_meta)<-chem_meta$SampleID
chem_meta=chem_meta[rownames(t.otu_counts),] ## reorder metadata to have same rows as original OTU table * must have same # of rows
head(chem_meta)


chem_meta_scale<-scale(chem_meta)

class(met_name)

metadat_fit<-envfit(NMDS_ra, chem_meta_scale, permutations = 999, strata = NULL, 
       choices=c(1,2),labels=met_name, na.rm = FALSE)
# ^ envfit function fits environmental vectors or factors onto an ordination. 
# The projections of points onto vectors have maximum correlation with corresponding environmental variables, and the factors show the averages of factor levels.

class(metadat_fit)

plot(NMDS_ra)  # messing with ordination plots and adding environmental factors
plot(metadat_fit, col = "red")

dev.off()

NMDS_env<-merge(NMDS,metadata,by="SampleID")
NMDS_env

#### Plotting NMDS with ggplot ####

ggplot(NMDS_env, aes(x=MDS1, y=MDS2, col=ele)) +geom_point(size=3,alpha=0.8) +theme_bw()+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))+labs(title = "Playing", col="Elevation")

dev.off()  


#### Mantel test ####
# comparing standardized chemical data (Euclidean distance) + Bray-Curtis dissimilarity on fungal composition

head(chem_meta)
chem_meta$SampleID<-NULL
chem_dist<-vegdist(scale(chem_meta), "euclidean") ## may need to transform these data prior to this step b/c comparing elements at varying concentraitons

head(chem_dist)
head(otu.RA.d)

mantel(otu.RA.d, chem_dist, method="pearson", permutations=999, na.rm = FALSE)
mantel(otu.RA.d, chem_dist, method="spearman", permutations=999, na.rm = FALSE)


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

otus_counts_taxa<-subset(otus_counts_taxa, Kingdom!="Unknown") ## keep only bacteria and archaean -- drop Unknowns 

otus_counts_taxa<-subset(otus_counts_taxa, Class!="Chloroplast") ## keep only bacteria -- exclude Chloroplast sequences
'Chloroplast' %in% otus_counts_taxa # check if Chloroplast counts are still in df, should be false because they've been removed

otu_tax_ID<-subset(otus_counts_taxa, select=c(OTUID, Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(otu_tax_ID)
otu_counts<-subset(otus_counts_taxa, select=-c(Kingdom, Phylum, Class, Order, Family, Genus, Species)) # subset only part of the metadata we need
head(otu_counts) ## these are RAW otu counts!

rownames(otu_counts)<-otu_counts$OTUID
otu_counts$OTUID<-NULL 
t.otu_counts<-data.frame(t(otu_counts)) ## a work around to get a table of samples as rows and OTU as columns for when you melt the data
#t.otu_counts$SampleID<-rownames(t.otu_counts) 
head(t.otu_counts) # should see long list of OTU IDs

### **** USE t.otu_counts FOR ANALYSES -- rows are Samples, columns are OTU IDs
# in vegan: rows are samples or site IDs, and columns are OTUs/ASVs/Species

### Import metadata ####
metadata<-as.data.frame(read_excel("SierraMapEle2SH.xls"))

head(metadata)
tail(metadata)
metadata<-subset(metadata, select=-c(BarcodeSequence, LinkerPrimerSequence, SiteCode, RepNum, SiteRep, SampType, TubeID, Description)) # subset only part of the metadata we need

rownames(metadata)<-metadata$SampleID
tail(metadata)

metadata=metadata[rownames(t.otu_counts),] ## reorder metadata to have same rows as original OTU table
# ^ this indexing method will only work if the two dfs have the same # of rows AND the same row names!

#### Alpha Diversity + Species Richness ####
# ## H <- diversity() -- Shannon entropy function
## exp(H) -- Shannon diversity (exp of entropy)

Shan_ent.16s<-vegan::diversity(t.otu_counts, index="shannon") # Shannon entropy
Shan_div.16s<- exp(Shan_ent.16s) # Shannon Diversity aka Hill number 1
div_a<-data.frame(Shan_ent.16s,Shan_div.16s)
class(div_a)
div_a$SampleID<-rownames(div_a)

S<-specnumber(t.otu_counts) # finds # of species per sample using RAW count data; if MARGIN = 2 it finds frequencies of species
S # tells you number of OTUs in every sample
S.freq<-specnumber(t.otu_counts, MARGIN = 2) # # finds how many times each ASV appeared in samples (frequency)
S.freq # # of OTUs that appear across ALL sites/samples

#### More Alpha Diversity Code - from Marko's EEOB 230 class ####

## We can also describe the data with diversity measures, such as:
SR <- rowSums(t.otu_counts > 0) ## Species richness
SR
H <- diversity(t.otu_counts) ## Shannon entropy
H
Div1 <- exp(H) ## Shannon's diversity (number of abundant species)
Div1
Div2 <- diversity(t.otu_counts, "inv") ## Simpson diversity (number of dominant species)
Div2
Eve <- H/log(SR) ## Pielou evenness
Eve
E10 <- Div1/SR ## Shannon evenness (Hill's ratio)
E10
E20 <- Div2/SR ## Simpson evenness (Hill's ratio)
E20
div_b <- data.frame(SR, H, Div1, Div2, E10, E20, Eve) ## create a dataframe for the above measures
head(div_b)
dim(div_b)
## you can also use this for Richness
specnumber(t.otu_counts)

#### Merge Alpha div w/ count/taxa data + metadata ####
head(div_a)

alpha.div.metadat <- merge(div_a,metadata, by="SampleID")
head(alpha.div.metadat)
class(alpha.div.metadat)

#alpha.div.metadat$Month2 <- factor(alpha.div.metadat$Month, levels = c("July","August","October")) ## reeorder month column so that it will be listed in chronological order in x axis

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


#### Alpha Diversity Visualization ####
head(alpha.div.metadat)

# shannon entropy by year
alpha.ent.raw.yr<-ggplot(alpha.div.metadat, aes(x=factor(Year), y=Shan_ent.16s, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Sampling Year", x="Year", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.yr,filename = "bacteria_shannon_entropy_raw_by_year_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon diversity by year
alpha.div.raw.yr<-ggplot(alpha.div.metadat, aes(x=factor(Year), y=Shan_div.16s, fill=factor(Year))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sampling Year", x="Year", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.yr,filename = "bacterial_shannon_diversity_raw_by_year_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon entropy by month
alpha.ent.raw.month<-ggplot(alpha.div.metadat, aes(x=Month2, y=Shan_ent.16s, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Sampling Month", x="Month", y="Shannon Entropy", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.month,filename = "bacteria_shannon_entropy_raw_by_month_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon diversity by month
alpha.div.raw.month<-ggplot(alpha.div.metadat, aes(x=Month2, y=Shan_div.16s, fill=Month2)) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Sampling Month", x="Month", y="Shannon Diversity", fill="Month")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.month,filename = "bacteria_shannon_diversity_raw_by_month_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation
alpha.ent.raw.elev<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_ent.16s, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Elevation", x="Elevation", y="Shannon Entropy", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.elev,filename = "bacteria_shannon_entropy_raw_by_elevation_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation
alpha.div.raw.elev<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_div.16s, fill=factor(Elevation))) +geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Elevation", x="Elevation", y="Shannon Diversity", fill="Elevation")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.elev,filename = "bacteria_shannon_diversity_raw_by_elevation_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year
alpha.ent.raw.elev.yr<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_ent.16s, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.elev.yr,filename = "bacteria_shan_entropy_raw_by_elevation.year_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year
alpha.div.raw.elev.yr<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_div.16s, fill=factor(Year)))+geom_boxplot(color="black")+scale_x_discrete()+theme_bw()+scale_fill_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", fill="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.elev.yr,filename = "bacteria_shan_diversity_raw_by_elevation.year_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon entropy across elevation + year (no box fill, color outlines only)
alpha.ent.raw.elev.yr.2<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_ent.16s, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Entropy by Elevation and Year", x="Elevation", y="Shannon Entropy", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.ent.raw.elev.yr.2,filename = "bacteria_shannon_entropy_raw_by_elevation.year.2_9.13.2020.pdf", width=8, height=6, dpi=600)

# shannon diversity across elevation + year (no box fill, color outlines only)
alpha.div.raw.elev.yr.2<-ggplot(alpha.div.metadat, aes(x=factor(Elevation), y=Shan_div.16s, col=factor(Year)))+geom_boxplot()+scale_x_discrete()+theme_bw()+scale_colour_manual(values = saturation(wes1, 1))+theme_classic()+
  labs(title = "Bacterial Shannon Diversity by Elevation and Year", x="Elevation", y="Shannon Diversity", col="Year")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5, size=15))
ggsave(alpha.div.raw.elev.yr.2,filename = "bacteria_shannon_diversity_raw_by_elevation.year.2_9.13.2020.pdf", width=8, height=6, dpi=600)

dev.off()

#### Visualizing Alpha Div - Marko EEOB 230 ####

head(metadata)
## let's make a few plots

site<-as.factor(metadata$Site)

elev<-as.factor(metadata$Elevation)

month<-factor(metadata$Month, levels = c("July","August","October")) 
# if you don't specify levels here, will order them alphabetically and not temporally

boxplot(SR~month, xlab="Month", ylab="Species Richness")
boxplot(SR~elev, xlab="Elevation", ylab="Species Richness")
boxplot(SR~site, xlab="Site", ylab="Species Richness")

boxplot(Div1~month, xlab="Month", ylab="Shannon's diversity")
boxplot(Div1~elev, xlab="Elevation", ylab="Shannon's diversity")
boxplot(Div1~metadata$Elevation, xlab="Elevation", ylab="Shannon's diversity") # same as line above, just different way to do it
boxplot(Div1~site, xlab="Site", ylab="Shannon's diversity")

boxplot(Div2~elev, xlab="Elevation", ylab="Simpson diversity")

boxplot(Eve~elev, xlab="Elevation", ylab="Pielou evenness")

dev.off()

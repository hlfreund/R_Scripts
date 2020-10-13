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

#### Relative Abundance at all taxa levels ####

head(otus_counts_taxa)
m.otus_counts_taxa<-melt(otus_counts_taxa)
head(m.otus_counts_taxa)
colnames(m.otus_counts_taxa)[which(names(m.otus_counts_taxa) == "variable")] <- "SampleID"
colnames(m.otus_counts_taxa)[which(names(m.otus_counts_taxa) == "value")] <- "Count"

### * below we use the dcast() function to "cast" the data into a wide format based on given elements (column names), taking sum of "Count"
### * decostand(df, method="total) is the function (with argument total) used to get relative abundance of OTU table

## phylum ....
phyla_counts <- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ### 
head(phyla_counts) # counts by phyla per sample
rownames(phyla_counts)<-phyla_counts$SampleID
phyla_counts$SampleID<-NULL
phyla_RelAb<-data.frame(decostand(phyla_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(phyla_RelAb) # sanity check to make sure the transformation worked!
phyla_RelAb$SampleID<-rownames(phyla_RelAb)
head(phyla_RelAb)
phyla_RelAb_melt<-melt(phyla_RelAb)
head(phyla_RelAb_melt)
colnames(phyla_RelAb_melt)[which(names(phyla_RelAb_melt) == "variable")] <- "Phylum"
colnames(phyla_RelAb_melt)[which(names(phyla_RelAb_melt) == "value")] <- "Count"
head(phyla_RelAb_melt) ## relative abundance based on sum of counts by phyla!

## class ...

class_counts <- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Class, value.var="Count", fun.aggregate=sum)) ### 
head(class_counts) # counts by phyla per sample
rownames(class_counts)<-class_counts$SampleID
class_counts$SampleID<-NULL
class_RelAb<-data.frame(decostand(class_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(class_RelAb) # sanity check to make sure the transformation worked!
class_RelAb$SampleID<-rownames(class_RelAb)
head(class_RelAb)
class_RelAb_melt<-melt(class_RelAb)
head(class_RelAb_melt)
colnames(class_RelAb_melt)[which(names(class_RelAb_melt) == "variable")] <- "Class"
colnames(class_RelAb_melt)[which(names(class_RelAb_melt) == "value")] <- "Count"
head(class_RelAb_melt)
class_RelAb_melt<-merge(class_RelAb_melt,metadata,by="SampleID")
head(class_RelAb_melt)

## order ...
order_counts<- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Order, value.var="Count", fun.aggregate=sum)) ### 
head(order_counts) # counts by phyla per sample
rownames(order_counts)<-order_counts$SampleID
order_counts$SampleID<-NULL
order_RelAb<-data.frame(decostand(order_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(order_RelAb) # sanity check to make sure the transformation worked!
order_RelAb$SampleID<-rownames(order_RelAb)
head(order_RelAb)
order_RelAb_melt<-melt(order_RelAb)
head(order_RelAb_melt)
colnames(order_RelAb_melt)[which(names(order_RelAb_melt) == "variable")] <- "Order"
colnames(order_RelAb_melt)[which(names(order_RelAb_melt) == "value")] <- "Count"
head(order_RelAb_melt)
order_RelAb_melt<-merge(order_RelAb_melt,metadata,by="SampleID")
head(order_RelAb_melt)

## family ...
fam_counts<- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Family, value.var="Count", fun.aggregate=sum)) ### 
head(fam_counts) # counts by phyla per sample
rownames(fam_counts)<-fam_counts$SampleID
fam_counts$SampleID<-NULL
fam_RelAb<-data.frame(decostand(fam_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(fam_RelAb) # sanity check to make sure the transformation worked!
fam_RelAb$SampleID<-rownames(fam_RelAb)
head(fam_RelAb)
fam_RelAb_melt<-melt(fam_RelAb)
head(fam_RelAb_melt)
colnames(fam_RelAb_melt)[which(names(fam_RelAb_melt) == "variable")] <- "Family"
colnames(fam_RelAb_melt)[which(names(fam_RelAb_melt) == "value")] <- "Count"
head(fam_RelAb_melt)
fam_RelAb_melt<-merge(fam_RelAb_melt,metadata,by="SampleID")
head(fam_RelAb_melt)

## genus ...
genus_counts<- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ### 
head(genus_counts) # counts by phyla per sample
rownames(genus_counts)<-genus_counts$SampleID
genus_counts$SampleID<-NULL
gen_RelAb<-data.frame(decostand(genus_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(gen_RelAb) # sanity check to make sure the transformation worked!
gen_RelAb$SampleID<-rownames(gen_RelAb)
head(gen_RelAb)
gen_RelAb_melt<-melt(gen_RelAb)
head(gen_RelAb_melt)
colnames(gen_RelAb_melt)[which(names(gen_RelAb_melt) == "variable")] <- "Genus"
colnames(gen_RelAb_melt)[which(names(gen_RelAb_melt) == "value")] <- "Count"
head(gen_RelAb_melt)
gen_RelAb_melt<-merge(gen_RelAb_melt,metadata,by="SampleID")
head(gen_RelAb_melt)

## species ...
spec_counts<- as.data.frame(dcast(m.otus_counts_taxa, SampleID~Species, value.var="Count", fun.aggregate=sum)) ### 
head(spec_counts) # counts by phyla per sample
rownames(spec_counts)<-spec_counts$SampleID
spec_counts$SampleID<-NULL
spec_RelAb<-data.frame(decostand(spec_counts, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(spec_RelAb) # sanity check to make sure the transformation worked!
spec_RelAb$SampleID<-rownames(spec_RelAb)
head(spec_RelAb)
spec_RelAb_melt<-melt(spec_RelAb)
head(spec_RelAb_melt)
colnames(spec_RelAb_melt)[which(names(spec_RelAb_melt) == "variable")] <- "Genus.Species"
colnames(spec_RelAb_melt)[which(names(spec_RelAb_melt) == "value")] <- "Count"
head(spec_RelAb_melt)
spec_RelAb_melt<-merge(spec_RelAb_melt,metadata,by="SampleID")
head(spec_RelAb_melt)

#### Stacked Bar Plots of Relative Abundances ####

#phyla_RelAb_melt2<-subset(phyla_RelAb_melt, (Count)>(1/100)) ## DROP BACTERIAL PHYLA that are less than 1% abundant!!!!!!1

bac_phyla_its_RelAb<-ggplot(phyla_RelAb_melt, aes(x=factor(SampleID), y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title = "Bacteria Phyla Per Sample", x="Sample ID", y="Relative Abundance", col="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=90,vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))

ggsave(bac_phyla_its_RelAb,filename = "16s_phyla_Sierra_RelativeAbundance_wArchaea_9.13.2020.pdf", width=15, height=10, dpi=600)
### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)

bac_class_its_RelAb<-ggplot(class_RelAb_melt, aes(x=factor(SampleID), y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title = "Bacterial Classes Per Sample", x="Sample ID", y="Relative Abundance", col="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=90,vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))

ggsave(bac_class_its_RelAb,filename = "16S_class_Sierra_RelativeAbundance_wArchaea_9.13.2020.pdf", width=20, height=10, dpi=600)
### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)

#### Prepping data + metadata for analyses by elevation, year, month, etc ####
taxa_all_phy_class<-subset(m.otus_counts_taxa, select=c(Phylum, Class, SampleID, Count)) # subset only part of the metadata we need
head(taxa_all_phy_class)
metadata2<-subset(metadata, select=c(SampleID, Year, Month, Elevation, Site)) # subset only part of the metadata we need
head(metadata2)
tax_dat<-merge(taxa_all_phy_class, metadata2, by="SampleID") ## excludes unknown phyla and classes

head(tax_dat)

#### RELATIVE ABUNDANCE BY ELEVATION ####
## phylum ....

phyla_ele <- as.data.frame(dcast(tax_dat,Elevation~Phylum, value.var="Count", fun.aggregate=sum)) ### 
head(phyla_ele) # counts by phyla per sample
rownames(phyla_ele)<-phyla_ele$Elevation
phyla_ele$Elevation<-NULL
phyla_ele_RelAb<-data.frame(decostand(phyla_ele, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == elevations in this case)
rowSums(phyla_ele_RelAb) # sanity check to make sure the transformation worked!
phyla_ele_RelAb$Elevation<-rownames(phyla_ele_RelAb)
head(phyla_ele_RelAb)
phyla_ele_RA_melt<-melt(phyla_ele_RelAb)
head(phyla_ele_RA_melt)
colnames(phyla_ele_RA_melt)[which(names(phyla_ele_RA_melt) == "variable")] <- "Phylum"
colnames(phyla_ele_RA_melt)[which(names(phyla_ele_RA_melt) == "value")] <- "Count"
head(phyla_ele_RA_melt) ## relative abundance based on sum of counts by phyla!

## class ...

class_ele <- as.data.frame(dcast(tax_dat, Elevation~Class, value.var="Count", fun.aggregate=sum)) ### 
head(class_ele) # counts by phyla per sample
rownames(class_ele)<-class_ele$Elevation
class_ele$Elevation<-NULL
class_ele_RelAb<-data.frame(decostand(class_ele, method="total", MARGIN=1, na.rm=TRUE))  # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(class_ele_RelAb) # sanity check to make sure the transformation worked!
class_ele_RelAb$Elevation<-rownames(class_ele_RelAb)
head(class_ele_RelAb)
class_ele_RA_melt<-melt(class_ele_RelAb)
head(class_ele_RA_melt)
colnames(class_ele_RA_melt)[which(names(class_ele_RA_melt) == "variable")] <- "Class"
colnames(class_ele_RA_melt)[which(names(class_ele_RA_melt) == "value")] <- "Count"
head(class_ele_RA_melt)

#### Stacked Bar Plots of Relative Abundances BY ELEVATION ####
head(phyla_ele_RA_melt)
phyla_ele_RA_melt$Elevation2 <- factor(phyla_ele_RA_melt$Elevation, levels = c("400","1100","2000","2700")) 

bacterial_phyla_its_RelAb<-ggplot(phyla_ele_RA_melt, aes(x=Elevation2, y=Count, fill=Phylum)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title = "Bacterial Phyla Across Elevation", x="Elevation", y="Relative Abundance", col="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))

ggsave(bacterial_phyla_its_RelAb,filename = "16S_phyla_Sierra_RelativeAbundance_Elevation_wArchaea_9.13.2020.pdf", width=15, height=10, dpi=600)
### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)
class_ele_RA_melt$Elevation2 <- factor(class_ele_RA_melt$Elevation, levels = c("400","1100","2000","2700")) 

bacterial_class_its_RelAb<-ggplot(class_ele_RA_melt, aes(x=Elevation2, y=Count, fill=Class)) + geom_bar(stat="identity",colour="black")+theme_bw()+theme_classic()+scale_x_discrete()+labs(title = "Bacterial Classes Across Elevation", x="Elevation", y="Relative Abundance", col="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11),plot.title = element_text(hjust=0.5))

ggsave(bacterial_class_its_RelAb,filename = "16S_class_Sierra_RelativeAbundance_Elevation_wArchaea_9.13.2020.pdf", width=20, height=10, dpi=600)
### ^^^^ plotting with relative abundance data (decostand "total" transformation where all counts are divided by total hits across taxa per sample)

#
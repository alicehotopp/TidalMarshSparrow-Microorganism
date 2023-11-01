# Date updated: 11/1/23

# Script description: Final analyses on ITS sequences from tidal marsh sparrow feather samples for the paper "Plumage microorganism communities of tidal marsh sparrows". 

# SET UP ####

## Library packages ####
library(dada2); packageVersion('dada2')
library(GenomeInfoDBData) # Usually have to just manually check this in the Packages tab
library(BiocManager)
library(colorspace)
library(Matrix)
library(ShortRead)
library(phyloseq); packageVersion("phyloseq")
library(vegan)
library(ape); packageVersion("ape")
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(DESeq2)
library(permute)
library(randomForest)
library(rfPermute)
library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(DECIPHER)
library(phangorn)
library(Rcpp)
library(devtools)
library(DivNet)
library(ANCOMBC)
library(DT)
library(microbiome)
library(tidyverse)
library(biomformat)
library(VennDiagram)
library(microbiomeMarker)
library(RColorBrewer); library(extrafont); library(ggplot2);  library(cowplot); library(gridExtra); library(wesanderson); library(ggsci)

# Reading in unrarefied phyloseqs ####

## All feathers with >1,000 reads
Phy_ITS_Feathers_1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/Phy_ITS_Feathers_1000.rds")
## Maine feather phyloseq ####
ITS_ME_F_1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_ME_F_1000.rds")

## Spring tide feather phyloseq #####
ITS_Tide_1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_Tide_1000.rds") 

## Sediment phyloseq ####
ITS_S_1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_S_1000.rds")

## Feather and sediment phyloseq ####

ITS_SF_1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_SF_1000.rds")
# 8625 taxa and 81 samples

# Rarefied phyloseqs ####

## All feathers, 1000 reads, reseqs merged, rarefied ####
### Rarefaction curve ####
library(vegan)
otutab.f <- otu_table(Phy_ITS_Feathers_1000)
class(otutab.f) <- "matrix"
otutab.f <- t(otutab.f)  # transpose observations to rows          
rare.F <- rarecurve(otutab.f, step=1, lwd=2, ylab="OTU",  label=F) 

# Calculating the number of OTUs per sample in the OTU matrix
otutab.f.count <- colSums(otutab.f != 0)
sort(otutab.f.count) # Ranges from 3 OTUs to 1137 OTUs

mean(sample_sums(Phy_ITS_Feathers_1000)) # 12011.69 mean reads

sort(sample_sums(Phy_ITS_Feathers_1000))
# Lowest number of reads: 2881-26813,  1064 reads

Phy_ITS_Feathers_1000.rar <- rarefy_even_depth(Phy_ITS_Feathers_1000,
                                                sample.size=1064, 
                                                replace=FALSE, 
                                                trimOTUs=TRUE, 
                                                rngseed=711, 
                                                verbose=TRUE)
# 2267 OTUs were removed

# 5151 taxa and 65 samples

saveRDS(Phy_ITS_Feathers_1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/Phy_ITS_Feathers_1000.rar.rds")

Phy_ITS_Feathers_1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/Phy_ITS_Feathers_1000.rar.rds")

## Maine, 1000 reads, reseqs merged, rarefied ####
### Rarefaction curve ####
library(vegan)
otutab.f.ME <- otu_table(ITS_ME_F_1000)
class(otutab.f.ME) <- "matrix"
otutab.f.ME <- t(otutab.f.ME)  # transpose observations to rows          
rare.F.ME <- rarecurve(otutab.f.ME, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 2000)) 

sort(sample_sums(ITS_ME_F_1000))
# Lowest number of reads: 2881-26813 1064 reads

ITS_ME_F_1000.rar <- rarefy_even_depth(ITS_ME_F_1000,
                                       sample.size=1064, 
                                       replace=FALSE, 
                                       trimOTUs=TRUE, 
                                       rngseed=711, 
                                       verbose=TRUE)
# 1828 OTUs removed

# 3578 taxa and 34 samples 
sum(sample_sums(ITS_ME_F_1000.rar)) # 36176
mean(sample_sums(ITS_ME_F_1000.rar)) # 1064

saveRDS(ITS_ME_F_1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_ME_F_1000.rar.rds")

ITS_ME_F_1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_ME_F_1000.rar.rds")

## Spring tide feathers, reseqs merged, rarefied #####

### Rarefaction curve ####
library(vegan)
otutab.f.tide <- otu_table(ITS_Tide_1000)
class(otutab.f.tide) <- "matrix"
otutab.f.tide <- t(otutab.f.tide)  # transpose observations to rows          
rare.F.tide <- rarecurve(otutab.f.tide, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 2000)) 

sort(sample_sums(ITS_Tide_1000))
# lowest reads: 1064 , sample 2881-26813

ITS_Tide_1000.rar <- rarefy_even_depth(ITS_Tide_1000,
                                       sample.size=1064, 
                                       replace=FALSE, 
                                       trimOTUs=TRUE, 
                                       rngseed=711, 
                                       verbose=TRUE)

# 620 OTUs removed

saveRDS(ITS_Tide_1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_Tide_1000.rar.rds")

ITS_Tide_1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_Tide_1000.rar.rds")


## Sediment phyloseq ####

### Rarefaction curve ####
library(vegan)
otutab.s <- otu_table(ITS_S_1000)
class(otutab.s) <- "matrix"
otutab.s <- t(otutab.s)  # transpose observations to rows          
rare.s <- rarecurve(otutab.s, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 2000))

sort(sample_sums(ITS_S_1000))
# sample ITS_141 has 1407 reads

ITS_S_1000.rar <- rarefy_even_depth(ITS_S_1000,
                                    sample.size=1407, 
                                    replace=FALSE, 
                                    trimOTUs=TRUE, 
                                    rngseed=711, 
                                    verbose=TRUE)

# 636 OTUs removed

# 1207 taxa and 16 samples 

sum(sample_sums(ITS_S_1000.rar)) # 22512

saveRDS(ITS_S_1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_S_1000.rar.rds")

ITS_S_1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_S_1000.rar.rds")


## Sediment and feather samples, 1000 reads, rsqs merged ####

### Rarefaction curve ####
library(vegan)
otutab.sf <- otu_table(ITS_SF_1000)
class(otutab.sf) <- "matrix"
otutab.sf <- t(otutab.sf)  # transpose observations to rows          
rare.sf <- rarecurve(otutab.sf, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 2000))


sort(sample_sums(ITS_SF_1000))
#  2881-26813, 1064 reads

ITS_SF_1000.rar <- rarefy_even_depth(ITS_SF_1000,
                                     sample.size=1064, 
                                     replace=FALSE, 
                                     trimOTUs=TRUE, 
                                     rngseed=711, 
                                     verbose=TRUE)

# 2797 OTUs removed

# 5828 taxa and 81 samples

saveRDS(ITS_SF_1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_SF_1000.rar.rds")

ITS_SF_1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_SF_1000.rar.rds")

#BASIC INFO ####

## Maine feathers - looking at number of individuals ####

sample_data(ITS_ME_F_1000.rar)$BandNum
unique(sample_data(ITS_ME_F_1000.rar)$BandNum)
# 34 samples, 30 individuals


sum(sample_sums((ITS_ME_F_1000))) # 562434 unrarefied reads
mean(sample_sums((ITS_ME_F_1000))) # average reads per sample: 16542.18

sum(sample_sums((ITS_ME_F_1000.rar))) # 36176 reads (rarefied), 3578 taxa
mean(sample_sums((ITS_ME_F_1000.rar))) # average reads per sample: 1,064

## Spring tide feathers 
sample_data(ITS_Tide_1000.rar)$Individual_ID
unique(sample_data(ITS_Tide_1000.rar)$BandNum)
# 34 samples, 34 individuals
table(sample_data(ITS_Tide_1000.rar)$Species)
# NESP SALS SESP 
#   5   19   10 

sum(sample_sums((ITS_Tide_1000))) # 222967 reads (unrarefied)
mean(sample_sums((ITS_Tide_1000))) # 6557.853 mean reads/sample

sum(sample_sums((ITS_Tide_1000.rar))) # 36176 reads (rarefied) 
mean(sample_sums((ITS_Tide_1000.rar))) # 1064 mean reads/sample

## Sediment samples
sum(sample_sums((ITS_S_1000))) # 198586 reads (unrarefied)
mean(sample_sums((ITS_S_1000))) # 12411.62 mean reads/sample

sum(sample_sums((ITS_S_1000.rar))) # 22512 reads
mean(sample_sums((ITS_S_1000.rar))) # 1407 mean reads/sample

## Maine feathers - calculating the percentages that each phylum comprises ######

sort(table(tax_table(ITS_ME_F_1000.rar)[, 2]))

# 9 phyla, 3597 SVs

# Ascomycota 

2228/3597
0.6194051

# Basidiomycota
917/3597
0.2526551

# Rozellomycota
23/3597
0.006394217

# Chytridiomycota
20/3597
0.005560189

# Glomeromycota 
10/3797
0.002633658

# Mucoromycota
5/3797
0.001316829

# Basidiobolomycota
2/3597
0.00055601895

# Mortierellomycota 
1/3597
0.0002780095

# Aphelidiomycota
1/3597
0.0002780095

# Low abundance phyla
23 + 20 + 10 + 5 + 2 + 1 + 1
62/3597 # 0.01723659 1.7%

# Number of assigned SVs
2228 + 917 + 23 + 20 + 10 + 5 + 2 + 1 + 1
3207 # total

# Total number of SVs that were unassigned to taxa at the level of phylum:

3207/3597
0.8915763 # 89.2% Assigned SVs

1 - 0.8915763 
0.1084237 # 10.8% unassigned SVs

Bar_NAs <- plot_bar(ITS_ME_F_1000.rar, "Phylum")

sort(table(tax_table(ITS_ME_F_1000.rar)[, 3]))
sort(table(tax_table(ITS_ME_F_1000.rar)[, 4]))
sort(table(tax_table(ITS_ME_F_1000.rar)[, 5]))


## Spring tide feathers - calculating the percentages that each phylum comprises ######

sort(table(tax_table(ITS_Tide_1000.rar)[, 2]))

# 2332 SVs, 8 phyla

# Ascomycota 

1503/2332
0.6445111

# Basidiomycota 

503/2332
0.2156947

# Chytridiomycota

24/2332
0.0102916

# Glomeromycota
18/2332
0.007718696

# Rozellomycota

16/2332
0.006861063

# Mucoromycota
4/2332
0.001715266

# Basidiobolomycota
4/2332
0.001715266

# Mortierellomycota
2/2332
0.0008576329

# Low abundance phyla
24 + 18 + 16 + 4 + 4 + 2
68/2332
# 0.02915952

1503 + 503 + 68 
2074 # Number of SVs assigned to a phylum
2074/2332 
0.8893654 # 88.9% assigned to a phylum
1 - 0.8893654
0.1106346 # 11.1% unassigned to a phylum

# 34.2 % of SVs unassigned at the level of phyla

## All feathers - exploring keratinotlyic taxa ####

K_degrade <- c("Aspergillus", "Trychophyton", "Trichophyton", "Engyodontium", "Doratomyces", "Paecilomyces", "Onygena", "Pseudogymnoascus", "Cadophora", "Neosetophoma", "Auxarthron", "Gynmoascus", "Microsporum", "Scopulariopsis", "Sepedonium", "Myriodontium", "Aphanoascus", "Arthroderma", "Chrysosporium", "Coccidoides", "Gymnoascoideus", "Talaromyces", "Myrothecium", "Tritirachium", "Trichoderma", "Candida", "Geotrichum", "Fusarium", "Penicillium", "Oidiodendron", "Verticillium", "Keratinomyces", "Alternaria", "Trichurus", "Curvularia", "Cladosporium", "Geomyces", "Gleomastis", "Monodictys", "Myrothecium", "Stachybotrys", "Urocladium", "Neosetophoma", "Pseudogymnoascus", "Aphanoascus", "Acrodontium", "Botryotricum", "Chaetomium", "Chrysosporium", "Gliocladium", "Keratinophyton", "Malbranhea", "Microsporum", "Phytophthora")

ITS_K.F <- subset_taxa(Phy_ITS_Feathers_1000.rar, Genus %in% K_degrade)
ITS_K.F <- prune_samples(sample_sums(ITS_K.F)>0, ITS_K.F)
# 160 taxa and 58 samples

saveRDS(ITS_K.F, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.F.rds")
ITS_K.F <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.F.rds")

# saving csv of the taxonomy table
ITS_K.F_tax <- as.data.frame(tax_table(ITS_K.F))
write.csv(ITS_K.F_tax, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.F_tax.csv")

sort(table(tax_table(ITS_K.F)[, 6]))

# 16 genera
#  Alternaria: 51; Fusarium: 24; Cladosporium: 18; Candida: 18; Penicillium: 14; Talaromyces: 8; Aspergillus: 5; Stachybotrys: 3; Oidoidendron: 3; Nesosetophoma: 3; Myrothecium: 3; Curvulariua: 3; Cadophora: 3; Trichoderma: 2; Engyodontium: 1; Acrodontium: 1

# SVs in aspergillus (can cause respiratory infection and death in avian hosts)
aspergillus <- subset_taxa(Phy_ITS_Feathers_1000.rar, Genus == "Aspergillus")
aspergillus <- prune_samples(sample_sums(aspergillus)>0,aspergillus)
aspergillus <- prune_taxa(taxa_sums(aspergillus)>0,aspergillus)
# 5 taxa, 6 samples

sum(sample_sums(aspergillus))
# 11 reads

## Sediment samples  - calculating the percentages that each phylum comprises ######

sort(table(tax_table(ITS_S_1000.rar)[, 2]))

# 1207 SVs, 9 phyla

# Ascomycota
639/1207
0.5294118

# Basidiomycota
145/1207
0.1201326

# Rozellomycota
41/1207
0.03396852

# Glomeromycota
37/1207
0.03065452

# Chytridiomycota
21/1207
0.01739851

# Neocallimastigomycota
2/1207
0.001657001

# Mucoromycota
2/1207
0.001657001

# Basidiobolomycota
2/1207
0.001657001

# Aphelidiomycota
1/1207
0.0008285004

# percentages of the phyla with not many SVs:
41 + 37 + 21 + 2 + 2 + 2 + 1
106/1207
# 0.08782104

# Number of assigned SVs

639 + 145 + 41 + 37 + 21 + 2 + 2 + 2 + 1 # 890 SVs assigned to taxonomy

890/1207 # 0.7373654 assigned

# Number of unassigned SVs
1207-890
317 # unassigned SVs

317/1207 
0.2626346 # 26.3 % unassigned



## Maine feather taxonomy plots ####
ITS_ME_F_1000.rar.rel <- transform_sample_counts(ITS_ME_F_1000.rar, function(OTU) OTU/sum(OTU))

ITS_ME_F_1000.rar.rel.glom <- tax_glom(ITS_ME_F_1000.rar.rel, taxrank = "Phylum") # Using tax_glom gets rid of the NAs

tax_ME_F_1000_phy <- plot_bar(ITS_ME_F_1000.rar.rel.glom, fill = "Phylum") 
# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/tax_ME_F_1000_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_ME_F_1000_phy) 
dev.off()

tax_ME_F_1000_phy_NAs <- plot_bar(ITS_ME_F_1000.rar.rel, fill = "Phylum") 

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/tax_ME_F_1000_phy_NAs.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_ME_F_1000_phy_NAs) 
dev.off()

## Tide feather taxonomy plots ####
ITS_Tide_1000.rar.rel <- transform_sample_counts(ITS_Tide_1000.rar, function(OTU) OTU/sum(OTU))

ITS_Tide_1000.rar.rel.glom <- tax_glom(ITS_Tide_1000.rar.rel, taxrank = "Phylum") # Using tax_glom gets rid of the NAs

tax_Tide_1000_phy <- plot_bar(ITS_Tide_1000.rar.rel.glom, fill = "Phylum") + facet_grid(~sample_Species, space = "free", scales = "free")
# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/tax_Tide_1000_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_Tide_1000_phy) 
dev.off()

## Feather and sediment taxonomy plots ####
ITS_SF_1000.rar.rel <- transform_sample_counts(ITS_SF_1000.rar, function(OTU) OTU/sum(OTU))
ITS_SF_1000.rar.rel.glom <- tax_glom(ITS_SF_1000.rar.rel, taxrank = "Phylum") # Using tax_glom gets rid of the NAs

tax_SF_1000_phy <- plot_bar(ITS_SF_1000.rar.rel.glom, fill = "Phylum") + facet_grid(~Sample_Type, scales = "free", space = "free") 
# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/tax_SF_1000.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_SF_1000_phy) 
dev.off()

# AIM 1: ACROSS HOST SPECIES COMPARISONS ####

## ALHPA DIV - Host species ####

### Violin plots ####

host.F.cycle.rar <- plot_richness(ITS_Tide_1000.rar, 
                                  x="Species", 
                                  measures=c("Observed","Shannon"), 
                                  title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Species)) + 
  geom_boxplot(width = 0.1, aes(group=Species)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/viol.tide.F.host.rar.png", width = 23, height = 21, units = 'cm', res = 300)
grid.arrange(host.F.cycle.rar) 
dev.off()

### Test for alpha div variation across host species ####

ITS_Tide_1000.rar.rich <- estimate_richness(ITS_Tide_1000.rar, measure=c("Observed", "Shannon")) 
ITS_Tide_1000.rar.rich.even <- ITS_Tide_1000.rar.rich$Shannon/log(ITS_Tide_1000.rar.rich$Observed)
ITS_Tide_1000.rar.rich.sd = as(sample_data(ITS_Tide_1000.rar), "matrix")
ITS_Tide_1000.rar.rich.sd = as.data.frame(ITS_Tide_1000.rar.rich.sd)
ITS_Tide_1000.rar.rich.df <- cbind(ITS_Tide_1000.rar.rich, ITS_Tide_1000.rar.rich.even, ITS_Tide_1000.rar.rich.sd)
write.csv(ITS_Tide_1000.rar.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_Tide_1000.rar.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals, 3) run kruskall-wallis if non-normal

# Effect of host species on observed diversity
a_Tide_F_species.obs <- aov(Observed ~ Species, data = ITS_Tide_1000.rar.rich.df)
plot(a_Tide_F_species.obs)
summary(a_Tide_F_species.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Species      2   9710    4855   1.459  0.248
# Residuals   31 103150    3327       
# No effect of host species on observed diversity

res.species.ITS.obs <- a_Tide_F_species.obs$residuals
shapiro.test(res.species.ITS.obs)
# Residuals normally distributed

# Effect of host species on shannon diversity
a_Tide_F_species.shan <- aov(Shannon ~ Species, data = ITS_Tide_1000.rar.rich.df)
plot(a_Tide_F_species.shan)
summary(a_Tide_F_species.shan)

#.             Df Sum Sq Mean Sq F value Pr(>F)
# Species      2   3.73   1.863   1.112  0.342
# Residuals   31  51.94   1.675 

res.species.ITS.shan <- a_Tide_F_species.shan$residuals
shapiro.test(res.species.ITS.shan)
# ANOVA residuals not normally distributed, running kruskall wallis

kruskal.test(Shannon ~ Species, data = ITS_Tide_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 3.9769, df = 2, p-value = 0.1369
# Effect of host species not significant on shannon div

# BETA DIV: Across host species (spring tide) ####
### PCOA ####

# Host species - jaccard
ord_tide_F.rar.its <- ordinate(ITS_Tide_1000.rar, #calculate similarities
                               method ="PCoA", #ordination type
                               "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.tide.species <- plot_ordination(ITS_Tide_1000.rar, ord_tide_F.rar.its, type="samples", color="Species") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.tide.species.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.tide.species) 
dev.off()

# PCOA - bray curtis

ord_tide_F.rar.its.b <- ordinate(ITS_Tide_1000.rar, #calculate similarities
                                 method ="PCoA", #ordination type
                                 "bray") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.tide.species.b <- plot_ordination(ITS_Tide_1000.rar, ord_tide_F.rar.its.b, type="samples", color="Species") + geom_point(size = 3)

# Export figure
require(gridExtra)
png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.tide.species.bray.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.tide.species.b) 
dev.off()


### HOMOGENEITY OF VARIANCE - JACCARD ####

meta.df.tide.feath.its <- as(sample_data(ITS_Tide_1000.rar), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.tide.F.its <- distance(ITS_Tide_1000.rar, method="jaccard")

# checking betadispersion by host species

disper.host.its <- betadisper(dis.tide.F.its, meta.df.tide.feath.its$Species)
disper.host.its
permutest(disper.host.its)

#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)   
# Groups     2 0.0082732 0.0041366 10.876    999  0.003 **
# Residuals 31 0.0117901 0.0003803    
# Host species have different dispersions

# ANOVA of the dispersions
anova(disper.host.its)
#           Df    Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     2 0.0082732 0.0041366  10.876 0.0002638 ***
# Residuals 31 0.0117901 0.0003803   

# Tukey's Honest Significant Differences to determine which group(s) are more variable
host.tukey.its <- TukeyHSD(disper.host.its)
plot(host.tukey.its)

#              diff         lwr         upr     p adj
#SALS-NESP  0.04561497  0.02148996 0.069739984 0.0001670
#SESP-NESP  0.03382318  0.00753356 0.060112807 0.0094026
#SESP-SALS -0.01179179 -0.03054372 0.006960142 0.2832344

## PERMANOVA of effect of host species on community composition, jaccard ####

perma.tide.host <- adonis2(distance(ITS_Tide_1000.rar, method="jaccard") ~ Species,
                           data = meta.df.tide.feath.its)

#         Df SumOfSqs    R2      F Pr(>F)  
# Species   2   1.0781 0.067 1.1131  0.041 *
# Residual 31  15.0130 0.933                
# Total    33  16.0911 1.000  
# Host species  - sig diff community comp, but hard to interpret results because of the differences in dispersion

### HOMOGENEITY OF VARIANCE - Bray Curtis dissimilarity #####

## Bray Curtis disatance between samples, for testing dispersion
dis.tide.F.its.b <- distance(ITS_Tide_1000.rar, method="bray")

# checking betadispersion by host species

disper.host.its.b <- betadisper(dis.tide.F.its.b, meta.df.tide.feath.its$Species)
disper.host.its.b
permutest(disper.host.its.b)

#           Df   Sum Sq   Mean Sq   F      N.Perm Pr(>F)
# Groups     2 0.005506 0.0027528 2.2607    999  0.117
# Residuals 31 0.037747 0.0012177
# Not sig diff in group dispersions

### PERMANOVA - effect of species on community comp, bray curtis dissimilarity ####

perma.tide.host.b <- adonis2(distance(ITS_Tide_1000.rar, method="bray") ~ Species,
                             data = meta.df.tide.feath.its)

#           Df SumOfSqs   R2      F    Pr(>F)  
# Species   2   1.1411 0.07239 1.2096  0.014 *
# Residual 31  14.6226 0.92761                
# Total    33  15.7637 1.00000 

# Pairwise Adonis to assess which species are different from each other 

library(pairwiseAdonis)
pair.host.its <- pairwise.adonis(dis.tide.F.its.b, meta.df.tide.feath.its$Species)
pair.host.its

#    pairs.       Df SumsOfSqs  F.Model    R2     p.value p.adjusted sig
# 1 SALS vs SESP  1 0.5685217 1.211682 0.04294966   0.083      0.249    
# 2 SALS vs NESP  1 0.5841810 1.246778 0.05363230   0.046      0.138    
# 3 SESP vs NESP  1 0.5527385 1.146286 0.08103091   0.066      0.198    

## DESEQ - HOST SPECIES ####

# Note on DESEQ comparison order: "By default, R will choose a reference level for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels." https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_Tide_1000, ~ Species)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
# Using the contrast argument to pull out pairwise comparisons - NESP and SALS
res.NESP.SALS = results(diagdds, contrast =c("Species", "NESP", "SALS"))
res.NESP.SALS = res.NESP.SALS[order(res.NESP.SALS$padj, na.last=NA), ]
alpha = 0.01
sigtab.NESP.SALS = res.NESP.SALS[(res.NESP.SALS$padj < alpha), ]
sigtab.NESP.SALS = cbind(as(sigtab.NESP.SALS, "data.frame"), as(tax_table(ITS_Tide_1000)[rownames(sigtab.NESP.SALS), ], "matrix")) 

head(sigtab.NESP.SALS)
dim(sigtab.NESP.SALS)
# 11 differentially abundant taxa between NESP and SALS. 

write.csv(sigtab.NESP.SALS, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.NESP.SALS.csv")

# Using the contrast argument to pull out pairwise comparisons - NESP and SESP
res.NESP.SESP = results(diagdds, contrast =c("Species", "NESP", "SESP"))
res.NESP.SESP = res.NESP.SESP[order(res.NESP.SESP$padj, na.last=NA), ]
alpha = 0.01
sigtab.NESP.SESP = res.NESP.SALS[(res.NESP.SESP$padj < alpha), ]
sigtab.NESP.SESP = cbind(as(sigtab.NESP.SESP, "data.frame"), as(tax_table(ITS_Tide_1000)[rownames(sigtab.NESP.SESP), ], "matrix")) 

head(sigtab.NESP.SESP)
dim(sigtab.NESP.SESP)
# 6 differentially abundant taxa between NESP and SESP. 

write.csv(sigtab.NESP.SESP, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.NESP.SESP.csv")

# Using the contrast argument to pull out pairwise comparisons - SALS and SESP
res.SALS.SESP = results(diagdds, contrast =c("Species", "SALS", "SESP"))
res.SALS.SESP = res.SALS.SESP[order(res.SALS.SESP$padj, na.last=NA), ]
alpha = 0.01
sigtab.SALS.SESP = res.SALS.SESP[(res.SALS.SESP$padj < alpha), ]
sigtab.SALS.SESP = cbind(as(sigtab.SALS.SESP, "data.frame"), as(tax_table(ITS_Tide_1000)[rownames(sigtab.SALS.SESP), ], "matrix")) 

head(sigtab.SALS.SESP)
dim(sigtab.SALS.SESP)
# 7 differentially abundant taxa between SALS and SESP. 

write.csv(sigtab.SALS.SESP, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.SALS.SESP.csv")

# KERATINOLYTC taxa in spring tide feathers ####
ITS_K.tide <- subset_taxa(ITS_Tide_1000.rar, Genus %in% K_degrade)

sort(table(tax_table(ITS_K.tide)[, 6]))
# 14 Genera
# Alternaria (12 SVs), Penicillium (9 SVs), Fusarium (9 SVs), Candida (9 SVs), Talaromyces (7 SVs), Cladosporium (6 SVs), Neosetophoma (2 SVs), Curvularia (2 SVs), Cadophora (2 SVs), Trichoderma (1 SV), Myrothecium (1 SV), Engyodontium (1 SV), Chaetomium (1 SV), Aspergillus (1 SV)

ITS_K.tide <- prune_taxa(taxa_sums(ITS_K.tide)>0, ITS_K.tide)
ITS_K.tide <- prune_samples(sample_sums(ITS_K.tide)>0, ITS_K.tide)
# 62 taxa and 28 samples

saveRDS(ITS_K.tide, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.tide.rds")

ITS_K.tide <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.tide.rds")

## ALHPA DIVERSITY - Keratinolytic SVs across host species ####

### Taxonomy plots ####

tax_tide.K_species <- plot_bar(ITS_K.tide, fill = "Genus") + facet_grid(~sample_Species, scales = "free", space = "free")

# Export figure
require(gridExtra)
png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/tax_tide.K_species.png", width = 21, height = 21, units = 'cm', res = 300)
grid.arrange(tax_tide.K_species) 
dev.off()


### Violin plots comparing host species  #####
viol.tide.K <- plot_richness(ITS_K.tide , 
                             x="Species", 
                             measures=c("Observed","Shannon"), 
                             title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Species)) + 
  geom_boxplot(width = 0.1, aes(group=Species)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)
png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/viol.tide.K.png", width = 21, height = 21, units = 'cm', res = 300)
grid.arrange(viol.tide.K) 
dev.off()

### Testing alpha diversity ####

ITS_K_tide.F.rich <- estimate_richness(ITS_K.tide, measure=c("Observed", "Shannon")) 
ITS_K_tide.F.rich.even <- ITS_K_tide.F.rich$Shannon/log(ITS_K_tide.F.rich$Observed)
ITS_K_tide.F.rich.sd = as(sample_data(ITS_K.tide), "matrix")
ITS_K_tide.F.rich.sd = as.data.frame(ITS_K_tide.F.rich.sd)
ITS_K_tide.F.rich.df <- cbind(ITS_K_tide.F.rich, ITS_K_tide.F.rich.even, ITS_K_tide.F.rich.sd)
write.csv(ITS_K_tide.F.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K_tide.F.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals, 3) run kruskall-wallis if not sig

# Effect of host species on observed diversity of keratinolytic SVs
a_K_tide.species.obs <- aov(Observed ~ Species, data = ITS_K_tide.F.rich.df)
plot(a_K_tide.species.obs)
summary(a_K_tide.species.obs)

#               Df Sum Sq Mean Sq F value Pr(>F)
# Species       2  16.18   8.092   1.419  0.261
# Residuals    25 142.53   5.701 

res.K_tide.species.obs  <- a_K_tide.species.obs$residuals
shapiro.test(res.K_tide.species.obs)
# W = 0.85133, p-value = 0.0009947 Not normal, run kruskall wallis

kruskal.test(Observed ~ Species, data = ITS_K_tide.F.rich.df)
# Kruskal-Wallis chi-squared = 5.8867, df = 2, p-value = 0.05269
# Not a sig difference in host species observed diversity of keratinolytic fungi

# Effect of host species on shannon diversity of keratinolytic SVs
a_K_tide.species.shan <- aov(Shannon ~ Species, data = ITS_K_tide.F.rich.df)
plot(a_K_tide.species.shan)
summary(a_K_tide.species.shan)

#               Df Sum Sq Mean Sq F value Pr(>F)  
# Species       2  1.363  0.6816   2.955 0.0705 .
# Residuals    25  5.766  0.2307 

res.K_tide.species.shan  <- a_K_tide.species.shan$residuals
shapiro.test(res.K_tide.species.shan)
# W = 0.97535, p-value = 0.7285 Normal, report results of ANOVA

## BETA DIV - KERATINOLYTIC HOST SPECIES ####
### PCOA ####

# Host species - jaccard
ord_tide_F.its.K <- ordinate(ITS_K.tide, #calculate similarities
                               method ="PCoA", #ordination type
                               "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.tide.species.K <- plot_ordination(ITS_K.tide, ord_tide_F.its.K, type="samples", color="Species") + geom_point(size = 3)

# Export figure
require(gridExtra)
png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.tide.species.K.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.tide.species.K) 
dev.off()

### DISPERSION & COMP ####
#### HOMOGENEITY OF VARIANCE - JACCARD ####
# Get metadata
meta.df.tide.feath.k.its <- as(sample_data(ITS_K.tide), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.tide.K.its <- distance(ITS_K.tide, method="jaccard")

# checking betadispersion by host species

disper.host.k.its <- betadisper(dis.tide.K.its, meta.df.tide.feath.k.its$Species)
disper.host.k.its
permutest(disper.host.k.its)

#            Df    Sum Sq   Mean Sq   F.   N.Perm Pr(>F)    
# Groups     2 0.0116291 0.0058145 17.238    999  0.001 ***
# Residuals 25 0.0084326 0.0003373  


# Tukey's Honest Significant Differences to determine which group(s) are more variable
host.tukey.k.its <- TukeyHSD(disper.host.k.its)
plot(disper.host.k.its)

#                diff         lwr         upr     p adj
# SALS-NESP  0.05941903  0.03399696 0.084841102 0.0000132
# SESP-NESP  0.04198204  0.01330904 0.070655031 0.0033728
# SESP-SALS -0.01743700 -0.03798113 0.003107138 0.1072607

# Pairs with NESP significantly different - due to small NESP sample size?


#### HOMOGENEITY OF VARIANCE - BRAY CURTIS ####

# Get metadata
meta.df.tide.feath.k.its <- as(sample_data(ITS_K.tide), "data.frame")

## Bray Curtis disatance between samples, for testing dispersion
dis.tide.K.its.b <- distance(ITS_K.tide, method="bray")

# checking betadispersion by host species

disper.host.k.its.b <- betadisper(dis.tide.K.its.b, meta.df.tide.feath.k.its$Species)
disper.host.k.its.b
permutest(disper.host.k.its.b)

#            Df    Sum Sq   Mean Sq   F     N.Perm Pr(>F)  
# Groups     2 0.0083745 0.0041872 5.1418    999  0.014 *
# Residuals 25 0.0203590 0.0008144 

# Tukey's Honest Significant Differences to determine which group(s) are more variable
host.tukey.k.its.b <- TukeyHSD(disper.host.k.its.b)
plot(disper.host.k.its.b)

#                diff         lwr         upr     p adj
# SALS-NESP  0.050850667  0.011349712 0.09035162 0.0098744
# SESP-NESP  0.041689625 -0.002862634 0.08624188 0.0697883
# SESP-SALS -0.009161042 -0.041082634 0.02276055 0.7570837

#### PERMANOVA - JACCARD ####

perma.tide.host.k <- adonis2(distance(ITS_K.tide, method="jaccard") ~ Species,
                           data = meta.df.tide.feath.k.its)

#           Df SumOfSqs   R2      F     Pr(>F)
# Species   2   0.9637 0.07334 0.9893   0.48
# Residual 25  12.1770 0.92666              
# Total    27  13.1407 1.00000   
# Host species do not differ in community composition of keratinolytic SVs

#### PERMANOVA - BRAY CURTIS ####
perma.tide.host.k.b <- adonis2(distance(ITS_K.tide, method="bray") ~ Species,
                             data = meta.df.tide.feath.k.its)

#           Df SumOfSqs      R2      F Pr(>F)
# Species   2   0.9858 0.07597 1.0278  0.367
# Residual 25  11.9901 0.92403              
# Total    27  12.9760 1.00000   
# Host species do not differ in community composition of keratinolytic SVs


# AIM 2: VARIATION WITHIN NESP ####

## ALPHA DIV: Across Maine NESP variables ####

### Test for alpha div variation across NESP ME feather variables ####

ITS_ME_F_1000.rar.rich <- estimate_richness(ITS_ME_F_1000.rar, measure=c("Observed", "Shannon")) 
ITS_ME_F_1000.rar.rich.even <- ITS_ME_F_1000.rar.rich$Shannon/log(ITS_ME_F_1000.rar.rich$Observed)
ITS_ME_F_1000.rar.rich.sd = as(sample_data(ITS_ME_F_1000.rar), "matrix")
ITS_ME_F_1000.rar.rich.sd = as.data.frame(ITS_ME_F_1000.rar.rich.sd)
ITS_ME_F_1000.rar.rich.df <- cbind(ITS_ME_F_1000.rar.rich, ITS_ME_F_1000.rar.rich.even, ITS_ME_F_1000.rar.rich.sd)
write.csv(ITS_ME_F_1000.rar.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_ME_F_1000.rar.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals, 3) run kruskall-wallis if not sig

# ANOVAs to test for differences in normally distributed alpha diversity

# Observed by marsh
a_ME_F_marsh.obs <- aov(Observed ~ Marsh, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_marsh.obs)
summary(a_ME_F_marsh.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)
# Marsh        3   9846    3282   0.529  0.666
# Residuals   30 186023    6201   

res.marsh.ITS.obs <- a_ME_F_marsh.obs$residuals
shapiro.test(res.marsh.ITS.obs)
# W = 0.93301, p-value = 0.0383 - Normal, report ANOVA results

# Shannon by marsh
a_ME_F_marsh.shan <- aov(Shannon ~ Marsh, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_marsh.shan)
summary(a_ME_F_marsh.shan)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Marsh        3   1.07   0.357   0.304  0.822
# Residuals   30  35.25   1.175 
res.marsh.ITS.shan <- a_ME_F_marsh.shan$residuals
shapiro.test(res.marsh.ITS.shan)
# W = 0.85962, p-value = 0.0004626 Not normally distributed, run kruskall-wallis

# Kruskall-wallis for testing variation in Shannon div by marsh
kruskal.test(Shannon ~ Marsh, data = ITS_ME_F_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.9097, df = 3, p-value = 0.5914

# Observed by month
a_ME_F_month.obs <- aov(Observed ~ Month, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_month.obs)
summary(a_ME_F_month.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)
# Month        3  20507    6836   1.169  0.338
# Residuals   30 175362    5845 

res.month.ITS.obs <- a_ME_F_month.obs$residuals
shapiro.test(res.month.ITS.obs)
# W = 0.93728, p-value = 0.05118 - Normal, report ANOVA results

# Shannon by month
a_ME_F_month.shan <- aov(Shannon ~ Month, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_month.shan)
summary(a_ME_F_month.shan)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Month        3   0.16   0.054   0.045  0.987
# Residuals   30  36.16   1.206   

res.month.ITS.shan <- a_ME_F_month.shan$residuals
shapiro.test(res.month.ITS.shan)
# W = 0.88021, p-value = 0.001433 - Not normal, run Kruskall-wallis
kruskal.test(Shannon ~ Month, data = ITS_ME_F_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 0.22689, df = 3, p-value = 0.9731

# No sig variation in Observed or Shannon div by month

# Observed by sex
a_ME_F_sex.obs <- aov(Observed ~ Sex, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_sex.obs)
summary(a_ME_F_sex.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)  
# Sex          1  19594   19594   3.557 0.0684 .
# Residuals   32 176276    5509   

res.sex.ITS.obs <- a_ME_F_sex.obs$residuals
shapiro.test(res.sex.ITS.obs)
# W = 0.9424, p-value = 0.0727 - Normal, report ANOVA results

# Shannon by sex
a_ME_F_sex.shan <- aov(Shannon ~ Sex, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_sex.shan)
summary(a_ME_F_sex.shan)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Sex          1   2.82   2.819   2.692  0.111
# Residuals   32  33.51   1.047   

res.sex.ITS.shan <- a_ME_F_sex.shan$residuals
shapiro.test(res.sex.ITS.shan)
# W = 0.88654, p-value = 0.002062 - Not normal, run Kruskall wallis
kruskal.test(Shannon ~ Sex, data = ITS_ME_F_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 2.4994, df = 1, p-value = 0.1139

# No sig variation in Observed or Shannon div by sex

# Observed by tidal range (meters)
a_ME_F_rng.obs <- aov(Observed ~ TidalRange_m, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_rng.obs)
summary(a_ME_F_rng.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange_m  4  19216    4804   0.789  0.542
# Residuals    29 176654    6092     

res.rng.ITS.obs <- a_ME_F_rng.obs$residuals
shapiro.test(res.rng.ITS.obs)
# W = 0.95479, p-value = 0.1707 - Normal, report ANOVA results

# Shannon by tidal range
a_ME_F_rng.shan <- aov(Shannon ~ TidalRange_m, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_rng.shan)
summary(a_ME_F_rng.shan)
#              Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange_m  4   1.24  0.3104   0.257  0.903
# Residuals    29  35.08  1.2098  

res.rng.ITS.shan <- a_ME_F_rng.shan$residuals
shapiro.test(res.rng.ITS.shan)
# W = 0.87881, p-value = 0.001323 - Not normal, run Kruskall-Wallis
kruskal.test(Shannon ~ TidalRange_m, data = ITS_ME_F_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.4272, df = 4, p-value = 0.8395
# No sig variation in observed or alpha diversity by tidal range in meters

# Observed by tidal cycle - spring, neap, or intermediate
a_ME_F_cycle.obs <- aov(Observed ~ TideCycle, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_cycle.obs)
summary(a_ME_F_cycle.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)  
# TideCycle    2  27610   13805   2.543 0.0949 .
# Residuals   31 168259    5428     

res.cycle.ITS.obs <- a_ME_F_cycle.obs$residuals
shapiro.test(res.cycle.ITS.obs)
# W = 0.94818, p-value = 0.1083 - Normal, report ANOVA results

# Shannon by tidal cycle - spring, neap, or intermediate
a_ME_F_cycle.shan <- aov(Shannon ~ TideCycle, data = ITS_ME_F_1000.rar.rich.df)
plot(a_ME_F_cycle.shan)
summary(a_ME_F_cycle.shan)
#               Df Sum Sq Mean Sq F value Pr(>F)
# TideCycle    2   1.07  0.5362   0.472  0.628
# Residuals   31  35.25  1.1372     

res.cycle.ITS.shan <- a_ME_F_cycle.shan$residuals
shapiro.test(res.cycle.ITS.shan)
# W = 0.85664, p-value = 0.0003954 - Not normal, run kruskall wallis
kruskal.test(Shannon ~ TideCycle, data = ITS_ME_F_1000.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.6683, df = 2, p-value = 0.4342

# No sig diff in alpha diversity by tidal cycle (spring, neap, or intermediate)


## BETA DIV: Across Maine NESP variables ####

### DISPERSION & COMMUNITY COMPOSITION ####

library(vegan)

#### HOMOGENITY OF VARIANCE - Jaccard #####

# data frame of the sample data
meta.df.Me.feath.its <- as(sample_data(ITS_ME_F_1000.rar), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.ME.F.its <- distance(ITS_ME_F_1000.rar, method="jaccard")

# sex
disper.sex.its <- betadisper(dis.ME.F.its, meta.df.Me.feath.its$Sex)
disper.sex.its
permutest(disper.sex.its)
#            Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)    
# Groups     1 0.0046640 0.0046640 26.24    999  0.001 ***
# Residuals 32 0.0056879 0.0001777  
# dispersion diff by sex

# ANOVA of the dispersions
anova(disper.sex.its)

# Response: Distances
#           Df    Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     1 0.0046640 0.0046640   26.24 1.392e-05 ***
# Residuals 32 0.0056879 0.0001777  

# Tukey's Honest Significant Differences to determine which group(s) are more variable
sex.tukey.its <- TukeyHSD(disper.sex.its)
plot(sex.tukey.its)

#            diff        lwr        upr    p adj
# CP-BP 0.02654776 0.01599112 0.03710439 1.39e-05


# month
disper.month.its <- betadisper(dis.ME.F.its, meta.df.Me.feath.its$Month)
disper.month.its
permutest(disper.month.its)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     2 0.038314 0.0191569 96.851    999  0.001 ***
# Residuals 31 0.006132 0.0001978  
# dispersion diff by month

# ANOVA of the dispersions
anova(disper.month.its)
#           Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     2 0.038314 0.0191569  96.851 4.637e-14 ***
# Residuals 31 0.006132 0.0001978       

# Tukey's Honest Significant Differences to determine which group(s) are more variable

month.tukey.its <- TukeyHSD(disper.month.its)
plot(month.tukey.its)
#                 diff         lwr         upr   p adj
# July-August  0.03863264  0.02376361  0.05350168 1.2e-06
# June-August -0.07521086 -0.09909705 -0.05132466 0.0e+00
# June-July   -0.11384350 -0.13504040 -0.09264661 0.0e+00
# Looks like July has the most different dispersion, which makes sense because it has the highest sample size. 

# marsh
disper.marsh.its <- betadisper(dis.ME.F.its, meta.df.Me.feath.its$Marsh)
disper.marsh.its
permutest(disper.marsh.its)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     3 0.075419 0.0251396 73.858    999  0.001 ***
# Residuals 30 0.010211 0.0003404  
# dispersion diff by marsh

# ANOVA of the dispersions
anova(disper.marsh.its)
#           Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     3 0.075419 0.0251396  73.858 5.919e-14 ***
# Residuals 30 0.010211 0.0003404 

# Tukey's Honest Significant Differences to determine which group(s) are more variable

marsh.tukey <- TukeyHSD(disper.marsh.its)
par(mar=c(2,2,2,2))
plot(marsh.tukey)

#           diff           lwr         upr     p adj
# HA-BB  0.16416042  0.1249439450  0.20337689 0.0000000
# NG-BB  0.18462789  0.1473349765  0.22192080 0.0000000
# PL-BB  0.10391148  0.0604665494  0.14735640 0.0000020
# NG-HA  0.02046747  0.0001677805  0.04076716 0.0475528
# PL-HA -0.06024894 -0.0903948837 -0.03010300 0.0000389
# PL-NG -0.08071641 -0.1083136458 -0.05311918 0.0000000

## All sig diff - closest dispersions are NG and HA. Biggest sample sizes.

# tidal range (m)
disper.rng.its <- betadisper(dis.ME.F.its, meta.df.Me.feath.its$TidalRange_m)
disper.rng.its
permutest(disper.rng.its)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     4 0.028547 0.0071368 19.961    999  0.001 ***
# Residuals 29 0.010368 0.0003575   
# dispersion diff by tidal range

# ANOVA of the dispersions
anova(disper.rng.its)
#           Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     4 0.028547 0.0071368  19.961 5.456e-08 ***
# Residuals 29 0.010368 0.0003575       

# Tukey's Honest Significant Differences to determine which group(s) are more variable
rng.tukey.its <- TukeyHSD(disper.rng.its)
plot(rng.tukey.its)

#                         diff          lwr           upr     p adj
# 3 to 3.5-2.5 to 3  0.0260450161 -0.002923357  0.0550133891 0.0940622
# 3.5 to 4-2.5 to 3  0.0263382489 -0.001143563  0.0538200605 0.0652766
# 4 to 4.5-2.5 to 3 -0.0245963792 -0.060075245  0.0108824871 0.2842397
# 4.5 to 5-2.5 to 3 -0.0671898308 -0.106054982 -0.0283246800 0.0002150
# 3.5 to 4-3 to 3.5  0.0002932328 -0.023943447  0.0245299125 0.9999996
# 4 to 4.5-3 to 3.5 -0.0506413953 -0.083670422 -0.0176123682 0.0010087
# 4.5 to 5-3 to 3.5 -0.0932348469 -0.129877262 -0.0565924313 0.0000004
# 4 to 4.5-3.5 to 4 -0.0509346280 -0.082667891 -0.0192013653 0.0005733
# 4.5 to 5-3.5 to 4 -0.0935280796 -0.129006946 -0.0580492133 0.0000002
# 4.5 to 5-4 to 4.5 -0.0425934516 -0.084572612 -0.0006142909 0.0454404

# Tidal cycle

disper.cycle.its <- betadisper(dis.ME.F.its, meta.df.Me.feath.its$TideCycle)
disper.cycle.its
permutest(disper.cycle.its )

#           Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)    
# Groups     2 0.0103159 0.0051580 18.48    999  0.001 ***
# Residuals 31 0.0086526 0.0002791 

# ANOVA of the dispersions
anova(disper.cycle.its)
#.           Df    Sum Sq   Mean Sq F value    Pr(>F)    
# Groups     2 0.0103159 0.0051580   18.48 5.203e-06 ***
# Residuals 31 0.0086526 0.0002791 

# Tukey's Honest Significant Differences to determine which group(s) are more variable
cycle.tukey.its <- TukeyHSD(disper.cycle.its)
plot(cycle.tukey.its)

#                          diff         lwr          upr     p adj
# Neap-Intermediate   -0.01386366 -0.02936682  0.001639494 0.0868102
# Spring-Intermediate -0.05158298 -0.07250185 -0.030664116 0.0000030
# Spring-Neap         -0.03771932 -0.05960626 -0.015832375 0.0005315


#### HOMOGENITY OF VARIANCE - Bray #####

meta.df.Me.feath.its <- as(sample_data(ITS_ME_F_1000.rar), "data.frame") # dataframe of the sample data
dis.ME.F.its.b <- distance(ITS_ME_F_1000.rar, method="bray") # distance 

# sex
disper.sex.its.b <- betadisper(dis.ME.F.its.b, meta.df.Me.feath.its$Sex)
disper.sex.its.b
permutest(disper.sex.its.b)

#             Df    Sum Sq   Mean Sq   F N.Perm Pr(>F)   
# Groups     1 0.0046476 0.0046476 7.3    999   0.01 **
# Residuals 32 0.0203730 0.0006367  

# marsh
disper.marsh.its.b <- betadisper(dis.ME.F.its.b, meta.df.Me.feath.its$Marsh)
disper.marsh.its.b
permutest(disper.marsh.its.b)

#            Df   Sum Sq   Mean Sq   F.   N.Perm Pr(>F)    
# Groups     3 0.075862 0.0252875 20.981    999  0.001 ***
# Residuals 30 0.036158 0.0012053

# month
disper.month.its.b <- betadisper(dis.ME.F.its.b, meta.df.Me.feath.its$Month)
disper.month.its.b
permutest(disper.month.its.b)

#            Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     2 0.038941 0.019470 28.506    999  0.001 ***
# Residuals 31 0.021174 0.000683 

# tide cycle
disper.cycle.its.b <- betadisper(dis.ME.F.its.b, meta.df.Me.feath.its$TideCycle)
disper.cycle.its.b
permutest(disper.cycle.its.b)

#           Df   Sum Sq   Mean Sq   F     N.Perm Pr(>F)   
# Groups     2 0.010650 0.0053250 6.0682    999  0.005 **
# Residuals 31 0.027203 0.0008775 

# tidal range 
disper.rng.its.b <- betadisper(dis.ME.F.its.b, meta.df.Me.feath.its$TidalRange)
disper.rng.its.b
permutest(disper.rng.its.b)

#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     3 0.0080248 0.0026749 2.6747    999  0.061 .
# Residuals 30 0.0300028 0.0010001    
# Not sig different 

#### PERMANOVA - Jaccard #####
# All variables had groups with different dispersions, so these results should be interpreted with caution

# effect of sex, jaccard

perma.ME.f.sex <- adonis2(distance(ITS_ME_F_1000.rar, method="jaccard") ~ Sex,
                          data = meta.df.Me.feath.its)

#          Df SumOfSqs      R2      F Pr(>F)
# Sex       1   0.4644 0.02896 0.9545  0.716
# Residual 32  15.5699 0.97104              
# Total    33  16.0343 1.00000 
# not sig

# effect of month, jaccard

perma.ME.f.month <- adonis2(distance(ITS_ME_F_1000.rar, method="jaccard") ~ Month,
                            data = meta.df.Me.feath.its)

#          Df SumOfSqs      R2      F Pr(>F)  
# Month     2   1.0404 0.06489 1.0755  0.053 .
# Residual 31  14.9939 0.93511                
# Total    33  16.0343 1.00000  
# No sig effect of month on community composition

# effect of marsh, jaccard

perma.ME.f.marsh <- adonis2(distance(ITS_ME_F_1000.rar, method="jaccard") ~ Marsh,
                            data = meta.df.Me.feath.its)

#          Df SumOfSqs      R2      F Pr(>F)  
# Marsh     3   1.5314 0.09551 1.0559  0.073 .
# Residual 30  14.5029 0.90449                
# Total    33  16.0343 1.00000 
# Marsh not sig

# effect of tidal range, jaccard

perma.ME.f.rng <- adonis2(distance(ITS_ME_F_1000.rar, method="jaccard") ~ TidalRange,
                          data = meta.df.Me.feath.its)

#             Df SumOfSqs   R2      F    Pr(>F)
# TidalRange  3   1.4898 0.09291 1.0243  0.229
# Residual   30  14.5446 0.90709              
# Total      33  16.0343 1.00000   
# tidal range not sig

# effect of tide cycle, jaccard

perma.ME.f.cyc <- adonis2(distance(ITS_ME_F_1000.rar, method="jaccard") ~ TideCycle,
                          data = meta.df.Me.feath.its)

#            Df SumOfSqs   R2      F   Pr(>F)  
# TideCycle  2    1.059 0.06605 1.0961  0.034 *
# Residual  31   14.975 0.93395                
# Total     33   16.034 1.00000  
# tide cycle sig, but dispersions different, so results hard to interpret

#### PERMANOVA - Bray Curtis #####

# Only running for tidal range because this was the only group with not sig different group dispersion
perma.ME.f.rng.b <- adonis2(distance(ITS_ME_F_1000.rar, method="bray") ~ TidalRange,
                            data = meta.df.Me.feath.its)

#            Df SumOfSqs    R2     F   Pr(>F)
# TidalRange  3   1.4888 0.09519 1.052  0.164
# Residual   30  14.1517 0.90481             
# Total      33  15.6404 1.00000  

### DESEQ Differential abundance - Maine variables ####

# Interpreting the DESEQ results sigtab: "Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula"
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

##### DESEQ - Sex #####
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_ME_F_1000, ~ Sex)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab.sex = res[(res$padj < alpha), ]
sigtab.sex = cbind(as(sigtab, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab), ], "matrix")) 
head(sigtab.sex)
dim(sigtab.sex)
# 14 differentially abundant taxa between sexes

write.csv(sigtab.sex, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.sex.csv")


# Phylum order - this orders the phylums based on logchange for plotting. 
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order - this orders the genera based on logchange for plotting
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

## if the Genus is empty, replace with the Family
sigtab$Genus[is.na(sigtab$Genus)] <- sigtab$Family[is.na(sigtab$Genus)]

# bind Genus and Species together
sigtab$Genus.species <- paste(sigtab$Genus, sigtab$Species)


library("ggplot2")

## graph differential abundance
ggplot(sigtab, aes(y=Genus.species, x=log2FoldChange, color=Genus.species)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) 

#### DESEQ - Tide cycle ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_ME_F_1000, ~ TideCycle)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
# Spring and Neap pairwise comparison
res.Spring.Neap = results(diagdds, contrast = c("TideCycle", "Spring", "Neap"))
res.Spring.Neap = res.Spring.Neap[order(res.Spring.Neap$padj, na.last=NA), ]
alpha = 0.01
sigtab.Spring.Neap = res.Spring.Neap[(res.Spring.Neap$padj < alpha), ]
sigtab.Spring.Neap = cbind(as(sigtab.Spring.Neap, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Spring.Neap), ], "matrix")) 

head(sigtab.Spring.Neap)
dim(sigtab.Spring.Neap)
# 4 differentially abundant taxa between Spring and Neap tide cycles

write.csv(sigtab.Spring.Neap, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.Spring.Neap.csv")

# Spring and Intermediate pairwise comparison
res.Spring.Int = results(diagdds, contrast = c("TideCycle", "Spring", "Intermediate"))
res.Spring.Int = res.Spring.Int[order(res.Spring.Int$padj, na.last=NA), ]
alpha = 0.01
sigtab.Spring.Int = res.Spring.Int[(res.Spring.Int$padj < alpha), ]
sigtab.Spring.Int = cbind(as(sigtab.Spring.Int, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Spring.Int), ], "matrix")) 

head(sigtab.Spring.Int)
dim(sigtab.Spring.Int)
# 2 differentially abundant taxa between Spring and Intermediate tide cycles

write.csv(sigtab.Spring.Int, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.Spring.Int.csv")

# Neap and Intermediate pairwise comparison
res.Neap.Int = results(diagdds, contrast = c("TideCycle", "Neap", "Intermediate"))
res.Neap.Int = res.Neap.Int[order(res.Neap.Int$padj, na.last=NA), ]
alpha = 0.01
sigtab.Neap.Int = res.Neap.Int[(res.Neap.Int$padj < alpha), ]
sigtab.Neap.Int = cbind(as(sigtab.Neap.Int, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Neap.Int), ], "matrix")) 
head(sigtab.Neap.Int)
dim(sigtab.Neap.Int)
# 5 differentially abundant taxa between Neap and Intermediate tide cycles

write.csv(sigtab.Neap.Int, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.Neap.Int.csv")

##### DESEQ - Marsh ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_ME_F_1000, ~ Marsh)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
# doing for each marsh comparison
# HA versus PL
res.HA.PL = results(diagdds, contrast =c("Marsh", "HA", "PL"))
res.HA.PL = res.HA.PL[order(res.HA.PL$padj, na.last=NA), ]
alpha = 0.01
sigtab.HA.PL = res.HA.PL[(res.HA.PL$padj < alpha), ]
sigtab.HA.PL = cbind(as(sigtab.HA.PL, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.HA.PL), ], "matrix")) #

head(sigtab.HA.PL)
dim(sigtab.HA.PL)
# 10 differentially abundant SVs

write.csv(sigtab.HA.PL, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.HA.PL.csv")

# HA versus NG
res.HA.NG = results(diagdds, contrast =c("Marsh", "HA", "NG"))
res.HA.NG = res.HA.NG[order(res.HA.NG$padj, na.last=NA), ]
alpha = 0.01
sigtab.HA.NG = res.HA.NG[(res.HA.NG$padj < alpha), ]
sigtab.HA.NG = cbind(as(sigtab.HA.NG, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.HA.NG), ], "matrix")) 

head(sigtab.HA.NG)
dim(sigtab.HA.NG)
# 5 differentially abundant SVs

write.csv(sigtab.HA.NG, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.HA.NG.csv")

# HA versus BB
res.HA.BB = results(diagdds, contrast =c("Marsh", "HA", "BB"))
res.HA.BB = res.HA.BB[order(res.HA.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.HA.BB = res.HA.BB[(res.HA.BB$padj < alpha), ]
sigtab.HA.BB = cbind(as(sigtab.HA.BB, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.HA.BB), ], "matrix")) 

head(sigtab.HA.BB)
dim(sigtab.HA.BB)
# 7 Differentially abundant SVs

write.csv(sigtab.HA.BB, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.HA.BB.csv")

## NG versus PL
res.NG.PL = results(diagdds, contrast =c("Marsh", "NG", "PL"))
res.NG.PL = res.NG.PL[order(res.NG.PL$padj, na.last=NA), ]
alpha = 0.01
sigtab.NG.PL= res.NG.PL[(res.NG.PL$padj < alpha), ]
sigtab.NG.PL = cbind(as(sigtab.NG.PL, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.NG.PL), ], "matrix")) 

head(sigtab.NG.PL)
dim(sigtab.NG.PL)
# 9 differentially abundant SVs

write.csv(sigtab.NG.PL, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.NG.PL.csv")

# Plotting differentially abundant taxa
# Phylum order
x = tapply(sigtab.NG.PL$log2FoldChange, sigtab.NG.PL$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.NG.PL$Phylum = factor(as.character(sigtab.NG.PL$Phylum), levels=names(x))
# Ascomycota

# Genus order - this orders the genera based on logchange for plotting
x = tapply(sigtab.NG.PL$log2FoldChange, sigtab.NG.PL$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.NG.PL$Genus = factor(as.character(sigtab.NG.PL$Genus), levels=names(x))

## if the Genus is empty, replace with the Family
sigtab.NG.PL$Genus[is.na(sigtab.NG.PL$Genus)] <- sigtab.NG.PL$Family[is.na(sigtab.NG.PL$Genus)]

# bind Genus and Species together
sigtab.NG.PL$Genus.species <- paste(sigtab.NG.PL$Genus, sigtab.NG.PL$Species)

library("ggplot2")

## graph differential abundance
DESEQ.NG.PL <- ggplot(sigtab.NG.PL, aes(y=Phylum, x=log2FoldChange, color=Phylum)) + #have to use color for Phylum becauase no taxonomy assigned below that level
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) 

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/DESEQ.NG.PL.png", width = 40, height = 30, units = 'cm', res = 300)
grid.arrange(DESEQ.NG.PL) 
dev.off()

## NG versus BB
res.NG.BB = results(diagdds, contrast =c("Marsh", "NG", "BB"))
res.NG.BB = res.NG.BB[order(res.NG.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.NG.BB = res.NG.BB[(res.NG.BB$padj < alpha), ]
sigtab.NG.BB = cbind(as(sigtab.NG.BB, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.NG.BB), ], "matrix")) 

head(sigtab.NG.BB)
dim(sigtab.NG.BB)
# 16 differentially abundant SVs 

write.csv(sigtab.NG.BB, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.NG.BB.csv")

## PL versus BB
res.PL.BB = results(diagdds, contrast =c("Marsh", "PL", "BB"))
res.PL.BB = res.PL.BB[order(res.PL.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.PL.BB = res.PL.BB[(res.PL.BB$padj < alpha), ]
sigtab.PL.BB = cbind(as(sigtab.PL.BB, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.PL.BB), ], "matrix")) 

head(sigtab.PL.BB)
dim(sigtab.PL.BB)
# 0 differentially abundant SVs between PL and BB


#### DESEQ - Month #####

# Jun and August
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_ME_F_1000, ~ Month)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res.Jun.Aug = results(diagdds, contrast = c("Month", "June", "August"))
res.Jun.Aug = res.Jun.Aug[order(res.Jun.Aug$padj, na.last=NA), ]
alpha = 0.01
sigtab.Jun.Aug = res.Jun.Aug[(res.Jun.Aug$padj < alpha), ]
sigtab.Jun.Aug = cbind(as(sigtab.Jun.Aug, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Jun.Aug), ], "matrix")) 

head(sigtab.Jun.Aug)
dim(sigtab.Jun.Aug)
# 0 differentially abundant taxa between Jun and Aug

# Jun and July
# calculate significance for those abundance calculations
res.Jun.Jul = results(diagdds, contrast = c("Month", "June", "July"))
res.Jun.Jul = res.Jun.Jul[order(res.Jun.Jul$padj, na.last=NA), ]
alpha = 0.01
sigtab.Jun.Jul = res.Jun.Jul[(res.Jun.Jul$padj < alpha), ]
sigtab.Jun.Jul = cbind(as(sigtab.Jun.Jul, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Jun.Jul), ], "matrix")) 

head(sigtab.Jun.Jul)
dim(sigtab.Jun.Jul)
# 14 differentially abundant taxa between Jun and July

write.csv(sigtab.Jun.Jul, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.Jun.Jul.csv")

# July and August
# calculate significance for those abundance calculations
res.Aug.Jul = results(diagdds, contrast = c("Month", "August", "July"))
res.Aug.Jul = res.Aug.Jul[order(res.Aug.Jul$padj, na.last=NA), ]
alpha = 0.01
sigtab.Aug.Jul = res.Aug.Jul[(res.Aug.Jul$padj < alpha), ]
sigtab.Aug.Jul = cbind(as(sigtab.Aug.Jul, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.Aug.Jul), ], "matrix")) 

head(sigtab.Aug.Jul)
dim(sigtab.Aug.Jul)
# 0 differentially abundant taxa between August and July

#### DESEQ - Tidal Range ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_ME_F_1000, ~ TidalRange_m)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
# 2.5 to 3 versus 3 to 3.5
res.2.5_to_3.3_to_3.5 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "3_to_3.5"))
res.2.5_to_3.3_to_3.5 = res.2.5_to_3.3_to_3.5[order(res.2.5_to_3.3_to_3.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.3_to_3.5 = res.2.5_to_3.3_to_3.5[(res.2.5_to_3.3_to_3.5$padj < alpha), ]
sigtab.2.5_to_3.3_to_3.5 = cbind(as(sigtab.2.5_to_3.3_to_3.5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.2.5_to_3.3_to_3.5), ], "matrix")) 

head(sigtab.2.5_to_3.3_to_3.5)

dim(sigtab.2.5_to_3.3_to_3.5)
# 11 differentially abundant SVs

write.csv(sigtab.2.5_to_3.3_to_3.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.2.5_to_3.3_to_3.5.csv")

# 2.5 to 3 versus 3.5 to 4
res.2.5_to_3.3.5_to_4 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "3.5_to_4"))
res.2.5_to_3.3.5_to_4 = res.2.5_to_3.3.5_to_4[order(res.2.5_to_3.3.5_to_4$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.3.5_to_4 = res.2.5_to_3.3.5_to_4[(res.2.5_to_3.3.5_to_4$padj < alpha), ]
sigtab.2.5_to_3.3.5_to_4 = cbind(as(sigtab.2.5_to_3.3.5_to_4, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.2.5_to_3.3.5_to_4), ], "matrix")) 

head(sigtab.2.5_to_3.3.5_to_4)

dim(sigtab.2.5_to_3.3.5_to_4)
# 15 differentially abundant SVs

write.csv(sigtab.2.5_to_3.3.5_to_4, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.2.5_to_3.3.5_to_4.csv")

# 2.5 to 3 versus 4 to 4.5
res.2.5_to_3.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "4_to_4.5"))
res.2.5_to_3.4_to_4.5 = res.2.5_to_3.4_to_4.5[order(res.2.5_to_3.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.4_to_4.5 = res.2.5_to_3.3.5_to_4[(res.2.5_to_3.4_to_4.5$padj < alpha), ]
sigtab.2.5_to_3.4_to_4.5 = cbind(as(sigtab.2.5_to_3.4_to_4.5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.2.5_to_3.4_to_4.5), ], "matrix")) 

head(sigtab.2.5_to_3.4_to_4.5)

dim(sigtab.2.5_to_3.4_to_4.5)
# 7 differentially abundant SVs

write.csv(sigtab.2.5_to_3.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.2.5_to_3.4_to_4.5.csv")

# 2.5 to 3 versus 4.5 to 5
res.2.5_to_3.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "4.5_to_5"))
res.2.5_to_3.4.5_to_5 = res.2.5_to_3.4.5_to_5[order(res.2.5_to_3.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.4.5_to_5 = res.2.5_to_3.4.5_to_5[(res.2.5_to_3.4.5_to_5$padj < alpha), ]
sigtab.2.5_to_3.4.5_to_5 = cbind(as(sigtab.2.5_to_3.4.5_to_5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.2.5_to_3.4.5_to_5), ], "matrix")) 

head(sigtab.2.5_to_3.4.5_to_5)

dim(sigtab.2.5_to_3.4.5_to_5)
# 0 differentially abundant SVs

# 3 to 3.5 versus 3.5 to 4
res.3_to_3.5.3.5_to_4 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "3.5_to_4"))
res.3_to_3.5.3.5_to_4 = res.3_to_3.5.3.5_to_4[order(res.3_to_3.5.3.5_to_4$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.3.5_to_4 = res.3_to_3.5.3.5_to_4[(res.3_to_3.5.3.5_to_4$padj < alpha), ]
sigtab.3_to_3.5.3.5_to_4 = cbind(as(sigtab.3_to_3.5.3.5_to_4, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.3_to_3.5.3.5_to_4), ], "matrix")) 

head(sigtab.3_to_3.5.3.5_to_4)

dim(sigtab.3_to_3.5.3.5_to_4)
# 16 differentially abundant SVs

write.csv(sigtab.3_to_3.5.3.5_to_4, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.3_to_3.5.3.5_to_4.csv")

# 3 to 3.5 versus 4 to 4.5
res.3_to_3.5.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "4_to_4.5"))
res.3_to_3.5.4_to_4.5 = res.3_to_3.5.4_to_4.5[order(res.3_to_3.5.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.4_to_4.5 = res.3_to_3.5.4_to_4.5[(res.3_to_3.5.4_to_4.5$padj < alpha), ]
sigtab.3_to_3.5.4_to_4.5 = cbind(as(sigtab.3_to_3.5.4_to_4.5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.3_to_3.5.4_to_4.5), ], "matrix")) 

head(sigtab.3_to_3.5.4_to_4.5)

dim(sigtab.3_to_3.5.4_to_4.5)
# 19 differentially abundant SVs

write.csv(sigtab.3_to_3.5.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.3_to_3.5.4_to_4.5.csv")

# 3 to 3.5 versus 4.5 to 5
res.3_to_3.5.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "4.5_to_5"))
res.3_to_3.5.4.5_to_5 = res.3_to_3.5.4.5_to_5[order(res.3_to_3.5.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.4.5_to_5 = res.3_to_3.5.4.5_to_5[(res.3_to_3.5.4.5_to_5$padj < alpha), ]
sigtab.3_to_3.5.4.5_to_5 = cbind(as(sigtab.3_to_3.5.4.5_to_5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.3_to_3.5.4.5_to_5), ], "matrix")) 

head(sigtab.3_to_3.5.4.5_to_5)

dim(sigtab.3_to_3.5.4.5_to_5)
# 3 differentially abundant SVs

write.csv(sigtab.3_to_3.5.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.3_to_3.5.4.5_to_5.csv")

# 3.5 to 4 versus 4 to 4.5
res.3.5_to_4.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "3.5_to_4", "4_to_4.5"))
res.3.5_to_4.4_to_4.5 = res.3.5_to_4.4_to_4.5[order(res.3.5_to_4.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3.5_to_4.4_to_4.5 = res.3.5_to_4.4_to_4.5[(res.3.5_to_4.4_to_4.5$padj < alpha), ]
sigtab.3.5_to_4.4_to_4.5 = cbind(as(sigtab.3.5_to_4.4_to_4.5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.3.5_to_4.4_to_4.5), ], "matrix")) 

head(sigtab.3.5_to_4.4_to_4.5)

dim(sigtab.3.5_to_4.4_to_4.5)
# 18 differentially abundant SVs

write.csv(sigtab.3.5_to_4.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.3.5_to_4.4_to_4.5.csv")

# 3.5 to 4 versus 4.5 to 5
res.3.5_to_4.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "3.5_to_4", "4.5_to_5"))
res.3.5_to_4.4.5_to_5 = res.3.5_to_4.4.5_to_5[order(res.3.5_to_4.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3.5_to_4.4.5_to_5 = res.3.5_to_4.4.5_to_5[(res.3.5_to_4.4.5_to_5$padj < alpha), ]
sigtab.3.5_to_4.4.5_to_5 = cbind(as(sigtab.3.5_to_4.4.5_to_5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.3.5_to_4.4.5_to_5), ], "matrix")) 

head(sigtab.3.5_to_4.4.5_to_5)

dim(sigtab.3.5_to_4.4.5_to_5)
# 13 differentially abundant SVs

write.csv(sigtab.3.5_to_4.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.3.5_to_4.4.5_to_5.csv")

# 4 to 4.5 versus 4.5 to 5
res.4_to_4.5.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "4_to_4.5", "4.5_to_5"))
res.4_to_4.5.4.5_to_5 = res.4_to_4.5.4.5_to_5[order(res.4_to_4.5.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.4_to_4.5.4.5_to_5 = res.4_to_4.5.4.5_to_5[(res.4_to_4.5.4.5_to_5$padj < alpha), ]
sigtab.4_to_4.5.4.5_to_5 = cbind(as(sigtab.4_to_4.5.4.5_to_5, "data.frame"), as(tax_table(ITS_ME_F_1000)[rownames(sigtab.4_to_4.5.4.5_to_5), ], "matrix")) 

head(sigtab.4_to_4.5.4.5_to_5)

dim(sigtab.4_to_4.5.4.5_to_5)
# 12 differentially abundant SVs

write.csv(sigtab.4_to_4.5.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.4_to_4.5.4.5_to_5.csv")


## KERATINOLYTIC taxa in ME feathers #####

K_degrade <- c("Aspergillus", "Trychophyton", "Trichophyton", "Engyodontium", "Doratomyces", "Paecilomyces", "Onygena", "Pseudogymnoascus", "Cadophora", "Neosetophoma", "Auxarthron", "Gynmoascus", "Microsporum", "Scopulariopsis", "Sepedonium", "Myriodontium", "Aphanoascus", "Arthroderma", "Chrysosporium", "Coccidoides", "Gymnoascoideus", "Talaromyces", "Myrothecium", "Tritirachium", "Trichoderma", "Candida", "Geotrichum", "Fusarium", "Penicillium", "Oidiodendron", "Verticillium", "Keratinomyces", "Alternaria", "Trichurus", "Curvularia", "Cladosporium", "Geomyces", "Gleomastis", "Monodictys", "Myrothecium", "Stachybotrys", "Urocladium", "Neosetophoma", "Pseudogymnoascus", "Aphanoascus", "Acrodontium", "Botryotricum", "Chaetomium", "Chrysosporium", "Gliocladium", "Keratinophyton", "Malbranhea", "Microsporum", "Phytophthora")

ITS_K_ME.F <- subset_taxa(ITS_ME_F_1000.rar, Genus %in% K_degrade)
# 109 taxa and 33

# Looking at taxa in the keratinolytic phyloseq
sort(table(tax_table(ITS_K_ME.F)[, 2])) # all in Phylum ascymycota
sort(table(tax_table(ITS_K_ME.F)[, 3]))
sort(table(tax_table(ITS_K_ME.F)[, 4]))
sort(table(tax_table(ITS_K_ME.F)[, 5]))
sort(table(tax_table(ITS_K_ME.F)[, 6]))

# 16 genera
# Alternaria (42 SVs), Cladosporium (12 SVs), Fusarium (11 SVs), Candida (11 SVs), Penicillium (9 SVs), Aspergillus (6 SVs), Talaromyces (4 SVs), Stachybotrys (3 SVs), Oidodendron (2 SVs), Curvularia (2 SVs), Trichoderma (1 SV), Neosetophoma (1 SV), Engyodontium (1 SV), Cadophora (1 SV), Acrodontium (1 SV)

ITS_K_ME.F <- prune_samples(sample_sums(ITS_K_ME.F)>0, ITS_K_ME.F)
ITS_K_ME.F <- prune_taxa(taxa_sums(ITS_K_ME.F)> 0, ITS_K_ME.F)
# 109 taxa and 33 samples

saveRDS(ITS_K_ME.F, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K_ME.F.rds")

ITS_K_ME.F <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K_ME.F.rds")

### ALPHA DIV: Keratinolytic taxa by Maine variables ######

#### Testing differences in alpha diversity - Keratinolytic taxa in ME NESP feathers #####

ITS_K_ME.F.rich <- estimate_richness(ITS_K_ME.F, measure=c("Observed", "Shannon")) 
ITS_K_ME.F.rich.even <- ITS_K_ME.F.rich$Shannon/log(ITS_K_ME.F.rich$Observed)
ITS_K_ME.F.rich.sd = as(sample_data(ITS_K_ME.F), "matrix")
ITS_K_ME.F.rich.sd = as.data.frame(ITS_K_ME.F.rich.sd)
ITS_K_ME.F.rich.df <- cbind(ITS_K_ME.F.rich, ITS_K_ME.F.rich.even, ITS_K_ME.F.rich.sd)
write.csv(ITS_K_ME.F.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K_ME.F.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals, 3) run kruskall-wallis if not sig

# Effect of sex on observed diversity of keratinolytic SVs
a_K_ME.F_sex.obs <- aov(Observed ~ Sex, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_sex.obs)
summary(a_K_ME.F_sex.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Sex          1   10.9   10.88   0.676  0.417
# Residuals   31  499.4   16.11  

res.K_ME.F_sex.obs  <- a_K_ME.F_sex.obs$residuals
shapiro.test(res.K_ME.F_sex.obs)
# W = 0.74465, p-value = 3.339e-06 # Not normal - run kruskall wallis

kruskal.test(Observed ~ Sex, data = ITS_K_ME.F.rich.df)
# Kruskal-Wallis chi-squared = 0.0041121, df = 1, p-value = 0.9489
# No difference in alpha diversity of keratinolytic SVs across sex

# Effect of sex on shannon diversity
a_K_ME.F_sex.shan <- aov(Shannon ~ Sex, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_sex.shan)
summary(a_K_ME.F_sex.shan)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Sex          1   0.13  0.1297   0.316  0.578
# Residuals   31  12.72  0.4105       
# No effect of sex on shannon diveristy

res.K_ME.F_sex.shan <- a_K_ME.F_sex.shan$residuals
shapiro.test(res.K_ME.F_sex.shan)
# W = 0.97021, p-value = 0.486 residuals normally distributed, report ANOVA results

# Effect of marsh on observed diversity
a_K_ME.F_marsh.obs <- aov(Observed ~ Marsh, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_marsh.obs)
summary(a_K_ME.F_marsh.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Marsh        3   19.0    6.34   0.374  0.772
# Residuals   29  491.2   16.94 

res.K_ME.F_marsh.obs  <- a_K_ME.F_marsh.obs$residuals
shapiro.test(res.K_ME.F_marsh.obs)
# W = 0.76953, p-value = 8.913e-06 # Not normal - run kruskall wallis

kruskal.test(Observed ~ Marsh, data = ITS_K_ME.F.rich.df)
# Kruskal-Wallis chi-squared = 0.27737, df = 3, p-value = 0.9642
# No difference in alpha diversity of keratinolytic SVs across marshes

# Effect of marsh on shannon diversity
a_K_ME.F_marsh.shan <- aov(Shannon ~ Marsh, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_marsh.shan)
summary(a_K_ME.F_marsh.shan)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Marsh        3   0.39  0.1300   0.303  0.823
# Residuals   29  12.46  0.4298   

res.K_ME.F_marsh.shan <- a_K_ME.F_marsh.shan$residuals
shapiro.test(res.K_ME.F_marsh.shan)
# W = 0.96111, p-value = 0.2777 Normal, report ANOVA results

# Effect of month on observed diversity
a_K_ME.F_month.obs <- aov(Observed ~ Month, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_month.obs)
summary(a_K_ME.F_month.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Month        2   30.4   15.19    0.95  0.398
# Residuals   30  479.9   15.99   

res.K_ME.F_month.obs  <- a_K_ME.F_month.obs$residuals
shapiro.test(res.K_ME.F_month.obs)
# W = 0.72666, p-value = 1.698e-06 # Not normal - run kruskall wallis

kruskal.test(Observed ~ Month, data = ITS_K_ME.F.rich.df)
# Kruskal-Wallis chi-squared = 4.2988, df = 2, p-value = 0.1166
# No difference in alpha diversity of keratinolytic SVs across months

# Effect of month on shannon diversity
a_K_ME.F_month.shan <- aov(Shannon ~ Month, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_month.shan)
summary(a_K_ME.F_month.shan)

#              Df Sum Sq Mean Sq F value Pr(>F)
# Month        2  0.861  0.4304   1.077  0.354
# Residuals   30 11.993  0.3998       
# No effect of tide cycle on shannon diveristy

res.K_ME.F_month.shan <- a_K_ME.F_month.shan$residuals
shapiro.test(res.K_ME.F_month.shan)
# W = 0.97107, p-value = 0.5103 Normal, report ANOVA results

# Effect of tide cycle on observed diversity
a_K_ME.F_cycle.obs <- aov(Observed ~ TideCycle, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_cycle.obs)
summary(a_K_ME.F_cycle.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# TideCycle    2   18.7   9.325   0.569  0.572
# Residuals   30  491.6  16.386   

res.K_ME.F_cycle.obs  <- a_K_ME.F_cycle.obs $residuals
shapiro.test(res.K_ME.F_cycle.obs)
# W = 0.71597, p-value = 1.15e-06 # Not normal - run kruskall wallis

kruskal.test(Observed ~ TideCycle, data = ITS_K_ME.F.rich.df)
# Kruskal-Wallis chi-squared = 2.5228, df = 2, p-value = 0.2833
# No difference in alpha diversity of keratinolytic SVs across tide cycle

# Effect of tide cycle on shannon diversity
a_K_ME.F_cycle.shan <- aov(Shannon ~ TideCycle, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_cycle.shan)
summary(a_K_ME.F_cycle.shan)

#              Df Sum Sq Mean Sq F value Pr(>F)
# TideCycle    2  0.909  0.4547   1.142  0.333
# Residuals   30 11.945  0.3982       
# No effect of tide cycle on shannon diveristy

res.K_ME.F_cycle.shan <- a_K_ME.F_cycle.shan$residuals
shapiro.test(res.K_ME.F_cycle.shan)
# W = 0.9717, p-value = 0.5285 Normal, report ANOVA results

# Effect of tidal range on observed diversity of keratinolytic SVs
a_K_ME.F_rng.obs <- aov(Observed ~ TidalRange, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_rng.obs)
summary(a_K_ME.F_rng.obs)

#              Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange   3   25.5   8.505   0.509  0.679
# Residuals   29  484.7  16.715    

res.K_ME.F_rng.obs  <- a_K_ME.F_rng.obs$residuals
shapiro.test(res.K_ME.F_rng.obs)
# W = 0.76075, p-value = 6.262e-06 - Not normal, run kruskall-wallace

kruskal.test(Observed ~ TidalRange, data = ITS_K_ME.F.rich.df)
# Kruskal-Wallis chi-squared = 2.4952, df = 3, p-value = 0.4762
# No effect of tidal range on observed diversity

# Effect of tidal range on shannon diversity of keratinolytic SVs
a_K_ME.F_rng.shan <- aov(Shannon ~ TidalRange, data = ITS_K_ME.F.rich.df)
plot(a_K_ME.F_rng.shan)
summary(a_K_ME.F_rng.shan)

#              Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange   3  0.304  0.1013   0.234  0.872
# Residuals   29 12.550  0.4328   

res.K_ME.F_rng.shan  <- a_K_ME.F_rng.shan$residuals
shapiro.test(res.K_ME.F_rng.shan)
# W = 0.96418, p-value = 0.3377. Results are normal, can report ANOVA results

### BETA DIV: KERATINOLYTIC SVS IN MAINE NESP ####

#### HOMOGENEITY OF VARIANCE AND COMMUNITY COMPOSITION - Keratinolyitc SVs in Maine feathers ####

meta.df.ME.K.its <- as(sample_data(ITS_K_ME.F), "data.frame")

##### Jaccard distance between samples, for testing dispersion #####
dis.ME.K.its <- distance(ITS_K_ME.F, method="jaccard")

# checking betadispersion by month
disper.K.its.month <- betadisper(dis.ME.K.its, meta.df.ME.K.its$Month)
disper.K.its.month
permutest(disper.K.its.month)
#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.041761 0.0208804 2.7329    999  0.085 .
# Residuals 30 0.229216 0.0076405  
# Month doesn't have sig diff dispersions

# checking betadispersion by marsh
disper.K.its.marsh <- betadisper(dis.ME.K.its, meta.df.ME.K.its$Month)
disper.K.its.marsh
permutest(disper.K.its.marsh)
#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.041761 0.0208804 2.7329    999  0.075 .
# Residuals 30 0.229216 0.0076405  
# Marsh doesn't have sig diff dispersions

# checking betadispersion by sex
disper.K.its.sex <- betadisper(dis.ME.K.its, meta.df.ME.K.its$Sex)
disper.K.its.sex
permutest(disper.K.its.sex)
#             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.06039 0.060393 5.1487    999  0.031 *
# Residuals 31 0.36363 0.011730     
# Sex has sig diff dispersions

# checking betadispersion by tide cycle
disper.K.its.cycle <- betadisper(dis.ME.K.its, meta.df.ME.K.its$TideCycle)
disper.K.its.cylce
permutest(disper.K.its.cycle)
#             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.04849 0.024245 2.0023    999  0.157
# Residuals 30 0.36326 0.012109      
# Tide cycle doesn't have sig diff dispersions

# checking betadispersion by tidal range
disper.K.its.range <- betadisper(dis.ME.K.its, meta.df.ME.K.its$TidalRange)
disper.K.its.range
permutest(disper.K.its.range)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     3 0.009753 0.0032511 0.4284    999  0.733
# Residuals 29 0.220064 0.0075884        
# Tidal range doesn't have sig diff dispersions

##### Bray Curtis disatance between samples, for testing dispersion ####
dis.ME.K.its.b <- distance(ITS_K_ME.F, method="bray")

# checking betadispersion by month
disper.K.its.month.b <- betadisper(dis.ME.K.its.b, meta.df.ME.K.its$Month)
disper.K.its.month.b
permutest(disper.K.its.month.b)
#             Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
# Groups     2 0.04187 0.020935 1.427    999  0.255
# Residuals 30 0.44012 0.014671  
# Month doesn't have sig diff dispersions

# checking betadispersion by marsh
disper.K.its.marsh.b <- betadisper(dis.ME.K.its.b, meta.df.ME.K.its$Month)
disper.K.its.marsh.b
permutest(disper.K.its.marsh.b)
#              Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
# Groups       2 0.04187 0.020935 1.427    999  0.252
# Residuals    30 0.44012 0.014671 
# Marsh doesn't have sig diff dispersions

# checking betadispersion by sex
disper.K.its.sex.b <- betadisper(dis.ME.K.its.b, meta.df.ME.K.its$Sex)
disper.K.its.sex.b
permutest(disper.K.its.sex.b)
#           Df  Sum Sq  Mean Sq    F    N.Perm Pr(>F)  
# Groups     1 0.09168 0.091678 4.0221    999  0.039 *
# Residuals 31 0.70661 0.022794     
# Sex has sig diff dispersions

# checking betadispersion by tide cycle
disper.K.its.cycle.b <- betadisper(dis.ME.K.its.b, meta.df.ME.K.its$TideCycle)
disper.K.its.cylce.b
permutest(disper.K.its.cycle.b)
#             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.08298 0.041488 1.5925    999  0.237
# Residuals  30 0.78155 0.026052       
# Tide cycle doesn't have sig diff dispersions

# checking betadispersion by tidal range
disper.K.its.range.b <- betadisper(dis.ME.K.its.b, meta.df.ME.K.its$TidalRange)
disper.K.its.range.b
permutest(disper.K.its.range.b)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     3 0.01440 0.0047986 0.3006    999  0.817
# Residuals 29 0.46297 0.0159646       
# Tidal range doesn't have sig diff dispersions

### PERMANOVA ####

#### Jaccard dissimilarity #####
# Marsh
perma.ME.K.marsh <- adonis2(distance(ITS_K_ME.F, method="jaccard") ~ Marsh,
                             data = meta.df.ME.K.its)

#           Df SumOfSqs    R2     F   Pr(>F)
# Marsh     3   1.3722 0.09171 0.976  0.532
# Residual 29  13.5899 0.90829             
# Total    32  14.9621 1.00000 

# Month
perma.ME.K.month <- adonis2(distance(ITS_K_ME.F, method="jaccard") ~ Month,
data = meta.df.ME.K.its)

#.         Df SumOfSqs    R2      F    Pr(>F)
# Month     2   0.8502 0.05682 0.9037  0.692
# Residual 30  14.1119 0.94318              
# Total    32  14.9621 1.00000 

# Tide Cycle
perma.ME.K.cycle <- adonis2(distance(ITS_K_ME.F, method="jaccard") ~ TideCycle,
                            data = meta.df.ME.K.its)

#.           Df SumOfSqs    R2      F Pr(>F)
# TideCycle  2   1.0889 0.07278 1.1774  0.133
# Residual  30  13.8732 0.92722              
# Total     32  14.9621 1.00000   

# Tidal Range
perma.ME.K.range <- adonis2(distance(ITS_K_ME.F, method="jaccard") ~ TidalRange,
                            data = meta.df.ME.K.its)

#.           Df SumOfSqs      R2      F Pr(>F)
# TidalRange  3   1.2809 0.08561 0.9051  0.743
# Residual   29  13.6812 0.91439              
# Total      32  14.9621 1.00000 

# Sex
perma.ME.K.sex <- adonis2(distance(ITS_K_ME.F, method="jaccard") ~ Sex,
                            data = meta.df.ME.K.its)

#.          Df SumOfSqs   R2      F    Pr(>F)
# Sex       1   0.4883 0.03263 1.0458  0.303
# Residual 31  14.4738 0.96737              
# Total    32  14.9621 1.00000   


#### Bray-Curtis dissimilarity #####
# Marsh
perma.ME.K.marsh.b <- adonis2(distance(ITS_K_ME.F, method="bray") ~ Marsh,
                            data = meta.df.ME.K.its)

#           Df SumOfSqs   R2     F    Pr(>F)
# Marsh     3   1.3347 0.09137 0.972  0.549
# Residual 29  13.2734 0.90863             
# Total    32  14.6081 1.00000   

# Month
perma.ME.K.month.b <- adonis2(distance(ITS_K_ME.F, method="bray") ~ Month,
                            data = meta.df.ME.K.its)

#.         Df SumOfSqs      R2    F  Pr(>F)
# Month     2   0.7929 0.05428 0.8609  0.717
# Residual 30  13.8152 0.94572              
# Total    32  14.6081 1.00000   

# Tide Cycle
perma.ME.K.cycle.b <- adonis2(distance(ITS_K_ME.F, method="bray") ~ TideCycle,
                            data = meta.df.ME.K.its)

#.           Df SumOfSqs  R2      F     Pr(>F)
# TideCycle  2   1.1283 0.07724 1.2555   0.12
# Residual  30  13.4798 0.92276              
# Total     32  14.6081 1.00000   

# Tidal Range
perma.ME.K.range.b <- adonis2(distance(ITS_K_ME.F, method="bray") ~ TidalRange,
                            data = meta.df.ME.K.its)

#.           Df SumOfSqs      R2      F Pr(>F)
# TidalRange  3   1.2187 0.08343 0.8799  0.734
# Residual   29  13.3893 0.91657              
# Total      32  14.6081 1.00000 

# Sex
perma.ME.K.sex.b <- adonis2(distance(ITS_K_ME.F, method="bray") ~ Sex,
                          data = meta.df.ME.K.its)

#          Df SumOfSqs      R2      F Pr(>F)
# Sex       1   0.4376 0.02996 0.9573  0.456
# Residual 31  14.1705 0.97004              
# Total    32  14.6081 1.00000  


# AIM 3 - SEDIMENT AND FEATHER COMPARISON #####

## ALPHA DIVERSITY ####

### Violin plots comparing sediment and feather samples #####
viol.SF.rar <- plot_richness(ITS_SF_1000.rar, 
                             x="Sample_Type", 
                             measures=c("Observed","Shannon"), 
                             title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Sample_Type)) + 
  geom_boxplot(width = 0.1, aes(group=Sample_Type)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/viol.SF.rar.png", width = 21, height = 21, units = 'cm', res = 300)
grid.arrange(viol.SF.rar) 
dev.off()

### Test for alpha div variation ####

SF.rar.rich <- estimate_richness(ITS_SF_1000.rar, measure=c("Observed", "Shannon")) 
SF.rar.rich.even <- SF.rar.rich$Shannon/log(SF.rar.rich$Observed)
SF.rar.rich.sd = as(sample_data(ITS_SF_1000.rar), "matrix")
SF.rar.rich.sd = as.data.frame(SF.rar.rich.sd)
SF.rar.rich.df <- cbind(SF.rar.rich, SF.rar.rich.even, SF.rar.rich.sd)
write.csv(SF.rar.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/SF.rar.rich.df.csv")


# Observed diversity by sample type
a_SvF.obs <- aov(Observed ~ Sample_Type, data = SF.rar.rich.df)
plot(a_SvF.obs)
summary(a_SvF.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)  
# Sample_Type  1  15409   15409    3.11 0.0817 .
# Residuals   79 391397    4954                   

res.SvF.ITS.obs <- a_SvF.obs$residuals
shapiro.test(res.SvF.ITS.obs)
# W = 0.95671, p-value = 0.007976
# non-normal residuals

# Kruskall-wallis
kruskal.test(Observed ~ Sample_Type, data = SF.rar.rich.df)
# chi-squared = 3.0415, df = 1, p-value = 0.08116
# Observed diversity not significantly different in sample types

# Shannon diversity by sample type
a_SvF.shan <- aov(Shannon ~ Sample_Type, data = SF.rar.rich.df)
plot(a_SvF.shan)
summary(a_SvF.shan)

#             Df Sum Sq Mean Sq F value Pr(>F)
# Sample_Type  1   1.33   1.327    0.97  0.328
# Residuals   79 108.04   1.368                   

res.SvF.ITS.shan <- a_SvF.shan$residuals
shapiro.test(res.SvF.ITS.shan)
# W = 0.88819, p-value = 3.55e-06
# non-normal residuals

# Kruskall-wallis
kruskal.test(Shannon ~ Sample_Type, data = SF.rar.rich.df)
# Kruskal-Wallis chi-squared = 3.3371, df = 1, p-value = 0.06773
# Shannon diversity not significantly different in sample types


## BETA DIVERSITY #####
### PCOA - SEDIMENT v FEATHERS ####

# Sample type - jaccard
ord_SF.rar <- ordinate(ITS_SF_1000.rar, #calculate similarities
                       method ="PCoA", #ordination type
                       "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF <- plot_ordination(ITS_SF_1000.rar, ord_SF.rar, type="samples", color="Sample_Type") + geom_point(size = 3)
# Don't look too dissimilar

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.SF.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF) 
dev.off()

# Bray Curtis
ord_SF.rar.b <- ordinate(ITS_SF_1000.rar, #calculate similarities
                       method ="PCoA", #ordination type
                       "bray") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF.b <- plot_ordination(ITS_SF_1000.rar, ord_SF.rar.b, type="samples", color="Sample_Type") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.SF.bray.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF) 
dev.off()


### DISPERSION & COMMUNITY COMP ####

#  Get sample data
meta.df.SF.its <- as(sample_data(ITS_SF_1000.rar), "data.frame")

#### HOMOGENEITY OF VARIANCE - Jaccard  ####
dis.SF.its <- distance(ITS_SF_1000.rar, method="jaccard")

# checking betadispersion by sample type

disper.SF.its <- betadisper(dis.SF.its, meta.df.SF.its$Sample_Type)
disper.SF.its
permutest(disper.SF.its)

#             Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0034206 0.0034206 36.035    999  0.001 ***
# Residuals 79 0.0074990 0.0000949   
# Sample types have different dispersions - can't rely on the PERMANOVA results

# PERMANOVA

perma.SF <- adonis2(distance(ITS_SF_1000.rar, method="jaccard") ~ Sample_Type,
                    data = meta.df.SF.its)

#             Df SumOfSqs      R2     F Pr(>F)  
# Sample_Type  1    0.589 0.01504 1.206  0.012 *
# Residual    79   38.605 0.98496               
# Total       80   39.194 1.00000   
# Sig difference in community comp, but paired with different dispersions

#### HOMOGENEITY OF VARIANCE - Bray-Curtis  ####

# get metadata
meta.df.SF.its <- as(sample_data(ITS_SF_1000.rar), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.SF.its.b <- distance(ITS_SF_1000.rar, method="bray")

# checking betadispersion by sample type

disper.SF.its.b <- betadisper(dis.SF.its.b, meta.df.SF.its$Sample_Type)
disper.SF.its.b
permutest(disper.SF.its.b)

#            Df   Sum Sq   Mean Sq   F    N.Perm Pr(>F)    
# Groups     1 0.003258 0.0032580 11.117    999  0.001 ***
# Residuals 79 0.023152 0.0002931   
# Sample types have different dispersions - can't rely on the PERMANOVA results

# PERMANOVA

perma.SF <- adonis2(distance(ITS_SF_1000.rar, method="bray") ~ Sample_Type,
                    data = meta.df.SF.its)

#             Df SumOfSqs    R2      F    Pr(>F)   
# Sample_Type  1    0.650 0.01687 1.3552  0.003 **
# Residual    79   37.888 0.98313                 
# Total       80   38.538 1.00000    


### DESEQ ####
library("DESeq2")
packageVersion("DESeq2")

# DESEQ Differential abundance - ME soil and ME feathers
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ITS_SF_1000, ~ Sample_Type) 

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab.FvS = res[(res$padj < alpha), ]
sigtab.FvS = cbind(as(sigtab, "data.frame"), as(tax_table(ITS_SF_1000)[rownames(sigtab), ], "matrix")) 

head(sigtab.FvS)
dim(sigtab.FvS)
# 14 differentially abundant taxa between feathers and sediment

write.csv(sigtab.FvS, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/DESEQ.sigtab.FvS.csv")

### Subsetting keratinolytic taxa in sediment ####

ITS_K.Sediment <- subset_taxa(ITS_S_1000.rar, Genus %in% K_degrade)
ITS_K.Sediment <- prune_samples(sample_sums(ITS_K.Sediment)>0, ITS_K.Sediment)
# 22 taxa and 9 samples - in about half of the sediment samples

saveRDS(ITS_K.Sediment, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.Sediment.rds")

ITS_K.Sediment <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_K.Sediment.rds")

sort(table(tax_table(ITS_K.Sediment)[, 6]))
# 10 Genera
# Alternaria (8 SVs), Talaromyces (4 SVs), Trichoderma (2 SVs), Penicillium (2 SVs), Oidiodendron (1 SV), Fusarium (1 SV), Cladosporium (1 SV), Candida (1 SV), Cadophora (1 SV), Aspergillus (1 SV)

## KERATINOLYTIC taxa in the sediment and feathers phyloseq ####

K_degrade <- c("Aspergillus", "Trychophyton", "Trichophyton", "Engyodontium", "Doratomyces", "Paecilomyces", "Onygena", "Pseudogymnoascus", "Cadophora", "Neosetophoma", "Auxarthron", "Gynmoascus", "Microsporum", "Scopulariopsis", "Sepedonium", "Myriodontium", "Aphanoascus", "Arthroderma", "Chrysosporium", "Coccidoides", "Gymnoascoideus", "Talaromyces", "Myrothecium", "Tritirachium", "Trichoderma", "Candida", "Geotrichum", "Fusarium", "Penicillium", "Oidiodendron", "Verticillium", "Keratinomyces", "Alternaria", "Trichurus", "Curvularia", "Cladosporium", "Geomyces", "Gleomastis", "Monodictys", "Myrothecium", "Stachybotrys", "Urocladium", "Neosetophoma", "Pseudogymnoascus", "Aphanoascus", "Acrodontium", "Botryotricum", "Chaetomium", "Chrysosporium", "Gliocladium", "Keratinophyton", "Malbranhea", "Microsporum", "Phytophthora")

ITS_SF_k <- subset_taxa(ITS_SF_1000.rar, Genus %in% K_degrade)
ITS_SF_k <- prune_taxa(taxa_sums(ITS_SF_k)>0, ITS_SF_k)
ITS_SF_k <- prune_samples(sample_sums(ITS_SF_k)>0, ITS_SF_k)
# 171 taxa and 67 samples

saveRDS(ITS_SF_k, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_SF_k.rds")

ITS_SF_k <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/ITS_SF_k.rds")


### Alpha div ####

#### Violin plots comparing sediment and feather samples #####
viol.SF.K <- plot_richness(ITS_SF_k, 
                             x="Sample_Type", 
                             measures=c("Observed","Shannon"), 
                             title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Sample_Type)) + 
  geom_boxplot(width = 0.1, aes(group=Sample_Type)) + 
  ylab("Observed Fungal Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/viol.SF.K.png", width = 21, height = 21, units = 'cm', res = 300)
grid.arrange(viol.SF.K) 
dev.off()


#### Test for alpha div variation ####

SF.K.rich <- estimate_richness(ITS_SF_k, measure=c("Observed", "Shannon")) 
SF.K.rich.even <- SF.K.rich$Shannon/log(SF.K.rich$Observed)
SF.K.rich.sd = as(sample_data(ITS_SF_k), "matrix")
SF.K.rich.sd = as.data.frame(SF.K.rich.sd)
SF.K.rich.df <- cbind(SF.K.rich, SF.K.rich.even, SF.K.rich.sd)
write.csv(SF.K.rich.df, file = "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/SF.K.rich.df.csv")

# Observed diversity by sample type
a_SvF.K.obs <- aov(Observed ~ Sample_Type, data = SF.K.rich.df)
plot(a_SvF.K.obs)
summary(a_SvF.K.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)
# Sample_Type  1   10.1   10.09   0.658   0.42
# Residuals   65  996.9   15.34                  

res.SvF.K.ITS.obs <- a_SvF.K.obs$residuals
shapiro.test(res.SvF.K.ITS.obs)
#W = 0.69684, p-value = 1.795e-10 non-normal residuals, run Kruskal wallis

kruskal.test(Observed ~ Sample_Type, data = SF.K.rich.df)
# Kruskal-Wallis chi-squared = 0.48922, df = 1, p-value = 0.4843Sample type not significanlty different

# Shannon diversity by sample type
a_SvF.K.shan <- aov(Shannon ~ Sample_Type, data = SF.K.rich.df)
plot(a_SvF.K.shan)
summary(a_SvF.K.shan)

#             Df Sum Sq Mean Sq F value Pr(>F)
# Sample_Type  1   1.33   1.327    0.97  0.328
# Residuals   79 108.04   1.368                     

res.SvF.K.ITS.shan <- a_SvF.K.shan$residuals
shapiro.test(res.SvF.K.ITS.shan)
# W = 0.93338, p-value = 0.001423 non-normal residuals, run Kruskal wallis

kruskal.test(Shannon ~ Sample_Type, data = SF.K.rich.df)
# Kruskal-Wallis chi-squared = 0.04147, df = 1, p-value = 0.8386 Sample type not significantly different

### BETA DIVERSITY ####

#### PCOA #####
# Sample type - jaccard
ord_ITS_SF_k <- ordinate(ITS_SF_k, #calculate similarities
                               method ="PCoA", #ordination type
                               "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF.k <- plot_ordination(ITS_SF_k, ord_ITS_SF_k, type="samples", color="Sample_Type") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.SF.k.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF.k) 
dev.off()

# Sample type - Bray curtis
ord_ITS_SF_k.b <- ordinate(ITS_SF_k, #calculate similarities
                         method ="PCoA", #ordination type
                         "bray") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF.k.b <- plot_ordination(ITS_SF_k, ord_ITS_SF_k.b, type="samples", color="Sample_Type") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/PCOA.SF.k.b.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF.k.b) 
dev.off()

### DISPERSION & COMP ####
#### HOMOGENEITY OF VARIANCE - JACCARD ####
# Get metadata
meta.df.SF.K.ITS <- as(sample_data(ITS_SF_k), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.SF.K.its <- distance(ITS_SF_k, method="jaccard")

# checking betadispersion by sample type

disper.SF.k.its <- betadisper(dis.SF.K.its, meta.df.SF.K.ITS$Sample_Type)
disper.SF.k.its
permutest(disper.SF.k.its)

#            Df   Sum Sq   Mean Sq    F   N.Perm Pr(>F)  
# Groups     1 0.006683 0.0066828 4.7276    999  0.031 *
# Residuals 65 0.091882 0.0014136 

#### HOMOGENEITY OF VARIANCE - BRAY CURTIS ####

# Get metadata
meta.df.SF.K.ITS <- as(sample_data(ITS_SF_k), "data.frame")

## Bray Curtis disatance between samples, for testing dispersion
dis.SF.K.its.b <- distance(ITS_SF_k, method="bray")

# checking betadispersion by sample type

disper.SF.k.its.b <- betadisper(dis.SF.K.its.b, meta.df.SF.K.ITS$Sample_Type)
disper.SF.k.its.b
permutest(disper.SF.k.its.b)

#           Df   Sum Sq   Mean Sq   F.    N.Perm Pr(>F)
# Groups     1 0.006754 0.0067543 2.7348    999   0.11
# Residuals 65 0.160534 0.0024698  
# Sample types don't differ significantly by bray curtis dissimilarity


#### PERMANOVA - JACCARD ####

perma.SF.k <- adonis2(distance(ITS_SF_k, method="jaccard") ~ Sample_Type,
                             data = meta.df.SF.K.ITS)

#            Df SumOfSqs      R2      F Pr(>F)
#Sample_Type  1    0.566 0.01777 1.1762  0.161
#Residual    65   31.269 0.98223              
#Total       66   31.835 1.00000      
# Sample types do not differ in community composition of keratinolytic SVs

#### PERMANOVA - BRAY CURTIS ####
perma.SF.k <- adonis2(distance(ITS_SF_k, method="bray") ~ Sample_Type,
                      data = meta.df.SF.K.ITS)

#           Df SumOfSqs      R2     F   Pr(>F)
# Sample_Type  1   0.5924 0.01888 1.251  0.147
# Residual    65  30.7814 0.98112             
# Total       66  31.3738 1.00000  
# Sample types do not differ in community composition of keratinolytic SVs

# CORE AND TRANSIENT TAXA ####

## All feathers ####
# Using unrarefied phyloseqs

# keep only taxa with positive sums
core_biomeF <- prune_taxa(taxa_sums(Phy_ITS_Feathers_1000) > 0, Phy_ITS_Feathers_1000)
# 7392 taxa and 65 samples 

#calcuate compositional version of the data (relative abundances)
core_biomeF.rel <- microbiome::transform(core_biomeF, "compositional")

#This returns the taxa that exceed the given prevalence and detection thresholds.
core_taxaF.60 <- core_members(core_biomeF.rel, detection = 1/10000, prevalence = 60/100) # 0 taxa at 60% prevalence
core_taxaF.50 <- core_members(core_biomeF.rel, detection = 1/10000, prevalence = 50/100) # 0 taxa at 50% prevalence
core_taxaF.40 <- core_members(core_biomeF.rel, detection = 1/10000, prevalence = 40/100) # 0 taxa at 40% prevalence
core_taxaF.30 <- core_members(core_biomeF.rel, detection = 1/10000, prevalence = 30/100) # 3 taxa at 30% prevalence
core_taxaF.20 <- core_members(core_biomeF.rel, detection = 1/10000, prevalence = 20/100) # 21 taxa at 20% prevalence
# No fungal taxa present in 60%, 50%, or 40% of samples, 3 SVs present in 30% of samples

# phyloseq of those in 20% of samples
phy_core_taxaF.20 <- core(core_biomeF.rel, detection = 1/10000, prevalence = .2)
phy_core_taxaF.20 <- prune_samples(sample_sums(phy_core_taxaF.20)>0, phy_core_taxaF.20)
phy_core_taxaF.20 <- prune_taxa(taxa_sums(phy_core_taxaF.20)>0, phy_core_taxaF.20)
# 21 taxa in 57 samples samples

write.csv(tax_table(phy_core_taxaF.20), "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/TaxTab_CoreF.20.csv")

# phyloseq of those in 30% of samples
phy_core_taxaF.30 <- core(core_biomeF.rel, detection = 1/10000, prevalence = .3)
phy_core_taxaF.30 <- prune_samples(sample_sums(phy_core_taxaF.30)>0, phy_core_taxaF.30)
phy_core_taxaF.30 <- prune_taxa(taxa_sums(phy_core_taxaF.30)>0, phy_core_taxaF.30)
# 3 taxa and 39 samples
# Feathers don't seem to have a large "core" fungal community

write.csv(tax_table(phy_core_taxaF.30), "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/TaxTab_CoreF.30.csv")

# rare feather taxa
# <20% prevalence threshold
phylo.rareF.2 <- rare(core_biomeF.rel, detection=1/10000, prevalence=20/100) # 0.1% detection threshold, <20% of samples 
phylo.rareF.2 <- prune_samples(sample_sums(phylo.rareF.2)>0, phylo.rareF.2)
phylo.rareF.2 <- prune_taxa(taxa_sums(phylo.rareF.2)>0, phylo.rareF.2)
# 7371 taxa and 65 samples


# <10% prevalence threshold
phylo.rareF.1 <- rare(core_biomeF.rel, detection=1/10000, prevalence=10/100) # 0.1% detection threshold, <20% of samples 
phylo.rareF.1 <- prune_samples(sample_sums(phylo.rareF.1)>0, phylo.rareF.1)
phylo.rareF.1 <- prune_taxa(taxa_sums(phylo.rareF.1)>0, phylo.rareF.1)
# 7254 taxa and 65 samples

## All sediment ####
# Using unrarefied phyloseqs

# keep only taxa with positive sums
core_biomeS <- prune_taxa(taxa_sums(ITS_S_1000) > 0, ITS_S_1000)
# 1843 taxa and 16 samples 

#calcuate compositional version of the data (relative abundances)
core_biomeS.rel <- microbiome::transform(core_biomeS, "compositional")

#This returns the taxa that exceed the given prevalence and detection thresholds.
core_taxaS.60 <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 60/100) # 0 taxa at 60% prevalence
core_taxaS.50 <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 50/100) # 0 taxa at 50% prevalence
core_taxaS.40 <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 40/100) # 0 taxa at 40% prevalence
core_taxaS.30 <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 30/100) # 14 taxa at 30% prevalence
core_taxaS.20 <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 20/100) # 33 taxa at 20% prevalence

phy_core_taxaS.20 <- core(core_biomeS.rel, detection = 1/10000, prevalence = .2)
phy_core_taxaS.20 <- prune_samples(sample_sums(phy_core_taxaS.20)>0, phy_core_taxaS.20)
phy_core_taxaS.20 <- prune_taxa(taxa_sums(phy_core_taxaS.20)>0, phy_core_taxaS.20)
# 33 taxa and 14 samples

write.csv(tax_table(phy_core_taxaS.20), "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/TaxTab_CoreS.20.csv")

phy_core_taxaS.30 <- core(core_biomeS.rel, detection = 1/10000, prevalence = .3)
phy_core_taxaS.30 <- prune_samples(sample_sums(phy_core_taxaS.30)>0, phy_core_taxaS.30)
phy_core_taxaS.30 <- prune_taxa(taxa_sums(phy_core_taxaS.30)>0, phy_core_taxaS.30)
# 14 taxa and 14 samples

write.csv(tax_table(phy_core_taxaS.30), "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_SavedFiles/TaxTab_CoreS.30.csv")

# rare sediment taxa 
# <20% prevalence threshold
phylo.rareS.2 <- rare(core_biomeS.rel, detection=1/10000, prevalence=20/100) # 0.1% detection threshold, <20% of samples 
phylo.rareS.2 <- prune_samples(sample_sums(phylo.rareS.2)>0, phylo.rareS.2)
phylo.rareS.2 <- prune_taxa(taxa_sums(phylo.rareS.2)>0, phylo.rareS.2)
# 1810 taxa and 16 samples

# <10% prevalence threshold
phylo.rareS.1 <- rare(core_biomeS.rel, detection=1/10000, prevalence=10/100) # 0.1% detection threshold, <10% of samples 
phylo.rareS.1 <- prune_samples(sample_sums(phylo.rareS.1)>0, phylo.rareS.1)
phylo.rareS.1 <- prune_taxa(taxa_sums(phylo.rareS.1)>0, phylo.rareS.1)
# 1572 taxa and 16 samples

## Phyloseq of the "core" feather and sediment communities ####
# present in 30% of samples
phy_core_taxaSF.30 <- merge_phyloseq(phy_core_taxaS.30, phy_core_taxaF.30)
# 17 taxa and 53 samples

# relative abundance taxonomy plot
phy_core_taxaSF.30.rel <- transform_sample_counts(phy_core_taxaSF.30, function(OTU) OTU/sum(OTU))

# phylum level
phy_core_taxaSF.30.rel.glom.phy <- tax_glom(phy_core_taxaSF.30.rel, taxrank = "Phylum")

tax_SF.30.phylum <- plot_bar(phy_core_taxaSF.30.rel.glom.phy, fill = "Phylum") + facet_grid(~Sample_Type, space = "free", scales = "free") 
# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/ITS_Figs/Tax_SF.30.phylum.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_SF.30.phylum) 
dev.off()

### VENN DIAGRAM ####

library(ggVennDiagram)

##### Accessing the taxonomy tables of feather and sediment phyloseqs (unrarefied) #####
tax_F <- as.data.frame(tax_table(Phy_ITS_Feathers_1000))
tax_S <- as.data.frame(tax_table(ITS_S_1000))

# ASV level
Venn.SF.ASV <- list(Feather = rownames(tax_F), Sediment = rownames(tax_S))

ggVennDiagram(Venn.SF.ASV)

# Unique to feathers: 6782 (79%)
# Unique to sediment: 1233 (14%)
# Shared between feathers and sediment: 610 (7%)

VennDiagram::get.venn.partitions(Venn.SF.ASV)

##### Shared taxa between feathers and sediment - rarefied phyloseqs ####
library(ggVennDiagram)

# Accessing the taxonomy tables of feather and sediment phyloseqs (rarefied)
tax_F.rar <- as.data.frame(tax_table(Phy_ITS_Feathers_1000.rar))
tax_S.rar <- as.data.frame(tax_table(ITS_S_1000.rar))

# ASV level
Venn.SF.ASV.rar <- list(Feather = rownames(tax_F.rar), Sediment = rownames(tax_S.rar))

ggVennDiagram(Venn.SF.ASV.rar)

# Unique to feathers: 4710 (80%)
# Unique to sediment: 792 (13%)
# Shared between feathers and sediment: 415 (7%)

####Shared taxa between core feather and sediment communities ####
library(ggVennDiagram)

# Accessing the taxonomy tables of feather and sediment core phyloseqs
tax_F.core <- as.data.frame(tax_table(phy_core_taxaF.30))
tax_S.core <- as.data.frame(tax_table(phy_core_taxaS.30))

# ASV level
Venn.SF.ASV.core <- list(Feather = rownames(tax_F.core), Sediment = rownames(tax_S.core))

ggVennDiagram(Venn.SF.ASV.core)

# Unique to feathers: 3 (18%)
# Unique to sediment: 14 (82%)
# Shared between feathers and sediment: 0 

VennDiagram::get.venn.partitions(Venn.SF.ASV)

#### Keratinolytic taxa shared between feather and sediment communities ####
library(ggVennDiagram)

# Accessing the taxonomy tables of feather and sediment core phyloseqs
tax_F.core.k <- as.data.frame(tax_table(ITS_K.F))
tax_S.core.k <- as.data.frame(tax_table(ITS_K.Sediment))

# ASV level
Venn.SF.ASV.core.k <- list(Feather = rownames(tax_F.core.k), Sediment = rownames(tax_S.core.k))

ggVennDiagram(Venn.SF.ASV.core.k)

# Unique to feathers: 148 (87%)
# Unique to sediment: 10 (6%)
# Shared between feathers and sediment: 12 (7%)

VennDiagram::get.venn.partitions(Venn.SF.ASV.core.k)




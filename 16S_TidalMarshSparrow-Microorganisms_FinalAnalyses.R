# Date updated: 11/1/23

# Script description: Final analyses on ITS sequences from tidal marsh sparrow feather samples for the paper "Plumage microorganism communities of tidal marsh sparrows". 

# SET UP
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
library(xlsx); library(RColorBrewer); library(extrafont); library(ggplot2);  library(cowplot); library(gridExtra); library(wesanderson); library(ggsci); library(vegan)


## Load unrarefied phyloseqs ####

### Tide feather phyloseq, runs 3&4 ####
tide_R3R4.F.1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide_R3R4.F.1000.rds")
# 1714 taxa and 20 samples
sum(sample_sums(tide_R3R4.F.1000)) # 137495 reads
mean(sample_sums(tide_R3R4.F.1000)) # 6874.75/sample

### Maine NESP phyloseq, runs 3&4 ####
ME_34_F <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34_F.rds")
# 2560 taxa and 32 samples
sum(sample_sums(ME_34_F)) # 218218 reads
mean(sample_sums(ME_34_F)) # 6819.312 reads/sample

### Maine sediment samples ####
ME_34_S <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34_S.rds")
# 6226 taxa and 12 samples
sum(sample_sums(ME_34_S)) # 420826 reads
mean(sample_sums(ME_34_S)) # 35068.83 reads per sample

### Phyloseq with Maine reseqs merged by BandNum, recaps as independent samples ####

ME_34.F.rsqmerg <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34.F.rsqmerg.rds")

# Will remove samples with <1,000 reads (2881-26804 has only 634, this is one of the two samples dominated by Cyanobacteria)

ME_34.F.rsqmerg.1000 <- prune_samples(sample_sums(ME_34.F.rsqmerg)>1000, ME_34.F.rsqmerg)
ME_34.F.rsqmerg.1000 <- prune_taxa(taxa_sums(ME_34.F.rsqmerg.1000)>0, ME_34.F.rsqmerg.1000)

# 2559 taxa and 28 samples. 

sum(sample_sums(ME_34.F.rsqmerg.1000)) # 217584 reads. 
mean(sample_sums(ME_34.F.rsqmerg.1000)) # 7770.86 reads/sample

saveRDS(ME_34.F.rsqmerg.1000, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34.F.rsqmerg.1000.rds")

ME_34.F.rsqmerg.1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34.F.rsqmerg.1000.rds")

### ALL FEATHERS Maine NESP (rsqs combined, recaps separate) and the Spring tide samples merged, all >1000 reads ####

all.f.1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/all.f.1000.rds")


### Spring tide phyloseq with the resequenced samples merged ####

tide_F.rsqmerg.1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide_F.rsqmerg.1000.rds")

## ALL SAMPLES S&F phyloseq ####

phy_16S_all <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_16S_all.rds")

## Rarefied phyloseqs ####

### Maine feathers, reseqs combined and recaps separate ####

#### Rarefaction curve ####
library(vegan)
otutab.f.16s.me <- otu_table(ME_34.F.rsqmerg.1000)
class(otutab.f.16s.me) <- "matrix"
otutab.f.16s.me <- t(otutab.f.16s.me)  # transpose observations to rows          
rare.f.16s.me <- rarecurve(otutab.f.16s.me, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 1000))

sort(sample_sums(ME_34.F.rsqmerg.1000))
# lowest reads: 1325, sample 2781-84902


ME_34.F.rsqmerg.1000.rar <- rarefy_even_depth(ME_34.F.rsqmerg.1000, 
                                          sample.size=1325, 
                                          replace=FALSE, 
                                          trimOTUs=TRUE, 
                                          rngseed=711, 
                                          verbose=TRUE)

# 631 SVs removed
# 1928 taxa and 28 samples

sum(sample_sums(ME_34.F.rsqmerg.1000.rar)) # 37,100 reads
mean(sample_sums(ME_34.F.rsqmerg.1000.rar)) # 1325 reads/sample

saveRDS(ME_34.F.rsqmerg.1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34.F.rsqmerg.1000.rds")

ME_34.F.rsqmerg.1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34.F.rsqmerg.1000.rds")

### Spring tide phyloseq, >1000 reads ####

#### Rarefaction curve ####
library(vegan)
otutab.f.16s.tide <- otu_table(tide_F.rsqmerg.1000)
class(otutab.f.16s.tide) <- "matrix"
otutab.f.16s.tide <- t(otutab.f.16s.tide)  # transpose observations to rows          
rare.f.16s.tide <- rarecurve(otutab.f.16s.tide, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 1000))

sort(sample_sums(tide_F.rsqmerg.1000))

# lowest reads is 2811-83706 with 1515 reads

tide_F.rsqmerg.1000.rar <- rarefy_even_depth(tide_F.rsqmerg.1000, 
                                    sample.size=1515, 
                                    replace=FALSE, 
                                    trimOTUs=TRUE, 
                                    rngseed=711, 
                                    verbose=TRUE)
# 368 SVs removed

saveRDS(tide_F.rsqmerg.1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide_F.rsqmerg.1000.rar.rds")


tide_F.rsqmerg.1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide_F.rsqmerg.1000.rar.rds")

# 1411 taxa and 20 samples

sum(sample_sums(tide_F.rsqmerg.1000.rar)) # 30300 reads
mean(sample_sums(tide_F.rsqmerg.1000.rar)) # 1515

### ALL Feathers - Maine and spring tide feather samples ####

sort(sample_sums(all.f.1000))

# lowest reads is 2781-84902 with 1325 reads

all.f.1000.rar <- rarefy_even_depth(all.f.1000, 
                                 sample.size=1325, 
                                 replace=FALSE, 
                                 trimOTUs=TRUE, 
                                 rngseed=711, 
                                 verbose=TRUE)
# 717 SVs removed

# 2241 taxa and 40 samples
sum(sample_sums(all.f.1000.rar)) # 53000 reads

saveRDS(all.f.1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/all.f.1000.rar.rds")

all.f.1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/all.f.1000.rar.rds")


### Maine sediment ####
#### Rarefaction curve ####
library(vegan)
otutab.f.16s.s <- otu_table(ME_34_S)
class(otutab.f.16s.s) <- "matrix"
otutab.f.16s.s <- t(otutab.f.16s.s)  # transpose observations to rows          
rare.f.16s.s <- rarecurve(otutab.f.16s.s, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 1000))

sort(sample_sums(ME_34_S))
# lowest read count: 5973, sample BB_Soil_072821

ME_34_S.rar <- rarefy_even_depth(ME_34_S, 
                                 sample.size=5973, 
                                 replace=FALSE, 
                                 trimOTUs=TRUE, 
                                 rngseed=711, 
                                 verbose=TRUE)
# 694 OTUs removed

sample_data(ME_34_S.rar)$Sample_Type <- "Sediment"


# 5532 taxa and 12 samples
sum(sample_sums(ME_34_S.rar)) # 71676 reads
mean(sample_sums(ME_34_S.rar)) # 5973 reads/sample

saveRDS(ME_34_S.rar, '/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34_S.rar.rds')

ME_34_S.rar <- readRDS('/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_34_S.rar.rds')

### Maine feathers and Maine sediment combined  ####
ME_SF.rsqmerg.1000 <- merge_phyloseq(ME_34.F.rsqmerg.1000, ME_34_S)

saveRDS(ME_SF.rsqmerg.1000, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rds")

ME_SF.rsqmerg.1000 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rds")
# 8198 taxa and 40 samples

#### Maine sediment and feather sample rarefaction curve ####
library(vegan)
otutab.16s.sf <- otu_table(ME_SF.rsqmerg.1000)
class(otutab.16s.sf) <- "matrix"
otutab.16s.sf <- t(otutab.16s.sf)  # transpose observations to rows          
rare.16s.sf <- rarecurve(otutab.16s.sf, step=1, lwd=2, ylab="OTU",  label=F, xlim=c(0, 1000))


sum(sample_sums(ME_SF.rsqmerg.1000)) #638410 reads

# rarefy
sort(sample_sums(ME_SF.rsqmerg.1000))
# lowest read count: 1325, sample 2781-84902

ME_SF.rsqmerg.1000.rar <- rarefy_even_depth(ME_SF.rsqmerg.1000, 
                                 sample.size=1325, 
                                 replace=FALSE, 
                                 trimOTUs=TRUE, 
                                 rngseed=711, 
                                 verbose=TRUE)

# 2830 OTUs removed, no longer present after random subsampling

saveRDS(ME_SF.rsqmerg.1000.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rar.rds")

ME_SF.rsqmerg.1000.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rar.rds")
# 5346 taxa and 39 samples

sum(sample_sums(ME_SF.rsqmerg.1000.rar)) # 53000 reads

## All samples, rarefied #####

sort(sample_sums(phy_16S_all))
#  1325 lowest read count

phy_16S_all.rar <- rarefy_even_depth(phy_16S_all, 
                                            sample.size=1325, 
                                            replace=FALSE, 
                                            trimOTUs=TRUE, 
                                            rngseed=711, 
                                            verbose=TRUE)

# 2932OTUs removed, no longer present after random subsampling

saveRDS(phy_16S_all.rar, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_16S_all.rar.rds")

phy_16S_all.rar <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_16S_all.rar.rds")

# BASIC INFO ####

## Maine feathers - calculating the percentages that each phylum comprises ######

sort(table(tax_table(ME_34.F.rsqmerg.1000.rar)[, 2]))

# 1928 total taxa, 22 phyla

# proteobacteria
661/1928
# 0.343

# Bacteroidota
460/1928
# 0.239

# Actinobacteriota
189/1928
# 0.098

# Firmicutes
115/1928 #0.0596

# Planctomycetota
108/1928 # 0.056

# Verrumicroviota
91/1928 #0.0472

# top 3 phyla cumulative
0.343+0.239+0.098
#  0.68

# top 5 phyla cumulative
0.343+0.239+0.098+0.0596+ 0.056
# 0.796

## Maine feathers taxonomy plot ####
ME_34.F.rsqmerg.1000.rar.rel <- transform_sample_counts(ME_34.F.rsqmerg.1000.rar, function(OTU) OTU/sum(OTU))

ME_34.F.rsqmerg.1000.rar.rel.glom <- tax_glom(ME_34.F.rsqmerg.1000.rar.rel, taxrank = "Phylum")

tax_ME.F_phy <- plot_bar(ME_34.F.rsqmerg.1000.rar.rel.glom, fill = "Phylum") 
# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_ME.F_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_ME.F_phy) 
dev.off()

## ALL feathers (Maine and spring tide) top phyla #####

sort(table(tax_table(all.f.1000.rar)[, 2]))
# 2253 taxa total

# Proteobacteria
804/2241
0.3587684

#Bacteroidota
512/2241
0.2284694

# Actinobacteriota
212/2241
0.09460062

# Firmicutes
134/2241
0.05979473

# Planctomycetota
121/2241
0.05399375

## Tax plot - all feathers ####

all.f.1000.rar.rel <- transform_sample_counts(all.f.1000.rar, function(OTU) OTU/sum(OTU))
all.f.1000.rar.glom <- tax_glom(all.f.1000.rar.rel, taxrank = "Phylum")

tax_All.F_phy.sp <- plot_bar(all.f.1000.rar.glom, fill = "Phylum") + facet_grid(~sample_Species, scales = "free", space = "free")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_All.F_phy.sp.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_All.F_phy.sp) 
dev.off()


## ALL feathers (Maine and spring tide) - Keratinolytic genera ####

Keratinolytic <- c("Bacillus", "Pseudomonas", "Enterococcus", "Staphylococcus", "Streptomyces", "Kocuria", "Arthrobacter", "Fervidobacterium", "Janthinobacterium", "Alcaligenes", "Flavobacterium", "Chryseobacterium", "Stenotrophomonas", "Xanthomonas", "Microbacterium", "Paenibacillus", "Terrabacter", "Caldioprobacter", "Thermoanaeobacter", "Actinomadura", "Amylocolatopsis", "Brevibacillus", "Fervidobacterium", "Terrabcter", "Serratia", "Nocardiosis", "Saccharomonospara", "Thermoactinomyces", "Saccarothrix")

all.f.k <- subset_taxa(all.f.1000.rar, Genus %in% Keratinolytic)
all.f.k <- prune_samples(sample_sums(all.f.k)>0, all.f.k)
all.f.k <- prune_taxa(taxa_sums(all.f.k)> 0, all.f.k)
# 91 SVs in 39 samples (all samples)

all.f.k.tax <- as.data.frame(tax_table(all.f.k))

write.csv(all.f.k.tax, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/all.f.k.tax.csv")

## Spring tide feathers - calculating the percentages that each phylum comprises ######

sort(table(tax_table(tide_F.rsqmerg.1000.rar)[, 2]))
# 1411 SVs

# Proteobacteria
534/1411
0.378455

# Bacteroidota
313/1411
0.2218285

# Actinobacteriota
148/1411
0.1048901

# Firmicutes
87/1411
0.0616584

# Planctomycetota
70/1411
0.04961021

0.04961021 + 0.0616584 + 0.1048901 + 0.2218285 + 0.378455
0.8164422


## Tax plot spring tide feathers #####
tide_F.rsqmerg.1000.rar.rel <- transform_sample_counts(tide_F.rsqmerg.1000.rar, function(OTU) OTU/sum(OTU))

tide_F.rsqmerg.1000.rar.rel.glom <- tax_glom(tide_F.rsqmerg.1000.rar.rel, taxrank = "Phylum")

tax_tide_phy <- plot_bar(tide_F.rsqmerg.1000.rar.rel.glom, fill = "Phylum", title = "B)") + facet_grid(~sample_Species, scales = "free", space = "free")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_tide_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_tide_phy) 
dev.off()

## Maine sediment - calculating the percentages that each phylum comprises ######

sort(table(tax_table(ME_34_S.rar)[, 2]))

# 5532 SVs total

# Proteobacteria 
1263/5532
0.228308

#Bacteriodota
903/5532
0.1632321

#Planctomycetota
706/5532
0.1276211

#Verrucomicrobiota  
617/5532
0.1115329

# Myxococcota
342/5532
0.06182213

# % the top 5 phyla comprise

0.06182213+0.1115329+0.1276211+0.1632321+0.228308
0.6925162 

## Sediment taxonomy plot ####

ME_34_S.rar.rel <- transform_sample_counts(ME_34_S.rar, function(OTU) OTU/sum(OTU))
ME_34_S.rar.rel.glom <- tax_glom(ME_34_S.rar.rel, taxrank = "Phylum")
tax_S_phy <- plot_bar(ME_34_S.rar.rel.glom, fill = "Phylum") 


# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_S_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_S_phy) 
dev.off()

## Tax plot - sediment versus feathers ####

ME_SF.rsqmerg.1000.rar.rel <- transform_sample_counts(ME_SF.rsqmerg.1000.rar, function(OTU) OTU/sum(OTU))

ME_SF.rsqmerg.1000.rar.rel.glom <- tax_glom(ME_SF.rsqmerg.1000.rar.rel, taxrank = "Phylum")

tax_SF_phy <- plot_bar(ME_SF.rsqmerg.1000.rar.rel.glom, fill = "Phylum", title = "A)") + facet_grid(~Sample_Type, scales = "free", space = "free")

## Tax plot ALL SAMPLES #####

phy_16S_all.rar.rel <- transform_sample_counts(phy_16S_all.rar, function(OTU) OTU/sum(OTU))
phy_16S_all.rar.rel.glom <- tax_glom(phy_16S_all.rar.rel, taxrank = "Phylum")
tax_16Sall_phy <- plot_bar(phy_16S_all.rar.rel.glom, fill = "Phylum") + facet_grid(Sample_Type~sample_Species, scales = "free", space = "free")

tax_AllSparrowsVSediment <- plot_bar(phy_16S_all.rar.rel.glom, fill = "Phylum", title = "A)") + facet_grid(~factor(sample_Species, levels = c("Nelson's", "Saltmarsh", "Seaside", "Sediment")), scales = "free", space = "free") + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), legend.key.size = unit(.4, 'cm'))

# Export figure
png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_AllSparrowsVSediment_phy.png", width = 35, height = 25, units = 'cm', res = 300)
grid.arrange(tax_16Sall_phy_fncy2) 
dev.off()

# AIM 1: HOST SPECIES COMPARISON #####

## ALPHA DIV: Across host species (spring tide) ######
### Violin plots ####

# Host species
viol.tide.species <- plot_richness(tide_F.rsqmerg.1000.rar, 
                                   x="Species", 
                                   measures=c("Observed","Shannon"), 
                                   title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Species)) +
  geom_boxplot(width = 0.1, aes(group=Species)) + 
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  theme(text=element_text(family="Times New Roman", face="bold", size=12))

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/viol.tide.species.png", width = 23, height = 21, units = 'cm', res = 300)
grid.arrange(viol.tide.species) 
dev.off()

### Test for alpha div variation ####

tide_F.rar.rich <- estimate_richness(tide_F.rsqmerg.1000.rar, measure=c("Observed", "Shannon")) 
tide_F.rar.rich.even <- tide_F.rar.rich$Shannon/log(tide_F.rar.rich$Observed)
tide_F.rar.rich.sd = as(sample_data(tide_F.rsqmerg.1000.rar), "matrix")
tide_F.rar.rich.sd = as.data.frame(tide_F.rar.rich.sd)
tide_F.rar.rich.df <- cbind(tide_F.rar.rich, tide_F.rar.rich.even, tide_F.rar.rich.sd)
write.csv(tide_F.rar.rich.df, file = "tide_F.rar.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals

# ANOVA of effect of host species on observed diversity
a_tide_F_sp.obs <- aov(Observed ~ Species, data = tide_F.rar.rich.df)
plot(a_tide_F_sp.obs)
sum.a_tide_F_sp.obs <- summary(a_tide_F_sp.obs)

#             Df Sum Sq Mean Sq F value Pr(>F)  
# Species      2  44890   22445   4.147 0.0341 *
# Residuals   17  92014    5413    
# Significant effect of host species on observed diversity

res.tide.sp.obs <- a_tide_F_sp.obs$residuals
shapiro.test(res.tide.sp.obs)
# W = 0.91915, p-value = 0.09544. Normal, can report results of the ANOVA

# tukey to test for significant differences between pairs of species
species.tukey.16s <- TukeyHSD(a_tide_F_sp.obs)
plot(species.tukey.16s) # NESP sig lower than SALS, NESP not sig diff from SESP

# shannon diversity by species
a_tide_F_sp.shan <- aov(Shannon ~ Species, data = tide_F.rar.rich.df)
plot(a_tide_F_sp.shan)
summary(a_tide_F_sp.shan)
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Species      2  2.683  1.3416   3.063 0.0731 .
# Residuals   17  7.447  0.4381 

res.tide.sp.shan <- a_tide_F_sp.shan$residuals
shapiro.test(res.tide.sp.shan)
# W = 0.93826, p-value = 0.2223. Normal, can report results of the ANOVA


## BETA DIV: Across host species (spring tide) ####
### PCOA ####

# Host species
ord_tide_F.rar <- ordinate(tide_F.rsqmerg.1000.rar, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.tide.species <- plot_ordination(tide_F.rsqmerg.1000.rar, ord_tide_F.rar, type="samples", color="Species") + geom_point(size = 3) +
  theme(text=element_text(family="Times New Roman", size=12)) + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/PCOA.tide.species.png", width = 20, height = 20, units = 'cm', res = 300)
grid.arrange(PCOA.tide.species) 
dev.off()


### DISPERSION AND COMMUNITY COMPOSITION ####
#### HOMOGENEITY OF VARIANCE - Jaccard ####
meta.df.tide.feath <- as(sample_data(tide_F.rsqmerg.1000.rar), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.tide.F <- phyloseq::distance(tide_F.rsqmerg.1000.rar, method="jaccard")

# checking betadispersion by host species

disper.host <- betadisper(dis.tide.F, meta.df.tide.feath$Species)
disper.host
permutest(disper.host)

#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.00994 0.0049701 0.4236    999  0.649
# Residuals 17 0.19945 0.0117322  
# Host species don't have different dispersions, can trust that variation in PERMANOVA is due to differences in centroids and not dispersions

#### PERMANOVA - Jaccard ####
# effect of host species on community composition, jaccard dissimilarity

perma.tide.host <- adonis2(phyloseq::distance(tide_F.rsqmerg.1000.rar, method="jaccard") ~ Species,
                           data = meta.df.tide.feath)

#adonis2(formula = distance(tide_F.rsqmerg.1000.rar, method = "jaccard") ~ Species, data = meta.df.tide.feath)
#         Df SumOfSqs      R2      F   Pr(>F)  
# Species   2   0.8991 0.14126 1.3983  0.035 *
# Residual 17   5.4655 0.85874                
# Total    19   6.3645 1.00000  

# Host - sig effect on community composition

# Pairwise Adonis

library(pairwiseAdonis)

pair.host <- pairwise.adonis(dis.tide.F, meta.df.tide.feath$Species)
pair.host # Nelson's have sig different community composition from saltmarsh sparrows, but not from seaside or saltmarsh from seaside

#### HOMOGENEITY OF VARIANCE - Bray-Curtis ####
meta.df.tide.feath <- as(sample_data(tide_F.rsqmerg.1000.rar), "data.frame") # sample data

## Jaccard disatance between samples, for testing dispersion
dis.tide.F.b <- phyloseq::distance(tide_F.rsqmerg.1000.rar, method="bray")

# checking betadispersion by host species

disper.host.b <- betadisper(dis.tide.F.b, meta.df.tide.feath$Species)
disper.host.b
permutest(disper.host.b)

#           Df  Sum Sq  Mean Sq     F.  N.Perm Pr(>F)
# Groups     2 0.02451 0.012254 0.5378    999  0.613
# Residuals 17 0.38734 0.022785 
# Host species do not differ in group dispersion


#### PERMANOVA - bray ####
# effect of host species on community composition, bray-curtis dissimilarity

perma.tide.host.b <- adonis2(phyloseq::distance(tide_F.rsqmerg.1000.rar, method="bray") ~ Species,
                           data = meta.df.tide.feath)

#adonis2(formula = distance(tide_F.rsqmerg.1000.rar, method = "jaccard") ~ Species, data = meta.df.tide.feath)
#           Df SumOfSqs   R2      F    Pr(>F)  
# Species   2   0.7486 0.15427 1.5505  0.067 .
# Residual 17   4.1040 0.84573                
# Total    19   4.8527 1.00000    

# Host doesn't have a significant effect on community composition calculated using bray-curtis dissimilarity

### DIFFERENTIALLY ABUNDANT TAXA - DESEQ ####

#### DESEQ Differential abundance - NESP/SALS ####
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(tide_F.rsqmerg.1000, ~ Species)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
# Pulling out pairwise comparisions with the contrast argument - NESP and SALS
res.NESP.SALS = results(diagdds, contrast = c("Species", "NESP", "SALS"))
res.NESP.SALS = res.NESP.SALS[order(res.NESP.SALS$padj, na.last=NA), ]
alpha = 0.01
sigtab.NESP.SALS = res.NESP.SALS[(res.NESP.SALS$padj < alpha), ]
sigtab.NESP.SALS = cbind(as(sigtab.NESP.SALS, "data.frame"), as(tax_table(tide_F.rsqmerg.1000)[rownames(sigtab.NESP.SALS), ], "matrix")) 

head(sigtab.NESP.SALS)
dim(sigtab.NESP.SALS)
# 3 SV differentially abundant between NESP and SALS

write.csv(sigtab.NESP.SALS, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.NESP.SALS.csv")

# Plotting
# Phylum order
x = tapply(sigtab.NESP.SALS$log2FoldChange, sigtab.NESP.SALS$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.NESP.SALS$Phylum = factor(as.character(sigtab.NESP.SALS$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab.NESP.SALS$log2FoldChange, sigtab.NESP.SALS$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.NESP.SALS$Genus = factor(as.character(sigtab.NESP.SALS$Genus), levels=names(x))

## if the Genus is empty, replace with the Family
sigtab.NESP.SALS$Genus[is.na(sigtab.NESP.SALS$Genus)] <- sigtab.NESP.SALS$Family[is.na(sigtab.NESP.SALS$Genus)]
# introduced NAs

library("ggplot2")

## graph differential abundance
DESEQ.NESPvSALS <- ggplot(sigtab.NESP.SALS, aes(y=Genus, x=log2FoldChange, color=Genus)) + #play with aestetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) 

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_Figs/DESEQ.NESPvSALS.png", width = 40, height = 30, units = 'cm', res = 300)
grid.arrange(DESEQ.NESPvSALS) 
dev.off()

# calculate significance for those abundance calculations
# Pulling out pairwise comparisions with the contrast argument - NESP and SESP
res.NESP.SESP = results(diagdds, contrast = c("Species", "NESP", "SESP"))
res.NESP.SESP = res.NESP.SESP[order(res.NESP.SESP$padj, na.last=NA), ]
alpha = 0.01
sigtab.NESP.SESP = res.NESP.SALS[(res.NESP.SESP$padj < alpha), ]
sigtab.NESP.SESP = cbind(as(sigtab.NESP.SESP, "data.frame"), as(tax_table(tide_F.rsqmerg.1000)[rownames(sigtab.NESP.SESP), ], "matrix")) 

head(sigtab.NESP.SESP)
dim(sigtab.NESP.SESP)
# 2 SV differentially abundant between NESP and SESP

write.csv(sigtab.NESP.SESP, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.NESP.SESP.csv")


# calculate significance for those abundance calculations
# Pulling out pairwise comparisions with the contrast argument - SALS and SESP
res.SALS.SESP = results(diagdds, contrast = c("Species", "SALS", "SESP"))
res.SALS.SESP = res.SALS.SESP[order(res.SALS.SESP$padj, na.last=NA), ]
alpha = 0.01
sigtab.SALS.SESP = res.SALS.SESP[(res.SALS.SESP$padj < alpha), ]
sigtab.SALS.SESP = cbind(as(sigtab.SALS.SESP, "data.frame"), as(tax_table(tide_F.rsqmerg.1000)[rownames(sigtab.SALS.SESP), ], "matrix")) 

head(sigtab.SALS.SESP)
dim(sigtab.SALS.SESP)
# 1 differentially abundant SV between SALS and SESP

write.csv(sigtab.SALS.SESP, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.SALS.SESP.csv")


## KERATINOLYTIC GENERA - Host Species ####

Keratinolytic <- c("Bacillus", "Pseudomonas", "Enterococcus", "Staphylococcus", "Streptomyces", "Kocuria", "Arthrobacter", "Fervidobacterium", "Janthinobacterium", "Alcaligenes", "Flavobacterium", "Chryseobacterium", "Stenotrophomonas", "Xanthomonas", "Microbacterium", "Paenibacillus", "Terrabacter", "Caldioprobacter", "Thermoanaeobacter", "Actinomadura", "Amylocolatopsis", "Brevibacillus", "Fervidobacterium", "Terrabcter", "Serratia", "Nocardiosis", "Saccharomonospara", "Thermoactinomyces", "Saccarothrix")

tide.rar.k <- subset_taxa(tide_F.rsqmerg.1000.rar, Genus %in% Keratinolytic)
tide.rar.k <- prune_samples(sample_sums(tide.rar.k)>0, tide.rar.k)
tide.rar.k <- prune_taxa(taxa_sums(tide.rar.k)> 0, tide.rar.k)
# 59 taxa and 20 samples

saveRDS(tide.rar.k, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide.rar.k")

tide.rar.k <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tide.rar.k")

### VIOLIN - keratinolytic taxa by host species ####

viol.tide.species.k <- plot_richness(tide.rar.k, 
                                     x="Species",
                                     measures=c("Observed","Shannon"), 
                                     title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Species)) + 
  geom_boxplot(width = 0.1, aes(group=Species)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/viol.tide.species.k.png", width = 23, height = 21, units = 'cm', res = 300)
grid.arrange(viol.tide.species.k) 
dev.off()

### Testing for alpha div variation ####

tide.rar.k.rich <- estimate_richness(tide.rar.k, measure=c("Observed", "Shannon")) 
tide.rar.k.rich.even <- tide.rar.k.rich$Shannon/log(tide.rar.k.rich$Observed)
tide.rar.k.rich.sd = as(sample_data(tide.rar.k), "matrix")
tide.rar.k.rich.sd = as.data.frame(tide.rar.k.rich.sd)
tide.rar.k.rich.df <- cbind(tide.rar.k.rich, tide.rar.k.rich.even, tide.rar.k.rich.sd)
write.csv(tide.rar.k.rich.df, file = "tide.rar.k.rich.df.csv")

# ANOVA of effect of host species on observed diversity
a_tide_sp.obs.k <- aov(Observed ~ Species, data = tide.rar.k.rich.df)
plot(a_tide_sp.obs.k)
summary(a_tide_sp.obs.k)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Species      2   84.4   42.19   1.785  0.198
# Residuals   17  401.8   23.64    
res.tide.sp.obs.k <- a_tide_sp.obs.k$residuals
shapiro.test(res.tide.sp.obs.k)
# W = 0.97128, p-value = 0.7817, normal, report ANOVA

# ANOVA of effect of host species on shannon diversity
a_tide_sp.shan.k <- aov(Shannon ~ Species, data = tide.rar.k.rich.df)
plot(a_tide_sp.shan.k)
summary(a_tide_sp.shan.k)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Species      2  0.909  0.4544   1.876  0.184
# Residuals   17  4.118  0.2423     
res.tide.sp.shan.k <- a_tide_sp.shan.k$residuals
shapiro.test(res.tide.sp.shan.k)
# W = 0.95446, p-value = 0.43997, normal, report ANOVA


## BETA DIV - Keratinolytic taxa by host species ####

### DISPERSION AND COMMUNITY COMPOSITION ######
#### HOMOEGENEITY OF VARIANCE - Jaccard ####
meta.df.tide.k <- as(sample_data(tide.rar.k), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.tide.F.k <- distance(tide.rar.k, method="jaccard")

# checking betadispersion by host species

disper.host.k <- betadisper(dis.tide.F.k, meta.df.tide.k$Species)
disper.host.k
permutest(disper.host.k)

#           Df  Sum Sq  Mean Sq   F    N.Perm Pr(>F)
#Groups     2 0.01922 0.009611 0.2552    999   0.77
#Residuals 17 0.64026 0.037662        
# Host species don't have different dispersions, can trust that variation in PERMANOVA is due to differences in centroids and not dispersions

#### HOMOEGENEITY OF VARIANCE - Bray ####
meta.df.tide.k <- as(sample_data(tide.rar.k), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.tide.F.k.b <- distance(tide.rar.k, method="bray")

# checking betadispersion by host species

disper.host.k.b <- betadisper(dis.tide.F.k.b, meta.df.tide.k$Species)
disper.host.k.b
permutest(disper.host.k.b)

#           Df  Sum Sq  Mean Sq    F     N.Perm Pr(>F)
# Groups     2 0.03047 0.015234 0.2859    999  0.729
# Residuals 17 0.90570 0.053277   


#### PERMANOVA - keratinolytic taxa by host species ####
# effect of host species, jaccard

perma.tide.host.k <- adonis2(distance(tide.rar.k, method="jaccard") ~ Species,
                             data = meta.df.tide.k)

#        Df SumOfSqs      R2      F Pr(>F)  
#Species   2   0.8337 0.14981 1.4978  0.059 .
#Residual 17   4.7316 0.85019                
#Total    19   5.5653 1.00000    

#### PERMANOVA - keratinolytic taxa by host species ####
# effect of host species, bray

perma.tide.host.k <- adonis2(distance(tide.rar.k, method="bray") ~ Species,
                             data = meta.df.tide.k)

#        Df SumOfSqs     R2      F Pr(>F)
# Species   2   0.6233 0.1496 1.4953   0.12
# Residual 17   3.5432 0.8504              
# Total    19   4.1665 1.0000  

# AIM 2: WITHIN-SPECIES VARIATION ####

## ALPHA DIV: Across Maine NESP variables ####

### Test for alpha div variation ####

ME_F.rar.rich <- estimate_richness(ME_34.F.rsqmerg.1000.rar, measure=c("Observed", "Shannon")) 
ME_F.rar.rich.even <- ME_F.rar.rich$Shannon/log(ME_F.rar.rich$Observed)
ME_F.rar.rich.sd = as(sample_data(ME_34.F.rsqmerg.1000.rar), "matrix")
ME_F.rar.rich.sd = as.data.frame(ME_F.rar.rich.sd)
ME_F.rar.rich.df <- cbind(ME_F.rar.rich, ME_F.rar.rich.even, ME_F.rar.rich.sd)
write.csv(ME_F.rar.rich.df, file = "ME_F.rar.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals

# ANOVA of effect of sex on observed diversity
a_ME_F_sex.obs <- aov(Observed ~ Sex, data = ME_F.rar.rich.df)
plot(a_ME_F_sex.obs)
summary(a_ME_F_sex.obs)
res.ME.sex.obs <- a_ME_F_sex.obs$residuals
shapiro.test(res.ME.sex.obs)
# W = 0.90691, p-value = 0.01667 residuals non-normal, so need to run kruskall-wallace

kruskal.test(Observed ~ Sex, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 0.52978, df = 1, p-value = 0.4667
# sex not sig on observed

a_ME_F_sex.shan <- aov(Shannon ~ Sex, data = ME_F.rar.rich.df)
plot(a_ME_F_sex.shan)
summary(a_ME_F_sex.shan)
res.ME.sex.shan <- a_ME_F_sex.shan$residuals
shapiro.test(res.ME.sex.shan)
# W = 0.77366, p-value = 3.791e-05 - non-normal, run Kruskall-Wallace

#Sex
kruskal.test(Shannon ~ Sex, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 0.20063, df = 1, p-value = 0.6542
# Sex not sig effect on shannon div


# Marsh
# ANOVA of effect of sex on observed diversity
a_ME_F_ma.obs <- aov(Observed ~ Marsh, data = ME_F.rar.rich.df)
plot(a_ME_F_ma.obs)
res.ME.ma.obs <- a_ME_F_ma.obs$residuals
shapiro.test(res.ME.ma.obs)
# W = 0.90426, p-value = 0.01442, non-normal residuals, use kruskall-wallis

kruskal.test(Observed ~ Marsh, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 3.9395, df = 3, p-value = 0.2681
# Marsh not sig on observed

#Marsh
a_ME_F_ma.shan <- aov(Shannon ~ Marsh, data = ME_F.rar.rich.df)
plot(a_ME_F_ma.shan)
res.ME.ma.shan <- a_ME_F_ma.shan$residuals
shapiro.test(res.ME.ma.shan)
# W = 0.85511, p-value = 0.001188. Non-normal, use Kruskall wallis
kruskal.test(Shannon ~ Marsh, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 3.8259, df = 3, p-value = 0.2809
# Marsh not sig effect on shannon div


#Month
a_ME_F_mo.obs <- aov(Observed ~ Month, data = ME_F.rar.rich.df)
plot(a_ME_F_mo.obs)
res.ME.mo.obs <- a_ME_F_mo.obs$residuals
shapiro.test(res.ME.mo.obs)
# W = 0.91283, p-value = 0.02314 - non-normal residuals, use kruskall wallis

kruskal.test(Observed ~ Month, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 2.4883, df = 2, p-value = 0.2882
# Month not sig on Observed

a_ME_F_mo.shan <- aov(Shannon ~ Month, data = ME_F.rar.rich.df)
plot(a_ME_F_mo.shan)
res.ME.mo.shan <- a_ME_F_mo.shan$residuals
shapiro.test(res.ME.mo.shan)
# W = 0.87978, p-value = 0.003971 - non-normal, use kruskall wallis

kruskal.test(Shannon ~ Month, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.621, df = 2, p-value = 0.4446
# month not sig effect on shannon div


#TideCycle - spring, neap, or intermediate
a_ME_F_cyc.obs <- aov(Observed ~ TideCycle, data = ME_F.rar.rich.df)
plot(a_ME_F_cyc.obs)
res.ME.cyc.obs <- a_ME_F_cyc.obs$residuals
shapiro.test(res.ME.cyc.obs)
# W = 0.9172, p-value = 0.02959, non-normal, use kruskal

kruskal.test(Observed ~ TideCycle, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.0244, df = 2, p-value = 0.5992
# not sig effect of tidal cycle on observed diversity

a_ME_F_cyc.shan <- aov(Shannon ~ TideCycle, data = ME_F.rar.rich.df)
plot(a_ME_F_cyc.shan)
res.ME.cyc.shan <- a_ME_F_cyc.shan$residuals
shapiro.test(res.ME.cyc.shan)
# W = 0.79739, p-value = 9.604e-05 non-normal, use kruskal
kruskal.test(Shannon ~ TideCycle, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 1.0417, df = 2, p-value = 0.594
# no sig effect of tide cycle on shannon diversity

# Tidal range - meters
a_ME_F_rng.m.obs <- aov(Observed ~ TidalRange_m, data = ME_F.rar.rich.df)
plot(a_ME_F_rng.m.obs)
res.ME.rng.m.obs <- a_ME_F_rng.m.obs$residuals
shapiro.test(res.ME.rng.m.obs) # W = 0.9505, p-value = 0.204
# Normal residuals, use ANOVA results 
summary(a_ME_F_rng.m.obs)
#               Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange_m  4  44141   11035   1.229  0.326
# Residuals    23 206470    8977 
# Tidal range(m) doesn't significantly effect observed diversity
# did get a warning about unbalanced samples, so going to report kruskall wallis results
kruskal.test(Observed ~ TidalRange_m, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 5.9255, df = 4, p-value = 0.2048

a_ME_F_rng.m.shan <- aov(Shannon ~ TidalRange_m, data = ME_F.rar.rich.df)
plot(a_ME_F_rng.m.shan)
res.ME.rng.m.shan <- a_ME_F_rng.m.shan$residuals
shapiro.test(res.ME.rng.m.shan)
# W = 0.9189, p-value = 0.03258 non-normal, run Kruskal wallis
kruskal.test(Shannon ~ TidalRange_m, data = ME_F.rar.rich.df)
# Kruskal-Wallis chi-squared = 4.2038, df = 4, p-value = 0.3791
# Tidal range(m) not sig effect on shannon div


## BETA DIV: Across Maine NESP variables ####

### PERMANOVAs ####

library(vegan)

# Effect of sex, month, marsh, tidal range, tide cycle on community composition in Maine NESP feather samples

# dataframe of the sample data
meta.df.Me.feath.16s <- as(sample_data(ME_34.F.rsqmerg.1000.rar), "data.frame")

#### HOMOGENEITY OF VARIANCE - Jaccard #####
## Jaccard disatance between samples, for testing dispersion
dis.ME.F.16s <- phyloseq::distance(ME_34.F.rsqmerg.1000.rar, method="jaccard")

# checking betadispersion of each variable

# sex
disper.sex.16s <- betadisper(dis.ME.F.16s, meta.df.Me.feath.16s$Sex)
disper.sex.16s
permutest(disper.sex.16s)
#           Df   Sum Sq   Mean Sq   F     N.Perm Pr(>F)  
#Groups     1 0.017802 0.0178021 3.5453    999  0.083 .
#Residuals 26 0.130553 0.0050213 

# ANOVA of the dispersions
anova(disper.sex.16s)
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
#Groups     1 0.017802 0.0178021  3.5453 0.07095 .
#Residuals 26 0.130553 0.0050213 
# ANOVA results are non-significant. Will go with the aggregation of results pointing towards non-significant difference in variances between sexes

# homogeneity of variance, month
disper.month.16s <- betadisper(dis.ME.F.16s, meta.df.Me.feath.16s$Month)
disper.month.16s
permutest(disper.month.16s)
plot(disper.month.16s)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
#Groups     2 0.029947 0.0149736 3.1538    999  0.073 .
#Residuals 25 0.118696 0.0047478 
# month doesn't have different dispersions

anova(disper.month.16s)
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
#Groups     2 0.029947 0.0149736  3.1538 0.06007 .
#Residuals 25 0.118696 0.0047478 

# Homogeneity of variance, marsh
disper.marsh.16s <- betadisper(dis.ME.F.16s, meta.df.Me.feath.16s$Marsh)
disper.marsh.16s
permutest(disper.marsh.16s)
plot(disper.marsh.16s)
#           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)    
# Groups     3 0.34874 0.116247 29.08    999  0.001 ***
# Residuals 24 0.09594 0.003997 

# ANOVA of the dispersions
anova(disper.marsh.16s)
#Response: Distances
#          Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups     3 0.34874 0.116247   29.08 3.668e-08 ***
# Residuals 24 0.09594 0.003997

# Tukey's Honest Significant Differences to determine which group(s) are more variable

marsh.tukey.16s <- TukeyHSD(disper.marsh.16s)
plot(marsh.tukey.16s)

#           diff         lwr         upr     p adj
#HA-BB  0.115495389 -0.02434735  0.25533812 0.1314137
#NG-BB  0.106964647 -0.02303627  0.23696556 0.1335336
#PL-BB -0.480392157 -0.69400546 -0.26677885 0.0000117
#NG-HA -0.008530742 -0.08622115  0.06915967 0.9900979
#PL-HA -0.595887545 -0.78234452 -0.40943057 0.0000000
#PL-NG -0.587356803 -0.76655070 -0.40816291 0.0000000

# Tidal Range
disper.rng <- betadisper(dis.ME.F.16s, meta.df.Me.feath.16s$TidalRange_m)
disper.rng
permutest(disper.rng)

#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     4 0.021653 0.0054132 0.8456    999  0.483
# Residuals 23 0.147241 0.0064018 

anova(disper.rng)
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     4 0.021653 0.0054132  0.8456 0.5108
# Residuals 23 0.147241 0.0064018    
# Not significanlty different group dispersions

# Tide Cycle
disper.cycle.16s <- betadisper(dis.ME.F.16s, meta.df.Me.feath.16s$TideCycle)
disper.cycle.16s
permutest(disper.cycle.16s)

#             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     2 0.047131 0.0235656 4.4919    999   0.02 *
# Residuals 25 0.131157 0.0052463       

# ANOVA of the dispersions
anova(disper.cycle.16s)
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     2 0.047131 0.0235656  4.4919 0.02155 *
# Residuals 25 0.131157 0.0052463  
# Tide cycle (spring, neap, or intermediate) has different dispersions

#### HOMOGENETIY OF VARIANCE - Bray Curtis #####

# dataframe of the sample data
meta.df.Me.feath.16s <- as(sample_data(ME_34.F.rsqmerg.1000.rar), "data.frame")

## Bray-curtis dissimilariy between samples, for testing dispersion
dis.ME.F.16s.b <- phyloseq::distance(ME_34.F.rsqmerg.1000.rar, method="bray")

# checking betadispersion of each variable

# sex
disper.sex.16s.b <- betadisper(dis.ME.F.16s.b, meta.df.Me.feath.16s$Sex)
disper.sex.16s.b
permutest(disper.sex.16s.b)

#            Df   Sum Sq  Mean Sq   F     N.Perm Pr(>F)
# Groups     1 0.023803 0.023803 2.0146    999  0.168
# Residuals 26 0.307198 0.011815          

# Month
disper.month.16s.b <- betadisper(dis.ME.F.16s.b, meta.df.Me.feath.16s$Month)
disper.month.16s.b
permutest(disper.month.16s.b)

#            Df   Sum Sq  Mean Sq   F   N.Perm Pr(>F)  
# Groups     2 0.057584 0.028792 2.701    999   0.09 .
# Residuals 25 0.266495 0.010660  

# Marsh
disper.marsh.16s.b <- betadisper(dis.ME.F.16s.b, meta.df.Me.feath.16s$Marsh)
disper.marsh.16s.b
permutest(disper.marsh.16s.b)

#            Df  Sum Sq  Mean Sq  F    N.Perm Pr(>F)   
# Groups     3 0.27741 0.092471 9.712    999  0.002 **
# Residuals 24 0.22851 0.009521

marsh.tukey.16s.b <- TukeyHSD(disper.marsh.16s.b)
plot(marsh.tukey.16s.b)

# Tidal Range
disper.rng.16s.b <- betadisper(dis.ME.F.16s.b, meta.df.Me.feath.16s$TidalRange_m)
disper.rng.16s.b
permutest(disper.rng.16s.b)

#.           Df  Sum Sq  Mean Sq   F    N.Perm Pr(>F)
# Groups     4 0.04595 0.011488 0.808    999  0.526
# Residuals 23 0.32700 0.014217      

# Tide Cycle
disper.cyc.16s.b <- betadisper(dis.ME.F.16s.b, meta.df.Me.feath.16s$TideCycle)
disper.cyc.16s.b
permutest(disper.cyc.16s.b)

#.          Df  Sum Sq  Mean Sq   F      N.Perm Pr(>F)  
# Groups     2 0.09148 0.045740 3.9081    999  0.031 *
# Residuals 25 0.29260 0.011704 


#### COMMUNITY COMPOSITION - Jaccard ######

# Sex
perma.ME.f.sex <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="jaccard") ~ Sex,
                        data = meta.df.Me.feath.16s)

#          Df SumOfSqs      R2   F   Pr(>F)
#Sex       1   0.4212 0.03985 1.0791  0.272
#Residual 26  10.1476 0.96015              
#Total    27  10.5687 1.00000    

# Marsh
perma.ME.f.marsh <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="jaccard") ~ Marsh,
                          data = meta.df.Me.feath.16s)

#          Df SumOfSqs    R2    F     Pr(>F)  
#Marsh     3   1.3301 0.12585 1.1518  0.093 .
#Residual 24   9.2386 0.87415                
#Total    27  10.5687 1.00000 

# Month
perma.ME.f.month <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="jaccard") ~ Month,
                            data = meta.df.Me.feath.16s)


#        Df SumOfSqs      R2    F Pr(>F)
#Month     2   0.8975 0.08492 1.16  0.132
#Residual 25   9.6713 0.91508            
#Total    27  10.5687 1.00000  


# Tidal Range
perma.ME.f.range <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="jaccard") ~ TidalRange_m,
                         data = meta.df.Me.feath.16s)


#              Df   SumOfSqs   R2     F    Pr(>F)  
# TidalRange_m  4   1.7645 0.16696 1.1524   0.08 .
# Residual     23   8.8042 0.83304                
# Total        27  10.5687 1.00000  
# Not sig

# Tidal Cycle
perma.ME.f.cyc <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="jaccard") ~ TideCycle,
                            data = meta.df.Me.feath.16s)

#           Df SumOfSqs    R2     F    Pr(>F)  
# TideCycle  2   0.9482 0.08972 1.232  0.073 .
# Residual  25   9.6205 0.91028               
# Total     27  10.5687 1.00000 
# Not sig

#### COMMUNITY COMPOSITION - BRAY CURTIS #####

# Sex
perma.ME.f.sex.b <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="bray") ~ Sex,
                                  data = meta.df.Me.feath.16s)

#          Df SumOfSqs   R2      F     Pr(>F)
# Sex       1   0.3709 0.04207 1.1419  0.261
# Residual 26   8.4453 0.95793              
# Total    27   8.8162 1.00000 

# Month
perma.ME.f.month.b <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="bray") ~ Month,
                            data = meta.df.Me.feath.16s)

#          Df SumOfSqs      R2    F    Pr(>F)
# Month     2   0.7584 0.08602 1.1765  0.168
# Residual 25   8.0578 0.91398              
# Total    27   8.8162 1.00000   

# Marsh
perma.ME.f.marsh.b <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="bray") ~ Marsh,
                              data = meta.df.Me.feath.16s)

#          Df SumOfSqs     R2    F    Pr(>F)  
# Marsh     3   1.2049 0.13666 1.2664  0.096 .
# Residual 24   7.6114 0.86334                
# Total    27   8.8162 1.00000  

# Tidal Range
perma.ME.f.rng.b <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="bray") ~ TidalRange_m,
                              data = meta.df.Me.feath.16s)

#               Df SumOfSqs    R2     F.    Pr(>F)
# TidalRange_m  4   1.5176 0.17214 1.1956  0.119
# Residual     23   7.2986 0.82786              
# Total        27   8.8162 1.00000 

# Tidal Cycle
perma.ME.f.cyc.b <- adonis2(distance(ME_34.F.rsqmerg.1000.rar, method="bray") ~ TideCycle,
                            data = meta.df.Me.feath.16s)

#           Df SumOfSqs     R2     F.    Pr(>F)
# TideCycle  2   0.7911 0.08973 1.2323  0.142
# Residual  25   8.0251 0.91027              
# Total     27   8.8162 1.00000 


#### DESEQ - sex #####
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_34.F.rsqmerg.1000, ~ Sex)

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
sigtab.sex = cbind(as(sigtab, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab), ], "matrix"))

head(sigtab.sex)
dim(sigtab.sex)
# 0 differentially abundant taxa between sexes 

#### DESEQ - Tide Cycle ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_34.F.rsqmerg.1000, ~ TideCycle)

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
sigtab.Spring.Neap = cbind(as(sigtab.Spring.Neap, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Spring.Neap), ], "matrix")) 

head(sigtab.Spring.Neap)
dim(sigtab.Spring.Neap)
# 4 differentially abundant taxa between Spring and Neap tide cycles

write.csv(sigtab.Spring.Neap, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.Spring.Neap.csv")


# Spring and Intermediate pairwise comparison
res.Spring.Int = results(diagdds, contrast = c("TideCycle", "Spring", "Intermediate"))
res.Spring.Int = res.Spring.Int[order(res.Spring.Int$padj, na.last=NA), ]
alpha = 0.01
sigtab.Spring.Int = res.Spring.Int[(res.Spring.Int$padj < alpha), ]
sigtab.Spring.Int = cbind(as(sigtab.Spring.Int, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Spring.Int), ], "matrix")) 

head(sigtab.Spring.Int)
dim(sigtab.Spring.Int)
# 4 differentially abundant taxa between Spring and Intermediate tide cycles

write.csv(sigtab.Spring.Int, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.Spring.Int.csv")

# Neap and Intermediate pairwise comparison
res.Neap.Int = results(diagdds, contrast = c("TideCycle", "Neap", "Intermediate"))
res.Neap.Int = res.Neap.Int[order(res.Neap.Int$padj, na.last=NA), ]
alpha = 0.01
sigtab.Neap.Int = res.Neap.Int[(res.Neap.Int$padj < alpha), ]
sigtab.Neap.Int = cbind(as(sigtab.Neap.Int, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Neap.Int), ], "matrix")) 
head(sigtab.Neap.Int)
dim(sigtab.Neap.Int)
# 4 differentially abundant taxa between Neap and Intermediate tide cycles

write.csv(sigtab.Neap.Int, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.Neap.Int.csv")

##### DESEQ - Marsh ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_34.F.rsqmerg.1000, ~ Marsh)

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
sigtab.HA.PL = cbind(as(sigtab.HA.PL, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.HA.PL), ], "matrix")) 

head(sigtab.HA.PL)
dim(sigtab.HA.PL)
# 0 differentially abundant SVs

# HA versus NG
res.HA.NG = results(diagdds, contrast =c("Marsh", "HA", "NG"))
res.HA.NG = res.HA.NG[order(res.HA.NG$padj, na.last=NA), ]
alpha = 0.01
sigtab.HA.NG = res.HA.NG[(res.HA.NG$padj < alpha), ]
sigtab.HA.NG = cbind(as(sigtab.HA.NG, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.HA.NG), ], "matrix")) 

head(sigtab.HA.NG)
dim(sigtab.HA.NG)
# 6 differentially abundant SVs

write.csv(sigtab.HA.NG, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.HA.NG.csv")

# HA versus BB
res.HA.BB = results(diagdds, contrast =c("Marsh", "HA", "BB"))
res.HA.BB = res.HA.BB[order(res.HA.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.HA.BB = res.HA.BB[(res.HA.BB$padj < alpha), ]
sigtab.HA.BB = cbind(as(sigtab.HA.BB, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.HA.BB), ], "matrix")) 

head(sigtab.HA.BB)
dim(sigtab.HA.BB)
# 11 Differentially abundant SVs

write.csv(sigtab.HA.BB, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.HA.BB.csv")


## NG versus PL
res.NG.PL = results(diagdds, contrast =c("Marsh", "NG", "PL"))
res.NG.PL = res.NG.PL[order(res.NG.PL$padj, na.last=NA), ]
alpha = 0.01
sigtab.NG.PL= res.NG.PL[(res.NG.PL$padj < alpha), ]
sigtab.NG.PL = cbind(as(sigtab.NG.PL, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.NG.PL), ], "matrix")) 

head(sigtab.NG.PL)
dim(sigtab.NG.PL)
# 0 differentially abundant SVs

## NG versus BB
res.NG.BB = results(diagdds, contrast =c("Marsh", "NG", "BB"))
res.NG.BB = res.NG.BB[order(res.NG.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.NG.BB = res.NG.BB[(res.NG.BB$padj < alpha), ]
sigtab.NG.BB = cbind(as(sigtab.NG.BB, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.NG.BB), ], "matrix")) 

head(sigtab.NG.BB)
dim(sigtab.NG.BB)
# 12 differentially abundant SVs 

write.csv(sigtab.NG.BB, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.NG.BB.csv")

## PL versus BB
res.PL.BB = results(diagdds, contrast =c("Marsh", "PL", "BB"))
res.PL.BB = res.PL.BB[order(res.PL.BB$padj, na.last=NA), ]
alpha = 0.01
sigtab.PL.BB = res.PL.BB[(res.PL.BB$padj < alpha), ]
sigtab.PL.BB = cbind(as(sigtab.PL.BB, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.PL.BB), ], "matrix")) 

head(sigtab.PL.BB)
dim(sigtab.PL.BB)
# 0 differentially abundant SVs between PL and BB

#### DESEQ - Month #####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_34.F.rsqmerg.1000, ~ Month)

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res.Jun.Aug = results(diagdds, contrast = c("Month", "June", "Aug"))
res.Jun.Aug = res.Jun.Aug[order(res.Jun.Aug$padj, na.last=NA), ]
alpha = 0.01
sigtab.Jun.Aug = res.Jun.Aug[(res.Jun.Aug$padj < alpha), ]
sigtab.Jun.Aug = cbind(as(sigtab.Jun.Aug, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Jun.Aug), ], "matrix")) 

head(sigtab.Jun.Aug)
dim(sigtab.Jun.Aug)
# 0 differentially abundant taxa between Jun and Aug

# Jun and July
# calculate significance for those abundance calculations
res.Jun.Jul = results(diagdds, contrast = c("Month", "June", "July"))
res.Jun.Jul = res.Jun.Jul[order(res.Jun.Jul$padj, na.last=NA), ]
alpha = 0.01
sigtab.Jun.Jul = res.Jun.Jul[(res.Jun.Jul$padj < alpha), ]
sigtab.Jun.Jul = cbind(as(sigtab.Jun.Jul, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Jun.Jul), ], "matrix")) 

head(sigtab.Jun.Jul)
dim(sigtab.Jun.Jul)
# 0 differentially abundant taxa between Jun and July

# July and August
# calculate significance for those abundance calculations
res.Aug.Jul = results(diagdds, contrast = c("Month", "Aug", "July"))
res.Aug.Jul = res.Aug.Jul[order(res.Aug.Jul$padj, na.last=NA), ]
alpha = 0.01
sigtab.Aug.Jul = res.Aug.Jul[(res.Aug.Jul$padj < alpha), ]
sigtab.Aug.Jul = cbind(as(sigtab.Aug.Jul, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.Aug.Jul), ], "matrix")) 

head(sigtab.Aug.Jul)
dim(sigtab.Aug.Jul)
# 0 differentially abundant taxa between August and July

#### DESEQ - Tidal Range ####

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_34.F.rsqmerg.1000, ~ TidalRange_m)

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
sigtab.2.5_to_3.3_to_3.5 = cbind(as(sigtab.2.5_to_3.3_to_3.5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.2.5_to_3.3_to_3.5), ], "matrix")) 

head(sigtab.2.5_to_3.3_to_3.5)
dim(sigtab.2.5_to_3.3_to_3.5)
# 10 differentially abundant SVs

write.csv(sigtab.2.5_to_3.3_to_3.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.2.5_to_3.3_to_3.5.csv")

# 2.5 to 3 versus 3.5 to 4
res.2.5_to_3.3.5_to_4 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "3.5_to_4"))
res.2.5_to_3.3.5_to_4 = res.2.5_to_3.3.5_to_4[order(res.2.5_to_3.3.5_to_4$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.3.5_to_4 = res.2.5_to_3.3.5_to_4[(res.2.5_to_3.3.5_to_4$padj < alpha), ]
sigtab.2.5_to_3.3.5_to_4 = cbind(as(sigtab.2.5_to_3.3.5_to_4, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.2.5_to_3.3.5_to_4), ], "matrix")) 

head(sigtab.2.5_to_3.3.5_to_4)
dim(sigtab.2.5_to_3.3.5_to_4)
# 4 differentially abundant SVs

write.csv(sigtab.2.5_to_3.3.5_to_4, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.2.5_to_3.3.5_to_4.csv")


# 2.5 to 3 versus 4 to 4.5
res.2.5_to_3.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "4_to_4.5"))
res.2.5_to_3.4_to_4.5 = res.2.5_to_3.4_to_4.5[order(res.2.5_to_3.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.4_to_4.5 = res.2.5_to_3.4_to_4.5[(res.2.5_to_3.4_to_4.5$padj < alpha), ]
sigtab.2.5_to_3.4_to_4.5 = cbind(as(sigtab.2.5_to_3.4_to_4.5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.2.5_to_3.4_to_4.5), ], "matrix")) 

head(sigtab.2.5_to_3.4_to_4.5)
dim(sigtab.2.5_to_3.4_to_4.5)
# 9 differentially abundant SVs

write.csv(sigtab.2.5_to_3.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.2.5_to_3.4_to_4.5.csv")

# 2.5 to 3 versus 4.5 to 5
res.2.5_to_3.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "2.5_to_3", "4.5_to_5"))
res.2.5_to_3.4.5_to_5 = res.2.5_to_3.4.5_to_5[order(res.2.5_to_3.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.2.5_to_3.4.5_to_5 = res.2.5_to_3.4.5_to_5[(res.2.5_to_3.4.5_to_5$padj < alpha), ]
sigtab.2.5_to_3.4.5_to_5 = cbind(as(sigtab.2.5_to_3.4.5_to_5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.2.5_to_3.4.5_to_5), ], "matrix")) 

head(sigtab.2.5_to_3.4.5_to_5)
dim(sigtab.2.5_to_3.4.5_to_5)
# 5 differentially abundant SVs

write.csv(sigtab.2.5_to_3.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.2.5_to_3.4.5_to_5.csv")

# 3 to 3.5 versus 3.5 to 4
res.3_to_3.5.3.5_to_4 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "3.5_to_4"))
res.3_to_3.5.3.5_to_4 = res.3_to_3.5.3.5_to_4[order(res.3_to_3.5.3.5_to_4$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.3.5_to_4 = res.3_to_3.5.3.5_to_4[(res.3_to_3.5.3.5_to_4$padj < alpha), ]
sigtab.3_to_3.5.3.5_to_4 = cbind(as(sigtab.3_to_3.5.3.5_to_4, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.3_to_3.5.3.5_to_4), ], "matrix")) 

head(sigtab.3_to_3.5.3.5_to_4)
dim(sigtab.3_to_3.5.3.5_to_4)
# 8 differentially abundant SVs

write.csv(sigtab.3_to_3.5.3.5_to_4, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.3_to_3.5.3.5_to_4.csv")

# 3 to 3.5 versus 4 to 4.5
res.3_to_3.5.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "4_to_4.5"))
res.3_to_3.5.4_to_4.5 = res.3_to_3.5.4_to_4.5[order(res.3_to_3.5.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.4_to_4.5 = res.3_to_3.5.4_to_4.5[(res.3_to_3.5.4_to_4.5$padj < alpha), ]
sigtab.3_to_3.5.4_to_4.5 = cbind(as(sigtab.3_to_3.5.4_to_4.5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.3_to_3.5.4_to_4.5), ], "matrix"))

head(sigtab.3_to_3.5.4_to_4.5)
dim(sigtab.3_to_3.5.4_to_4.5)
# 9 differentially abundant SVs

write.csv(sigtab.3_to_3.5.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.3_to_3.5.4_to_4.5.csv")


# 3 to 3.5 versus 4.5 to 5
res.3_to_3.5.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "3_to_3.5", "4.5_to_5"))
res.3_to_3.5.4.5_to_5 = res.3_to_3.5.4.5_to_5[order(res.3_to_3.5.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3_to_3.5.4.5_to_5 = res.3_to_3.5.4.5_to_5[(res.3_to_3.5.4.5_to_5$padj < alpha), ]
sigtab.3_to_3.5.4.5_to_5 = cbind(as(sigtab.3_to_3.5.4.5_to_5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.3_to_3.5.4.5_to_5), ], "matrix")) 

head(sigtab.3_to_3.5.4.5_to_5)
dim(sigtab.3_to_3.5.4.5_to_5)
# 11 differentially abundant SVs

write.csv(sigtab.3_to_3.5.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.3_to_3.5.4.5_to_5.csv")

# 3.5 to 4 versus 4 to 4.5
res.3.5_to_4.4_to_4.5 = results(diagdds, contrast = c("TidalRange_m", "3.5_to_4", "4_to_4.5"))
res.3.5_to_4.4_to_4.5 = res.3.5_to_4.4_to_4.5[order(res.3.5_to_4.4_to_4.5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3.5_to_4.4_to_4.5 = res.3.5_to_4.4_to_4.5[(res.3.5_to_4.4_to_4.5$padj < alpha), ]
sigtab.3.5_to_4.4_to_4.5 = cbind(as(sigtab.3.5_to_4.4_to_4.5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.3.5_to_4.4_to_4.5), ], "matrix")) 

head(sigtab.3.5_to_4.4_to_4.5)
dim(sigtab.3.5_to_4.4_to_4.5)
# 10 differentially abundant SVs

write.csv(sigtab.3.5_to_4.4_to_4.5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.3.5_to_4.4_to_4.5.csv")

# 3.5 to 4 versus 4.5 to 5
res.3.5_to_4.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "3.5_to_4", "4.5_to_5"))
res.3.5_to_4.4.5_to_5 = res.3.5_to_4.4.5_to_5[order(res.3.5_to_4.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.3.5_to_4.4.5_to_5 = res.3.5_to_4.4.5_to_5[(res.3.5_to_4.4.5_to_5$padj < alpha), ]
sigtab.3.5_to_4.4.5_to_5 = cbind(as(sigtab.3.5_to_4.4.5_to_5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.3.5_to_4.4.5_to_5), ], "matrix")) 

head(sigtab.3.5_to_4.4.5_to_5)
dim(sigtab.3.5_to_4.4.5_to_5)
# 4 differentially abundant SVs

write.csv(sigtab.3.5_to_4.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.3.5_to_4.4.5_to_5.csv")

# 4 to 4.5 versus 4.5 to 5
res.4_to_4.5.4.5_to_5 = results(diagdds, contrast = c("TidalRange_m", "4_to_4.5", "4.5_to_5"))
res.4_to_4.5.4.5_to_5 = res.4_to_4.5.4.5_to_5[order(res.4_to_4.5.4.5_to_5$padj, na.last=NA), ]
alpha = 0.01
sigtab.4_to_4.5.4.5_to_5 = res.4_to_4.5.4.5_to_5[(res.4_to_4.5.4.5_to_5$padj < alpha), ]
sigtab.4_to_4.5.4.5_to_5 = cbind(as(sigtab.4_to_4.5.4.5_to_5, "data.frame"), as(tax_table(ME_34.F.rsqmerg.1000)[rownames(sigtab.4_to_4.5.4.5_to_5), ], "matrix"))

head(sigtab.4_to_4.5.4.5_to_5)
dim(sigtab.4_to_4.5.4.5_to_5)
# 9 differentially abundant SVs

write.csv(sigtab.4_to_4.5.4.5_to_5, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.4_to_4.5.4.5_to_5.csv")


# KERATINOLYTIC GENERA - Maine NESP #####

Keratinolytic <- c("Bacillus", "Pseudomonas", "Enterococcus", "Staphylococcus", "Streptomyces", "Kocuria", "Arthrobacter", "Fervidobacterium", "Janthinobacterium", "Alcaligenes", "Flavobacterium", "Chryseobacterium", "Stenotrophomonas", "Xanthomonas", "Microbacterium", "Paenibacillus", "Terrabacter", "Caldioprobacter", "Thermoanaeobacter", "Actinomadura", "Amylocolatopsis", "Brevibacillus", "Fervidobacterium", "Terrabcter", "Serratia", "Nocardiosis", "Saccharomonospara", "Thermoactinomyces", "Saccarothrix")

ME.F.rar.k <- subset_taxa(ME_34.F.rsqmerg.1000.rar, Genus %in% Keratinolytic)
ME.F.rar.k <- prune_samples(sample_sums(ME.F.rar.k)>0, ME.F.rar.k)
ME.F.rar.k <- prune_taxa(taxa_sums(ME.F.rar.k)> 0, ME.F.rar.k)
# 76 taxa and 27 samples

saveRDS(ME.F.rar.k, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME.F.rar.k")

ME.F.rar.k <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME.F.rar.k")

##### ALPHA DIV: Keratinolytic taxa by Maine variables #####

###### Testing for alpha diversity variation by Maine variables #####

ME.F.rar.k.rich <- estimate_richness(ME.F.rar.k, measure=c("Observed", "Shannon")) 
ME.F.rar.k.rich.even <- ME.F.rar.k.rich$Shannon/log(ME.F.rar.k.rich$Observed)
ME.F.rar.k.rich.sd = as(sample_data(ME.F.rar.k), "matrix")
ME.F.rar.k.rich.sd = as.data.frame(ME.F.rar.k.rich.sd)
ME.F.rar.k.rich.df <- cbind(ME.F.rar.k.rich, ME.F.rar.k.rich.even, ME.F.rar.k.rich.sd)
write.csv(ME.F.rar.k.rich.df, file = "ME.F.rar.k.rich.df.csv")

# Steps: 1) Run ANOVA, 2) use shapiro-wilk to test for normality of residuals

# ANOVA of effect of sex on observed diversity
a_ME_F_sex.obs.k <- aov(Observed ~ Sex, data = ME.F.rar.k.rich.df)
plot(a_ME_F_sex.obs.k)
summary(a_ME_F_sex.obs.k)
# Df Sum Sq Mean Sq F value Pr(>F)
# Sex          1   46.8   46.80   2.134  0.156
# Residuals   25  548.2   21.93   
res.ME.sex.obs.k <- a_ME_F_sex.obs.k$residuals
shapiro.test(res.ME.sex.obs.k)
# W = 0.91007, p-value = 0.02291. non-normal, test with kruskal wallis
kruskal.test(Observed ~ Sex, data = ME.F.rar.k.rich.df)  
# Kruskal-Wallis chi-squared = 2.6872, df = 1, p-value = 0.1012
# non-significant effect of sex on observed diversity 

# ANOVA of effect of sex on shannon diversity
a_ME_F_sex.shan.k <- aov(Shannon ~ Sex, data = ME.F.rar.k.rich.df)
plot(a_ME_F_sex.shan.k)
summary(a_ME_F_sex.shan.k)
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Sex          1  0.754  0.7540   3.189 0.0863 .
# Residuals   25  5.911  0.2364  

res.ME.sex.shan.k <- a_ME_F_sex.shan.k$residuals
shapiro.test(res.ME.sex.shan.k)
# W = 0.97666, p-value = 0.7804. Normal, can use ANOVA results

# ANOVA of effect of month on observed diversity
a_ME_F_mo.obs.k <- aov(Observed ~ Month, data = ME.F.rar.k.rich.df)
plot(a_ME_F_mo.obs.k)
summary(a_ME_F_mo.obs.k)
#             Df Sum Sq Mean Sq F value Pr(>F)
# Month        2   71.4   35.71   1.637  0.216
# Residuals   24  523.5   21.81  
# No significant effect of month on observed diversity of keratinolytic bacteria
res.ME.mo.obs.k <- a_ME_F_mo.obs.k$residuals
shapiro.test(res.ME.mo.obs.k)
# W = 0.93983, p-value = 0.1207. normal, report ANOVA results

# ANOVA of effect of month on shannon diversity
a_ME_F_mo.shan.k <- aov(Shannon ~ Month, data = ME.F.rar.k.rich.df)
plot(a_ME_F_mo.shan.k)
summary(a_ME_F_mo.shan.k)
#             Df Sum Sq Mean Sq F value Pr(>F)  
# Month        2  1.263  0.6315   2.806 0.0804 .
# Residuals   24  5.402  0.2251  
# No significant effect of month on shannon diversity of keratinoltyic bacteria
res.ME.mo.shan.k <- a_ME_F_mo.shan.k$residuals
shapiro.test(res.ME.mo.shan.k)
# W = 0.9797, p-value = 0.8554. Normal, can use ANOVA results

# ANOVA of effect of marsh on observed diversity
a_ME_F_ma.obs.k <- aov(Observed ~ Marsh, data = ME.F.rar.k.rich.df)
plot(a_ME_F_ma.obs.k)
summary(a_ME_F_ma.obs.k)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Marsh        3   29.1    9.71   0.395  0.758
#Residuals   23  565.8   24.60   
# No significant effect of marsh on observed diversity of keratinolytic bacteria
res.ME.ma.obs.k <- a_ME_F_ma.obs.k$residuals
shapiro.test(res.ME.ma.obs.k)
# W = 0.9304, p-value = 0.0707. normal, report ANOVA results

# ANOVA of effect of marsh on shannon diversity
a_ME_F_ma.shan.k <- aov(Shannon ~ Marsh, data = ME.F.rar.k.rich.df)
plot(a_ME_F_ma.shan.k)
summary(a_ME_F_ma.shan.k)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Marsh        3  0.319  0.1062   0.385  0.765
#Residuals   23  6.346  0.2759    
# No significant effect of month on shannon diversity of keratinoltyic bacteria
res.ME.ma.shan.k <- a_ME_F_ma.shan.k$residuals
shapiro.test(res.ME.ma.shan.k)
# W = 0.97118, p-value = 0.6332. Normal, can use ANOVA results

# ANOVA of effect of tidal range (m) on observed diversity
a_ME_F_rng.obs.k <- aov(Observed ~ TidalRange_m, data = ME.F.rar.k.rich.df)
plot(a_ME_F_rng.obs.k)
summary(a_ME_F_rng.obs.k)
res.ME.rng.obs.k <- a_ME_F_rng.obs.k$residuals
shapiro.test(res.ME.rng.obs.k)
# W = 0.95476, p-value = 0.2791 Normal, can use ANOVA resultss
#                Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange_m  4  122.1   30.52    1.42   0.26
# Residuals    22  472.9   21.49  

# No significant effect of tidal range on observed diversity of keratinolytic bacteria

# ANOVA of effect of tide range on shannon diversity
a_ME_F_rng.shan.k <- aov(Shannon ~ TidalRange_m, data = ME.F.rar.k.rich.df)
plot(a_ME_F_rng.shan.k)
summary(a_ME_F_rng.shan.k)
#             Df Sum Sq Mean Sq F value Pr(>F)
# TidalRange_m  4  1.442  0.3605   1.519  0.231
# Residuals    22  5.223  0.2374    
# No significant effect of tidal range on shannon diversity of keratinoltyic bacteria
res.ME.rng.shan.k <- a_ME_F_rng.shan.k$residuals
shapiro.test(res.ME.rng.shan.k)
# W = 0.98456, p-value = 0.9476 Normal, can use ANOVA results


## BETA DIV: Keratinolytic taxa in Maine NESP ####

###### Homogeneity of variance - Jaccard #####
meta.df.ME.F.rar.k <- as(sample_data(ME.F.rar.k), "data.frame")

## Jaccard disatance between samples, for testing dispersion
dis.ME.F.k <- distance(ME.F.rar.k, method="jaccard")

# checking betadispersion of each variable

# sex
disper.sex.k <- betadisper(dis.ME.F.k, meta.df.ME.F.rar.k$Sex)
disper.sex.k
permutest(disper.sex.k)

#          Df  Sum Sq  Mean Sq    F   N.Perm Pr(>F)
#Groups     1 0.02638 0.026380 2.2719    999  0.152
#Residuals 25 0.29028 0.011611  
# not sig diff dispersions

# month
disper.month.k <- betadisper(dis.ME.F.k, meta.df.ME.F.rar.k$Month)
disper.month.k
permutest(disper.month.k)

#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     2 0.12530 0.062648 8.3525    999  0.003 **
#Residuals 24 0.18001 0.007500   
# Month has significantly different dispersion

month.tukey.16s.k <- TukeyHSD(disper.month.k)
#                diff         lwr         upr     p adj
# July-Aug   0.10493946 -0.03113271  0.24101163 0.1531731
# June-Aug  -0.04321544 -0.18963660  0.10320573 0.7441700
# June-July -0.14815490 -0.24180619 -0.05450361 0.00166783
plot(month.tukey.16s.k)


# marsh
disper.marsh.k <- betadisper(dis.ME.F.k, meta.df.ME.F.rar.k$Marsh)
disper.marsh.k
permutest(disper.marsh.k)

#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.30960 0.103198 10.444    999  0.001 ***
#Residuals 23 0.22726 0.009881     
# Marshes have diff dispersions

marsh.tukey.16s.k <- TukeyHSD(disper.marsh.k)

#           diff        lwr        upr     p adj
# HA-BB  0.04003463 -0.1845645  0.2646338 0.9597995
# NG-BB  0.06900979 -0.1360202  0.2740398 0.7884309
# PL-BB -0.50000000 -0.8368987 -0.1631013 0.0022714
# NG-HA  0.02897515 -0.1006972  0.1586475 0.9251649
# PL-HA -0.54003463 -0.8371514 -0.2429179 0.0002380
#PL-NG -0.56900979 -0.8516242 -0.2863954 0.0000637

plot(marsh.tukey.16s.k)

# tide cycle
disper.cycle.k <- betadisper(dis.ME.F.k, meta.df.ME.F.rar.k$TideCycle)
disper.cycle.k
permutest(disper.cycle.k)

#           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.019111 0.019111 1.7139    999  0.199
#Residuals 25 0.278763 0.011151                     
# Tide cycle doesn't have diff dispersions

# tidal range
disper.rng.k <- betadisper(dis.ME.F.k, meta.df.ME.F.rar.k$TidalRange_m)
disper.rng.k
permutest(disper.rng.k)

#           Df  Sum Sq   Mean Sq   F    N.Perm Pr(>F)  
# Groups     4 0.12765 0.031913 3.8012    999  0.013 *
# Residuals 22 0.18470 0.008396     
# Tidal range does have different dispersions

rng.tukey.16s.k <- TukeyHSD(disper.rng.k)
#                        diff         lwr         upr     p adj
# 3 to 3.5-2.5 to 3  0.23698832  0.01501834 0.458958305 0.0325621
# 3.5 to 4-2.5 to 3  0.19833309 -0.01658842 0.413254606 0.0799885
# 4 to 4.5-2.5 to 3  0.17274066 -0.07542932 0.420910649 0.2698634
# 4.5 to 5-2.5 to 3  0.10015608 -0.11476543 0.315077596 0.6445478
# 3.5 to 4-3 to 3.5 -0.03865523 -0.18547457 0.108164115 0.9333816
# 4 to 4.5-3 to 3.5 -0.06424766 -0.25647930 0.127983986 0.8562507
# 4.5 to 5-3 to 3.5 -0.13683224 -0.28365158 0.009987106 0.0757421
# 4 to 4.5-3.5 to 4 -0.02559243 -0.20964022 0.158455357 0.9934739
# 4.5 to 5-3.5 to 4 -0.09817701 -0.23410531 0.037751289 0.2381059
# 4.5 to 5-4 to 4.5 -0.07258458 -0.25663237 0.111463207 0.7677684

plot(rng.tukey.16s.k)


###### PERMANOVAs - Jaccard #####

# effect of sex, jaccard
perma.ME.F.k.sex <- adonis2(distance(ME.F.rar.k, method="jaccard") ~ Sex,
                            data = meta.df.ME.F.rar.k)

#        Df SumOfSqs      R2      F   Pr(>F)
# Sex       1   0.3911 0.04106 1.0705  0.355
# Residual 25   9.1336 0.95894              
# Total    26   9.5247 1.00000    

# effect of marsh, jaccard
perma.ME.F.k.marsh <- adonis2(distance(ME.F.rar.k, method="jaccard") ~ Marsh,
                           data = meta.df.ME.F.rar.k)

#          f SumOfSqs      R2     F Pr(>F)
#Marsh     3   1.2583 0.13211 1.167  0.193
#Residual 23   8.2663 0.86789             
#Total    26   9.5247 1.00000  

          
# effect of month, jaccard
perma.ME.F.k.month <- adonis2(distance(ME.F.rar.k, method="jaccard") ~ Month,
                            data = meta.df.ME.F.rar.k)

#          Df SumOfSqs      R2      F Pr(>F)
#Month     2   0.8836 0.09277 1.2271  0.155
#Residual 24   8.6410 0.90723              
#Total    26   9.5247 1.00000  

# effect of tide cycle, jaccard
perma.ME.F.k.cycle <- adonis2(distance(ME.F.rar.k, method="jaccard") ~ TideCycle,
                              data = meta.df.ME.F.rar.k)

#           Df SumOfSqs    R2    F    Pr(>F)
# TideCycle  2   0.7324 0.0769 0.9996  0.463
# Residual  24   8.7923 0.9231              
# Total     26   9.5247 1.0000    

# effect of tidal range, jaccard
perma.ME.F.k.rng <- adonis2(distance(ME.F.rar.k, method="jaccard") ~ TidalRange,
                              data = meta.df.ME.F.rar.k)

#            Df SumOfSqs      R2      F Pr(>F)
#TidalRange  4   1.6594 0.17422 1.1603  0.148
#Residual   22   7.8653 0.82578              
#Total      26   9.5247 1.00000 

###### Homogeneity of variance - Bray #####
meta.df.ME.F.rar.k <- as(sample_data(ME.F.rar.k), "data.frame")

## Bray disatance between samples, for testing dispersion
dis.ME.F.k.b <- distance(ME.F.rar.k, method="bray")

# Checking betadispersion of each variable

# sex
disper.sex.k.b <- betadisper(dis.ME.F.k.b, meta.df.ME.F.rar.k$Sex)
disper.sex.k.b
permutest(disper.sex.k.b)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.03603 0.036033 1.5383    999  0.246
# Residuals 25 0.58560 0.023424

# month
disper.month.k.b <- betadisper(dis.ME.F.k.b, meta.df.ME.F.rar.k$Month)
disper.month.k.b
permutest(disper.month.k.b)
#            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
# Groups     2 0.21779 0.10890 7.2742    999  0.002 **
# Residuals 24 0.35929 0.01497   

# marsh
disper.marsh.k.b <- betadisper(dis.ME.F.k.b, meta.df.ME.F.rar.k$Marsh)
disper.marsh.k.b
permutest(disper.marsh.k.b)
#            Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
# Groups     3 0.23717 0.079056 3.8787    999  0.016 *
# Residuals 23 0.46879 0.020382   

# range
disper.rng.k.b <- betadisper(dis.ME.F.k.b, meta.df.ME.F.rar.k$TidalRange_m)
disper.rng.k.b
permutest(disper.rng.k.b)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups     4 0.22486 0.056216 3.2265    999  0.031 *
# Residuals 22 0.38332 0.017423

# cycle
disper.cyc.k.b <- betadisper(dis.ME.F.k.b, meta.df.ME.F.rar.k$TideCycle)
disper.cyc.k.b
permutest(disper.cyc.k.b)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.06690 0.033452 1.5301    999  0.222
# Residuals 24 0.52471 0.021863      

### PERMANOVAs - Bray Curtis ####

# effect of sex, Bray
perma.ME.F.k.sex.b <- adonis2(distance(ME.F.rar.k, method="bray") ~ Sex,
                              data = meta.df.ME.F.rar.k)

#         Df SumOfSqs      R2     F   Pr(>F)
# Sex       1   0.3421 0.04354 1.1381  0.293
# Residual 25   7.5139 0.95646              
# Total    26   7.8560 1.00000    

# effect of month, Bray
perma.ME.F.k.month.b <- adonis2(distance(ME.F.rar.k, method="bray") ~ Month,
                              data = meta.df.ME.F.rar.k)

#         Df SumOfSqs      R2      F Pr(>F)
# Month     2   0.7355 0.09363 1.2396  0.201
# Residual 24   7.1205 0.90637              
# Total    26   7.8560 1.00000  

# effect of marsh, Bray
perma.ME.F.k.marsh.b <- adonis2(distance(ME.F.rar.k, method="bray") ~ Marsh,
                                data = meta.df.ME.F.rar.k)

#          Df SumOfSqs    R2      F    Pr(>F)
# Marsh     3   1.1497 0.14635 1.3144  0.153
# Residual 23   6.7063 0.85365              
# Total    26   7.8560 1.00000 

# effect of tidal range, Bray
perma.ME.F.k.rng.b <- adonis2(distance(ME.F.rar.k, method="bray") ~ TidalRange_m,
                            data = meta.df.ME.F.rar.k)

#              Df SumOfSqs      R2      F Pr(>F)  
# TidalRange_m  4   1.5786 0.20095 1.3832  0.078 .
# Residual     22   6.2773 0.79905                
# Total        26   7.8560 1.00000 

# effect of tide cycle, Bray
perma.ME.F.k.cyc.b <- adonis2(distance(ME.F.rar.k, method="bray") ~ TideCycle,
                              data = meta.df.ME.F.rar.k)

#            Df SumOfSqs    R2     F  Pr(>F)
# TideCycle  2   0.5628 0.07164 0.926  0.553
# Residual  24   7.2932 0.92836             
# Total     26   7.8560 1.00000   

# AIM 3: SEDIMENT AND FEATHERS COMPARISON - Maine-collected sediment and Maine NESP feathers #####
## ALPHA DIV: Maine NESP v Maine sediment ####

### Violin plot #####

viol.SF.rar <- plot_richness(ME_SF.rsqmerg.1000.rar, 
                                       x="Sample_Type", 
                                       measures=c("Observed","Shannon"), 
                                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Sample_Type)) + 
  geom_boxplot(width = 0.1, aes(group=Sample_Type)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/viol.SF.rar.png", width = 21, height = 21, units = 'cm', res = 300)
grid.arrange(viol.SF.rar) 
dev.off()

### Test for alpha div variation ####

SF.rar.rich <- estimate_richness(ME_SF.rsqmerg.1000.rar, measure=c("Observed", "Shannon")) 
SF.rar.rich.even <- SF.rar.rich$Shannon/log(SF.rar.rich$Observed)
SF.rar.rich.sd = as(sample_data(ME_SF.rsqmerg.1000.rar), "matrix")
SF.rar.rich.sd = as.data.frame(SF.rar.rich.sd)
SF.rar.rich.df <- cbind(SF.rar.rich, SF.rar.rich.even, SF.rar.rich.sd)
write.csv(SF.rar.rich.df, file = "SF.rar.rich.df.csv")

#check the distribution of your data, which will change your stats approach for testing for variation
shapiro.test(SF.rar.rich.df$Shannon)
# W = 0.87425, p-value = 0.0004367 NOT NORMAL
shapiro.test(SF.rar.rich.df$Observed)
# W = 0.87527, p-value = 0.0004642 NOT NORMAL

# Kruskal-Wallis tests for effect of sample type on observed and shannon's 

kruskal.test(Observed ~ Sample_Type, data = SF.rar.rich.df)
# Kruskal-Wallis chi-squared = 21.116, df = 1, p-value = 4.322e-06
# Highly sig

kruskal.test(Shannon ~ Sample_Type, data = SF.rar.rich.df)
# Kruskal-Wallis chi-squared = 24.001, df = 1, p-value = 9.629e-07
# Highly sig

## BETA DIVERSITY - FEATHERS AND SEDIMENT ####
#### PCOA ####

# Jaccard dissimilarity
# Sample type 
ord_SF.rar <- ordinate(ME_SF.rsqmerg.1000.rar, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF <- plot_ordination(ME_SF.rsqmerg.1000.rar, ord_SF.rar, type="samples", color="Sample_Type") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/PCOA.SF.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF) 
dev.off()

# Bray Curtis dissimilarity
# Sample type 
ord_SF.rar.b <- ordinate(ME_SF.rsqmerg.1000.rar, #calculate similarities
                       method ="PCoA", #ordination type
                       "bray") #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

PCOA.SF.b <- plot_ordination(ME_SF.rsqmerg.1000.rar, ord_SF.rar.b, type="samples", color="Sample_Type") + geom_point(size = 3)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/PCOA.SF.bray.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF.b) 
dev.off()


### DISPERSION AND COMMUNITY COMPOSITION - FEATHERDS AND SEDIMENT ####

##### HOMOGENEIYTY OF VARIANCE - Jaccard ####

meta.df.SF <- as(sample_data(ME_SF.rsqmerg.1000.rar), "data.frame") # get metadata

## Jaccard disatance between samples, for testing dispersion
dis.SF <- distance(ME_SF.rsqmerg.1000.rar, method="jaccard")

# checking betadispersion by sample type

disper.sediment <- betadisper(dis.SF, meta.df.SF$Sample_Type)
disper.sediment
permutest(disper.sediment)

#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.001238 0.0012381 0.3725    999  0.545
#Residuals 37 0.122969 0.0033235    
# Sample types don't have different dispersions, can trust that variation in PERMANOVA is due to differences in centroids and not dispersions

##### HOMOGENEIYTY OF VARIANCE - Bray ####

meta.df.SF <- as(sample_data(ME_SF.rsqmerg.1000.rar), "data.frame") # get metadata

## Jaccard disatance between samples, for testing dispersion
dis.SF.b <- distance(ME_SF.rsqmerg.1000.rar, method="jaccard")

# checking betadispersion by sample type

disper.sediment.b <- betadisper(dis.SF.b, meta.df.SF$Sample_Type)
disper.sediment.b
permutest(disper.sediment.b)

#            Df   Sum Sq   Mean Sq   F     N.Perm Pr(>F)
# Groups     1 0.001368 0.0013679 0.4004    999  0.527
# Residuals 38 0.129827 0.0034165

##### COMMUNITY COMPOSITION - Jaccard ####
# Effect of sample type on community composition

perma.SF <- adonis2(distance(ME_SF.rsqmerg.1000.rar, method="jaccard") ~ Sample_Type,
                           data = meta.df.SF)

#              Df SumOfSqs   R2      F    Pr(>F)    
# Sample_Type  1   1.8313 0.10961 4.5546  0.001 ***
# Residual    37  14.8770 0.89039                  
# Total       38  16.7083 1.00000  

##### COMMUNITY COMPOSITION - Bray ####
# Effect of sample type on community composition

perma.SF.b <- adonis2(distance(ME_SF.rsqmerg.1000.rar, method="bray") ~ Sample_Type,
                    data = meta.df.SF)

#             Df SumOfSqs    R2     F    Pr(>F)    
# Sample_Type  1    2.716 0.1742 8.0158  0.001 ***
# Residual    38   12.875 0.8258                  
# Total       39   15.591 1.0000    

### DIFFERENTIALLY ABUNDANT TAXA - DESEQ ####
# ME soil and ME feathers
library("DESeq2")
packageVersion("DESeq2")

# DESEQ Differential abundance - ME soil and ME feathers
# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(ME_SF.rsqmerg.1000, ~ Sample_Type) 

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
sigtab.FvS = cbind(as(sigtab.FvS, "data.frame"), as(tax_table(ME_SF.rsqmerg.1000)[rownames(sigtab.FvS), ], "matrix"))

head(sigtab.FvS)
dim(sigtab.FvS)
# 212 differentially abundant taxa between feathers and sediment

write.csv(sigtab.FvS, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_SavedFiles/DESEQ.sigtab.FvS.csv")

# Phylum order
x = tapply(sigtab.FvS$log2FoldChange, sigtab.FvS$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.FvS$Phylum = factor(as.character(sigtab.FvS$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab.FvS$log2FoldChange, sigtab.FvS$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.FvS$Genus = factor(as.character(sigtab.FvS$Genus), levels=names(x))

## if the Genus is empty, replace with the Family
sigtab.FvS$Genus[is.na(sigtab.FvS$Genus)] <- sigtab.FvS$Family[is.na(sigtab.FvS$Genus)]


library("ggplot2")

## graph differential abundance
DESEQ.SvF <- ggplot(sigtab.FvS, aes(y=Genus, x=log2FoldChange, color=Genus)) + #play with aestetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) +
  theme(text=element_text(family="Times New Roman", size=12)) + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 
  

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/CodeForManuscript/16S_Figs/DESEQ.SvF.png", width = 40, height = 30, units = 'cm', res = 300)
grid.arrange(DESEQ.SvF) 
dev.off()

# KERATINOLYTIC GENERA  - ME SEDIMENT ONLY ####

Keratinolytic <- c("Bacillus", "Pseudomonas", "Enterococcus", "Staphylococcus", "Streptomyces", "Kocuria", "Arthrobacter", "Fervidobacterium", "Janthinobacterium", "Alcaligenes", "Flavobacterium", "Chryseobacterium", "Stenotrophomonas", "Xanthomonas", "Microbacterium", "Paenibacillus", "Terrabacter", "Caldioprobacter", "Thermoanaeobacter", "Actinomadura", "Amylocolatopsis", "Brevibacillus", "Fervidobacterium", "Terrabcter", "Serratia", "Nocardiosis", "Saccharomonospara", "Thermoactinomyces", "Saccarothrix")

ME_S.k <- subset_taxa(ME_34_S.rar, Genus %in% Keratinolytic)
ME_S.k <- prune_taxa(taxa_sums(ME_S.k)> 0, ME_S.k)
ME_S.k <- prune_samples(sample_sums(ME_S.k)>0, ME_S.k)
# 29 taxa and 10 samples

saveRDS(ME_S.k, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_S.k.rds")

ME_S.rar.k <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_S.k.rds")

# KERATINOLYTIC GENERA  - ME SEDIMENT AND FEATHERS ####

Keratinolytic <- c("Bacillus", "Pseudomonas", "Enterococcus", "Staphylococcus", "Streptomyces", "Kocuria", "Arthrobacter", "Fervidobacterium", "Janthinobacterium", "Alcaligenes", "Flavobacterium", "Chryseobacterium", "Stenotrophomonas", "Xanthomonas", "Microbacterium", "Paenibacillus", "Terrabacter", "Caldioprobacter", "Thermoanaeobacter", "Actinomadura", "Amylocolatopsis", "Brevibacillus", "Fervidobacterium", "Terrabcter", "Serratia", "Nocardiosis", "Saccharomonospara", "Thermoactinomyces", "Saccarothrix")

ME_SF.rar.k <- subset_taxa(ME_SF.rsqmerg.1000.rar, Genus %in% Keratinolytic)
ME_SF.rar.k <- prune_taxa(taxa_sums(ME_SF.rar.k)> 0, ME_SF.rar.k)
ME_SF.rar.k <- prune_samples(sample_sums(ME_SF.rar.k)>0, ME_SF.rar.k)
# 93 taxa and 36 samples

saveRDS(ME_SF.rar.k, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rar.k")

ME_SF.rar.k <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/ME_SF.rsqmerg.1000.rar.k")

## ALPHA DIV -  keratinolytic taxa in feathers and sediment ####

### VIOLIN - keratinolytic taxa in sediment and feathers ####

viol.SF.k <- plot_richness(ME_SF.rar.k, 
                                     x="Sample_Type", 
                                     measures=c("Observed","Shannon"), 
                                     title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Sample_Type)) + 
  geom_boxplot(width = 0.1, aes(group=Sample_Type)) + 
  ylab("Observed Bacterial Richness (SVs)")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/viol.SF.k.png", width = 23, height = 21, units = 'cm', res = 300)
grid.arrange(viol.SF.k) 
dev.off()

### Testing for differences in alpha div ####

SF.rar.k.rich <- estimate_richness(ME_SF.rar.k, measure=c("Observed", "Shannon")) 
SF.rar.k.rich.even <- SF.rar.k.rich$Shannon/log(SF.rar.k.rich$Observed)
SF.rar.k.rich.sd = as(sample_data(ME_SF.rar.k), "matrix")
SF.rar.k.rich.sd = as.data.frame(SF.rar.k.rich.sd)
SF.rar.k.rich.df <- cbind(SF.rar.k.rich, SF.rar.k.rich.even, SF.rar.k.rich.sd)
write.csv(SF.rar.k.rich.df, file = "SF.rar.k.rich.df.csv")

#check the distribution of your data, which will change your stats approach for testing for variation
shapiro.test(SF.rar.k.rich.df$Shannon)
# W = 0.93957, p-value = 0.04925 NOT NORMAL
shapiro.test(SF.rar.k.rich.df$Observed)
# W = 0.9186, p-value = 0.01145 NOT NORMAL

# Kruskal-Wallis tests for the non-normally distributed observed and shannon div

kruskal.test(Observed ~ Sample_Type, data = SF.rar.k.rich.df)
# Kruskal-Wallis chi-squared = 14.811, df = 1, p-value = 0.0001189

kruskal.test(Shannon ~ Sample_Type, data = SF.rar.k.rich.df)
# Kruskal-Wallis chi-squared = 14.839, df = 1, p-value = 0.0001171

# Both Observed and Shannon div sig diff by sample type

### Feathers versus sediment violin plots - all taxa versus keratinolytic ####

viol.sf.all <- grid.arrange(viol.SF.rar, viol.SF.k)

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/viol.sf.all.png", width = 23, height = 21, units = 'cm', res = 300)
grid.arrange(viol.sf.all) 
dev.off()

## BETA DIVERSITY - FEATHERS V SEDIMENT, KERATINOLYTIC SVS ####

### HOMOGENEITY OF VARIANCE - JACCARD #####
meta.df.SF.k <- as(sample_data(ME_SF.rar.k), "data.frame") # get sample data

## Jaccard disatance between samples, for testing dispersion
dis.S.k <- distance(ME_SF.rar.k, method="jaccard")

# checking betadispersion by sample type
disper.sediment.k <- betadisper(dis.S.k, meta.df.SF.k$Sample_Type)
disper.sediment.k
permutest(disper.sediment.k)

#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.01352 0.013518 1.2395    999  0.266
#Residuals 34 0.37081 0.010906     
# Sample types don't have different dispersions, can trust that variation in PERMANOVA is due to differences in centroids and not dispersions

### HOMOGENEITY OF VARIANCE - BRAY-CURTIS #####
meta.df.SF.k <- as(sample_data(ME_SF.rar.k), "data.frame") # get sample data

## Bray-curtis disatance between samples, for testing dispersion
dis.S.k.b <- distance(ME_SF.rar.k, method="bray")

# checking betadispersion by sample type
disper.sediment.k.b <- betadisper(dis.S.k.b, meta.df.SF.k$Sample_Type)
disper.sediment.k.b
permutest(disper.sediment.k.b)

#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.05369 0.053692 2.5262    999  0.135
# Residuals 34 0.72264 0.021254     
# Sample types don't have different dispersions, can trust that variation in PERMANOVA is due to differences in centroids and not dispersions

### PERMANOVA - JACCARD #####
# effect of sample type, jaccard
perma.SF.k <- adonis2(distance(ME_SF.rar.k, method="jaccard") ~ Sample_Type,
                    data = meta.df.SF.k)

#             Df SumOfSqs      R2    F    Pr(>F)    
# Sample_Type  1   1.6299 0.11009 4.2062  0.001 ***
# Residual    34  13.1750 0.88991                  
# Total       35  14.8049 1.00000  

### PERMANOVA - BRAY-CURTIS #####
# effect of sample type, bray curtis
perma.SF.k.b <- adonis2(distance(ME_SF.rar.k, method="bray") ~ Sample_Type,
                      data = meta.df.SF.k)

#             Df SumOfSqs      R2    F    Pr(>F)    
# Sample_Type  1   2.1982 0.16281 6.6122  0.001 ***
# Residual    34  11.3033 0.83719                  
# Total       35  13.5015 1.00000  

## CORE AND TRANSIENT TAXA #####
library(microbiome)
### Maine sediment core and transient ####
#### core >60% prevalence threshold ####
# using the unrarefied phyloseq with the reseqs/recaps combined, >1000 reads

# keep only taxa with positive sums
core_biomeS <- prune_taxa(taxa_sums(ME_34_S) > 0, ME_34_S)

#calcuate compositional version of the data (relative abundances)
core_biomeS.rel <- microbiome::transform(core_biomeS, "compositional")

#This returns the taxa that exceed the given prevalence and detection thresholds.
core_taxa_S <- core_members(core_biomeS.rel, detection = 1/10000, prevalence = 60/100) # 12 taxa

#A full phyloseq object of the core microbiota is obtained as follows
phylo.coreS.6 <- core(core_biomeS.rel, detection = 1/10000, prevalence = .6)
phylo.coreS.6 <- prune_samples(sample_sums(phylo.coreS.6)>0, phylo.coreS.6)
phylo.coreS.6 <- prune_taxa(taxa_sums(phylo.coreS.6)>0, phylo.coreS.6)
# 89 SVs, 12 samples

saveRDS(phylo.coreS.6, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.coreS.6.rds")

phylo.coreS.6 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.coreS.6.rds")

#### transient <20% prevalence threshold ####

phylo.rareS.2 <- rare(core_biomeS.rel, detection=1/10000, prevalence=20/100) # 0.1% detection threshold, <20% of samples 
phylo.rareS.2 <- prune_samples(sample_sums(phylo.rareS.2)>0, phylo.rareS.2)
phylo.rareS.2 <- prune_taxa(taxa_sums(phylo.rareS.2)>0, phylo.rareS.2)
# 5153 taxa and 12 samples 

saveRDS(phylo.rareS.2, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.rareS.2.rds")

phylo.rareS.2 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.rareS.2.rds")

#### core >50% prevalence threshold ####

phylo.coreS.5 <- core(core_biomeS.rel, detection = 1/10000, prevalence = .5)
phylo.coreS.5 <- prune_samples(sample_sums(phylo.coreS.5)>0, phylo.coreS.5)
phylo.coreS.5 <- prune_taxa(taxa_sums(phylo.coreS.5)>0, phylo.coreS.5)
# 147 taxa and 12 samples

#### transient <10% prevalence threshold ####

phylo.rareS.1 <- rare(core_biomeS.rel, detection=1/10000, prevalence=10/100) # 0.1% detection threshold, <10% of samples 
phylo.rareS.1 <- prune_samples(sample_sums(phylo.rareS.1)>0, phylo.rareS.1)
phylo.rareS.1 <- prune_taxa(taxa_sums(phylo.rareS.1)>0, phylo.rareS.1)
# 4098 taxa and 12 samples

### Maine feather Core and transient ####
# using the unrarefied phyloseq with the reseqs/recaps combined, >1000 reads

# keep only taxa with positive sums
core_biomeF.Me.rsemerg <- prune_taxa(taxa_sums(ME_34.F.rsqmerg.1000) > 0, ME_34.F.rsqmerg.1000)

#calcuate compositional version of the data (relative abundances)
core_biomeF.Me.rsqmerg.rel <- microbiome::transform(core_biomeF.Me.rsemerg, "compositional")

#This returns the taxa that exceed the given prevalence and detection thresholds.
core_taxa_F.Me.rsqmerg <- core_members(core_biomeF.Me.rsqmerg.rel, detection = 1/10000, prevalence = 60/100) # 8 taxa

#A full phyloseq object of the core microbiota is obtained as follows
phylo.coreF.Me.rsqmerg.6 <- core(core_biomeF.Me.rsqmerg.rel, detection = 1/10000, prevalence = .6)
phylo.coreF.Me.rsqmerg.6 <- prune_samples(sample_sums(phylo.coreF.Me.rsqmerg.6)>0, phylo.coreF.Me.rsqmerg.6)
phylo.coreF.Me.rsqmerg.6 <- prune_taxa(taxa_sums(phylo.coreF.Me.rsqmerg.6)>0, phylo.coreF.Me.rsqmerg.6)
# 8 taxa and 24 samples

saveRDS(phylo.coreF.Me.rsqmerg.6, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.coreF.Me.rsqmerg.6.rds")

plot_bar(phylo.coreF.Me.rsqmerg.6, fill = "Genus")

sort(table(tax_table(phylo.coreF.Me.rsqmerg.6)[, 6]))
sort(table(tax_table(phylo.coreF.Me.rsqmerg.6)[, 7]))

write.csv(tax_table(phylo.coreF.Me.rsqmerg.6), "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/tax_feathcore_60.csv")

#### transient <20% prevalence threshold ####

phylo.rareF.Me.rsqmerg.2 <- rare(core_biomeF.Me.rsqmerg.rel, detection=1/10000, prevalence=20/100) # 0.1% detection threshold, <20% of samples 
phylo.rareF.Me.rsqmerg.2 <- prune_samples(sample_sums(phylo.rareF.Me.rsqmerg.2)>0, phylo.rareF.Me.rsqmerg.2)
phylo.rareF.Me.rsqmerg.2 <- prune_taxa(taxa_sums(phylo.rareF.Me.rsqmerg.2)>0, phylo.rareF.Me.rsqmerg.2)
# 2281 taxa and 27 samples 

saveRDS(phylo.rareF.Me.rsqmerg.2, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.rareF.Me.rsqmerg.2")

phylo.rareF.Me.rsqmerg.2 <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phylo.rareF.Me.rsqmerg.2")


#### core >50% prevalence threshold ####

phylo.coreF.Me.rsqmerg.5 <- core(core_biomeF.Me.rsqmerg.rel, detection = 1/10000, prevalence = .5)
phylo.coreF.Me.rsqmerg.5 <- prune_samples(sample_sums(phylo.coreF.Me.rsqmerg.5)>0, phylo.coreF.Me.rsqmerg.5)
phylo.coreF.Me.rsqmerg.5 <- prune_taxa(taxa_sums(phylo.coreF.Me.rsqmerg.5)>0, phylo.coreF.Me.rsqmerg.5)
# 33 taxa and 26 samples

#### transient <10% prevalence threshold ####

phylo.rareF.Me.rsqmerg.1 <- rare(core_biomeF.Me.rsqmerg.rel, detection=1/10000, prevalence=10/100) # 0.1% detection threshold, <10% of samples 
phylo.rareF.Me.rsqmerg.1 <- prune_samples(sample_sums(phylo.rareF.Me.rsqmerg.1)>0, phylo.rareF.Me.rsqmerg.1)
phylo.rareF.Me.rsqmerg.1 <- prune_taxa(taxa_sums(phylo.rareF.Me.rsqmerg.1)>0, phylo.rareF.Me.rsqmerg.1)
# 2038 taxa and 27 samples


### Exploring sediment core: >60% prevalence, Transient: <20% prevalence ####
##### Phyloseq of the sediment core, transient, and entire microbial communities ####

# Adding a column with community type to each of the phyloseqs
sample_data(phylo.coreS.6)$CommunityType <- "Core"
sample_data(phylo.rareS.2)$CommunityType <- "Peripheral"
sample_data(core_biomeS.rel)$CommunityType <- "Entire"

# EDITING SAMPLE NAMES
sample_names(phylo.coreS.6) <- paste("C", sample_names(phylo.coreS.6))
sample_names(phylo.rareS.2) <- paste("P", sample_names(phylo.rareS.2))
sample_names(core_biomeS.rel) <- paste("E", sample_names(core_biomeS.rel))

phy_S.60C.20T.E <- merge_phyloseq(phylo.coreS.6, phylo.rareS.2, core_biomeS.rel)
phy_S.60C.20T.E <- prune_samples(sample_sums(phy_S.60C.20T.E)>0, phy_S.60C.20T.E)
phy_S.60C.20T.E <- prune_taxa(taxa_sums(phy_S.60C.20T.E)>0, phy_S.60C.20T.E)
# 6226 taxa and 36 samples

saveRDS(phy_S.60C.20T.E, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_S.60C.20T.E.rds")

##### Soil taxonomy plot of 60% core v 20% transient v entire communities ####

phy_S.60C.20T.E.rel <- transform_sample_counts(phy_S.60C.20T.E, function(OTU) OTU/sum(OTU))
phy_S.60C.20T.E.rel.glom <- tax_glom(phy_S.60C.20T.E.rel , taxrank = "Phylum")
plot_bar(phy_S.60C.20T.E.rel.glom, fill = "Phylum")

tax_S.core6.trans2.entire <- plot_bar(phy_S.60C.20T.E.rel.glom, fill = "Phylum") +
  facet_grid(~CommunityType, scales = "free", space = "free")


# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_S.core6.trans2.entire.png", width = 35, height = 28, units = 'cm', res = 300)
grid.arrange(tax_S.core6.trans2.entire) 
dev.off()

##### Phyloseq sedmiment core and rare, entire removed #####
phy_S.60C.20T <- subset_samples(phy_S.60C.20T.E, CommunityType != "Entire")
phy_S.60C.20T <- prune_taxa(taxa_sums(phy_S.60C.20T)>0, phy_S.60C.20T)
# 5242 taxa and 24 samples

saveRDS(phy_S.60C.20T, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_S.60C.20T.rds")


##### Tax plot of transient feather and sediment microbiota #####

phy_transientSF <- merge_phyloseq(phylo.rareS.2, phylo.rareF.Me.rsqmerg.2)
# 7121 taxa and 39 samples

phy_transientSF.rel <- transform_sample_counts(phy_transientSF, function(OTU) OTU/sum(OTU))

phy_transientSF.rel.glom <- tax_glom(phy_transientSF.rel, taxrank = "Phylum")

plot_bar(phy_transientSF.rel.glom, fill = "Phylum") + facet_grid(~Sample_Type, scales = "free", space = "free")

##### PCOA of all (community types assigned) #####

# Jaccard, Presence/absence
ord.S.60C.20T.jac <- ordinate(phy_S.60C.20T.E, #calculate similarities
                               method ="PCoA", #ordination type
                               "jaccard") # uses presence/absence only

PCOA_S.core6.trans2.entire <- plot_ordination(phy_S.60C.20T.E, ord.S.60C.20T.jac, type="samples", color="CommunityType") + geom_point(size = 2.2)
# Axis 1: 14.3%, Axis 2: 6.8%


# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/PCOA_S.core6.trans2.entire.png", width = 23, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA_S.core6.trans2.entire) 
dev.off()


### Exploring feather core: >60% prevalence, Peripheral: <20% prevalence ####
##### Phyloseq of the NESP core, transient, and entire microbial communities ####

# Adding a column with community type to each of the phyloseqs
sample_data(phylo.coreF.Me.rsqmerg.6)$CommunityType <- "Core"
sample_data(phylo.rareF.Me.rsqmerg.2)$CommunityType <- "Peripheral"
sample_data(core_biomeF.Me.rsqmerg.rel)$CommunityType <- "Entire"

# EDITING SAMPLE NAMES
sample_names(phylo.coreF.Me.rsqmerg.6) <- paste("C", sample_names(phylo.coreF.Me.rsqmerg.6))
sample_names(phylo.rareF.Me.rsqmerg.2) <- paste("P", sample_names(phylo.rareF.Me.rsqmerg.2))
sample_names(core_biomeF.Me.rsqmerg.rel) <- paste("E", sample_names(core_biomeF.Me.rsqmerg.rel))

phy_ME.F.60C.20T.E.rsqmerg <- merge_phyloseq(phylo.coreF.Me.rsqmerg.6, phylo.rareF.Me.rsqmerg.2, core_biomeF.Me.rsqmerg.rel)
phy_ME.F.60C.20T.E.rsqmerg <- prune_samples(sample_sums(phy_ME.F.60C.20T.E.rsqmerg)>0, phy_ME.F.60C.20T.E.rsqmerg)
phy_ME.F.60C.20T.E.rsqmerg <- prune_taxa(taxa_sums(phy_ME.F.60C.20T.E.rsqmerg)>0, phy_ME.F.60C.20T.E.rsqmerg)
# 2474 taxa and 78 samples

saveRDS(phy_ME.F.60C.20T.E.rsqmerg, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_ME.F.60C.20T.E.rsqmerg.rds")

##### Taxonomy plot of 60% core v 20% transient v entire communities ####

phy_ME.F.60C.20T.E.rsqmerg.rel <- transform_sample_counts(phy_ME.F.60C.20T.E.rsqmerg, function(OTU) OTU/sum(OTU))
phy_ME.F.60C.20T.E.rsqmerg.rel.glom <- tax_glom(phy_ME.F.60C.20T.E.rsqmerg.rel, taxrank = "Phylum")
plot_bar(phy_ME.F.60C.20T.E.rsqmerg.rel.glom, fill = "Phylum")

tax_F.core6.trans2.entire <- plot_bar(phy_ME.F.60C.20T.E.rsqmerg.rel.glom, fill = "Phylum") +
  facet_grid(~CommunityType, scales = "free", space = "free")
  
##### Phyloseq of the feather core and peripheral commmunities, the entire removed #####
phy_ME.F.60C.20T.rsqmerg <- subset_samples(phy_ME.F.60C.20T.E.rsqmerg, CommunityType != "Entire")
phy_ME.F.60C.20T.rsqmerg <- prune_taxa(taxa_sums(phy_ME.F.60C.20T.rsqmerg)>0, phy_ME.F.60C.20T.rsqmerg)
# 1760 taxa and 55 samples

saveRDS(phy_ME.F.60C.20T.rsqmerg, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_ME.F.60C.20T.rsqmerg.rds")

phy_ME.F.60C.20T.rsqmerg.rel <- transform_sample_counts(phy_ME.F.60C.20T.rsqmerg, function(OTU) OTU/sum(OTU))
phy_ME.F.60C.20T.rsqmerg.rel.glom <- tax_glom(phy_ME.F.60C.20T.rsqmerg.rel, taxrank = "Phylum")
plot_bar(phy_ME.F.60C.20T.rsqmerg.rel.glom, fill = "Phylum")

tax_F.core6.trans2 <- plot_bar(phy_ME.F.60C.20T.rsqmerg.rel.glom, fill = "Phylum") +
  facet_grid(~CommunityType, scales = "free", space = "free")


# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_F.core6.trans2.entire.png", width = 35, height = 28, units = 'cm', res = 300)
grid.arrange(tax_F.core6.trans2.entire) 
dev.off()

## Comparison: feather and sediment core and transient microbes ####

#First, merge the feather and soil phyloseqs. Starting with the >60% prevalence for core and <20% prevalence for transient

phy_FS.60C.20T.E <- merge_phyloseq(phy_ME.F.60C.20T.E.rsqmerg, phy_S.60C.20T.E)
phy_FS.60C.20T.E <- prune_taxa(taxa_sums(phy_FS.60C.20T.E)>0, phy_FS.60C.20T.E)
phy_FS.60C.20T.E <- prune_samples(sample_sums(phy_FS.60C.20T.E)>0, phy_FS.60C.20T.E)
# 7696 taxa and 119 samples 

# Phyloseq of just the soil and feather core and transient, entire removed

phy_FS.60C.20T <- subset_samples(phy_FS.60C.20T.E, CommunityType != "Entire")
phy_FS.60C.20T <- prune_taxa(taxa_sums(phy_FS.60C.20T)>0, phy_FS.60C.20T)
phy_FS.60C.20T <- prune_samples(sample_sums(phy_FS.60C.20T)>0, phy_FS.60C.20T)
# 6734 taxa and 79 samples

saveRDS(phy_FS.60C.20T, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_FS.60C.20T.rds")

### PCOA  ####
# Jaccard, Presence/absence
ord.FS.60C.20T.jac <- ordinate(phy_FS.60C.20T.E, #calculate similarities
                               method ="PCoA", #ordination type
                               "jaccard", binary = "TRUE") # uses presence/absence only

PCOA.SF.Core6.Trans2.Entire <- plot_ordination(phy_FS.60C.20T.E, ord.FS.60C.20T.jac, type="samples", color="CommunityType", shape = "Sample_Type") + geom_point(size = 2.2)
# Axis 1: 14.3%, Axis 2: 8.6%

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/PCOA.SF.Core6.Trans2.Entire.png", width = 26, height = 23, units = 'cm', res = 300)
grid.arrange(PCOA.SF.Core6.Trans2.Entire) 
dev.off()

### Taxonomy plot ####

phy_FS.60C.20T.rel <- transform_sample_counts(phy_FS.60C.20T, function(OTU) OTU/sum(OTU))
phy_FS.60C.20T.rel.glom <- tax_glom(phy_FS.60C.20T.rel, taxrank = "Phylum")
plot_bar(phy_FS.60C.20T.rel.glom, fill = "Phylum")

tax_SF.core6.trans2 <- plot_bar(phy_FS.60C.20T.rel.glom, fill = "Phylum") +
  facet_grid(CommunityType ~ Sample_Type, scales = "free", space = "free")

# Export figure
require(gridExtra)

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_SF.core6.trans2.png", width = 47, height = 35, units = 'cm', res = 300)
grid.arrange(tax_SF.core6.trans2) 
dev.off()

# same taxonomy plot as above, removing the sample names

tax_SF.core6.trans2.fncy <- plot_bar(phy_FS.60C.20T.rel.glom, fill = "Phylum") +
  facet_grid(CommunityType ~ Sample_Type, scales = "free", space = "free") + 
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(text=element_text(family="Times New Roman", size=13))


png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_SF.core6.trans2.fncy.png", width = 40, height = 30, units = 'cm', res = 300)
grid.arrange(tax_SF.core6.trans2.fncy) 
dev.off()

# combining rows to save space

tax_SF.core6.trans2.fncy2 <- plot_bar(phy_FS.60C.20T.rel.glom, fill = "Phylum") +
  facet_wrap(~ Sample_Type, scales = "free") + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())

png("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/Figs_103122/tax_SF.core6.trans2.fncy.png", width = 40, height = 30, units = 'cm', res = 300)
grid.arrange(tax_SF.core6.trans2.fncy) 
dev.off()


### Phyloseq of just the feather and sediment core microbiota (60% prevalence) #####

phy_FS.60C <- subset_samples(phy_FS.60C.20T, CommunityType == "Core")
phy_FS.60C <- prune_taxa(taxa_sums(phy_FS.60C)>0, phy_FS.60C)
# 96 taxa and 36 samples

saveRDS(phy_FS.60C, "/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_FS.60C")

phy_FS.60C <- readRDS("/Users/alicehotopp/OneDrive - University of Maine System/2021/Feather_16S_ITS/Sequencing_Fall2021/Runs_Combined_Analysis/16S/phy_FS.60C")

#### Taxonomy plot of feather versus sediment core microbiota (60% prevalence)

phy_FS.60C.rel <- transform_sample_counts(phy_FS.60C, function(OTU) OTU/sum(OTU))
phy_FS.60C.rel.glom <- tax_glom(phy_FS.60C.rel, taxrank = "Phylum")
plot_bar(phy_FS.60C.rel.glom, fill = "Phylum")

tax_SF.core60 <- plot_bar(phy_FS.60C.rel.glom, fill = "Phylum") +
  facet_grid(~ Sample_Type, scales = "free", space = "free")

### Venn Diagram of shared taxa between sediment and feather transient ####

library(ggVennDiagram)

# Accessing the taxonomy tables of the transient feather and transient soil phyloseqs
tax_20.F.trans <- as.data.frame(tax_table(phylo.rareF.Me.rsqmerg.2))
tax_20.S.trans <- as.data.frame(tax_table(phylo.rareS.2))

# Replace Genus with Family if the Genus is NA

#feathers
tax_20.F.trans$Genus[is.na(tax_20.F.trans$Genus)] <- tax_20.F.trans$Family[is.na(tax_20.F.trans$Genus)]

# soil
tax_20.S.trans$Genus[is.na(tax_20.S.trans$Genus)] <- tax_20.S.trans$Family[is.na(tax_20.S.trans$Genus)]

# Adding a Genus.species column to the tax tables 

# feathers
tax_20.F.trans$Genus.species <- paste(tax_20.F.trans$Genus, tax_20.F.trans$Species)

# soil
tax_20.S.trans$Genus.species <- paste(tax_20.S.trans$Genus, tax_20.S.trans$Species)

# Finding the shared SVs between the transient soil and the transient feathers 
library(ggVennDiagram)

# Venn of the Genus.species combos shared
VennTrans <- list(Feather = tax_20.F.trans$Genus.species, Sediment = tax_20.S.trans$Genus.species)

ggVennDiagram(VennTrans)

# 269 (28%) Genus.species unique to transient feathers
# 349 (36%) Genus.species unique to transient soil
# 352 (36%) Genus.species SHARED between transient feathers and transient soil

# ASV level
Venn.Trans.ASV <- list(Feather = rownames(tax_20.F.trans), Sediment = rownames(tax_20.S.trans))

ggVennDiagram(Venn.Trans.ASV)

# 1968 (28%) SVs unique to feathers
# 4840 (68%) SVs unique to soil
# 313 (4%) SVs shared between feathers and soil


### Venn Diagram of shared taxa between sediment and feather core ####

library(ggVennDiagram)

# Accessing the taxonomy tables of the transient feather and soil core phyloseqs
tax_60.F.core <- as.data.frame(tax_table(phylo.coreF.Me.rsqmerg.6))
tax_60.S.core <- as.data.frame(tax_table(phylo.coreS.6))

# Replace Genus with Family if the Genus is NA

#feathers
tax_60.F.core$Genus[is.na(tax_60.F.core$Genus)] <- tax_60.F.core$Family[is.na(tax_60.F.core$Genus)]

# soil
tax_60.S.core$Genus[is.na(tax_60.S.core$Genus)] <- tax_60.S.core$Family[is.na(tax_60.S.core$Genus)]

# Adding a Genus.species column to the tax tables 

# feathers
tax_60.F.core$Genus.species <- paste(tax_60.F.core$Genus, tax_60.F.core$Species)

# soil
tax_60.S.core$Genus.species <- paste(tax_60.S.core$Genus, tax_60.S.core$Species)

# Finding the shared SVs between the transient soil and the transient feathers 
library(ggVennDiagram)

# Venn of the Genus.species combos shared
Venn.Core <- list(Feather = tax_60.F.core$Genus.species, Sediment = tax_60.S.core$Genus.species)

ggVennDiagram(Venn.Core)

# 7 (12%) Genus.species unique to core feathers
# 49 (86%) Genus.species unique to core soil
# 1 (2%) Genus.species SHARED between core feathers and transient soil

# ASV level
Venn.Core.ASV <- list(Feather = rownames(tax_60.F.core), Sediment = rownames(tax_60.S.core))

ggVennDiagram(Venn.Core.ASV)

# 11 (11%) SVs unique to feathers
# 88 (88%) SVs unique to soil
# 1 (1%) SVs shared between feathers and soil

VennDiagram::get.venn.partitions(Venn.Core.ASV)

##### Venn Diagram - Entire feather community versus Entire soil community #####

# Unrarefied phyloseqs

# Accessing the taxonomy tables of the entire feather and entire soil phyloseqs
tax_ME.F <- as.data.frame(tax_table(ME_34.F.rsqmerg.1000))
tax_ME.S <- as.data.frame(tax_table(ME_34_S))

# Replace Genus with Family if the Genus is NA

#feathers
tax_ME.F$Genus[is.na(tax_ME.F$Genus)] <- tax_ME.F$Family[is.na(tax_ME.F$Genus)]

# soil
tax_ME.S$Genus[is.na(tax_ME.S$Genus)] <- tax_ME.S$Family[is.na(tax_ME.S$Genus)]

# Adding a Genus.species column to the tax tables 

# feathers
tax_ME.F$Genus.species <- paste(tax_ME.F$Genus, tax_ME.F$Species)

# soil
tax_ME.S$Genus.species <- paste(tax_ME.S$Genus, tax_ME.S$Species)

# Finding the shared SVs between the transient feathers and the core soil
library(ggVennDiagram)

# Genus.species Venn results
VennAll <- list(Feather = tax_ME.F$Genus.species, Sediment = tax_ME.S$Genus.species)

ggVennDiagram(VennAll)

# 273 (27%) Genus.species unique to feathers
# 369 (36%) Genus.species shared
# 372 (37%) Genus.species unique to sediment

# ASV level
VennAll <- list(Feather = rownames(tax_ME.F), Sediment = rownames(tax_ME.S))

ggVennDiagram(VennAll)

# 1910 (23%) SVs unique to feathers
# 5662 (70%) SVs unique to soil
# 564 (7%) SVs shared between feathers and soil

##### Transient feather and soil versus entire feather and soil ####

# ASV level
VennAll <- list(FeatherWholeCommunity = rownames(tax_ME.F), SedimentWholeCommunity = rownames(tax_ME.S), FeatherTransient = rownames(tax_20.F.trans), SedimentTransient = rownames(tax_20.S.trans))

ggVennDiagram(VennAll)


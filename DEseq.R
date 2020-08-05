#######
#' Using DESeq to determine taxa responing to the different fertilization types 
#' Since DESeq only does pairwise comparison I have done it in 3 parts to cover 
#' all possible combinations (no Vs synthetic, no Vs organic, organic Vs 
#' synthetic)

#' First create phyloseq object:
# make phyloseq otu table and taxonomy
OTU = otu_table(as.matrix(otu_filtered), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_filtered))
# add map
map_16S$sampling_week<-as.factor(map_16S$sampling_week)
rownames(map_16S) <- map_16S$sequence_name
# make phyloseq map
phyloseq_map <- sample_data(map_16S)
# make phyloseq object
physeq <- merge_phyloseq(OTU,TAX,phyloseq_map)
physeq

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

#Nitrogen free vs standard fertilization miscanthus 16
physeq_mis = subset_samples(physeq, plant != "switchgrass")
diagdds = phyloseq_to_deseq2(physeq_mis, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('treatment', 'standard fertilization', 'nitrogen free'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_mis)[rownames(sigtab), ], "matrix"))


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")

write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgen_mis.csv")

#Nitrogen free vs standard fertilization switchgrass 2016
physeq_16 = subset_samples(physeq, year != "2017")
physeq_swg16 = subset_samples(physeq_16, plant != "miscanthus")
diagdds = phyloseq_to_deseq2(physeq_swg16, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('treatment', 'standard fertilization', 'nitrogen free'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_swg16)[rownames(sigtab), ], "matrix"))


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")


write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgen_swg16.csv")

#Nitrogen free vs standard fertilization switchgrass 2017
physeq_swg17 = subset_samples(physeq, year != "2016")
diagdds = phyloseq_to_deseq2(physeq_swg17, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('treatment', 'standard fertilization', 'nitrogen free'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_swg17)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")


write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgen_swg17.csv")


#Nitrogen free vs standard fertilization switchgrass 2016 2017
physeq_swg = subset_samples(physeq, plant != "miscanthus")
diagdds = phyloseq_to_deseq2(physeq_swg, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('treatment', 'standard fertilization', 'nitrogen free'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_swg)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")


write.csv(sigtabfami,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabfami_swg1617.csv")


#Nitrogen free vs standard fertilization
diagdds = phyloseq_to_deseq2(physeq, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('treatment', 'standard fertilization', 'nitrogen free'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")


write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgenfert.csv")

#Switchgrass vs miscanthus
diagdds = phyloseq_to_deseq2(physeq, ~ plant)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('plant', 'switchgrass', 'miscanthus'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")



write.csv(sigtab,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgenplant.csv")


#Switchgrass vs miscanthus under nitrogen free
physeq_nf = subset_samples(physeq, treatment != "standard fertilization")
diagdds = phyloseq_to_deseq2(physeq_nf, ~ plant)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('plant', 'switchgrass', 'miscanthus'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")



write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgenplant_undernf.csv")

#Switchgrass vs miscanthus under standard fertilization
physeq_sf = subset_samples(physeq, treatment != "nitrogen free")
diagdds = phyloseq_to_deseq2(physeq_sf, ~ plant)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('plant', 'switchgrass', 'miscanthus'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(legend.position = "top")



write.csv(sigtabgen,"C:/Users/jaffyhu/Desktop/Soil/new/sigtabgenplant_undersf.csv")



#Before indicator species analysis
write.csv(otu_filtered,"C:/Users/jaffyhu/Desktop/Soil/new/otu_table.csv")
write.csv(map_16S,"C:/Users/jaffyhu/Desktop/Soil/new/map_16S.csv")


#' Indicator species analysis - need processing power! 
#' (done on high performance computational cluster at MSU)
#' Here is the job script for hpcc (copy past everythingfrom here down until the next #' and save it as <somename>.sb)

#!/bin/bash --login

#SBATCH --time=15:00:00         # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24      # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=520G              # memory required per node - amount of memory (in bytes)
#SBATCH --job-name indispec     # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=indispec.log

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 R/3.5.1-X11-20180604

cd /mnt/research/ShadeLab/WorkingSpace/Stopnisek/core_microbiome/  #use path to your own directory where data and R script is located
Rscript Indicspec_16S.R


#' Now you need the R script named above as Indicspec_16S.R (but you can use any other name - if 
#' you do change, change also the names in the code).
#' Everything below is now the R script (Indicspec_16S.R) that the above job will look for and 
#' execute and run. Same as before, copy and past into a text editor and save as Indicspec_16S.R or any_other_name.R.
 
install.packages("indicspecies")
library("indicspecies")
set.seed(1)
otu16S = read.csv("C:/Users/jaffyhu/Desktop/Soil/new/otu_table.csv",header= T, row.names = 1)  #read OTU table, here named indval_input.csv
map = read.csv("C:/Users/jaffyhu/Desktop/Soil/new/map_16S.csv", header=TRUE)  #read map file, here named indval_map.csv
map$treatment <- factor(map$treatment, levels = c("nitrogen free","standard fertilization"))
group_treatment = droplevels(map$treatment)  # I have chosen to use treatment regimes to search against the indicator species

indicspec16S_treatment = multipatt(as.data.frame(t(otu16S)), group_treatment, control = how(nperm = 999))
sum_indicspec16S_treatment = summary(indicspec16S_treatment, indvalcomp=TRUE)
saveRDS(sum_indicspec16S_treatment,"C:/Users/jaffyhu/Desktop/Soil/new/HPCC_Indicator_files/sum_indicspec16S_treatment.rds")
write.csv(sum_indicspec16S_treatment,"C:/Users/jaffyhu/Desktop/Soil/new/HPCC Indicator files/indicspec16S_treatment.csv")

#Switchgrass vs miscanthus
map$plant <- factor(map$plant, levels = c("switchgrass","miscanthus"))
group_plant = droplevels(map$plant)  # I have chosen to use plant regimes to search against the indicator species

indicspec16S_plant = multipatt(as.data.frame(t(otu16S)), group_plant, control = how(nperm = 999))
sum_indicspec16S_plant = summary(indicspec16S_plant, indvalcomp=TRUE)
saveRDS(sum_indicspec16S_plant,"C:/Users/jaffyhu/Desktop/Soil/new/HPCC_Indicator_files/sum_indicspec16S_plant.rds")
write.csv(sum_indicspec16S_plant,"C:/Users/jaffyhu/Desktop/Soil/new/HPCC Indicator files/indicspec16S_plant.csv")

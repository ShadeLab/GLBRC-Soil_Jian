######################################################
# Identifying taxa differently abundant between switchgrass and miscanthus (DESeq2)
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

diagdds = phyloseq_to_deseq2(physeq, ~ plant)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, contrast=c('plant', 'switchgrass', 'miscanthus'), cooksCutoff = FALSE)
alpha = 0.001

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
write.csv(sigtab,"C:/Users/jaffyhu/Desktop/Soil/new/Netural_model/sigtab_plant.csv")


######################################################
#Occupancy abundance relationship
swg_otu <- otu_rare[,map_16S$plant=="switchgrass"]
mis_otu <- otu_rare[,map_16S$plant=="miscanthus"]
#selecting variety unique OTUs
swg_otu <- swg_otu[rowSums(swg_otu)>0,]
mis_otu <- mis_otu[rowSums(mis_otu)>0,]

swg_uniq <- swg_otu[!(rownames(swg_otu) %in% rownames(mis_otu)),]
mis_uniq <- mis_otu[!(rownames(mis_otu) %in% rownames(swg_otu)),]

#remove the bulk samples
otu_rare <- otu_rare[rowSums(otu_rare)>0,]
otu_PA <- 1*((otu_rare>0)==1)
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Occ <- rowSums(otu_PA)/ncol(otu_PA)

com_abund <- rowSums(otu_rare)/ncol(otu_rare)
abund <- decostand(otu_rare, method='total', MARGIN=2)
com_abund <- rowSums(abund)/ncol(abund)

df_occ <- data.frame(otu=names(Occ), occ=Occ) 
df_abun <- data.frame(otu=names(com_abund), abun=log10(com_abund))
occ_abun <- left_join(df_occ, df_abun)

occ_abun$plant_found <- 'shared'
occ_abun$plant_found[occ_abun$otu %in% rownames(swg_uniq)] <- 'switchgrass'
occ_abun$plant_found[occ_abun$otu %in% rownames(mis_uniq)] <- 'miscanthus'

lastValue <- function(x) tail(x[!is.na(x)], 1)
last_taxons<- apply(taxonomy_filtered, 1, lastValue)
taxonomy_filtered$last_taxon <- last_taxons
taxonomy_filtered$final_names <- paste(taxonomy_filtered$otu_id, taxonomy_filtered$last_taxon, sep='-')
taxonomy_filtered$otu <- taxonomy_filtered$OTU

size_occ<- data.frame(otu=as.factor(rownames(otu_PA)),otu_PA) %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[,c("treatment", "plant", "sampling_week", "sequence_name", "year", "time_numeric", "sampling_Rdate")], by='sequence_name') %>%
  left_join(taxonomy_filtered, by='otu') %>%
  group_by(otu, sampling_Rdate, Family, Genus, final_names) %>%
  mutate(sum_abun=sum(abun),
         n_rep=length(abun),
         rel_occ=sum_abun/n_rep,
         is_present=1*((sum_abun>0)==1))

#calculating presence across time
size_occ %>%
  group_by(otu, sampling_Rdate) %>%
  dplyr::summarize(n=sum(is_present),
                   presence= 1*((n>0)==1)) %>%
  group_by(otu) %>%
  dplyr::summarize(
    total_presence=sum(presence)
  ) -> tmp_occ

combined_occ_data <- left_join(tmp_occ, occ_abun)

combined_occ_data%>%
  #filter(total_presence==1) %>%
  group_by(total_presence,plant_found) %>%
  summarise(counts=length(otu)) %>%
  group_by(total_presence) %>%
  mutate(tot_counts=sum(counts),
         rel_contribution=counts/tot_counts)

##*********************************
#Sloan neutral model 
spp=t(otu_rare)
taxon=as.vector(size_occ$final_names)

#Models for the whole community
obs.np=sncm.fit(spp, taxon=F, stats=FALSE, pool=NULL)
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
sta.np.16S <- sta.np

above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

Predict <- obs.np
Predict$otu <- rownames(Predict)
UpPredictOTU <- unique(Predict$otu[ap==TRUE])
DownPedictedOTU <- unique(Predict$otu[bp==TRUE])
UpPredictOTU
DownPedictedOTU
combined_occ_data$Predict <- 'Neturel'
combined_occ_data$Predict[combined_occ_data$otu %in% UpPredictOTU]<- 'Up'
combined_occ_data$Predict[combined_occ_data$otu %in% DownPedictedOTU]<- 'Down'
########
# Fig Netural

ggplot(data=combined_occ_data, aes(x=abun, y=occ)) +
  theme_bw()+
  geom_point(pch=21,  size=3, aes(fill=plant_found)) +
  geom_line(color='black', data=obs.np, size=2, aes(y=obs.np$freq.pred, x=log10(obs.np$p))) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)))+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)))+
  scale_fill_manual(aes(breaks = plant_found), values=c('darkorange','white', 'black')) +
  xlim(-6,-1.5)+
  labs(x=paste("log10(mean abundance)\n (n=",sta.np.16S$Richness," OTUs)", sep=''), y=paste("Occupancy (n=",sta.np.16S$Samples," samples)",sep=''), fill='plant') +
  theme(legend.position="top",
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = guide_legend(override.aes = list(alpha=1)),
         fill = guide_legend(title=NULL))


core_plant <- combined_occ_data%>%
  filter(occ==1)

core_plant_up <- combined_occ_data%>%
           filter(occ==1)%>%
           filter(Predict=='Up')

write.csv(core_plant,"C:/Users/jaffyhu/Desktop/Soil/new/Core_plant.csv")
write.csv(core_plant_up,"C:/Users/jaffyhu/Desktop/Soil/new/Core_plant_up.csv")

combined_occ_data <- left_join(combined_occ_data,taxonomy_filtered, by='otu')
write.csv(combined_occ_data,"C:/Users/jaffyhu/Desktop/Soil/new/combined_occ_plant.csv")



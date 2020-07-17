# Total R Code for GLBRC soil analyses
### Start Common
### DIRECTIONS, SET WORKING DIRECTORY TO "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil>"
### Before working on Analysis, run the entire Common block of Code (all the way until "### End Common")
setwd("C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil")
install.packages("dplyr")
install.packages("tidyr")
install.packages("reshape2")
install.packages("RSQLite")
install.packages("stringr")



library(dplyr)
library(tidyr)
library(reshape2)
library(RSQLite)
library(stringr)


# Preping the environmental metadata
glbrc <- dbConnect(RSQLite::SQLite(), "InputFiles/GLBRC_bioenergy_db.db" )

#Content of the DB
dbListTables(glbrc)

#Content of selected tables in the DB
dbListFields(glbrc, 'plant')
dbListFields(glbrc, 'sequencing')
dbListFields(glbrc, 'soil')
dbListFields(glbrc, 'nucleic_acids')
dbListFields(glbrc, 'plot')
dbListFields(glbrc, 'sampling')

#get tables into R
glbrc_NA <- dbGetQuery(glbrc, "select * from nucleic_acids") 
glbrc_soil <- dbGetQuery(glbrc, 'select * from soil')
glbrc_plant <- dbGetQuery(glbrc, 'select * from plant')
glbrc_plot <- dbGetQuery(glbrc, 'select * from plot')
glbrc_sampling <- dbGetQuery(glbrc, 'select * from sampling')
glbrc_sequncing <- dbGetQuery(glbrc, 'select * from sequencing')
glbrc_sequncing <- glbrc_sequncing %>% 
  mutate(nucleic_acid_name = str_trim(glbrc_sequncing$nucleic_acid_name, side = "both"))

#joining tables to create complete map file
metadata <- full_join(glbrc_sampling, glbrc_plot, by='plotID')
metadata <- full_join(metadata, glbrc_soil, by='sampleID')
metadata <- full_join(metadata, glbrc_plant, by='sampleID')
metadata <- full_join(metadata, glbrc_NA, by='sampleID')
metadata <- full_join(metadata, glbrc_sequncing, by='nucleic_acid_name')

map_full <- metadata

# #creating numeric time column
map_full$sampling_date <- paste0(map_full$month,'-', map_full$day,'-',map_full$year)
map_full$sampling_date <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time <- as.POSIXct(map_full$sampling_date, format='%m-%d-%Y')
time_numeric <- as.numeric(time)
map_time <- cbind(map_full, time_numeric)

#adding weather data
weather <- read.csv("InputFiles/kbs_weather_09212017.csv", encoding = 'UTF-8')
dim(weather)
head(weather)
weather$sampling_date <- as.POSIXct(weather$date, format='%d.%m.%y')

#subsetting weather file for sample dates
sub_weather <- weather[weather$sampling_date %in% map_time$sampling_date,] 

#merging dataframes - map file and weather
map_complete <- full_join(map_time, weather)

#subset metadata for 2106 and 2017 soil itag
all_16S_soil1617 <- subset(map_complete, exclude_from_analysis == "N" & sequencing_type == 'Illumina 16S iTag' & source == 'soil')
dim(all_16S_soil1617)

# write out file
write.table(all_16S_soil1617, file = "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil/all_16S_soil1617.txt", sep = "\t", row.names=FALSE, quote=FALSE)

head(all_16S_soil1617)  

#finding duplicates
all_16S_soil1617$help_name = as.character(lapply(strsplit(as.character(all_16S_soil1617$nucleic_acid_name), split="D"), "[", 1))
unique(all_16S_soil1617$help_name)
n_occur <- data.frame(table(all_16S_soil1617$help_name))
n_occur[n_occur$Freq > 1,]
duplicate_df <- all_16S_soil1617[all_16S_soil1617$help_name %in% n_occur$Var1[n_occur$Freq > 1],]
list_dupli <- duplicate_df$sequence_name #list of duplicate samples

D1 <- duplicate_df[grep('D1', duplicate_df$sequence_name),]
D1$removing <- 'remove'
dim(D1)

map_full <- full_join(all_16S_soil1617, D1) 
map_full$removing

# write out file
write.table(map_full, file = "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil/map_full.txt", sep = "\t", row.names=FALSE, quote=FALSE)


# Read in OTU table
otu <- read.table("otu.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)

# Get the name of the samples we have data for
samples <- colnames(otu)

# put taxonomy into its own variable
taxonomy <- read.csv('InputFiles/taxonomy_combined_merged_trimmed_otus.csv', header = T, row.names = 1, na.strings= c("NA", " ", ""))

# subset the map to include only those samples we have sequence data
map_small <- map_full[map_full$sequence_name %in% samples,]

# Remove rows that have N/A for the sequence name
map_small<- map_small[complete.cases(map_small$sequence_name),]


# Subset the samples to those we want to analyse (IE remove the duplicates)
samples <- samples[samples %in% map_small$sequence_name]
#Subset the OTU table to only the samples we want to analyze(IE remove the duplicates)
otu_sub <- otu[,colnames(otu) %in% samples]
# Order the samples
otu_sub <- otu_sub[,order(colnames(otu_sub))]
# Order the samples of the map the same way
map_small <- map_small[order(map_small$sequence_name),]
# Check to make sure they all match with each other
colnames(otu_sub) == map_small$sequence_name

# Map file that has duplicates of samples removed
map_16S <- map_small
# OTU table that removes duplicates of single samples
otu <- otu_sub
otu_CM <- otu
#removing Eukaryota from the OTU table
tax_short <- taxonomy[!grepl("Mitochondria", taxonomy$Family),]
tax_short <- tax_short[!grepl("Chloroplast", tax_short$Class),]

otu <- otu[rownames(otu) %in% rownames(tax_short),]

taxonomy_full <- taxonomy
taxonomy <- taxonomy[rowSums(otu)>0,]
otu <- otu[rowSums(otu)>0,]
otu_soil <- otu[,map_16S$source=="soil"]

tax_filtered <- tax_short%>%
  mutate(otu = rownames(tax_short)) %>%
  filter(!is.na(Phylum))

silva_bact_only <- read.csv('InputFiles/silva_bacteria_only_glbrc.csv', header=T)


keep_otus <- silva_bact_only %>%
  filter(lca_tax_slv != 'Unclassified;') %>%
  separate(lca_tax_slv, into=c("Kingdom", "Phylum", "Class", 
                               "Order", "Family", "Genus", "Species"), sep=";", remove=F) %>%
  mutate(taxonomy = lca_tax_slv) %>%
  filter(Kingdom!= 'Eukaryota') %>%
  select(-lca_tax_slv) %>%
  bind_rows(tax_filtered) %>%
  select(-taxonomy)


otu_filtered <- otu[rownames(otu) %in% keep_otus$otu,]
write.table(otu_filtered, file = "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil/otu_filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(tax_filtered, file = "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil/tax_filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)

taxonomy_filtered <- taxonomy_full[rownames(taxonomy_full) %in% rownames(otu_filtered),]


#rarecurve of soil samples

library(vegan)
set.seed(13)
colSums(otu_filtered)
min(colSums(otu_filtered)) #19757
#otu_rare <- t(rrarefy(t(otu), min(colSums(otu))))
otu_rare<- t(rrarefy(t(otu_filtered), min(colSums(otu_filtered))))
otu_rare <- otu_rare[,colSums(otu_rare)>19756]
map_16S <- map_16S[map_16S$sequence_name%in%colnames(otu_rare),]


curve_colors <- rep("darkgreen", ncol(otu))
curve_colors[map_16S$plant=="switchgrass"] <- "burlywood"
curve_colors[map_16S$plant=="miscanthus"] <- "burlywood4"

rarecurve(t(otu), step=1000, sample=max(colSums(otu)), label=FALSE, col = curve_colors)


# ### Multiplot code taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### End Common

### Start alpha diversity Analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(Heatplus)
library("limma")
library(ggrepel)
library(codyn)
library(gridExtra)
library(grid)
library(egg)

#install Heatplus
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Heatplus")

#install limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

#########################
#### Alpha Diversity ####
#########################

library(ggplot2)
library(scales)

otu_rare.PA <- 1*(otu_rare>0)

s <- specnumber(otu_rare,MARGIN=2) # richness
h <- diversity(t(otu_rare), index = "shannon") # Shannon index
pielou=h/log(s) # Pielou's evenness

#Adding sampling week property
map_16S$sampling_week <- 0
map_16S$sampling_week[map_16S$sampling_date == '2017-04-24'] <- 1 
map_16S$sampling_week[map_16S$sampling_date == '2016-04-18'] <- 1

map_16S$sampling_week[map_16S$sampling_date == '2017-05-15'] <- 2 
map_16S$sampling_week[map_16S$sampling_date == '2016-05-09'] <- 2

map_16S$sampling_week[map_16S$sampling_date == '2016-05-31'] <- 3
map_16S$sampling_week[map_16S$sampling_date == '2017-06-05'] <- 3 

map_16S$sampling_week[map_16S$sampling_date == '2016-06-20'] <- 4
map_16S$sampling_week[map_16S$sampling_date == '2017-06-26'] <- 4

map_16S$sampling_week[map_16S$sampling_date == '2016-07-12'] <- 5 
map_16S$sampling_week[map_16S$sampling_date == '2017-07-17'] <- 5

map_16S$sampling_week[map_16S$sampling_date == '2016-08-01'] <- 6 
map_16S$sampling_week[map_16S$sampling_date == '2017-08-07'] <- 6 

map_16S$sampling_week[map_16S$sampling_date == '2016-08-22'] <- 7 
map_16S$sampling_week[map_16S$sampling_date == '2017-08-28'] <- 7 

map_16S$sampling_week[map_16S$sampling_date == '2016-09-12'] <- 8 
map_16S$sampling_week[map_16S$sampling_date == '2017-09-18'] <- 8 

map_16S$sampling_week[map_16S$sampling_date == '2016-10-03'] <- 9

map_16S$sampling_week[map_16S$sampling_date == '2016-11-07'] <- 10

map_16S$sampling_Rdate <- as.Date(map_16S$sampling_date)


### Settin up Contextual Data Maps

map_df <- data.frame(map_16S)
map.div <- map_df
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 
map.div.mis <- map.div[map.div$plant=="miscanthus",]
map.div.swg <- map.div[map.div$plant=="switchgrass",]
map.div.swg16 <- map.div.swg[map.div.swg$year=="2016",]
map.div.swg17 <- map.div.swg[map.div.swg$year=="2017",]

# Melt the map and alpha diversity variables
library(MASS)
library(reshape2)
library(reshape)

map.alpha.mis16 <- melt(map.div.mis, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))
map.alpha.swg16 <- melt(map.div.swg16, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))
map.alpha.swg17 <- melt(map.div.swg17, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))


misc.map.2016 <- map_16S[map_16S$plant=="miscanthus",]
switch.map <- map_16S[map_16S$plant=="switchgrass",]
switch.map.2016 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2016,]
switch.map.2017 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2017,]



cor.test(map.div.swg16$time_numeric, map.div.swg16$Richness) #t = -1.8558, df = 64, p-value = 0.06809
cor.test(map.div.swg16$time_numeric, map.div.swg16$Shannon) #t = -2.4316, df = 64, p-value = 0.01784
cor.test(map.div.swg16$time_numeric, map.div.swg16$Pielou) #t = -2.7737, df = 64, p-value = 0.007253

cor.test(map.div.swg17$time_numeric, map.div.swg17$Richness) #t = -5.192, df = 60, p-value = 2.608e-06
cor.test(map.div.swg17$time_numeric, map.div.swg17$Shannon) #t = -3.8252, df = 60, p-value = 0.0003136
cor.test(map.div.swg17$time_numeric, map.div.swg17$Pielou) #t = -2.3996, df = 60, p-value = 0.01954

cor.test(map.div.mis$time_numeric, map.div.mis$Richness) #t = -0.12533, df = 72, p-value = 0.9006
cor.test(map.div.mis$time_numeric, map.div.mis$Shannon) #t = -0.52522, df = 72, p-value = 0.601
cor.test(map.div.mis$time_numeric, map.div.mis$Pielou) #t = -0.81141, df = 72, p-value = 0.4198



otu_rare_2016 <- otu_rare[,map_16S$year=="2016"]
switch_rare_otu.2017 <- otu_rare[,map_16S$year=="2017"]

misc_rare_otu.2016 <- otu_rare[,map_16S $ plant=="miscanthus"]
switch_rare_otu.2016 <- otu_rare[,map_16S$plant=="switchgrass" & map_16S$year=="2016"]

### soil misc unique times
misc_unique_times <- unique(misc.map.2016$time_numeric)[order(unique(misc.map.2016$time_numeric))]

### soil switchgrass unique times
switch_unique_times <- unique(switch.map$time_numeric)[order(unique(switch.map$time_numeric))]
switch_unique_times.2016 <- unique(switch.map.2016$time_numeric)[order(unique(switch.map.2016$time_numeric))]
switch_unique_times.2017 <-  unique(switch.map.2017$time_numeric)[order(unique(switch.map.2017$time_numeric))]


# Look at species accumulation
misc_accumulation <- rep(1, length(misc_unique_times))
z <- NULL
for( i in 1:length(unique(misc.map.2016$time_numeric))){
  x <- misc_rare_otu.2016[,misc.map.2016$time_numeric==misc_unique_times[i]]
  y <- matrix(x, ncol=sum(misc.map.2016$time_numeric==misc_unique_times[i]))
  row.names(y) <- row.names(misc_rare_otu.2016)
  z <- c(z, row.names(misc_rare_otu.2016[rowSums(y)!=0,]))
  misc_accumulation[i] <- length(unique(z))
}

switch_accumulation.2016 <- rep(1, length(switch_unique_times.2016))
z <- NULL
for( i in 1:length(unique(switch.map.2016$time_numeric))){
  x <- switch_rare_otu.2016[,switch.map.2016$time_numeric==switch_unique_times.2016[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2016[i] <- length(unique(z))
}

switch_accumulation.2017 <- rep(1, length(switch_unique_times.2017))
z <- NULL
for( i in 1:length(unique(switch.map.2017$time_numeric))){
  x <- switch_rare_otu.2017[,switch.map.2017$time_numeric==switch_unique_times.2017[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2017[i] <- length(unique(z))
}


Spec_accum.2016 <- NULL
Spec_accum.2016$Species_Accumulation <- c(misc_accumulation, switch_accumulation.2016)
Spec_accum.2016$Date <- c(unique(misc.map.2016$sampling_Rdate)[order(unique(misc.map.2016$sampling_Rdate))], unique(switch.map.2016$sampling_Rdate)[order(unique(switch.map.2016$sampling_Rdate))])
Spec_accum.2016 <- as.data.frame(Spec_accum.2016)
Spec_accum.2016$Plant <- c(rep("Miscanthus", 10), rep("Switchgrass", 9))
Spec_accum.2016

Spec_accum.2017 <- NULL
Spec_accum.2017$Species_Accumulation <- switch_accumulation.2017
Spec_accum.2017$Date <- unique(switch.map.2017$sampling_Rdate)[order(unique(switch.map.2017$sampling_Rdate))]
Spec_accum.2017 <- as.data.frame(Spec_accum.2017)
Spec_accum.2017$Plant <- rep("Switchgrass", nrow(Spec_accum.2017)) 
Spec_accum.2017

NewDates_Switch2017 <- c(as.Date("2016-04-18"), as.Date("2016-05-09"), as.Date("2016-05-31"), as.Date("2016-06-20"),as.Date("2016-07-12"), as.Date("2016-08-01"), as.Date("2016-08-22"), as.Date("2016-09-12"))

Spec_accum_total <- NULL
Spec_accum_total$Species_Accumulation <- c(misc_accumulation, switch_accumulation.2016, switch_accumulation.2017)
Spec_accum_total$Date <- c(unique(misc.map.2016$sampling_Rdate)[order(unique(misc.map.2016$sampling_Rdate))], unique(switch.map.2016$sampling_Rdate)[order(unique(switch.map.2016$sampling_Rdate))], NewDates_Switch2017)
Spec_accum_total <- as.data.frame(Spec_accum_total)
Spec_accum_total$Plant <- c(rep("Miscanthus", length(misc_accumulation)), rep("Switchgrass", length(switch_accumulation.2016)), rep("Switchgrass", length(switch_accumulation.2017)))



# Code of soil richness accumulation through the season
switch.Richness.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),switch.Richness.2016)
summary(fit_model) #p=0.00389 **
TukeyHSD(fit_model)

misc.Richness <- map.alpha.mis16[map.alpha.mis16$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),misc.Richness)
summary(fit_model) #p=0.00462 **
TukeyHSD(fit_model)

switch.Richness.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Richness",]
fit_model <- aov(value~factor(sampling_Rdate),switch.Richness.2017)
summary(fit_model) #p=0.000122 ***
TukeyHSD(fit_model)

Richness.df <- rbind(switch.Richness.2016, switch.Richness.2017, misc.Richness)
Richness.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Richness.2016)), rep("Switchgrass 2017", nrow(switch.Richness.2017)), rep("Miscanthus 2016", nrow(misc.Richness))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Richness.df$FakeDate <- as.Date(gsub(x=Richness.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))



### Spec Accumulation Curves
Spec_accum_total$Year <- c(rep("2016",length(misc_accumulation)), rep("2016",length(switch_accumulation.2016)), rep("2017",length(switch_accumulation.2017)))
Spec_accum_total
Spec_accum_total$PlantYear <- paste(Spec_accum_total$Plant, Spec_accum_total$Year)

Fig1A <- ggplot(Spec_accum_total, aes(x=Date, y=Species_Accumulation)) + 
  geom_point(aes(color=Plant, shape=PlantYear), size=2) +
  geom_line(aes(color=Plant, linetype=Year))+
  scale_shape_manual(values=c(15,16,1))+
  scale_linetype_manual(values = c(1,2))+
  theme_bw()+
  labs(y="Total Observed Taxa") + 
  scale_color_manual(values=c("darkgreen","darkolivegreen3"))  + 
  scale_x_date(date_breaks = "1 month", labels=date_format("%b"))+
  guides(shape=FALSE)

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig1A_SpeciesAccumulation.eps", Fig1A, device = "eps", width=6, height=4, units = "in")



#accumulative observed taxa under different treatment, plant type and year.


### Settin up Contextual Data Maps

map_df <- data.frame(map_16S)
map.div <- map_df
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou 
map.div.soil <- map.div[map.div$source=="soil",]
map.div.mis <- map.div.soil[map.div.soil$plant=="miscanthus",]
map.div.misnf <- map.div.mis[map.div.mis$treatment=="nitrogen free",]
map.div.missf <- map.div.mis[map.div.mis$treatment=="standard fertilization",]
map.div.swg <- map.div.soil[map.div.soil$plant=="switchgrass",]
map.div.swg16 <- map.div.swg[map.div.swg$year=="2016",]
map.div.swg16nf <- map.div.swg16[map.div.swg16$treatment=="nitrogen free",]
map.div.swg16sf <- map.div.swg16[map.div.swg16$treatment=="standard fertilization",]
map.div.swg17 <- map.div.swg[map.div.swg$year=="2017",]
map.div.swg17nf <- map.div.swg17[map.div.swg17$treatment=="nitrogen free",]
map.div.swg17sf <- map.div.swg17[map.div.swg17$treatment=="standard fertilization",]

# Melt the map and alpha diversity variables
library(MASS)
library(reshape2)
library(reshape)

map.alpha.mis16 <- melt(map.div.mis, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))
map.alpha.swg16 <- melt(map.div.swg16, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))
map.alpha.swg17 <- melt(map.div.swg17, id.vars=c("sequence_name","treatment", "source", "plant", "time_numeric", "sampling_Rdate","year","sampling_week"), measure.vars=c("Richness", "Shannon", "Pielou"))



misc.map.2016 <- map_16S[map_16S$plant=="miscanthus",]
misc.map.2016nf <- misc.map.2016[misc.map.2016$treatment=="nitrogen free",]
misc.map.2016sf <- misc.map.2016[misc.map.2016$treatment=="standard fertilization",]

switch.map.2016 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2016,]
switch.map.2016nf <- switch.map.2016[switch.map.2016$treatment=="nitrogen free",]
switch.map.2016sf <- switch.map.2016[switch.map.2016$treatment=="standard fertilization",]

switch.map.2017 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2017,]
switch.map.2017nf <- switch.map.2017[switch.map.2017$treatment=="nitrogen free",]
switch.map.2017sf <- switch.map.2017[switch.map.2017$treatment=="standard fertilization",]


cor.test(map.div.swg16$time_numeric, map.div.swg16$Richness) #t = -1.8558, df = 64, p-value = 0.06809
cor.test(map.div.swg16$time_numeric, map.div.swg16$Shannon) #t = -2.4316, df = 64, p-value = 0.01784
cor.test(map.div.swg16$time_numeric, map.div.swg16$Pielou) #t = -2.7737, df = 64, p-value = 0.007253


cor.test(map.div.swg17$time_numeric, map.div.swg17$Richness) #t = -5.192, df = 60, p-value = 2.608e-06
cor.test(map.div.swg17$time_numeric, map.div.swg17$Shannon) #t = -3.8252, df = 60, p-value = 0.0003136
cor.test(map.div.swg17$time_numeric, map.div.swg17$Pielou) #t = -2.3996, df = 60, p-value = 0.01954


cor.test(map.div.mis$time_numeric, map.div.mis$Richness) #t = -0.12533, df = 72, p-value = 0.9006
cor.test(map.div.mis$time_numeric, map.div.mis$Shannon) #t = -0.52522, df = 72, p-value = 0.601
cor.test(map.div.mis$time_numeric, map.div.mis$Pielou) #t = -0.81141, df = 72, p-value = 0.4198



switch_rare_otu.2017 <- otu_rare[,map_16S$year==2017]
switch_rare_otu.2017nf <- switch_rare_otu.2017[,switch.map.2017$treatment=="nitrogen free"]
switch_rare_otu.2017sf <- switch_rare_otu.2017[,switch.map.2017$treatment=="standard fertilization"]

misc_rare_otu.2016 <- otu_rare[,map_16S$plant=="miscanthus"]
misc_rare_otu.2016nf <- misc_rare_otu.2016[,misc.map.2016$treatment=="nitrogen free"]
misc_rare_otu.2016sf <- misc_rare_otu.2016[,misc.map.2016$treatment=="standard fertilization"]

switch_rare_otu.2016 <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$year==2016]
switch_rare_otu.2016nf <- switch_rare_otu.2016[,switch.map.2016$treatment=="nitrogen free"]
switch_rare_otu.2016sf <- switch_rare_otu.2016[,switch.map.2016$treatment=="standard fertilization"]


### soil misc unique times
misc_unique_timesnf <- unique(misc.map.2016nf$time_numeric)[order(unique(misc.map.2016nf$time_numeric))]
misc_unique_timessf <- unique(misc.map.2016sf$time_numeric)[order(unique(misc.map.2016sf$time_numeric))]

### soil switchgrass unique times

switch_unique_times.2016nf <- unique(switch.map.2016nf$time_numeric)[order(unique(switch.map.2016nf$time_numeric))]
switch_unique_times.2016sf <- unique(switch.map.2016sf$time_numeric)[order(unique(switch.map.2016sf$time_numeric))]

switch_unique_times.2017nf <- unique(switch.map.2017nf$time_numeric)[order(unique(switch.map.2017nf$time_numeric))]
switch_unique_times.2017sf <- unique(switch.map.2017sf$time_numeric)[order(unique(switch.map.2017sf$time_numeric))]

# Look at species accumulation
nfmisc_accumulation <- rep(1, length(misc_unique_timesnf))
z <- NULL
for( i in 1:length(unique(misc.map.2016nf$time_numeric))){
  x <- misc_rare_otu.2016nf[,misc.map.2016nf$time_numeric==misc_unique_timesnf[i]]
  y <- matrix(x, ncol=sum(misc.map.2016nf$time_numeric==misc_unique_timesnf[i]))
  row.names(y) <- row.names(misc_rare_otu.2016nf)
  z <- c(z, row.names(misc_rare_otu.2016nf[rowSums(y)!=0,]))
  nfmisc_accumulation[i] <- length(unique(z))
}

sfmisc_accumulation <- rep(1, length(misc_unique_timessf))
z <- NULL
for( i in 1:length(unique(misc.map.2016sf$time_numeric))){
  x <- misc_rare_otu.2016sf[,misc.map.2016sf$time_numeric==misc_unique_timessf[i]]
  y <- matrix(x, ncol=sum(misc.map.2016sf$time_numeric==misc_unique_timessf[i]))
  row.names(y) <- row.names(misc_rare_otu.2016sf)
  z <- c(z, row.names(misc_rare_otu.2016sf[rowSums(y)!=0,]))
  sfmisc_accumulation[i] <- length(unique(z))
}

switch_accumulation.2016nf <- rep(1, length(switch_unique_times.2016nf))
z <- NULL
for( i in 1:length(unique(switch.map.2016nf$time_numeric))){
  x <- switch_rare_otu.2016nf[,switch.map.2016nf$time_numeric==switch_unique_times.2016nf[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2016nf[i] <- length(unique(z))
}

switch_accumulation.2016sf <- rep(1, length(switch_unique_times.2016sf))
z <- NULL
for( i in 1:length(unique(switch.map.2016sf$time_numeric))){
  x <- switch_rare_otu.2016sf[,switch.map.2016sf$time_numeric==switch_unique_times.2016sf[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2016sf[i] <- length(unique(z))
}

switch_accumulation.2017nf <- rep(1, length(switch_unique_times.2017nf))
z <- NULL
for( i in 1:length(unique(switch.map.2017nf$time_numeric))){
  x <- switch_rare_otu.2017nf[,switch.map.2017nf$time_numeric==switch_unique_times.2017nf[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2017nf[i] <- length(unique(z))
}

switch_accumulation.2017sf <- rep(1, length(switch_unique_times.2017sf))
z <- NULL
for( i in 1:length(unique(switch.map.2017sf$time_numeric))){
  x <- switch_rare_otu.2017sf[,switch.map.2017sf$time_numeric==switch_unique_times.2017sf[i]]
  z <- c(z, row.names(x[rowSums(x)!=0,]))
  switch_accumulation.2017sf[i] <- length(unique(z))
}


Spec_accum.2016 <- NULL
Spec_accum.2016$Species_Accumulation <- c(sfmisc_accumulation, nfmisc_accumulation, switch_accumulation.2016nf, switch_accumulation.2016sf)
Spec_accum.2016$Date <- c(unique(misc.map.2016nf$sampling_Rdate)[order(unique(misc.map.2016nf$sampling_Rdate))], unique(misc.map.2016sf$sampling_Rdate)[order(unique(misc.map.2016sf$sampling_Rdate))], unique(switch.map.2016nf$sampling_Rdate)[order(unique(switch.map.2016nf$sampling_Rdate))], unique(switch.map.2016sf$sampling_Rdate)[order(unique(switch.map.2016sf$sampling_Rdate))])
Spec_accum.2016 <- as.data.frame(Spec_accum.2016)
Spec_accum.2016$Plant <- c(rep("Miscanthus", 10), rep("Switchgrass", 9))
Spec_accum.2016

Spec_accum.2017 <- NULL
Spec_accum.2017$Species_Accumulation <- switch_accumulation.2017
Spec_accum.2017$Date <- unique(switch.map.2017$sampling_Rdate)[order(unique(switch.map.2017$sampling_Rdate))]
Spec_accum.2017 <- as.data.frame(Spec_accum.2017)
Spec_accum.2017$Plant <- rep("Switchgrass", nrow(Spec_accum.2017)) 


NewDates_Switch2017 <- c(as.Date("2016-04-18"), as.Date("2016-05-09"), as.Date("2016-05-31"), as.Date("2016-06-20"),as.Date("2016-07-12"), as.Date("2016-08-01"), as.Date("2016-08-22"), as.Date("2016-09-12"))

### Spec Accumulation Curves
Spec_accum_total <- NULL
Spec_accum_total$Species_Accumulation <- c(nfmisc_accumulation, sfmisc_accumulation, switch_accumulation.2016nf, switch_accumulation.2016sf, switch_accumulation.2017nf, switch_accumulation.2017sf)
Spec_accum_total$Date <- c(unique(misc.map.2016nf$sampling_Rdate)[order(unique(misc.map.2016nf$sampling_Rdate))], unique(misc.map.2016sf$sampling_Rdate)[order(unique(misc.map.2016sf$sampling_Rdate))], unique(switch.map.2016nf$sampling_Rdate)[order(unique(switch.map.2016nf$sampling_Rdate))], unique(switch.map.2016sf$sampling_Rdate)[order(unique(switch.map.2016sf$sampling_Rdate))], unique(switch.map.2017nf$sampling_Rdate)[order(unique(switch.map.2017nf$sampling_Rdate))], unique(switch.map.2017sf$sampling_Rdate)[order(unique(switch.map.2017sf$sampling_Rdate))])
Spec_accum_total <- as.data.frame(Spec_accum_total)
Spec_accum_total$FertilizationStatus <- c(rep("Nitrogen free", length(nfmisc_accumulation)), rep("Standard fertilization", length(sfmisc_accumulation)), rep("Nitrogen free", length(switch_accumulation.2016nf)), rep("Standard fertilization", length(switch_accumulation.2016sf)), rep("Nitrogen free", length(switch_accumulation.2017nf)), rep("Standard fertilization", length(switch_accumulation.2017sf)))
Spec_accum_total$Plant <- c(rep("Miscanthus", length(nfmisc_accumulation)), rep("Miscanthus", length(sfmisc_accumulation)), rep("Switchgrass", length(switch_accumulation.2016nf)), rep("Switchgrass", length(switch_accumulation.2016sf)), rep("Switchgrass", length(switch_accumulation.2017nf)), rep("Switchgrass", length(switch_accumulation.2017sf)))
Spec_accum_total$Year <- c(rep("2016",length(nfmisc_accumulation)), rep("2016",length(sfmisc_accumulation)), rep("2016",length(switch_accumulation.2016nf)), rep("2016",length(switch_accumulation.2016sf)), rep("2017",length(switch_accumulation.2017nf)), rep("2017",length(switch_accumulation.2017sf)))
Spec_accum_total$PlantFertilization <- c(rep("Miscanthus16NF", length(nfmisc_accumulation)), rep("Miscanthus16SF", length(sfmisc_accumulation)), rep("Switchgrass16NF", length(switch_accumulation.2016nf)), rep("Switchgrass16SF", length(switch_accumulation.2016sf)), rep("Switchgrass17NF", length(switch_accumulation.2017nf)), rep("Switchgrass17SF", length(switch_accumulation.2017sf)))
Spec_accum_total
Spec_accum_total$PlantYear <- paste(Spec_accum_total$Plant, Spec_accum_total$Year)

#Adding sampling week property
Spec_accum_total$sampling_week <- 0
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-04-24'] <- 1 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-04-18'] <- 1

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-05-15'] <- 2 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-05-09'] <- 2

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-05-31'] <- 3
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-06-05'] <- 3 

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-06-20'] <- 4
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-06-26'] <- 4

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-07-12'] <- 5 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-07-17'] <- 5

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-08-01'] <- 6 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-08-07'] <- 6 

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-08-22'] <- 7 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-08-28'] <- 7 

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-09-12'] <- 8 
Spec_accum_total$sampling_week[Spec_accum_total$Date == '2017-09-18'] <- 8 

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-10-03'] <- 9

Spec_accum_total$sampling_week[Spec_accum_total$Date == '2016-11-07'] <- 10



library(ggplot2)
library(scales)
Fig1B <- ggplot(Spec_accum_total, aes(x=sampling_week, y=Species_Accumulation)) + 
  geom_point(aes(color=PlantFertilization), size=2) +
  geom_line(aes(color=PlantFertilization, linetype=FertilizationStatus))+
  scale_shape_manual()+
  scale_linetype_manual(values = c(2,1))+
  theme_bw()+
  labs(y="Total Observed Taxa") +
  scale_color_manual(values=c("green", "green", "red", "red", "blue", "blue"))  +
  guides(shape=FALSE)

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig1B_SpeciesAccumulation.eps", Fig1B, device = "eps", width=6, height=4, units = "in")



library(ggpubr)


# Plot Pielou Evenness Through Season
switch.Pielou.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Pielou",]

switch.Pielou.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Pielou",]

misc.Pielou.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Pielou",]


Pielou.df <- rbind(switch.Pielou.2016, switch.Pielou.2017, misc.Pielou.2016)
Pielou.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Pielou.2016)), rep("Switchgrass 2017", nrow(switch.Pielou.2017)), rep("Miscanthus 2016", nrow(misc.Pielou.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Pielou.df$FakeDate <- as.Date(gsub(x=Pielou.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Pielou.df$treatment <- factor(Pielou.df$treatment, levels=c("nitrogen free", "standard fertilization"))


SupplementaryFigure2 <- ggplot(data = Pielou.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Pielou")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Pielou.eps", SupplementaryFigure2, device = "eps", width=6, height=4, units = "in")

# Plot Richness Through Season
switch.Richness.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Richness",]

switch.Richness.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Richness",]

misc.Richness.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Richness",]


Richness.df <- rbind(switch.Richness.2016, switch.Richness.2017, misc.Richness.2016)
Richness.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Richness.2016)), rep("Switchgrass 2017", nrow(switch.Richness.2017)), rep("Miscanthus 2016", nrow(misc.Richness.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Richness.df$FakeDate <- as.Date(gsub(x=Richness.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Richness.df$treatment <- factor(Richness.df$treatment, levels=c("nitrogen free", "standard fertilization"))


SupplementaryFigure2 <- ggplot(data = Richness.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Richness")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Richness.eps", SupplementaryFigure2, device = "eps", width=6, height=4, units = "in")

# Plot Shannon Through Season
switch.Shannon.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Shannon",]

switch.Shannon.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Shannon",]

misc.Shannon.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Shannon",]


Shannon.df <- rbind(switch.Shannon.2016, switch.Shannon.2017, misc.Shannon.2016)
Shannon.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Shannon.2016)), rep("Switchgrass 2017", nrow(switch.Shannon.2017)), rep("Miscanthus 2016", nrow(misc.Shannon.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Shannon.df$FakeDate <- as.Date(gsub(x=Shannon.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Shannon.df$treatment <- factor(Shannon.df$treatment, levels=c("nitrogen free", "standard fertilization"))


SupplementaryFigure2 <- ggplot(data = Shannon.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Shannon")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Shannon.eps", SupplementaryFigure2, device = "eps", width=6, height=4, units = "in")



library(ggpubr)

# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.misnf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness ##

### 1. Compare bacterial and archaeal richness among sampling_week using one-way ANOVA
Aov_richness_samplingweek <- lm(mis.map_aov$Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_richness_samplingweek
drop1(Aov_richness_samplingweek,~.,test="F") # type III SS and F Tests 0.345
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_richness_samplingweek)
RC_samplingweek_preds <- predict(Aov_richness_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.5649, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
library(car)
leveneTest(Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous
boxplot(mis.map_aov$Richness ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

library(moments)
skew_xts <-  skewness(RC_samplingweek_resids) # -0.308 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 2.69 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 9.0175, df = 9, p-value = 0.4357
ggboxplot(mis.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal richness among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)

#Standard fertilization miscanthus2016
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.missf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness ##

### 1. Compare bacterial and archaeal richness among sampling_week using one-way ANOVA
Aov_richness_samplingweek <- lm(mis.map_aov$Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_richness_samplingweek
drop1(Aov_richness_samplingweek,~.,test="F") # type III SS and F Tests 0.08423
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_richness_samplingweek)
RC_samplingweek_preds <- predict(Aov_richness_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.1465, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
library(car)
leveneTest(Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.03681
boxplot(mis.map_aov$Richness ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

library(moments)
skew_xts <-  skewness(RC_samplingweek_resids) # 0.349 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.245 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 14.487, df = 9, p-value = 0.106
ggboxplot(mis.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal richness among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)


#Switchgrass 2016 Nitrogen free
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16nf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness in switch2016##

### 1. Compare bacterial and archaeal richness among sampling_week using one-way ANOVA
Aov_richness_samplingweek <- lm(swg16.map_aov$Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_richness_samplingweek
drop1(Aov_richness_samplingweek,~.,test="F") # type III SS and F Tests 0.01843
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_richness_samplingweek)
RC_samplingweek_preds <- predict(Aov_richness_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.08507, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
library(car)
leveneTest(Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.5138
boxplot(swg16.map_aov$Richness ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

library(moments)
skew_xts <-  skewness(RC_samplingweek_resids) # -0.6058 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 19.114, df = 8, p-value = 0.01426
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal richness among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)

#Switchgrass 2016 Standard Fertilization
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16sf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness in switch2016##

### 1. Compare bacterial and archaeal richness among sampling_week using one-way ANOVA
Aov_richness_samplingweek <- lm(swg16.map_aov$Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_richness_samplingweek
drop1(Aov_richness_samplingweek,~.,test="F") # type III SS and F Tests 0.01843
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_richness_samplingweek)
RC_samplingweek_preds <- predict(Aov_richness_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.08507, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
library(car)
leveneTest(Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.5138
boxplot(swg16.map_aov$Richness ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

library(moments)
skew_xts <-  skewness(RC_samplingweek_resids) # -0.6058 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 19.114, df = 8, p-value = 0.01426
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal richness among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)



# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2017 
swg17.map_aov <- map.div.swg17sf
class(swg17.map_aov)
swg17.map_aov$sampling_week<-as.factor(swg17.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg17.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness in switch2017##

### 1. Compare bacterial and archaeal richness among sampling_week using one-way ANOVA
Aov_richness_samplingweek <- lm(swg17.map_aov$Richness ~ sampling_week, data=swg17.map_aov, na.action=na.exclude)
Aov_richness_samplingweek
drop1(Aov_richness_samplingweek,~.,test="F") # type III SS and F Tests 0.000149
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_richness_samplingweek)
RC_samplingweek_preds <- predict(Aov_richness_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.7774, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
library(car)
leveneTest(Richness ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous
boxplot(swg17.map_aov$Richness ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

library(moments)
skew_xts <-  skewness(RC_samplingweek_resids) # -0.061 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 23.436, df = 7, p-value = 0.001431
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal richness among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, swg17.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)





############################################################################################
######################### BACTERIAL COMMUNITIES COMPOSITION #####################
############################################################################################
BiocManager::install("phyloseq")
library(phyloseq)

# 1. BACTERIA COMPOSITION
# read bacterial taxonomy
taxonomy_filtered
rownames(taxonomy_filtered) <- rownames(otu_filtered)
# make phyloseq otu table and taxonomy
OTU = otu_table(as.matrix(otu_filtered), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_filtered))
# add map
map_16S$sampling_week<-as.factor(map_16S$sampling_week)
rownames(map_16S) <- map_16S$sequence_name
# make phyloseq map
phyloseq_map <- sample_data(map_16S)
# make phyloseq object
PHYL_16S <- merge_phyloseq(OTU,TAX,phyloseq_map)
PHYL_16S

PHYL_16S_phylum <- PHYL_16S %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)

PHYL_16S_phylum$PlantYearTreat <- paste(PHYL_16S_phylum$plant,PHYL_16S_phylum$Year, PHYL_16S_phylum$treatment, sep="")

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861","#ffffff", "#000000"
)

# Plot by sampling week
Fig3a <- ggplot(PHYL_16S_phylum, aes(x = sampling_week, y = Abundance, fill = Phylum)) + 
  facet_grid(PlantYearTreat~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition \n Bacterial Communities by Sampling Week")  
ggsave("C:/Users/jaffyhu/Desktop/Soil/Figure_3a.eps", Fig3a, device = "eps", width = 12, height = 7, units= "in", dpi = 600)

# Plot by treatment
Fig3b <- ggplot(PHYL_16S_phylum, aes(x = treatment, y = Abundance, fill = Phylum)) + 
  facet_grid(PlantYearTreat~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition \n Bacterial Communities by Sampling Week")  
ggsave("C:/Users/jaffyhu/Desktop/Soil/Figure_3b.eps", Fig3b, device = "eps", width = 12, height = 7, units= "in", dpi = 600)

################################
#### Beta Diversity Analyses####
################################
Date_to_Week <- data.frame(sampling_date=unique(map_16S$sampling_date)[order(unique(map_16S$sampling_date))], Week=c(1:10,1:8))

map_16S <- left_join(map_16S, Date_to_Week, by="sampling_date")
map.2016 <- map_16S[map_16S$Year==2016,]
swg.map <- map_16S[map_16S$plant=="switchgrass",]
swg.map.2016 <- map.2016[map.2016$plant=="switchgrass",]
mis.map.2016 <- map.2016[map.2016$plant=="miscanthus",]
swg.map.2017 <- map_16S[map_16S$Year==2017,]


otura.swg <- otu_rare[,map_16S$plant=="switchgrass"]
otura.2016 <- otu_rare[,map_16S$Year==2016]
otura.swg.2016 <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$Year==2016]
otura.mis.2016 <- otu_rare[,map_16S$plant=="miscanthus"&map_16S$Year==2016]
otura.swg.2017 <- otu_rare[,map_16S$Year==2017]

dist.otu <- vegdist(t(otu_rare), method="bray")
dist.otu.swg <- vegdist(t(otura.swg), method="bray")
dist.otu.2016 <- vegdist(t(otura.2016), method="bray")
dist.otu.2016swg <- vegdist(t(otura.swg.2016), method="bray")
dist.otu.2016mis <- vegdist(t(otura.mis.2016), method="bray")
dist.otu.2017swg <- vegdist(t(otura.swg.2017), method="bray")

pcoa.16S <- cmdscale(dist.otu, eig=TRUE) 
pcoa.swg <- cmdscale(dist.otu.swg, eig=TRUE)
pcoa.2016 <- cmdscale(dist.otu.2016, eig=TRUE)
pcoa.swg.2016 <- cmdscale(dist.otu.2016swg, eig=TRUE)
pcoa.mis.2016 <- cmdscale(dist.otu.2016mis, eig=TRUE)
pcoa.swg.2017 <- cmdscale(dist.otu.2017swg, eig=TRUE)

adonis(dist.otu~map_16S$time_numeric)  #0.001 ***
adonis(dist.otu.2016swg~swg.map.2016$time_numeric) #0.001 ***
adonis(dist.otu.2016mis~mis.map.2016$time_numeric) #0.001 ***
adonis(dist.otu.2017swg~swg.map.2017$time_numeric) #0.001 ***

adonis(dist.otu.2016~map.2016$treatment) #0.001 ***
adonis(dist.otu.2016swg~swg.map.2016$treatment) #0.001 ***
adonis(dist.otu.2016mis~mis.map.2016$treatment) #0.001 ***
adonis(dist.otu.2017swg~swg.map.2017$treatment) #0.001 ***

#Crop
adonis(dist.otu.2016~map.2016$plant) #0.001 ***

#Year
adonis(dist.otu.swg~swg.map$year) #0.001 ***

#######################################
### Setting up collapsed soil pcoa ###
#######################################

### Start by getting the dates for soil smaples
soil.dates <- unique(map_16S$sampling_date)[order(unique(map_16S$sampling_date))]
swg2016.dates <- unique(swg.map.2016$sampling_date)[order(unique(swg.map.2016$sampling_date))]
swg2017.dates <- unique(swg.map.2017$sampling_date)[order(unique(swg.map.2017$sampling_date))]
mis2016.dates <- unique(mis.map.2016$sampling_date)[order(unique(mis.map.2016$sampling_date))]

### Set up the points from the soil PCoA
soil.points <- pcoa.16S$points
swg2016.points <- pcoa.swg.2016$points
swg2017.points <- pcoa.swg.2017$points
mis2016.points <- pcoa.mis.2016$points

### determine the average location of each crop on the PCoA at each timepoint
swg2016.points.collapsed <- NULL
swg2016.samples.per.date <- NULL
for (i in 1:length(swg2016.dates)){
  x <- swg2016.points[swg.map.2016$sampling_date==swg2016.dates[i]&swg.map.2016$plant=="switchgrass",]
  swg2016.samples.per.date <- c(swg2016.samples.per.date, nrow(x))
  swg2016.points.collapsed <- rbind(swg2016.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
swg2016.points.collapsed <- as.data.frame(swg2016.points.collapsed)
colnames(swg2016.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
swg2016.points.collapsed$sampling_date <- swg2016.dates
swg2016.points.collapsed[,1:2] <- swg2016.points.collapsed[,1:2]/swg2016.samples.per.date

swg2017.points.collapsed <- NULL
swg2017.samples.per.date <- NULL
for (i in 1:length(swg2017.dates)){
  x <- swg2017.points[swg.map.2017$sampling_date==swg2017.dates[i]&swg.map.2017$plant=="switchgrass",]
  swg2017.samples.per.date <- c(swg2017.samples.per.date, nrow(x))
  swg2017.points.collapsed <- rbind(swg2017.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
swg2017.points.collapsed <- as.data.frame(swg2017.points.collapsed)
colnames(swg2017.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
swg2017.points.collapsed$sampling_date <- swg2017.dates
swg2017.points.collapsed[,1:2] <- swg2017.points.collapsed[,1:2]/swg2017.samples.per.date

mis2016.points.collapsed <- NULL
mis2016.samples.per.date <- NULL
for (i in 1:length(mis2016.dates)){
  x <- mis2016.points[mis.map.2016$sampling_date==mis2016.dates[i]&mis.map.2016$plant=="miscanthus",]
  mis2016.samples.per.date <- c(mis2016.samples.per.date, nrow(x))
  mis2016.points.collapsed <- rbind(mis2016.points.collapsed, c(colSums(x), sd(x[,1]), sd(x[,2])))
}
mis2016.points.collapsed <- as.data.frame(mis2016.points.collapsed)
colnames(mis2016.points.collapsed) <- c("Axis1", "Axis2", "sd_axis1", "sd_axis2")
mis2016.points.collapsed$sampling_date <- mis2016.dates
mis2016.points.collapsed[,1:2] <- mis2016.points.collapsed[,1:2]/mis2016.samples.per.date

soil.points.collapsed <- rbind(swg2016.points.collapsed, swg2017.points.collapsed, mis2016.points.collapsed)
soil.points.collapsed <- left_join(soil.points.collapsed, Date_to_Week, by="sampling_date")
soil.points.collapsed$plant <- c(rep("Switchgrass", 17), rep("Miscanthus", 10))
soil.points.collapsed$Week <- factor(x = soil.points.collapsed$Week, levels=c(1,2,3,4,5,6,7,8,9,10))
soil.points.collapsed$WeekNumeric <- as.numeric(soil.points.collapsed$Week)
soil.points.collapsed$Year <- c(rep(2016,9), rep(2017,8), rep(2016,10))
soil.points.collapsed$Year <- factor(soil.points.collapsed$Year, levels = c(2016,2017))
soil.points.collapsed$plant <- factor(soil.points.collapsed$plant, levels = c("Switchgrass", "Miscanthus"))
soil.points.collapsed <- soil.points.collapsed[order(soil.points.collapsed$sampling_date),]

soil2016.points.collapsed <- subset(soil.points.collapsed, soil.points.collapsed$Year==2016)
switchgrass2016.points.collapsed <- subset(soil2016.points.collapsed, soil2016.points.collapsed$plant=="Switchgrass")
miscanthus2016.point.collapsed <- subset(soil2016.points.collapsed, soil2016.points.collapsed$plant=="Miscanthus")
switchgrass2017.points.collapsed <- subset(soil.points.collapsed, soil.points.collapsed$Year==2017)

### Set up the weather maps for Weather Environmental fits
swg2016.weather <- swg.map.2016[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
swg2016.weather <- unique(swg2016.weather)

swg2017.weather <- swg.map.2017[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
swg2017.weather <- unique(swg2017.weather)

mis2016.weather <- mis.map.2016[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
mis2016.weather <- unique(mis2016.weather)

### Set up EnvFit for ggplot (remembering to use vegan:::ordiArrowMul to adjust arrow sizes)
soil.weather <- rbind(swg2016.weather, swg2017.weather, mis2016.weather)
soil.weather <- soil.weather[order(soil.weather$sampling_date),]
soil.weather$sampling_date == soil.points.collapsed$sampling_date
soil.weather$Week <- soil.points.collapsed$WeekNumeric

envfit.soil.weather <- envfit(soil.points.collapsed[,1:2], soil.weather[,c(1:11,14)])
envfit.soil.weather.df<-as.data.frame(scores(envfit.soil.weather, display = "vectors"))
colnames(envfit.soil.weather.df) <- c("Dim1", "Dim2")


soil.leaf.chemistry <- map_16S[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.soil.lc <- envfit(pcoa.16S$points, soil.leaf.chemistry, na.rm = TRUE)
envfit.soil.lc.df<-as.data.frame(scores(envfit.soil.lc, display = "vectors"))


soil.chemsitry <- map_16S[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.soil.sc <- envfit(pcoa.16S$points, soil.chemsitry)
envfit.soil.sc.df<-as.data.frame(scores(envfit.soil.sc, display = "vectors"))

envfit.soil.total <- rbind(envfit.soil.weather.df, envfit.soil.lc.df, envfit.soil.sc.df)

envfit.soil.total$r2 <- c(envfit.soil.weather$vectors$r, envfit.soil.lc$vectors$r, envfit.soil.sc$vectors$r)
envfit.soil.total$pval <- c(envfit.soil.weather$vectors$pvals, envfit.soil.lc$vectors$pvals, envfit.soil.sc$vectors$pvals)

scale_arrow <- function(arrows, data, at = c(0, 0), fill = 0.75) {
  u <- c(range(data[,1], range(data[,2])))
  u <- u - rep(at, each = 2)
  r <- c(range(arrows[, 1], na.rm = TRUE), range(arrows[, 2], na.rm = TRUE))
  rev <- sign(diff(u))[-2]
  
  if (rev[1] < 0) {
    u[1:2] <- u[2:1]
  }
  if (rev[2] < 0) {
    u[3:4] <- u[4:3]
  }
  u <- u/r
  u <- u[is.finite(u) & u > 0]
  invisible(fill * min(u))
}


arrow_scaling <- scale_arrow(arrows=envfit.soil.total[,1:2], data=soil.points.collapsed)

envfit.soil.total$Axis1 <- envfit.soil.total$Dim1 * arrow_scaling
envfit.soil.total$Axis2 <- envfit.soil.total$Dim2 * arrow_scaling
envfit.soil.total$Variable <- row.names(envfit.soil.total)

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/TableS2_EnvFit.txt", x=envfit.soil.total[,1:4], sep="\t", quote=FALSE)

### Subset to p <0.05 and r2>0.4
envfit.soil.sub <- subset (envfit.soil.total , envfit.soil.total$pval<0.05)
envfit.soil.sub <- subset(envfit.soil.sub, envfit.soil.sub$r2>0.4)




### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.soil <- pcoa.16S$eig[1]/sum(pcoa.16S$eig)
Ax2.soil <- pcoa.16S$eig[2]/sum(pcoa.16S$eig)

### Plot Plant PCoA
Fig2B1 <- ggplot(soil.points.collapsed, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, shape=Year)) +
  scale_color_manual(values=c("darkolivegreen3","darkgreen"))+
  scale_shape_manual(values = c(19,1))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=soil.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.soil,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.soil,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.soil.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.soil.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/Fig2B1_SoilPCoA.eps", Fig2B1, device = "eps", height = 6, width = 6, units = "in")

### Make Whole PCoA
pcoa.whole <- cmdscale(dist.otu, eig=TRUE)

### Set up dataframe for plotting
points.whole <- pcoa.whole$points
points.whole <- as.data.frame(points.whole)
colnames(points.whole) <- c("Axis1", "Axis2")
points.whole$plant <- factor(map_16S$plant, levels=c("switchgrass", "miscanthus"))
points.whole$Year <- factor(map_16S$Year, levels=c(2016,2017))
points.whole$sampling_date <- map_16S$sampling_date
points.whole <- left_join(points.whole, Date_to_Week, by="sampling_date")
points.whole$Week <- factor(points.whole$Week, levels=c(1:10))
points.whole$PlantYear <- paste(points.whole$plant,points.whole$Year, sep="")
points.whole$PlantYear <- factor(points.whole$PlantYear, levels=c("miscanthus2016", "switchgrass2016", "switchgrass2017"))

points.whole$Fert <- map_16S$treatment
points.whole$Fert <- factor(points.whole$Fert, levels = c("standard fertilization", "nitrogen free"))
points.whole$SampleType <- points.whole$Fert


### Determine % variation explained on each axis
Ax1.whole <- pcoa.whole$eig[1]/sum(pcoa.whole$eig)
Ax2.whole <- pcoa.whole$eig[2]/sum(pcoa.whole$eig)

### ggplot Whole PCoA
Fig2A <- ggplot(points.whole, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=PlantYear, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(15,0))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("red", "darkgreen","darkolivegreen3"))+
  scale_fill_manual(values=c("red", "darkgreen", "darkolivegreen3"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  geom_segment(data = envfit.soil.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.soil.sub, aes(label=Variable))+
  theme(legend.position = "none")+
  guides(fill=guide_legend(override.aes=list(colour=c(Miscanthus2016="darkgreen",Switchgrass2016="darkolivegreen3", Switchgrass2017="white"))), color=FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig2A_WholePCoA_NoLegend.eps", Fig2A, width = 6, height = 6, device = "eps", units="in")




### Plot Switch2016 PCoA
### Set up the weather maps for Weather Environmental fits
swg2016.weather <- swg.map.2016[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
swg2016.weather <- unique(swg2016.weather)

swg2016.weather <- swg2016.weather[order(swg2016.weather$sampling_date),]
swg2016.weather$sampling_date == switchgrass2016.points.collapsed$sampling_date
swg2016.weather$Week <- switchgrass2016.points.collapsed$WeekNumeric


envfit.swg2016.weather <- envfit(switchgrass2016.points.collapsed[,1:2], swg2016.weather[,c(1:11,14)])
envfit.swg2016.weather.df<-as.data.frame(scores(envfit.swg2016.weather, display = "vectors"))
colnames(envfit.swg2016.weather.df) <- c("Dim1", "Dim2")


swg2016.leaf.chemistry <- swg.map.2016[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.swg2016.lc <- envfit(pcoa.swg.2016$points, swg2016.leaf.chemistry, na.rm = TRUE)
envfit.swg2016.lc.df<-as.data.frame(scores(envfit.swg2016.lc, display = "vectors"))


swg2016.soil.chemsitry <- swg.map.2016[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.swg2016.sc <- envfit(pcoa.swg.2016$points, swg2016.soil.chemsitry)
envfit.swg2016.sc.df<-as.data.frame(scores(envfit.swg2016.sc, display = "vectors"))

envfit.swg2016.total <- rbind(envfit.swg2016.weather.df, envfit.swg2016.lc.df, envfit.swg2016.sc.df)

envfit.swg2016.total$r2 <- c(envfit.swg2016.weather$vectors$r, envfit.swg2016.lc$vectors$r, envfit.swg2016.sc$vectors$r)
envfit.swg2016.total$pval <- c(envfit.swg2016.weather$vectors$pvals, envfit.swg2016.lc$vectors$pvals, envfit.swg2016.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.swg2016.total[,1:2], data=swg2016.points.collapsed)

envfit.swg2016.total$Axis1 <- envfit.swg2016.total$Dim1 * arrow_scaling
envfit.swg2016.total$Axis2 <- envfit.swg2016.total$Dim2 * arrow_scaling
envfit.swg2016.total$Variable <- row.names(envfit.swg2016.total)

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/TableS3_swg2016EnvFit.txt", x=envfit.soil.total[,1:4], sep="\t", quote=FALSE)

### Subset to p <0.05 r2>0.4
envfit.swg2016.sub <- subset (envfit.swg2016.total , envfit.swg2016.total$pval<0.05)
envfit.swg2016.sub <- subset (envfit.swg2016.sub, envfit.swg2016.sub$r2>0.4)



### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.swg2016 <- pcoa.swg.2016$eig[1]/sum(pcoa.swg.2016$eig)
Ax2.swg2016 <- pcoa.swg.2016$eig[2]/sum(pcoa.swg.2016$eig)

### Plot Plant PCoA
Figswg2016 <- ggplot(switchgrass2016.points.collapsed, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, shape=Year)) +
  scale_color_manual(values=c("darkolivegreen3"))+
  scale_shape_manual(values = c(19))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=switchgrass2016.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=switchgrass2016.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=switchgrass2016.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=switchgrass2016.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.swg2016,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.swg2016,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.swg2016.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.swg2016.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/Fig_swg2016PCoA.eps", Figswg2016, device = "eps", height = 6, width = 6, units = "in")

### Make Whole swg2016 PCoA

### Set up dataframe for plotting
points.swg2016 <- pcoa.swg.2016$points
points.swg2016 <- as.data.frame(points.swg2016)
colnames(points.swg2016) <- c("Axis1", "Axis2")
points.swg2016$plant <- factor(swg.map.2016$plant, levels=c("switchgrass"))
points.swg2016$Year <- factor(swg.map.2016$Year, levels=c(2016))
points.swg2016$sampling_date <- swg.map.2016$sampling_date
points.swg2016<- left_join(points.swg2016, Date_to_Week, by="sampling_date")
points.swg2016$Week <- factor(points.swg2016$Week, levels=c(1:9))
points.swg2016$PlantYear <- paste(points.swg2016$plant,points.swg2016$Year, sep="")
points.swg2016$PlantYear <- factor(points.swg2016$PlantYear, levels=c("switchgrass2016"))

points.swg2016$Fert <- swg.map.2016$treatment
points.swg2016$Fert <- factor(points.swg2016$Fert, levels = c("standard fertilization", "nitrogen free"))
points.swg2016$SampleType <- points.swg2016$Fert

library(cowplot)

Fig2B <- ggplot(data = points.swg2016, aes(x=Axis1, y=Axis2), size=5)+
  coord_fixed()+
  geom_point(aes(shape=SampleType, size=Week, color=PlantYear))+
  scale_shape_manual(values = c(15,0))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("darkgreen"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "none")+
  guides(shape=guide_legend("SampleType"))+
  guides(color=guide_legend("PlantYear"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig2B.eps", Fig2B, width = 6, height = 6, device = "eps", units="in")





### Plot Miscanthus2016 PCoA
### Set up the weather maps for Weather Environmental fits
mis2016.weather <- mis.map.2016[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
mis2016.weather <- unique(mis2016.weather)

mis2016.weather <- mis2016.weather[order(mis2016.weather$sampling_date),]
mis2016.weather$sampling_date ==miscanthus2016.point.collapsed$sampling_date
mis2016.weather$Week <-miscanthus2016.point.collapsed$WeekNumeric


envfit.mis2016.weather <- envfit(miscanthus2016.point.collapsed[,1:2], mis2016.weather[,c(1:11,14)])
envfit.mis2016.weather.df<-as.data.frame(scores(envfit.mis2016.weather, display = "vectors"))
colnames(envfit.mis2016.weather.df) <- c("Dim1", "Dim2")


mis2016.leaf.chemistry <- mis.map.2016[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.mis2016.lc <- envfit(pcoa.mis.2016$points, mis2016.leaf.chemistry, na.rm = TRUE)
envfit.mis2016.lc.df<-as.data.frame(scores(envfit.mis2016.lc, display = "vectors"))


mis2016.soil.chemsitry <- mis.map.2016[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.mis2016.sc <- envfit(pcoa.mis.2016$points, mis2016.soil.chemsitry)
envfit.mis2016.sc.df<-as.data.frame(scores(envfit.mis2016.sc, display = "vectors"))

envfit.mis2016.total <- rbind(envfit.mis2016.weather.df, envfit.mis2016.lc.df, envfit.mis2016.sc.df)

envfit.mis2016.total$r2 <- c(envfit.mis2016.weather$vectors$r, envfit.mis2016.lc$vectors$r, envfit.mis2016.sc$vectors$r)
envfit.mis2016.total$pval <- c(envfit.mis2016.weather$vectors$pvals, envfit.mis2016.lc$vectors$pvals, envfit.mis2016.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.mis2016.total[,1:2], data=mis2016.points.collapsed)

envfit.mis2016.total$Axis1 <- envfit.mis2016.total$Dim1 * arrow_scaling
envfit.mis2016.total$Axis2 <- envfit.mis2016.total$Dim2 * arrow_scaling
envfit.mis2016.total$Variable <- row.names(envfit.mis2016.total)

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/TableS3_mis2016EnvFit.txt", x=envfit.soil.total[,1:4], sep="\t", quote=FALSE)

### Subset to p <0.05 r2>0.4
envfit.mis2016.sub <- subset (envfit.mis2016.total , envfit.mis2016.total$pval<0.05)
envfit.mis2016.sub <- subset (envfit.mis2016.sub , envfit.mis2016.sub$r2>0.4)


### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.mis2016 <- pcoa.mis.2016$eig[1]/sum(pcoa.mis.2016$eig)
Ax2.mis2016 <- pcoa.mis.2016$eig[2]/sum(pcoa.mis.2016$eig)

### Plot mis2016 PCoA
Figmis2016 <- ggplot(miscanthus2016.point.collapsed, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, shape=Year)) +
  scale_color_manual(values=c("darkgreen"))+
  scale_shape_manual(values = c(19))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=miscanthus2016.point.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=miscanthus2016.point.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=miscanthus2016.point.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=miscanthus2016.point.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.mis2016,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.mis2016,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.mis2016.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.mis2016.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/Fig_mis2016PCoA.eps", Figmis2016, device = "eps", height = 6, width = 6, units = "in")





### Make Whole mis2016 PCoA

### Set up dataframe for plotting
points.mis2016 <- pcoa.mis.2016$points
points.mis2016 <- as.data.frame(points.mis2016)
colnames(points.mis2016) <- c("Axis1", "Axis2")
points.mis2016$plant <- factor(mis.map.2016$plant, levels=c("miscanthus"))
points.mis2016$Year <- factor(mis.map.2016$Year, levels=c(2016))
points.mis2016$sampling_date <- mis.map.2016$sampling_date
points.mis2016<- left_join(points.mis2016, Date_to_Week, by="sampling_date")
points.mis2016$Week <- factor(points.mis2016$Week, levels=c(1:10))
points.mis2016$PlantYear <- paste(points.mis2016$plant,points.mis2016$Year, sep="")
points.mis2016$PlantYear <- factor(points.mis2016$PlantYear, levels=c("miscanthus2016"))

points.mis2016$Fert <- mis.map.2016$treatment
points.mis2016$Fert <- factor(points.mis2016$Fert, levels = c("standard fertilization", "nitrogen free"))
points.mis2016$SampleType <- points.mis2016$Fert

library(cowplot)

Fig2C <- ggplot(data = points.mis2016, aes(x=Axis1, y=Axis2), size=5)+
  coord_fixed()+
  geom_point(aes(shape=SampleType, size=Week, color=PlantYear))+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values = c(15,0))+
  scale_size_manual(values=Point_Sizes)+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  geom_segment(data = envfit.mis2016.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.mis2016.sub, aes(label=Variable))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig2C.eps", Fig2C, width = 6, height = 6, device = "eps", units="in")





### Plot Switch2017 PCoA
### Set up the weather maps for Weather Environmental fits
swg2017.weather <- swg.map.2017[,c("precipitation", "Air_temp_mean", "air_temp_max", "Air_Temp_Min", "Air_Pressure", "RH", "AH", "Wind_Speed_Mean", "Solar_Radiation", "PAR", "soil_temp_5_cm_bare_avg", "sampling_date", "time_numeric")]
swg2017.weather <- unique(swg2017.weather)

swg2017.weather <- swg2017.weather[order(swg2017.weather$sampling_date),]
swg2017.weather$sampling_date == switchgrass2017.points.collapsed$sampling_date
swg2017.weather$Week <- switchgrass2017.points.collapsed$WeekNumeric


envfit.swg2017.weather <- envfit(switchgrass2017.points.collapsed[,1:2], swg2017.weather[,c(1:11,14)])
envfit.swg2017.weather.df<-as.data.frame(scores(envfit.swg2017.weather, display = "vectors"))
colnames(envfit.swg2017.weather.df) <- c("Dim1", "Dim2")


swg2017.leaf.chemistry <- swg.map.2017[,c("LDMC_mg_per_g", "nitrogen_percent", "carbon_percent", "carbon_per_nitrogen", "height_mean_cm")]
envfit.swg2017.lc <- envfit(pcoa.swg.2017$points, swg2017.leaf.chemistry, na.rm = TRUE)
envfit.swg2017.lc.df<-as.data.frame(scores(envfit.swg2017.lc, display = "vectors"))


swg2017.soil.chemsitry <- swg.map.2017[,c("pH", "P_ppm", "K_ppm", "Ca_ppm", "Mg_ppm", "organic_matter", "NO3N_ppm", "NH4_ppm", "soil_moisture_percent", "soil_temp_10cm")]
envfit.swg2017.sc <- envfit(pcoa.swg.2017$points, swg2017.soil.chemsitry)
envfit.swg2017.sc.df<-as.data.frame(scores(envfit.swg2017.sc, display = "vectors"))

envfit.swg2017.total <- rbind(envfit.swg2017.weather.df, envfit.swg2017.lc.df, envfit.swg2017.sc.df)

envfit.swg2017.total$r2 <- c(envfit.swg2017.weather$vectors$r, envfit.swg2017.lc$vectors$r, envfit.swg2017.sc$vectors$r)
envfit.swg2017.total$pval <- c(envfit.swg2017.weather$vectors$pvals, envfit.swg2017.lc$vectors$pvals, envfit.swg2017.sc$vectors$pvals)

arrow_scaling <- scale_arrow(arrows=envfit.swg2017.total[,1:2], data=swg2017.points.collapsed)

envfit.swg2017.total$Axis1 <- envfit.swg2017.total$Dim1 * arrow_scaling
envfit.swg2017.total$Axis2 <- envfit.swg2017.total$Dim2 * arrow_scaling
envfit.swg2017.total$Variable <- row.names(envfit.swg2017.total)

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/TableS3_swg2017EnvFit.txt", x=envfit.soil.total[,1:4], sep="\t", quote=FALSE)

### Subset to p <0.05, r2>0.4
envfit.swg2017.sub <- subset (envfit.swg2017.total , envfit.swg2017.total$pval<0.05)
envfit.swg2017.sub <- subset(envfit.swg2017.sub, envfit.swg2017.sub$r2>0.4)



### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.swg2017 <- pcoa.swg.2017$eig[1]/sum(pcoa.swg.2017$eig)
Ax2.swg2017 <- pcoa.swg.2017$eig[2]/sum(pcoa.swg.2017$eig)

### Plot Plant PCoA
Figswg2017 <- ggplot(switchgrass2017.points.collapsed, aes(x=Axis1, y=Axis2))+
  geom_point(aes(size=Week, color=plant, shape=Year)) +
  scale_color_manual(values=c("darkolivegreen3"))+
  scale_shape_manual(values = c(1))+
  scale_size_manual(values=Point_Sizes)+
  coord_fixed()+
  geom_segment(data=switchgrass2017.points.collapsed, aes(x=Axis1,xend=Axis1+sd_axis1,y=Axis2,yend=Axis2,color=plant ))+
  geom_segment(data=switchgrass2017.points.collapsed, aes(x=Axis1,xend=Axis1-sd_axis1,y=Axis2,yend=Axis2,color=plant))+
  geom_segment(data=switchgrass2017.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2+sd_axis2,color=plant))+
  geom_segment(data=switchgrass2017.points.collapsed, aes(x=Axis1,xend=Axis1,y=Axis2,yend=Axis2-sd_axis2,color=plant))+
  xlab(label = paste(round(Ax1.swg2017,digits = 3)*100, "% Var. Explained", sep = ""))+
  ylab(label= paste(round(Ax2.swg2017,digits = 3)*100, "% Var. Explained", sep = ""))+
  geom_segment(data = envfit.swg2017.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.swg2017.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/Fig_swg2017PCoA.eps", Figswg2017, device = "eps", height = 6, width = 6, units = "in")

### Make Whole swg2017 PCoA

### Set up dataframe for plotting
points.swg2017 <- pcoa.swg.2017$points
points.swg2017 <- as.data.frame(points.swg2017)
colnames(points.swg2017) <- c("Axis1", "Axis2")
points.swg2017$plant <- factor(swg.map.2017$plant, levels=c("switchgrass"))
points.swg2017$Year <- factor(swg.map.2017$Year, levels=c(2017))
points.swg2017$sampling_date <- swg.map.2017$sampling_date
points.swg2017<- left_join(points.swg2017, Date_to_Week, by="sampling_date")
points.swg2017$Week <- factor(points.swg2017$Week, levels=c(1:8))
points.swg2017$PlantYear <- paste(points.swg2017$plant,points.swg2017$Year, sep="")
points.swg2017$PlantYear <- factor(points.swg2017$PlantYear, levels=c("switchgrass2017"))

points.swg2017$Fert <- swg.map.2017$treatment
points.swg2017$Fert <- factor(points.swg2017$Fert, levels = c("standard fertilization", "nitrogen free"))
points.swg2017$SampleType <- points.swg2017$Fert

library(cowplot)

Fig2D <- ggplot(data = points.swg2017, aes(x=Axis1, y=Axis2), size=5)+
  coord_fixed()+
  geom_point(aes(shape=SampleType, size=Week, color=PlantYear))+
  scale_color_manual(values=c("darkolivegreen3"))+
  scale_shape_manual(values = c(15,0))+
  scale_size_manual(values=Point_Sizes)+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  geom_segment(data = envfit.swg2017.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.swg2017.sub, aes(label=Variable))+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig2D.eps", Fig2D, width = 6, height = 6, device = "eps", units="in")



Fig2A_leg <- ggplot(points.whole, aes(x=Axis1, y=Axis2))+
  coord_fixed()+
  geom_point(aes(shape=SampleType, color=PlantYear, size=Week, fill=PlantYear))+
  scale_shape_manual(values = c(15,0))+
  scale_size_manual(values=Point_Sizes)+
  scale_color_manual(values=c("red", "darkgreen","darkolivegreen3"))+
  scale_fill_manual(values=c("red", "darkgreen","darkolivegreen3"))+
  xlab(label = paste(round(Ax1.whole, 3)*100, "% Var. Explained", sep=""))+
  ylab(label= paste(round(Ax2.whole, 3)*100, "% Var. Explained", sep=""))+
  theme(legend.position = "right", legend.box = "horizontal")+
  guides(fill=guide_legend(override.aes=list(shape=21, color=c(miscanthus2016="red", swtichgrass2016="darkgreen",switchgrass2017="darkolivegreen3"),fill=c(miscanthus2016="red", switchgrass2016="darkgreen",switchgrass2017="darkolivegreen3"), size=3)), color=FALSE, shape=guide_legend(override.aes = list(size=3)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
library(cowplot)
Fig2_Legend <- get_legend(Fig2A_leg)

setEPS()
postscript("C:/Users/jaffyhu/Desktop/Soil/new/Figure2.eps", width=10, height=10,pointsize=10, paper="special")
plot_grid(plot_grid(Fig2A, Fig2B, Fig2C, Fig2D, align = "h"), Fig2_Legend, ncol=1)
dev.off()

###########################
#' Venn diagram - Figure S7
###########################
#' Subsetting the data to the crop, year and treatment
otu_rare=otu_rare[rowSums(otu_rare)>0,]

misc_otu <- otu_rare[,map_16S$plant=="miscanthus"]
nfmisc_otu <- otu_rare[,map_16S$plant=="miscanthus" & (map_16S$treatment=='nitrogen free')]
sfmisc_otu <- otu_rare[,map_16S$plant=="miscanthus" & (map_16S$treatment=='standard fertilization')]

swit16_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2016)]
nfswit16_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2016)& (map_16S$treatment=='nitrogen free')]
sfswit16_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2016)& (map_16S$treatment=='standard fertilization')]

swit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2017)]
nfswit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2017) & (map_16S$treatment=='nitrogen free')]
sfswit17_otu <- otu_rare[,map_16S$plant=="switchgrass" & (map_16S$year==2017) & (map_16S$treatment=='standard fertilization')]

# make presence absence list from soil into 1 & 0
swit16_otu_venn <- 1*(rowSums(swit16_otu)>0)
swit17_otu_venn <- 1*(rowSums(swit17_otu)>0)
misc_otu_venn <- 1*(rowSums(misc_otu)>0)
nfswit16_otu_venn <- 1*(rowSums(nfswit16_otu)>0)
sfswit16_otu_venn <- 1*(rowSums(sfswit16_otu)>0)
nfswit17_otu_venn<- 1*(rowSums(nfswit17_otu)>0)
sfswit17_otu_venn<- 1*(rowSums(sfswit17_otu)>0)
nfmisc16_otu_venn <- 1*(rowSums(nfmisc_otu)>0)
sfmisc16_otu_venn <- 1*(rowSums(sfmisc_otu)>0)

#' plot Venn
venn_soil_data <- cbind(nfswit16_otu_venn, sfswit16_otu_venn, nfswit17_otu_venn, sfswit17_otu_venn, nfmisc16_otu_venn, sfmisc16_otu_venn)
colnames(venn_soil_data) <- c("Switchgrass 2016NF", "Switchgrass 2016SF", "Switchgrass 2017NF", "Switchgrass 2017SF", "Miscanthus 2016NF", "Miscanthus 2016SF")
venn_soil_data=venn_soil_data[rowSums(venn_soil_data)>0,]
v_soil=vennCounts(venn_soil_data)
v_soil_2=round(v_soil[,"Counts"]/sum(v_soil[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_soil, circle.col = c('darkolivegreen3', 'darkgreen','blue', 'sky blue', 'red', 'brown'), lwd=6, cex=1.2, scale=F)

# plot Venn switchgrass and treatment
venn_soil_data <- cbind(nfswit16_otu_venn, sfswit16_otu_venn, nfswit17_otu_venn, sfswit17_otu_venn)
colnames(venn_soil_data) <- c("Switchgrass16NF", "Switchgrass16SF", "Switchgrass17NF", "Switchgrass17SF")
venn_soil_data=venn_soil_data[rowSums(venn_soil_data)>0,]
v_soil=vennCounts(venn_soil_data)
v_soil_2=round(v_soil[,"Counts"]/sum(v_soil[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_soil, circle.col = c('darkolivegreen3', 'darkgreen','blue', 'sky blue'), lwd=6, cex=1.2, scale=F)

# plot Venn switchgrass and miscanthus 16, and treatment
venn_soil_data <- cbind(nfswit16_otu_venn, sfswit16_otu_venn, nfmisc16_otu_venn, sfmisc16_otu_venn)
colnames(venn_soil_data) <- c("Switchgrass 2016NF", "Switchgrass 2016SF", "Miscanthus 2016NF", "Miscanthus 2016SF")
venn_soil_data=venn_soil_data[rowSums(venn_soil_data)>0,]
v_soil=vennCounts(venn_soil_data)
v_soil_2=round(v_soil[,"Counts"]/sum(v_soil[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_soil, circle.col = c('darkolivegreen3', 'darkgreen','blue', 'sky blue', 'red', 'brown'), lwd=6, cex=1.2, scale=F)

# plot Venn switchgrass and miscanthus
venn_soil_data <- cbind(swit16_otu_venn, swit17_otu_venn, misc_otu_venn)
colnames(venn_soil_data) <- c("Switchgrass 2016", "Switchgrass 2017", "Miscanthus 2016")
venn_soil_data=venn_soil_data[rowSums(venn_soil_data)>0,]
v_soil=vennCounts(venn_soil_data)
v_soil_2=round(v_soil[,"Counts"]/sum(v_soil[,"Counts"]),2) #calculate percentage of each group
vennDiagram(v_soil, circle.col = c('darkolivegreen3', 'blue', 'red'), lwd=6, cex=1.2, scale=F)

############################################################################################
######################### BACTERIAL COMMUNITIES COMPOSITION #####################
############################################################################################
BiocManager::install("phyloseq")
library(phyloseq)


# BACTERIA COMPOSITION
# read bacterial taxonomy
taxonomy_filtered
rownames(taxonomy_filtered) <- rownames(otu_filtered)
# make phyloseq otu table and taxonomy
OTU = otu_table(as.matrix(otu_filtered), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxonomy_filtered))
# add map
map_16S$sampling_week<-as.factor(map_16S$sampling_week)
rownames(map_16S) <- map_16S$sequence_name
# make phyloseq map
phyloseq_map <- sample_data(map_16S)
# make phyloseq object
PHYL_16S <- merge_phyloseq(OTU,TAX,phyloseq_map)
PHYL_16S

PHYL_16S_phylum <- PHYL_16S %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)

PHYL_16S_phylum$PlantYearTreat <- paste(PHYL_16S_phylum$plant,PHYL_16S_phylum$Year, PHYL_16S_phylum$treatment, sep="")

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861","#ffffff", "#000000"
)

# Plot by sampling week
Fig5 <- ggplot(PHYL_16S_phylum, aes(x = sampling_week, y = Abundance, fill = Phylum)) + 
  facet_grid(PlantYearTreat~.) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition \n Bacterial Communities by Sampling Week")  
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Figure5.eps", Fig5, device = "eps", width = 12, height = 7, units= "in", dpi = 600)



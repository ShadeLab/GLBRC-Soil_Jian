# Total R Code for GLBRC soil analyses
### Start Common
### DIRECTIONS, SET WORKING DIRECTORY TO "C:/Users/jaffyhu/Documents/MobaXterm/home/PAPER_GradySorensenStopnisek_NatComm_2019/R/Soil>"
### Before working on Analysis, run the entire Common block of Code (all the way until "### End Common")
setwd("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC Soil")

install.packages(c("agricolae", "BiocManager", "car", "codyn", "cowplot", "dplyr", "egg", "FSA", "ggplot2", "ggpubr", "ggrepel", "gridExtra", "moments", "rcompanion", "reshape", "reshape2", "RSQLite", "scales", "stringr", "tibble", "tidyr", "vegan"))

library(agricolae)
library(BiocManager)
library(car)
library(codyn)
library(cowplot)
library(dplyr)
library(egg)
library(FSA)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(moments)
library(rcompanion)
library(reshape)
library(reshape2)
library(RSQLite)
library(scales)
library(stringr)
library(tibble)
library(tidyr)
library(vegan)




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

# Read in OTU table
otu <- read.table("InputFiles/otu.txt",sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)

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

taxonomy_filtered <- taxonomy_full[rownames(taxonomy_full) %in% rownames(otu_filtered),]


#rarecurve of soil samples

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


#########################
#### Alpha Diversity ####
#########################

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
Fig1A
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig1A_SpeciesAccumulation.eps", Fig1A, device = "eps", width=6, height=4, units = "in")


#Accumulative observed taxa under different treatment, plant type and year.

### Setting up Contextual Data Maps

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


misc.map.2016 <- map_16S[map_16S$plant=="miscanthus",]
misc.map.2016nf <- misc.map.2016[misc.map.2016$treatment=="nitrogen free",]
misc.map.2016sf <- misc.map.2016[misc.map.2016$treatment=="standard fertilization",]

switch.map.2016 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2016,]
switch.map.2016nf <- switch.map.2016[switch.map.2016$treatment=="nitrogen free",]
switch.map.2016sf <- switch.map.2016[switch.map.2016$treatment=="standard fertilization",]

switch.map.2017 <- map_16S[map_16S$plant=="switchgrass"&map_16S$year==2017,]
switch.map.2017nf <- switch.map.2017[switch.map.2017$treatment=="nitrogen free",]
switch.map.2017sf <- switch.map.2017[switch.map.2017$treatment=="standard fertilization",]


cor.test(map.div.swg16nf$time_numeric, map.div.swg16nf$Richness) #t = -0.47864, df = 29, p-value = 0.6358
cor.test(map.div.swg16sf$time_numeric, map.div.swg16sf$Richness) #t = -2.2814, df = 33, p-value = 0.02911
cor.test(map.div.swg16nf$time_numeric, map.div.swg16nf$Shannon) #t = -1.186, df = 29, p-value = 0.2453
cor.test(map.div.swg16sf$time_numeric, map.div.swg16sf$Shannon) #t = -2.815, df = 33, p-value = 0.008163
cor.test(map.div.swg16nf$time_numeric, map.div.swg16nf$Pielou) #t = -1.6564, df = 29, p-value = 0.1084
cor.test(map.div.swg16sf$time_numeric, map.div.swg16sf$Pielou) #t = -2.9745, df = 33, p-value = 0.005453


cor.test(map.div.swg17nf$time_numeric, map.div.swg17nf$Richness) #t = -3.1405, df = 30, p-value = 0.003773
cor.test(map.div.swg17sf$time_numeric, map.div.swg17sf$Richness) #t = -4.6105, df = 28, p-value = 8.041e-05
cor.test(map.div.swg17nf$time_numeric, map.div.swg17nf$Shannon) #t = -2.1666, df = 30, p-value = 0.03833
cor.test(map.div.swg17sf$time_numeric, map.div.swg17sf$Shannon) #t = -3.7512, df = 28, p-value = 0.0008154
cor.test(map.div.swg17nf$time_numeric, map.div.swg17nf$Pielou) #t = -1.2653, df = 30, p-value = 0.2155
cor.test(map.div.swg17sf$time_numeric, map.div.swg17sf$Pielou) #t = -2.6138, df = 28, p-value = 0.01425

cor.test(map.div.misnf$time_numeric, map.div.misnf$Richness) #t = -0.21517, df = 33, p-value = 0.831
cor.test(map.div.missf$time_numeric, map.div.missf$Richness) #t = -0.21312, df = 37, p-value = 0.8324
cor.test(map.div.misnf$time_numeric, map.div.misnf$Shannon) #t = -0.55241, df = 33, p-value = 0.5844
cor.test(map.div.missf$time_numeric, map.div.missf$Shannon) #t = -0.35131, df = 37, p-value = 0.7274
cor.test(map.div.misnf$time_numeric, map.div.misnf$Pielou) #t = -0.78981, df = 33, p-value = 0.4353
cor.test(map.div.missf$time_numeric, map.div.missf$Pielou) #t = -0.4514, df = 37, p-value = 0.6543



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


Fig1B <- ggplot(Spec_accum_total, aes(x=sampling_week, y=Species_Accumulation)) + 
  geom_point(aes(color=PlantFertilization), size=2) +
  geom_line(aes(color=PlantFertilization, linetype=FertilizationStatus))+
  scale_shape_manual()+
  scale_linetype_manual(values = c(2,1))+
  theme_bw()+
  labs(y="Total Observed Taxa") +
  scale_color_manual(values=c("green", "green", "red", "red", "blue", "blue"))  +
  guides(shape=FALSE)
Fig1B
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig1B_SpeciesAccumulation.eps", Fig1B, device = "eps", width=6, height=4, units = "in")


# Plot Pielou Evenness Through Season
switch.Pielou.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Pielou",]
switch.Pielou.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Pielou",]
misc.Pielou.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Pielou",]


Pielou.df <- rbind(switch.Pielou.2016, switch.Pielou.2017, misc.Pielou.2016)
Pielou.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Pielou.2016)), rep("Switchgrass 2017", nrow(switch.Pielou.2017)), rep("Miscanthus 2016", nrow(misc.Pielou.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Pielou.df$FakeDate <- as.Date(gsub(x=Pielou.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Pielou.df$treatment <- factor(Pielou.df$treatment, levels=c("nitrogen free", "standard fertilization"))


Pielou <- ggplot(data = Pielou.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Pielou")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")
Pielou
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Pielou.eps", Pielou, device = "eps", width=6, height=4, units = "in")

# Plot Richness Through Season
switch.Richness.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Richness",]
switch.Richness.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Richness",]
misc.Richness.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Richness",]


Richness.df <- rbind(switch.Richness.2016, switch.Richness.2017, misc.Richness.2016)
Richness.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Richness.2016)), rep("Switchgrass 2017", nrow(switch.Richness.2017)), rep("Miscanthus 2016", nrow(misc.Richness.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Richness.df$FakeDate <- as.Date(gsub(x=Richness.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Richness.df$treatment <- factor(Richness.df$treatment, levels=c("nitrogen free", "standard fertilization"))


Richness <- ggplot(data = Richness.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Richness")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")
Richness
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Richness.eps", Richness, device = "eps", width=6, height=4, units = "in")

# Plot Shannon Through Season
switch.Shannon.2016 <- map.alpha.swg16[map.alpha.swg16$variable=="Shannon",]
switch.Shannon.2017 <- map.alpha.swg17[map.alpha.swg17$variable=="Shannon",]
misc.Shannon.2016 <- map.alpha.mis16[map.alpha.mis16$variable=="Shannon",]


Shannon.df <- rbind(switch.Shannon.2016, switch.Shannon.2017, misc.Shannon.2016)
Shannon.df$Factor <- factor(c(rep("Switchgrass 2016", nrow(switch.Shannon.2016)), rep("Switchgrass 2017", nrow(switch.Shannon.2017)), rep("Miscanthus 2016", nrow(misc.Shannon.2016))), levels = c("Miscanthus 2016", "Switchgrass 2016", "Switchgrass 2017"))
Shannon.df$FakeDate <- as.Date(gsub(x=Shannon.df$sampling_Rdate, pattern = "2017", replacement = "2016" ))

Shannon.df$treatment <- factor(Shannon.df$treatment, levels=c("nitrogen free", "standard fertilization"))


Shannon <- ggplot(data = Shannon.df, aes(x=FakeDate, y=value))+
  geom_boxplot(mapping=aes(group=FakeDate, fill=Factor, color=as.character(year)), width=8, position="dodge")+
  facet_grid(rows=vars(treatment), cols=vars(Factor), scales = "free")+
  scale_fill_manual(values=c("darkgreen", "darkolivegreen3", "White"))+
  scale_color_manual(values=c("black","darkolivegreen3"))+
  scale_x_date(date_breaks = "1 month", labels=date_format("%b") ) +
  ylab(label = "Shannon")+
  xlab(label="Date")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.75), legend.position = "none")
Shannon
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/FigureS2_Shannon.eps", Shannon, device = "eps", width=6, height=4, units = "in")


# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.misnf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial richness ##

### 1. Compare bacterial richness among sampling_week using one-way ANOVA
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
leveneTest(Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.9844
boxplot(mis.map_aov$Richness ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.308 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 2.69 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 9.0175, df = 9, p-value = 0.4357
ggboxplot(mis.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial richness among samplingweeks under Nitrogen free

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "output/DT_RC_samplingweek.df.csv")
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

## Statistical test of bacterial richness ##

### 1. Compare bacterial richness among sampling_week using one-way ANOVA
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
leveneTest(Richness ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.03681
boxplot(mis.map_aov$Richness ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # 0.349 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.245 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 14.487, df = 9, p-value = 0.106
ggboxplot(mis.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial richness among samplingweeks under standard fertilization

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Richness~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "output/DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Richness_samplingweek<- HSD.test(Aov_richness_samplingweek, "sampling_week", group = F, console = TRUE)


#Switchgrass 2016 Nitrogen free
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week
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
leveneTest(Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.5138
boxplot(swg16.map_aov$Richness ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.6058 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 19.114, df = 8, p-value = 0.01426
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial richness among samplingweeks

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
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.1997, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.7983
boxplot(swg16.map_aov$Richness ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.3568 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 15.441, df = 8, p-value = 0.05111
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial richness among samplingweeks

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
leveneTest(Richness ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.5879
boxplot(swg17.map_aov$Richness ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.061 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal richness among sampling_week using Kruskal-Wallis Test
kruskal.test(Richness ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 18.441, df = 7, p-value = 0.01013
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Richness")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial richness among samplingweeks

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

#Shannon 
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.misnf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial  Shannon ##

### 1. Compare bacterial Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(mis.map_aov$Shannon ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.2357
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.6159, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.7408
boxplot(mis.map_aov$Shannon ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.308 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 2.69 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 10.985, df = 9, p-value = 0.2767
ggboxplot(mis.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial Shannon among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)

#Standard fertilization miscanthus2016
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.missf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial Shannon ##

### 1. Compare bacterial Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(mis.map_aov$Shannon ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.1224
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.1084, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.08686
boxplot(mis.map_aov$Shannon ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # 0.359 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.245 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 14.487, df = 9, p-value = 0.106
ggboxplot(mis.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial Shannon among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)


#Switchgrass 2016 Nitrogen free
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16nf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal Shannon in switch2016##

### 1. Compare bacterial and archaeal Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(swg16.map_aov$Shannon ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.01936
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.5204, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.72
boxplot(swg16.map_aov$Shannon ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.535 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 18.023, df = 8, p-value = 0.02105
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal Shannon among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)

#Switchgrass 2016 Standard Fertilization
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16sf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial Shannon in switch2016##

### 1. Compare bacterial and archaeal Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(swg16.map_aov$Shannon ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.004623
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.8957, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.6959
boxplot(swg16.map_aov$Shannon ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.039 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 18.307, df = 8, p-value = 0.01904
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial and archaeal Shannon among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)



# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2017 
swg17.map_aov <- map.div.swg17nf
class(swg17.map_aov)
swg17.map_aov$sampling_week<-as.factor(swg17.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg17.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal Shannon in switch2017##

### 1. Compare bacterial and archaeal Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(swg17.map_aov$Shannon ~ sampling_week, data=swg17.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.1394
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.7602, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.7397
boxplot(swg17.map_aov$Shannon ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.0268 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 10.188, df = 7, p-value = 0.1782
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal Shannon among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, swg17.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)


# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2017 
swg17.map_aov <- map.div.swg17sf
class(swg17.map_aov)
swg17.map_aov$sampling_week<-as.factor(swg17.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg17.map_aov, give.attr=F)

## Statistical test of bacterial and archaeal Shannon in switch2017##

### 1. Compare bacterial and archaeal Shannon among sampling_week using one-way ANOVA
Aov_Shannon_samplingweek <- lm(swg17.map_aov$Shannon ~ sampling_week, data=swg17.map_aov, na.action=na.exclude)
Aov_Shannon_samplingweek
drop1(Aov_Shannon_samplingweek,~.,test="F") # type III SS and F Tests 0.002297
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Shannon_samplingweek)
RC_samplingweek_preds <- predict(Aov_Shannon_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.5255, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.3296
boxplot(swg17.map_aov$Shannon ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.2889 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Shannon among sampling_week using Kruskal-Wallis Test
kruskal.test(Shannon ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 18.034, df = 7, p-value = 0.01182
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Shannon")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial and archaeal Shannon among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Shannon~sampling_week, swg17.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Shannon_samplingweek<- HSD.test(Aov_Shannon_samplingweek, "sampling_week", group = F, console = TRUE)

#Pielou
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.misnf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial Pielou ##

### 1. Compare bacterial Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(mis.map_aov$Pielou ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.2389
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.3014, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.3293
boxplot(mis.map_aov$Pielou ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) #0.0759 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 2.69 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 11.887, df = 9, p-value = 0.2197
ggboxplot(mis.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal Pielou among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)

#Standard fertilization miscanthus2016
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in misnf 
mis.map_aov <- map.div.missf
class(mis.map_aov)
mis.map_aov$sampling_week<-as.factor(mis.map_aov$sampling_week) # inform R that sampling_week is factor
str(mis.map_aov, give.attr=F)

## Statistical test of bacterial Pielou ##

### 1. Compare bacterial Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(mis.map_aov$Pielou ~ sampling_week, data=mis.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.08423
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.5027, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=mis.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.3878
boxplot(mis.map_aov$Pielou ~ sampling_week, data = mis.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # 0.31799 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.245 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = mis.map_aov) # Kruskal-Wallis chi-squared = 10.918, df = 9, p-value = 0.2814
ggboxplot(mis.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal Pielou among samplingweeks under NF

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, mis.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)


#Switchgrass 2016 Nitrogen free
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16nf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial Pielou in switch2016 nitrogen free

### 1. Compare bacterial Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(swg16.map_aov$Pielou ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.01679
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.8159, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.7672
boxplot(swg16.map_aov$Pielou ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.3736 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 19.114, df = 8, p-value = 0.01426
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial Pielou among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)

#Switchgrass 2016 Standard Fertilization
# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2016 
swg16.map_aov <- map.div.swg16sf
class(swg16.map_aov)
swg16.map_aov$sampling_week<-as.factor(swg16.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg16.map_aov, give.attr=F)

## Statistical test of bacterial  Pielou in switch2016 standard fertilization

### 1. Compare bacterial Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(swg16.map_aov$Pielou ~ sampling_week, data=swg16.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.001379
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.847, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=swg16.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.8272
boxplot(swg16.map_aov$Pielou ~ sampling_week, data = swg16.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.2102 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.7748 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = swg16.map_aov) # Kruskal-Wallis chi-squared = 19.744, df = 8, p-value = 0.01135
ggboxplot(swg16.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial Pielou among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, swg16.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)



# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2017 
swg17.map_aov <- map.div.swg17nf
class(swg17.map_aov)
swg17.map_aov$sampling_week<-as.factor(swg17.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg17.map_aov, give.attr=F)

## Statistical test of bacterial Pielou in switch2017##

### 1. Compare bacterial and archaeal Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(swg17.map_aov$Pielou ~ sampling_week, data=swg17.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.2739
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.06976, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.6471
boxplot(swg17.map_aov$Pielou ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # -0.061 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial and archaeal Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 8.4091, df = 7, p-value = 0.2979
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are no significant differences of bacterial and archaeal Pielou among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, swg17.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)


# ANOVA and Kruskal-Wallis test to compare alpha diversity among sampling_week, treatment in switch2017 
swg17.map_aov <- map.div.swg17sf
class(swg17.map_aov)
swg17.map_aov$sampling_week<-as.factor(swg17.map_aov$sampling_week) # inform R that sampling_week is factor
str(swg17.map_aov, give.attr=F)

## Statistical test of bacterial Pielou in switch2017##

### 1. Compare bacterial Pielou among sampling_week using one-way ANOVA
Aov_Pielou_samplingweek <- lm(swg17.map_aov$Pielou ~ sampling_week, data=swg17.map_aov, na.action=na.exclude)
Aov_Pielou_samplingweek
drop1(Aov_Pielou_samplingweek,~.,test="F") # type III SS and F Tests 0.01323
# testing assumptions
# Generate residual and predicted values
RC_samplingweek_resids <- residuals(Aov_Pielou_samplingweek)
RC_samplingweek_preds <- predict(Aov_Pielou_samplingweek)
# Look at a plot of residual vs. predicted values
plot(RC_samplingweek_resids ~ RC_samplingweek_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_samplingweek_resids) # ALERT!! p-value = 0.5755, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ sampling_week, data=swg17.map_aov, na.action=na.exclude) # GOOD variances among group are homogenous 0.2579
boxplot(swg17.map_aov$Pielou ~ sampling_week, data = swg17.map_aov) # there are no outliers
plot(density(RC_samplingweek_resids)) # density is not bad
qqnorm(RC_samplingweek_resids)
qqline(RC_samplingweek_resids) # I think the data normality is fine
hist(RC_samplingweek_resids)

skew_xts <-  skewness(RC_samplingweek_resids) # 0.2258 (within range of -2 to 2)

# use libraries in Fina's paper (2020) 
kurtosis(RC_samplingweek_resids,method = 'sample') # 3.079 (should be in the range of -7 to 7)

### 1. Compare bacterial Pielou among sampling_week using Kruskal-Wallis Test
kruskal.test(Pielou ~ sampling_week, data = swg17.map_aov) # Kruskal-Wallis chi-squared = 14.927, df = 7, p-value = 0.03695
ggboxplot(swg17.map_aov, x = "sampling_week", y = "Pielou")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial Pielou among samplingweeks

# Do Post Hoc Dunn's Test
DT_RC_samplingweek <- dunnTest(Pielou~sampling_week, swg17.map_aov, method = "bh", kw=TRUE)
print(DT_RC_samplingweek,dunn.test.results=TRUE)
DT_RC_samplingweek$res
DT_RC_samplingweek.df <- as.data.frame(DT_RC_samplingweek$res)
write.csv(DT_RC_samplingweek.df, file = "DT_RC_samplingweek.df.csv")
DT_RC_samplingweek_letter = cldList(P.adj ~ Comparison,
                                    data = DT_RC_samplingweek$res,
                                    threshold = 0.05)

# Do Tukey's HSD Post Hoc Test
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = TRUE, console = TRUE)
hsd_Pielou_samplingweek<- HSD.test(Aov_Pielou_samplingweek, "sampling_week", group = F, console = TRUE)


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

PHYL_16S_Family <- PHYL_16S %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at Family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.03) %>%                         # Filter out low abundance taxa
  arrange(Family)

PHYL_16S_Family$PlantYearTreat <- paste(PHYL_16S_Family$plant,PHYL_16S_Family$Year, PHYL_16S_Family$treatment, sep="")

Family_colors <- c(
  '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff')

p.soil <- ggplot(data=PHYL_16S_Family, aes(x=sampling_week, y=Abundance, fill=Family))
barplot.baccomp <- p.soil + geom_bar(aes(), stat="identity", position="fill") + 
  facet_grid(PlantYearTreat~.) +
  scale_fill_manual(values=Family_colors)+
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3))+
  labs(title="Bacteria",y= "Relative Abundance")+
  theme(plot.title = element_text(size = rel(1.5), face="bold"),
        axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        axis.text.x = element_text(hjust = 1),
        axis.title=element_text(size=10,face="bold"),
        legend.text=element_text(size = 8),
        legend.title = element_text(size=10),
        panel.grid = element_blank(), 
        panel.background = element_blank())
barplot.baccomp 
ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig_BactCompo.eps", barplot.baccomp , device = "eps", height = 12, width = 8, units = "in")


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

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Table_EnvFit_PlantYear.txt", x=envfit.soil.total[,1:4], sep="\t", quote=FALSE)

### Subset to p <0.05 and r2>0.4
envfit.soil.sub <- subset (envfit.soil.total , envfit.soil.total$pval<0.05)
envfit.soil.sub <- subset(envfit.soil.sub, envfit.soil.sub$r2>0.4)




### Set up some of the plotting specifics
Point_Sizes <- seq(from=2, to=6, length.out = 10)
Ax1.soil <- pcoa.16S$eig[1]/sum(pcoa.16S$eig)
Ax2.soil <- pcoa.16S$eig[2]/sum(pcoa.16S$eig)

### Plot PlantYear PCoA
FigAPlantYear <- ggplot(soil.points.collapsed, aes(x=Axis1, y=Axis2))+
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
FigAPlantYear
ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/FigA_PlantYearPCoA.eps", FigAPlantYear, device = "eps", height = 6, width = 6, units = "in")

### Make Plant treatment time PCoA
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
Fig2A
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

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/TableS3_swg2016EnvFit.txt", x=envfit.swg2016.total[,1:4], sep="\t", quote=FALSE)

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
Figswg2016
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
  geom_segment(data = envfit.swg2016.sub,
               aes(x = 0, xend = Axis1, y = 0, yend = Axis2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey")+
  geom_text(data=envfit.swg2016.sub, aes(label=Variable))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig2B
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig2B_swg2016Pcoa.eps", Fig2B, width = 6, height = 6, device = "eps", units="in")





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

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/TableS3_mis2016EnvFit.txt", x=envfit.mis2016.total[,1:4], sep="\t", quote=FALSE)

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

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig_mis2016PCoA.eps", Figmis2016, device = "eps", height = 6, width = 6, units = "in")





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
Fig2C 
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig2C_MisPcoa.eps", Fig2C, width = 6, height = 6, device = "eps", units="in")





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

write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/TableS3_swg2017EnvFit.txt", x=envfit.swg2017.total[,1:4], sep="\t", quote=FALSE)

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

ggsave(filename = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Fig_swg2017PCoA.eps", Figswg2017, device = "eps", height = 6, width = 6, units = "in")

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
Fig2D
ggsave("C:/Users/jaffyhu/Desktop/Soil/new/Fig2D_Swg16Pcoa.eps", Fig2D, width = 6, height = 6, device = "eps", units="in")



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
Fig2_Legend <- get_legend(Fig2A_leg)

setEPS()
postscript("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Figure2.eps", width=10, height=10,pointsize=10, paper="special")
plot_grid(plot_grid(Fig2A, Fig2B, Fig2C, Fig2D, align = "h"), Fig2_Legend, ncol=1)
dev.off()

###########################
#' Venn diagram 
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
#install limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)
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


#############################
#Occupancy abundance analysis ----------------------------------------
#############################

# *** Useful function --------------------------------------------------------
#' thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R
blank2na = function(x,na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    
    # the levels will be reset here
    x = factor(x)
    
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}

lastValue <- function(x) tail(x[!is.na(x)], 1)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}

#make new taxonomy dataframe, but add OTUs as new column
keep_otus <- tibble::rownames_to_column(taxonomy_filtered, "otu")

#in taxonomy table, add column naming the highest resolution taxonomy achieved for each OTU
rownames(keep_otus) <- keep_otus$otu
keep_otus[] = lapply(keep_otus, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
last_taxons<- apply(keep_otus, 1, lastValue)
keep_otus$last_taxon <- last_taxons
keep_otus$final_names <- paste(keep_otus$last_taxon, keep_otus$otu, sep=' - ')



#Core analysis 
#First, convert rarefied OTU table abundances to relative abundances
rel_otu_rare <- decostand(otu_rare, method="total", MARGIN=2)


#Using relabund table (above), give each OTU in each sample a row.
#Add relative abundance, metadata, and taxonomy
selected_otus <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% gather(sequence_name, abun, -otu) %>% 
  left_join(map_16S[, c('sequence_name','rep','time_numeric', 'treatment' ,'source', 'plant', 'month', 'sampling_date', 'year')], by = 'sequence_name') %>%
  left_join(keep_otus, by='otu')



#miscanthus 2016 Nitrogen Free
#Subset rarefied OTU table to only miscanthus 2016 Nitrogen Free
mis16nf_otu <- otu_rare[,map_16S$year=="2016"&map_16S$treatment=="nitrogen free"&map_16S$plant=="miscanthus"]

#subset mapping data to only miscanthus 2016 Standard Nitrogen free
map_mis16nf <- map_16S %>%
  filter(year == '2016'&treatment =='nitrogen free'&plant == 'miscanthus')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:10)
result <- NULL
source='miscanthus 2016nf'

for(i in weeks) {
  if(i %in% unique(map_mis16nf$sampling_week)) {
    name_sample <- map_mis16nf$sequence_name[map_mis16nf$sampling_week == i]
    timed_matrix <- mis16nf_otu[,colnames(mis16nf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s16 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s16$otu))
view(tempCore)
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_mis16nf.txt", tempCore, quote=FALSE)


#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied OTU table to PA matrix 
mis16nf_otu <- mis16nf_otu[rowSums(mis16nf_otu)>0,]
mis16nf_otu_PA <- 1*((mis16nf_otu>0)==1)
Occ_mis16nf <- rowSums(mis16nf_otu_PA)/ncol(mis16nf_otu_PA)
mis16nf_otu_rel <- decostand(mis16nf_otu, method="total", MARGIN=2)
Mean_abund_mis16nf <- apply(mis16nf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
mis16nf_df_occ <- data.frame(otu=names(Occ_mis16nf), occ=Occ_mis16nf) 
#create dataframe for abundance data generated above.
mis16nf_df_abun <- data.frame(otu=names(Mean_abund_mis16nf), abun=log10(Mean_abund_mis16nf))

#combine the two
mis16nf_occ_abun <- left_join(mis16nf_df_abun, mis16nf_df_occ, by='otu')

#label all otus as 'NotCore'
mis16nf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
mis16nf_occ_abun$unique[(mis16nf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#FigureMisNF (occupancy abundance plot)----------------
mis16nf_occ_abun_plot <- ggplot(data=mis16nf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(mis16nf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(mis16nf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "A")
mis16nf_occ_abun_plot

#mis16nf_core_abundance

selected_otus_mis16nf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_mis16nf %>%
  filter(treatment == 'nitrogen free' & plant == 'miscanthus' & year == 2016) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_mis16nf %>%
  filter(treatment == 'nitrogen free' & plant == 'miscanthus' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

mis16nf_core <- tmp3[,c(1,15,12)]
mis16nf_core <- as.data.frame(mis16nf_core)
mis16nf_core_wide <- spread(mis16nf_core, key='sampling_date', value='z_score')
mis16nf_core_wide[is.na(mis16nf_core_wide)] <- 0
rownames(mis16nf_core_wide) <- mis16nf_core_wide$otu
mis16nf_core_wide$otu <-  NULL
set.seed(20)
clusters_mis16nf <- hclust(dist(mis16nf_core_wide),'complete')
memb_mis16nf <- cutree(clusters_mis16nf, k=7)
mis16nf_dend <- plot(clusters_mis16nf, main=NULL)
rect.hclust(clusters_mis16nf, k=7)

mis16nf_clusters <- data.frame(memb_mis16nf)
mis16nf_clusters$otu <-  rownames(mis16nf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
mis16nf_clusters$order <- seq(1:2892)
test1 <-as.data.frame(clusters_mis16nf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, mis16nf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 1,2,3,5,6,7,4

mis16nf_clusters$stage <- 'A'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==2] <- 'B'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==3] <- 'C'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==4] <- 'D'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==5] <- 'E'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==6] <- 'F'
mis16nf_clusters$stage[mis16nf_clusters$memb_mis16nf==7] <- 'G'

#Figure --------------------------------------
mis16nf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'miscanthus', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(mis16nf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_mis16nf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_mis16nf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_mis16nf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_mis16nf), group=memb_mis16nf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "A")
mis16nf_core_abundance

mis16nf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'miscanthus', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(mis16nf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "A")
mis16nf_core_taxonomy



#miscanthus 2016 Standard Fertilization
#Subset rarefied OTU table to only miscanthus 2016 Standard Fertilization
mis16sf_otu <- otu_rare[,map_16S$year=="2016"&map_16S$treatment=="standard fertilization"&map_16S$plant=="miscanthus"]

#subset mapping data to only miscanthus 2016 Standard Fertilization
map_mis16sf <- map_16S %>%
  filter(year == '2016'&treatment =='standard fertilization'&plant == 'miscanthus')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:10)
result <- NULL
source='miscanthus 2016sf'

for(i in weeks) {
  if(i %in% unique(map_mis16sf$sampling_week)) {
    name_sample <- map_mis16sf$sequence_name[map_mis16sf$sampling_week == i]
    timed_matrix <- mis16sf_otu[,colnames(mis16sf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s16 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s16$otu))
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_mis16sf.txt", tempCore, quote=FALSE)



#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied OTU table to PA matrix 
mis16sf_otu <- mis16sf_otu[rowSums(mis16sf_otu)>0,]
mis16sf_otu_PA <- 1*((mis16sf_otu>0)==1)
Occ_mis16sf <- rowSums(mis16sf_otu_PA)/ncol(mis16sf_otu_PA)
mis16sf_otu_rel <- decostand(mis16sf_otu, method="total", MARGIN=2)
Mean_abund_mis16sf <- apply(mis16sf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
mis16sf_df_occ <- data.frame(otu=names(Occ_mis16sf), occ=Occ_mis16sf) 
#create dataframe for abundance data generated above.
mis16sf_df_abun <- data.frame(otu=names(Mean_abund_mis16sf), abun=log10(Mean_abund_mis16sf))

#combine the two
mis16sf_occ_abun <- left_join(mis16sf_df_abun, mis16sf_df_occ, by='otu')

#label all otus as 'NotCore'
mis16sf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
mis16sf_occ_abun$unique[(mis16sf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#Figure6B (occupancy abundance plot)----------------
mis16sf_occ_abun_plot <- ggplot(data=mis16sf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(mis16sf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(mis16sf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "B")
mis16sf_occ_abun_plot


#mis16sf_core_abundance

selected_otus_mis16sf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_mis16sf %>%
  filter(treatment == 'standard fertilization' & plant == 'miscanthus' & year == 2016) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_mis16sf %>%
  filter(treatment == 'standard fertilization' & plant == 'miscanthus' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

mis16sf_core <- tmp3[,c(1,15,12)]
mis16sf_core <- as.data.frame(mis16sf_core)
mis16sf_core_wide <- spread(mis16sf_core, key='sampling_date', value='z_score')
mis16sf_core_wide[is.na(mis16sf_core_wide)] <- 0
rownames(mis16sf_core_wide) <- mis16sf_core_wide$otu
mis16sf_core_wide$otu <-  NULL
set.seed(20)
clusters_mis16sf <- hclust(dist(mis16sf_core_wide),'complete')
memb_mis16sf <- cutree(clusters_mis16sf, k=7)
mis16sf_dend <- plot(clusters_mis16sf, main=NULL)
rect.hclust(clusters_mis16sf, k=7)

mis16sf_clusters <- data.frame(memb_mis16sf)
mis16sf_clusters$otu <-  rownames(mis16sf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
mis16sf_clusters$order <- seq(1:2599)
test1 <-as.data.frame(clusters_mis16sf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, mis16sf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 4,1,3,6,2,7,5

mis16sf_clusters$stage <- 'A'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==2] <- 'B'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==3] <- 'C'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==4] <- 'D'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==5] <- 'E'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==6] <- 'F'
mis16sf_clusters$stage[mis16sf_clusters$memb_mis16sf==7] <- 'G'

#Figure --------------------------------------
mis16sf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'miscanthus', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(mis16sf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_mis16sf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_mis16sf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_mis16sf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_mis16sf), group=memb_mis16sf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "B")
mis16sf_core_abundance

mis16sf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'miscanthus', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(mis16sf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "B")
mis16sf_core_taxonomy




#switchgrass 2016 nitrogen free
#Subset rarefied OTU table to only switchgrass nitrogen free
swg16nf_otu <- otu_rare[,map_16S$year=="2016"&map_16S$treatment=="nitrogen free"&map_16S$plant=="switchgrass"]

#subset mapping data to only switchgrass 2016 nitrogen free
map_swg16nf <- map_16S %>%
  filter(year == '2016'&treatment =='nitrogen free'&plant == 'switchgrass')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:9)
result <- NULL
source='switchgrass 2016nf'

for(i in weeks) {
  if(i %in% unique(map_swg16nf$sampling_week)) {
    name_sample <- map_swg16nf$sequence_name[map_swg16nf$sampling_week == i]
    timed_matrix <- swg16nf_otu[,colnames(swg16nf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s16 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s16$otu))
view(tempCore)
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_swg16nf.txt", tempCore, quote=FALSE)



#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied OTU table to PA matrix 
swg16nf_otu <- swg16nf_otu[rowSums(swg16nf_otu)>0,]
swg16nf_otu_PA <- 1*((swg16nf_otu>0)==1)
Occ_swg16nf <- rowSums(swg16nf_otu_PA)/ncol(swg16nf_otu_PA)
swg16nf_otu_rel <- decostand(swg16nf_otu, method="total", MARGIN=2)
Mean_abund_swg16nf <- apply(swg16nf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
swg16nf_df_occ <- data.frame(otu=names(Occ_swg16nf), occ=Occ_swg16nf) 
#create dataframe for abundance data generated above.
swg16nf_df_abun <- data.frame(otu=names(Mean_abund_swg16nf), abun=log10(Mean_abund_swg16nf))

#combine the two
swg16nf_occ_abun <- left_join(swg16nf_df_abun, swg16nf_df_occ, by='otu')

#label all otus as 'NotCore'
swg16nf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
swg16nf_occ_abun$unique[(swg16nf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#Figure6C (occupancy abundance plot)----------------
swg16nf_occ_abun_plot <- ggplot(data=swg16nf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(swg16nf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(swg16nf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "C")
swg16nf_occ_abun_plot


#swg16nf_core_abundance

selected_otus_swg16nf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_swg16nf %>%
  filter(treatment == 'nitrogen free' & plant == 'switchgrass' & year == 2016) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_swg16nf %>%
  filter(treatment == 'nitrogen free' & plant == 'switchgrass' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

swg16nf_core <- tmp3[,c(1,15,12)]
swg16nf_core <- as.data.frame(swg16nf_core)
swg16nf_core_wide <- spread(swg16nf_core, key='sampling_date', value='z_score')
swg16nf_core_wide[is.na(swg16nf_core_wide)] <- 0
rownames(swg16nf_core_wide) <- swg16nf_core_wide$otu
swg16nf_core_wide$otu <-  NULL
set.seed(20)
clusters_swg16nf <- hclust(dist(swg16nf_core_wide),'complete')
memb_swg16nf <- cutree(clusters_swg16nf, k=7)
swg16nf_dend <- plot(clusters_swg16nf, main=NULL)
rect.hclust(clusters_swg16nf, k=7)

swg16nf_clusters <- data.frame(memb_swg16nf)
swg16nf_clusters$otu <-  rownames(swg16nf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
swg16nf_clusters$order <- seq(1:2928)
test1 <-as.data.frame(clusters_swg16nf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, swg16nf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 3,5,6,1,2,7,4

swg16nf_clusters$stage <- 'A'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==2] <- 'B'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==3] <- 'C'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==4] <- 'D'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==5] <- 'E'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==6] <- 'F'
swg16nf_clusters$stage[swg16nf_clusters$memb_swg16nf==7] <- 'G'

#Figure --------------------------------------
swg16nf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'switchgrass', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg16nf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_swg16nf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_swg16nf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_swg16nf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_swg16nf), group=memb_swg16nf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "C")
swg16nf_core_abundance

swg16nf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'switchgrass', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg16nf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "C")
swg16nf_core_taxonomy



#switchgrass 2016 Standard Fertilization
#Subset rarefied OTU table to only switchgrass 2016 Standard Fertilization
swg16sf_otu <- otu_rare[,map_16S$year=="2016"&map_16S$treatment=="standard fertilization"&map_16S$plant=="switchgrass"]

#subset mapping data to only switchgrass 2016 Standard Fertilization
map_swg16sf <- map_16S %>%
  filter(year == '2016'&treatment =='standard fertilization'&plant == 'switchgrass')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:9)
result <- NULL
source='switchgrass 2016sf'

for(i in weeks) {
  if(i %in% unique(map_swg16sf$sampling_week)) {
    name_sample <- map_swg16sf$sequence_name[map_swg16sf$sampling_week == i]
    timed_matrix <- swg16sf_otu[,colnames(swg16sf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s16 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s16$otu))
view(tempCore)
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_swg16sf.txt", tempCore, quote=FALSE)



#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied OTU table to PA matrix 
swg16sf_otu <- swg16sf_otu[rowSums(swg16sf_otu)>0,]
swg16sf_otu_PA <- 1*((swg16sf_otu>0)==1)
Occ_swg16sf <- rowSums(swg16sf_otu_PA)/ncol(swg16sf_otu_PA)
swg16sf_otu_rel <- decostand(swg16sf_otu, method="total", MARGIN=2)
Mean_abund_swg16sf <- apply(swg16sf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
swg16sf_df_occ <- data.frame(otu=names(Occ_swg16sf), occ=Occ_swg16sf) 
#create dataframe for abundance data generated above.
swg16sf_df_abun <- data.frame(otu=names(Mean_abund_swg16sf), abun=log10(Mean_abund_swg16sf))

#combine the two
swg16sf_occ_abun <- left_join(swg16sf_df_abun, swg16sf_df_occ, by='otu')

#label all otus as 'NotCore'
swg16sf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
swg16sf_occ_abun$unique[(swg16sf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#Figure (occupancy abundance plot)----------------
swg16sf_occ_abun_plot <- ggplot(data=swg16sf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(swg16sf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(swg16sf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "C")
swg16sf_occ_abun_plot


#swg16sf_core_abundance

selected_otus_swg16sf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_swg16sf %>%
  filter(treatment == 'standard fertilization' & plant == 'switchgrass' & year == 2016) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_swg16sf %>%
  filter(treatment == 'standard fertilization' & plant == 'switchgrass' & year == 2016) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

swg16sf_core <- tmp3[,c(1,15,12)]
swg16sf_core <- as.data.frame(swg16sf_core)
swg16sf_core_wide <- spread(swg16sf_core, key='sampling_date', value='z_score')
swg16sf_core_wide[is.na(swg16sf_core_wide)] <- 0
rownames(swg16sf_core_wide) <- swg16sf_core_wide$otu
swg16sf_core_wide$otu <-  NULL
set.seed(20)
clusters_swg16sf <- hclust(dist(swg16sf_core_wide),'complete')
memb_swg16sf <- cutree(clusters_swg16sf, k=7)
swg16sf_dend <- plot(clusters_swg16sf, main=NULL)
rect.hclust(clusters_swg16sf, k=7)

swg16sf_clusters <- data.frame(memb_swg16sf)
swg16sf_clusters$otu <-  rownames(swg16sf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
swg16sf_clusters$order <- seq(1:2837)
test1 <-as.data.frame(clusters_swg16sf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, swg16sf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 3,4,7,1,6,2,5

swg16sf_clusters$stage <- 'A'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==2] <- 'B'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==3] <- 'C'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==4] <- 'D'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==5] <- 'E'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==6] <- 'F'
swg16sf_clusters$stage[swg16sf_clusters$memb_swg16sf==7] <- 'G'

#Figure --------------------------------------
swg16sf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'switchgrass', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg16sf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_swg16sf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_swg16sf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_swg16sf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_swg16sf), group=memb_swg16sf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "D")
swg16sf_core_abundance

swg16sf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2016, plant == 'switchgrass', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg16sf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "D")
swg16sf_core_taxonomy



#Subset rarefied OTU table to only switchgrass 2017 nitrogen free
swg17nf_otu <- otu_rare[,map_16S$year=="2017"&map_16S$treatment=="nitrogen free"]

#subset mapping data to only switchgrass 2017
map_swg17nf <- map_16S %>%
  filter(year == '2017'&treatment =='nitrogen free')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:8)
result <- NULL
source='switchgrass 2017nf'

for(i in weeks) {
  if(i %in% unique(map_swg17nf$sampling_week)) {
    name_sample <- map_swg17nf$sequence_name[map_swg17nf$sampling_week == i]
    timed_matrix <- swg17nf_otu[,colnames(swg17nf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s17 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s17$otu))
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_swg17nf.txt", tempCore, quote=FALSE)



#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied leaf OTU table to PA matrix 
swg17nf_otu <- swg17nf_otu[rowSums(swg17nf_otu)>0,]
swg17nf_otu_PA <- 1*((swg17nf_otu>0)==1)
Occ_swg17nf <- rowSums(swg17nf_otu_PA)/ncol(swg17nf_otu_PA)
swg17nf_otu_rel <- decostand(swg17nf_otu, method="total", MARGIN=2)
Mean_abund_swg17nf <- apply(swg17nf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
swg17nf_df_occ <- data.frame(otu=names(Occ_swg17nf), occ=Occ_swg17nf) 
#create dataframe for abundance data generated above.
swg17nf_df_abun <- data.frame(otu=names(Mean_abund_swg17nf), abun=log10(Mean_abund_swg17nf))

#combine the two
swg17nf_occ_abun <- left_join(swg17nf_df_abun, swg17nf_df_occ, by='otu')

#label all otus as 'NotCore'
swg17nf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
swg17nf_occ_abun$unique[(swg17nf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#Figure (occupancy abundance plot)----------------
swg17nf_occ_abun_plot <- ggplot(data=swg17nf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(swg17nf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(swg17nf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "E")
swg17nf_occ_abun_plot


#swg17nf_core_abundance

selected_otus_swg17nf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_swg17nf %>%
  filter(treatment == 'nitrogen free' & plant == 'switchgrass' & year == 2017) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_swg17nf %>%
  filter(treatment == 'nitrogen free' & plant == 'switchgrass' & year == 2017) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

swg17nf_core <- tmp3[,c(1,15,12)]
swg17nf_core <- as.data.frame(swg17nf_core)
swg17nf_core_wide <- spread(swg17nf_core, key='sampling_date', value='z_score')
swg17nf_core_wide[is.na(swg17nf_core_wide)] <- 0
rownames(swg17nf_core_wide) <- swg17nf_core_wide$otu
swg17nf_core_wide$otu <-  NULL
set.seed(20)
clusters_swg17nf <- hclust(dist(swg17nf_core_wide),'complete')
memb_swg17nf <- cutree(clusters_swg17nf, k=7)
swg17nf_dend <- plot(clusters_swg17nf, main=NULL)
rect.hclust(clusters_swg17nf, k=7)

swg17nf_clusters <- data.frame(memb_swg17nf)
swg17nf_clusters$otu <-  rownames(swg17nf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
swg17nf_clusters$order <- seq(1:2357)
test1 <-as.data.frame(clusters_swg17nf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, swg17nf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 1,7,6,5,4,3,2

swg17nf_clusters$stage <- 'A'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==2] <- 'B'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==3] <- 'C'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==4] <- 'D'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==5] <- 'E'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==6] <- 'F'
swg17nf_clusters$stage[swg17nf_clusters$memb_swg17nf==7] <- 'G'

#Figure --------------------------------------
swg17nf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2017, plant == 'switchgrass', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg17nf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_swg17nf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_swg17nf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_swg17nf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_swg17nf), group=memb_swg17nf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "E")
swg17nf_core_abundance

swg17nf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2017, plant == 'switchgrass', treatment == 'nitrogen free') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg17nf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "E")
swg17nf_core_taxonomy




#switchgrass 2017 Standard Fertilization
#Subset rarefied OTU table to only switchgrass 2017 Standard Fertilization
swg17sf_otu <- otu_rare[,map_16S$year=="2017"&map_16S$treatment=="standard fertilization"]

#subset mapping data to only switchgrass 2017
map_swg17sf <- map_16S %>%
  filter(year == '2017'&treatment =='standard fertilization')

#generate abundance/occupancy data for each otu for each week
weeks <- c(0:8)
result <- NULL
source='switchgrass 2017sf'

for(i in weeks) {
  if(i %in% unique(map_swg17sf$sampling_week)) {
    name_sample <- map_swg17sf$sequence_name[map_swg17sf$sampling_week == i]
    timed_matrix <- swg17sf_otu[,colnames(swg17sf_otu) %in% name_sample]
    timed_matrix <- timed_matrix[rowSums(timed_matrix)>0,]
    timed_matrix_PA <- 1*((timed_matrix>0)==1)
    timed_matrix_PA <- timed_matrix_PA[rowSums(timed_matrix_PA)>0,]
    Occ <- rowSums(timed_matrix_PA)/ncol(timed_matrix_PA)
    rel_abun <- decostand(timed_matrix, method="total", MARGIN=2)
    Mean_rel_abund <- apply(rel_abun, 1, mean)
    df_o <- data.frame(otu=names(Occ), occ=Occ) 
    df_a <- data.frame(otu=names(Mean_rel_abund), abun=log10(Mean_rel_abund))
    table <- left_join(df_a, df_o, by='otu') %>% mutate(week = i, source = source)
    result <- rbind(result, table)
  }
  else {
  }
}

#determine otus at occupancy of 1 for any week)
occ1_s17 <- result[result$occ==1,]
tempCore <- as.character(unique(occ1_s17$otu))
write.table(file = "C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/core_swg17sf.txt", tempCore, quote=FALSE)



#The above was just to determine which OTUs were at occupancy of 1 for a given week.
#now prepare the data for the actual occupancy/abundance plot.

#convert rarefied leaf OTU table to PA matrix 
swg17sf_otu <- swg17sf_otu[rowSums(swg17sf_otu)>0,]
swg17sf_otu_PA <- 1*((swg17sf_otu>0)==1)
Occ_swg17sf <- rowSums(swg17sf_otu_PA)/ncol(swg17sf_otu_PA)
swg17sf_otu_rel <- decostand(swg17sf_otu, method="total", MARGIN=2)
Mean_abund_swg17sf <- apply(swg17sf_otu_rel, 1, mean)

#create dataframe for occupancy data generated above.
swg17sf_df_occ <- data.frame(otu=names(Occ_swg17sf), occ=Occ_swg17sf) 
#create dataframe for abundance data generated above.
swg17sf_df_abun <- data.frame(otu=names(Mean_abund_swg17sf), abun=log10(Mean_abund_swg17sf))

#combine the two
swg17sf_occ_abun <- left_join(swg17sf_df_abun, swg17sf_df_occ, by='otu')

#label all otus as 'NotCore'
swg17sf_occ_abun$unique <- 'NotCore'
#change label  to 'core' for otus that were occ=1 at any time point.
swg17sf_occ_abun$unique[(swg17sf_occ_abun$otu %in% tempCore)] <- 'Core (occ=1 at any time pt)'

#Figure6F (occupancy abundance plot)----------------
swg17sf_occ_abun_plot <- ggplot(data=swg17sf_occ_abun, aes(x=abun, y=occ, fill=unique)) +
  theme_bw()+
  geom_point(size=3, pch=21, alpha=.8) +
  scale_fill_manual(breaks=unique, values=c('green4','white')) +
  labs(x=paste('log(mean relative abundance per OTU)\n (n=',nrow(swg17sf_occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(swg17sf_otu_PA),')', sep=''), fill=NULL) +
  theme(legend.position = 'none',
        legend.background = element_rect(fill=alpha(0.1)),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme(axis.title.x= element_text(size = 18))+
  labs(tag = "F")
swg17sf_occ_abun_plot



#swg17sf_core_abundance

selected_otus_swg17sf <- selected_otus[selected_otus$otu %in% tempCore,]
selected_otus_swg17sf %>%
  filter(treatment == 'standard fertilization' & plant == 'switchgrass' & year == 2017) %>%
  group_by(sampling_date, final_names, Phylum, Class, Order, Family, Genus) %>%
  dplyr::summarise(n=sum(abun>0)/length(abun),
                   all=length(abun),
                   rep_ab=mean(abun),
                   sd_rep=sd(abun)
  ) %>%
  filter(n>0) -> temp

#df with stats for the whole dataset per OTU
selected_otus_swg17sf %>%
  filter(treatment == 'standard fertilization' & plant == 'switchgrass' & year == 2017) %>%
  group_by(final_names, otu) %>%
  dplyr::summarise(
    all_ab=mean(abun),
    all_sd=sd(abun)
  ) -> temp2

#combining df and calculating the z-score
z_df <- left_join(temp, temp2)
z_df %>% arrange(otu) %>%
  mutate(
    z_score=(rep_ab-all_ab)/all_sd
  ) %>%
  arrange(Class, Family) -> tmp3

swg17sf_core <- tmp3[,c(1,15,12)]
swg17sf_core <- as.data.frame(swg17sf_core)
swg17sf_core_wide <- spread(swg17sf_core, key='sampling_date', value='z_score')
swg17sf_core_wide[is.na(swg17sf_core_wide)] <- 0
rownames(swg17sf_core_wide) <- swg17sf_core_wide$otu
swg17sf_core_wide$otu <-  NULL
set.seed(20)
clusters_swg17sf <- hclust(dist(swg17sf_core_wide),'complete')
memb_swg17sf <- cutree(clusters_swg17sf, k=7)
swg17sf_dend <- plot(clusters_swg17sf, main=NULL)
rect.hclust(clusters_swg17sf, k=7)

swg17sf_clusters <- data.frame(memb_swg17sf)
swg17sf_clusters$otu <-  rownames(swg17sf_clusters)
#determine the order of the clusters in the dendrogram so I can label the dendrogram.
swg17sf_clusters$order <- seq(1:2644)
test1 <-as.data.frame(clusters_swg17sf$order)
colnames(test1)<-c('orderinTree')
testing<-merge(test1, swg17sf_clusters, sort=FALSE, by.x="orderinTree",by.y="order", all=TRUE)
#left to right: 2,7,3,1,5,6,4

swg17sf_clusters$stage <- 'A'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==2] <- 'B'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==3] <- 'C'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==4] <- 'D'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==5] <- 'E'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==6] <- 'F'
swg17sf_clusters$stage[swg17sf_clusters$memb_swg17sf==7] <- 'G'

#Figure --------------------------------------
swg17sf_core_abundance <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare)  %>%
  gather(sequence_name, abun, -otu) %>%
  left_join(map_16S[, c('sequence_name','source', 'plant', 'sampling_date', 'year', 'sampling_week','rep', 'treatment')], by = 'sequence_name') %>%
  filter(year == 2017, plant == 'switchgrass', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg17sf_clusters, by='otu') %>%
  group_by(sampling_date) %>%
  mutate(sample_size = length(unique(sequence_name))) %>%
  group_by(sampling_week, sampling_date, memb_swg17sf, stage, sequence_name) %>%
  summarise(n_totalrelabun=sum(abun)) %>% ##get total abundance of genus for each sample/rep
  group_by(sampling_week, sampling_date, memb_swg17sf, stage) %>%
  summarise(n_relabun = mean(n_totalrelabun),
            n_sd = sd(n_totalrelabun)) %>% ##get mean, SD abundance of genus for each sample/rep
  filter(!is.na(memb_swg17sf)) %>%
  ggplot(aes(x=as.Date(sampling_date), y=n_relabun, color=as.factor(memb_swg17sf), group=memb_swg17sf))+
  geom_line(size=2,linetype = "dashed")+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=n_relabun-n_sd, ymax=n_relabun+n_sd), width=10, size=1.1,
                position=position_dodge(.1))+
  scale_color_manual(values = c('dodgerblue1','dodgerblue4','coral1','coral4','green2','green4',
                                'red', 'brown'))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.background = element_rect(fill = "transparent"),
                          legend.position = c(.4,.8), axis.text.x =element_blank(),
                          legend.key.size = unit(.8, "cm"), legend.title = element_text(size=12)) +
  labs(x=NULL, y=NULL, color='Core OTU Cluster')+
  theme(axis.title.y = element_text(size = 16),
        text = element_text(size=22))+
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='bottom')+
  
  ylab("Relative abundance")+
  ylim(-0.1, 0.5)+
  labs(tag = "F")
swg17sf_core_abundance

swg17sf_core_taxonomy <- data.frame(otu = as.factor(row.names(rel_otu_rare)), rel_otu_rare) %>% 
  gather(sequence_name, abun, -otu) %>%  
  left_join(map_16S[, c('sequence_name','treatment', 'plant', 'sampling_date', 'year', 'sampling_week','rep')], by = 'sequence_name') %>%
  filter(year == 2017, plant == 'switchgrass', treatment == 'standard fertilization') %>%
  filter(otu %in% tempCore) %>%
  left_join(swg17sf_clusters, by='otu') %>%
  left_join(tax_filtered, by='otu') %>% 
  group_by(Phylum, sampling_week)%>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sequence_name)),
            norm_ra=n_count/n_sample) %>%
  ggplot(aes(x=as.factor(sampling_week), y=norm_ra, color=Phylum, group=Phylum)) +
  geom_line(size=2) +
  theme_classic() + theme(strip.background = element_blank(), 
                          legend.position='none')+
  labs(x='Sampling week', y=NULL, color=NULL)+
  ylab("Relative abundance")+
  labs(tag = "F")
swg17sf_core_taxonomy


grid.arrange(ggarrange(mis16nf_occ_abun_plot,mis16sf_occ_abun_plot, 
                       swg16nf_occ_abun_plot,swg16sf_occ_abun_plot, 
                       swg17nf_occ_abun_plot, swg17sf_occ_abun_plot,
                       ncol=2, nrow=3))

grid.arrange(ggarrange(mis16nf_core_abundance,mis16sf_core_abundance, 
                       swg16nf_core_abundance,swg16sf_core_abundance, 
                       swg17nf_core_abundance, swg17sf_core_abundance,
                       ncol=2, nrow=3))

grid.arrange(ggarrange(mis16nf_core_taxonomy,mis16sf_core_taxonomy, 
                       swg16nf_core_taxonomy,swg16sf_core_taxonomy, 
                       swg17nf_core_taxonomy, swg17sf_core_taxonomy,
                       ncol=2, nrow=3))



###############################################################
######################### CORE MICROBIOTA #####################
###############################################################
#All samples together
# OCCUPANCY VS ABUNDANCE
# Occupancy
# 1. Bacteria
otu_PA <- 1*((otu_filtered>0)==1)
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Occ <- rowSums(otu_PA)/ncol(otu_PA)
class(Occ)
df.Occ <- as.data.frame(Occ)
head(df.Occ)
df.Occ=rownames_to_column(df.Occ, var = "OTU")
dim(df.Occ)

# Taxonomy
# 1. Bacteria

rownames(taxonomy_filtered) <- rownames(otu_filtered)
head(taxonomy_filtered)
dim(taxonomy_filtered)
taxonomy_filtered=rownames_to_column(taxonomy_filtered, var = "OTU")
dim(df.Occ)
head(df.Occ)
dim(taxonomy_filtered)
df.Occ.tax <- merge(df.Occ, taxonomy_filtered, by.x =c("OTU"), by.y = c("OTU"))
head(df.Occ.tax)
dim(df.Occ.tax) ### all OTU with occupancy and taxonomy!!!!!!!!!!!

# Cumulative Relative Abundance per OTU in all samples
# 1. Bacteria
otu_rel <- decostand(otu_filtered, method="total", MARGIN=2)
com_abund <- rowSums(otu_rel)
df.com_abund <- as.data.frame(com_abund)
head(df.com_abund)
df.com_abund$RelAbund=df.com_abund$com_abund/202
sum(df.com_abund$com_abund)
sum(df.com_abund$RelAbund)
df.com_abund$PercentRelAbund=df.com_abund$RelAbund*100
sum(df.com_abund$PercentRelAbund)
df.com_abund=rownames_to_column(df.com_abund, var = "OTU")
head(df.com_abund)
dim(df.com_abund) ### all OTU with CumulativeRelAbund, percent CumulativeRelAbund!!!!!!!!!!!


# Relative Abundance (or mean Relative Abundance) each OTU in each sample
# 1. Bacteria
otu <- otu[rowSums(otu_filtered)>0,]
otu.ID=rownames_to_column(otu,var = "OTU")
otu.melt=melt(otu.ID,variable.name = "Sample")
sum(otu.melt$value) #34044184
otu.melt$relabund = otu.melt$value/34044184
head(otu.melt)
sum(otu.melt$relabund)
otu.melt$percentrelabund=otu.melt$relabund*100
sum(otu.melt$percentrelabund)
head(otu.melt)
dim(otu.melt) ### all OTU with relabund, percent relabund for each sample!!!!!!!!!!!



# merge occupancy 1 and cumulative relative abundance 
# 1. Bacteria
Occ_RelAbund <- merge(df.Occ.tax, df.com_abund, by.x =c("OTU"), by.y = c("OTU"))
df.Occ.tax1 <- subset(df.Occ.tax, df.Occ.tax$Occ==1)
head(df.Occ.tax1)
dim(df.Occ.tax1)
Occ1_RelAbund <- merge(df.Occ.tax1, df.com_abund, by.x =c("OTU"), by.y = c("OTU"))
dim(Occ1_RelAbund)
head(Occ1_RelAbund)
sort_Occ1_RelAbund <- Occ1_RelAbund[order(Occ1_RelAbund$RelAbund, decreasing = TRUE),]
write.table(sort_Occ1_RelAbund, file = 'C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/sort_Occ1_RelAbund.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
sum(sort_Occ1_RelAbund$PercentRelAbund)


#class count= bacteria
class_count.bac <- sort_Occ1_RelAbund %>%
  group_by(Class,Phylum) %>%
  summarise(class_count=n())

sort_class_count <- class_count.bac[order(class_count.bac$class_count, decreasing = TRUE),]
library("readr")
write_csv(sort_class_count,"C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/sort_class_count.csv")
class.count=read.csv('C:/Users/jaffyhu/Desktop/Soil/new//GLBRC soil/output/sort_class_count.csv', header=TRUE)

# write.table(gc.count, file = 'gc.count.bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
# gc.count=read.table('gc.count.bac.txt', sep='\t', header=TRUE)
# gc.count$Genus_Class <- as.character(gc.count$Genus_Class)
# gc.count$Genus_Class <- factor(gc.count$Genus_Class, levels=unique(gc.count$Genus_Class))
class.count$Class <- as.character(class.count$Class)
class.count$Class <- factor(class.count$Class, levels=unique(class.count$Class))

Bac.taxa.num=ggplot(class.count, aes(x = Class, y = class_count, fill=Phylum))+ 
  geom_bar(position = "dodge",stat = "identity")+
  ylab("Number of taxa")+
  scale_y_continuous(expand = c(0,0.5))+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#7fffd4', '#8b8378', '#00008b', '#caff70', '#8b0a50', '#228b22', '#8b6914'))+
  theme_bw()+
  coord_flip()+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
Bac.taxa.num

#plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
#labs(x= "Genus (or Class)", y="Number of taxa")

par(mfrow=c(1,2))
# Plot Core Bacteria
# 1. Bacteria
dim(sort_Occ1_RelAbund)
top100.occ1 <- sort_Occ1_RelAbund[1:100,]
CoreBac <- ggplot(sort_Occ1_RelAbund,aes(x=fct_reorder(Class, RelAbund, .desc=T), y=PercentRelAbund, fill=Phylum))+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#7fffd4', '#8b8378', '#00008b', '#caff70', '#8b0a50', '#228b22', '#8b6914'))+
  labs(title = "Bacteria", y= "Relative Abundance (%)", x="Class")+
  theme_bw()+
  theme(plot.title = element_text(size=16, face="bold"),
        axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position = "right",
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
CoreBac

####### Plotting core microbiome members ##############
library(grid)
library(gtable)
grid.newpage()
gA <- ggplotGrob(CoreBac)
gB <- ggplotGrob(Bac.taxa.num)
grid::grid.newpage()
setEPS()
postscript("C:/Users/jaffyhu/Desktop/Soil/new//GLBRC soil/output/Core_micrbiome.eps", height = 7, width = 14)
grid::grid.draw(cbind(gA, gB))
dev.off()
graphics.off()


# Occupancy-Abundance Plot
# 1. Bacteria
color_top <- df.com_abund$RelAbund
color_top <- Occ
color_top[] <- 'black' 
Occ.1 <- Occ[Occ==1]
color_top[names(color_top) %in% names(Occ.1)] <- 'orange'
# Default plot in R
plot(log10(df.com_abund$RelAbund), Occ, col=color_top, pch=20, ylab='Occupancy', xlab='log(Mean of relative abundance)')
# Plot using ggplot2
occ_bacless1 <- subset(Occ_RelAbund, Occ != 1)
Occ.RelAbun.Bac <- ggplot()+
  geom_point(aes(x=log10(RelAbund), y=Occ), data=occ_bacless1, size=2.5, alpha=0.5, pch=21, colour="darkgrey")+
  labs(title="Bacteria",y= "Occupancy", x="Log10(mean of relative abundance)")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face="bold"),
        axis.text=element_text(size=9), 
        axis.title=element_text(size=15,face="bold"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size = 8, face="bold"),
        legend.title = element_blank())+
  geom_point(data = Occ1_RelAbund, aes(x=log10(RelAbund), y=Occ, colour=Phylum), size=2.5, alpha=0.5)+
  scale_colour_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#7fffd4', '#8b8378', '#00008b', '#caff70', '#8b0a50', '#228b22', '#8b6914'))
#b <- Occ.RelAbun.Bac+theme(legend.position = "none")
#ggarrange(b, f, ncol = 2, nrow = 1)
Occ.RelAbun.Bac

plot <- ggarrange(Occ.RelAbun.Bac, ncol = 1, nrow = 1, align = "hv")

ggsave("C:/Users/jaffyhu/Desktop/Soil/new/GLBRC soil/output/Occ_abunWhole.tiff", Occ.RelAbun.Bac, device = "tiff",
       width = 10, height = 7, 
       units= "in", dpi = 600)

#################################################################
### Look at betadispersion of crops through time and treatment###
#################################################################
#Switchgrass 2016 nitrogen free
otu.swg16nf <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$treatment=="nitrogen free"&map_16S$Year==2016]
map.swg16nf <- map_16S[map_16S$plant=="switchgrass"&map_16S$treatment=="nitrogen free"&map_16S$Year==2016,]

dist.swg16nf <- vegdist(t(otu.swg16nf), method="bray")

Dispersion.swg16nf <- betadisper(dist.swg16nf, map.swg16nf$sampling_week)


names(Dispersion.swg16nf$distances)==map.swg16nf$sequence_name
Dispersion.swg16nf.df <- data.frame(Distance_to_Median=Dispersion.swg16nf$distances, Date=map.swg16nf$sampling_week)

Disp.swg16nf <- ggplot(Dispersion.swg16nf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("Switchgrass 2016 nitrogen free ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.swg16nf
ggsave("C:/Users/jaffyhu/Desktop/soil/new/GLBRC soil/output/swg16nf_BetaDispersion.eps", Disp.swg16nf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.swg16nf.df, Distance_to_Median~Date))


#Switchgrass 2016 standard fertilization
otu.swg16sf <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$treatment=="standard fertilization"&map_16S$Year==2016]
map.swg16sf <- map_16S[map_16S$plant=="switchgrass"&map_16S$treatment=="standard fertilization"&map_16S$Year==2016,]

dist.swg16sf <- vegdist(t(otu.swg16sf), method="bray")

Dispersion.swg16sf <- betadisper(dist.swg16sf, map.swg16sf$sampling_week)


names(Dispersion.swg16sf$distances)==map.swg16sf$sequence_name
Dispersion.swg16sf.df <- data.frame(Distance_to_Median=Dispersion.swg16sf$distances, Date=map.swg16sf$sampling_week)

Disp.swg16sf <- ggplot(Dispersion.swg16sf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("Switchgrass 2016 standard fertilization ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.swg16sf
ggsave("C:/Users/jaffyhu/Desktop/soil/new//GLBRC soil/output/swg16sf_BetaDispersion.eps", Disp.swg16sf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.swg16sf.df, Distance_to_Median~Date))

#Switchgrass 2017 nitrogen free
otu.swg17nf <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$treatment=="nitrogen free"&map_16S$Year==2017]
map.swg17nf <- map_16S[map_16S$plant=="switchgrass"&map_16S$treatment=="nitrogen free"&map_16S$Year==2017,]

dist.swg17nf <- vegdist(t(otu.swg17nf), method="bray")

Dispersion.swg17nf <- betadisper(dist.swg17nf, map.swg17nf$sampling_week)


names(Dispersion.swg17nf$distances)==map.swg17nf$sequence_name
Dispersion.swg17nf.df <- data.frame(Distance_to_Median=Dispersion.swg17nf$distances, Date=map.swg17nf$sampling_week)

Disp.swg17nf <- ggplot(Dispersion.swg17nf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("Switchgrass 2017 nitrogen free ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.swg17nf
ggsave("C:/Users/jaffyhu/Desktop/soil/new/GLBRC soil/output/swg17nf_BetaDispersion.eps", Disp.swg17nf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.swg17nf.df, Distance_to_Median~Date))


#Switchgrass 2017 standard fertilization
otu.swg17sf <- otu_rare[,map_16S$plant=="switchgrass"&map_16S$treatment=="standard fertilization"&map_16S$Year==2017]
map.swg17sf <- map_16S[map_16S$plant=="switchgrass"&map_16S$treatment=="standard fertilization"&map_16S$Year==2017,]

dist.swg17sf <- vegdist(t(otu.swg17sf), method="bray")

Dispersion.swg17sf <- betadisper(dist.swg17sf, map.swg17sf$sampling_week)


names(Dispersion.swg17sf$distances)==map.swg17sf$sequence_name
Dispersion.swg17sf.df <- data.frame(Distance_to_Median=Dispersion.swg17sf$distances, Date=map.swg17sf$sampling_week)

Disp.swg17sf <- ggplot(Dispersion.swg17sf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("Switchgrass 2016 standard fertilization ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.swg17sf
ggsave("C:/Users/jaffyhu/Desktop/soil/new/GLBRC soil/output/swg17sf_BetaDispersion.eps", Disp.swg17sf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.swg17sf.df, Distance_to_Median~Date))


#Miscanthus 2016 nitrogen free
otu.mis16nf <- otu_rare[,map_16S$plant=="miscanthus"&map_16S$treatment=="nitrogen free"&map_16S$Year==2016]
map.mis16nf <- map_16S[map_16S$plant=="miscanthus"&map_16S$treatment=="nitrogen free"&map_16S$Year==2016,]

dist.mis16nf <- vegdist(t(otu.mis16nf), method="bray")

Dispersion.mis16nf <- betadisper(dist.mis16nf, map.mis16nf$sampling_week)


names(Dispersion.mis16nf$distances)==map.mis16nf$sequence_name
Dispersion.mis16nf.df <- data.frame(Distance_to_Median=Dispersion.mis16nf$distances, Date=map.mis16nf$sampling_week)

Disp.mis16nf <- ggplot(Dispersion.mis16nf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("miscanthus 2016 nitrogen free ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.mis16nf
ggsave("C:/Users/jaffyhu/Desktop/soil/new/GLBRC soil/output/mis16nf_BetaDispersion.eps", Disp.mis16nf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.mis16nf.df, Distance_to_Median~Date))


#miscanthus 2016 standard fertilization
otu.mis16sf <- otu_rare[,map_16S$plant=="miscanthus"&map_16S$treatment=="standard fertilization"&map_16S$Year==2016]
map.mis16sf <- map_16S[map_16S$plant=="miscanthus"&map_16S$treatment=="standard fertilization"&map_16S$Year==2016,]

dist.mis16sf <- vegdist(t(otu.mis16sf), method="bray")

Dispersion.mis16sf <- betadisper(dist.mis16sf, map.mis16sf$sampling_week)


names(Dispersion.mis16sf$distances)==map.mis16sf$sequence_name
Dispersion.mis16sf.df <- data.frame(Distance_to_Median=Dispersion.mis16sf$distances, Date=map.mis16sf$sampling_week)

Disp.mis16sf <- ggplot(Dispersion.mis16sf.df, aes(x=Date, y=Distance_to_Median, group=Date))+
  geom_boxplot()+
  geom_point()+
  ggtitle("miscanthus 2016 standard fertilization ")+
  theme(axis.text.x = element_text(angle = 60))+
  ylim(c(0,0.5))
Disp.mis16sf
ggsave("C:/Users/jaffyhu/Desktop/soil/new/GLBRC soil/output/mis16sf_BetaDispersion.eps", Disp.mis16sf, device = "eps", width = 4, height = 4, units = "in")

TukeyHSD(aov(data = Dispersion.mis16sf.df, Distance_to_Median~Date))

grid.arrange(ggarrange(Disp.mis16nf,Disp.mis16sf, 
                       Disp.swg16nf,Disp.swg16sf, 
                       Disp.swg17nf, Disp.swg17sf,
                       ncol=2, nrow=3))





# Data preparation for MENA analysis  ------
# Before continuing, you must run "RecreateGrady2019_OTUtable.R". It produces 3 files:
#1)otu_rare_16s.txt: the OTU table from Grady et al 2019 (rarefied to 1000 reads)
#2)core_16s.txt: the core phyllosphere OTUs from Grady et al. 2019
#3)tax_16s.txt: the taxonomy of the core phyllosphere OTUs from Grady et al. 2019


#16S: Upload the rarefied OTU table from Grady et al. 2019 code. 
otu_rare_16s <- read.table("otu_rare_16s.txt",header=TRUE, sep='\t',row.names=1)
#convert to rel abund table
#otu_rare_16s <- decostand(otu_rare_16s, method="total", MARGIN=2)


#upload switchgrass 2017 core OTU names and format as vector
core_16s <-read.table("core_16s.txt",sep='\t')
core_16s <- as.vector(core_16s$V1)
#subset 16s OTU table to only the samples that are also in the fungal data.
otu_rare_16s<- otu_rare_16s[,colnames(otu_rare_16s)%in%colnames(rare_core_fungal)] 

#subset the rarefied 16s OTU table for 16s core leaf otus
rare_core_16s <- as.data.frame(subset(otu_rare_16s, rownames(otu_rare_16s) %in% core_16s))
#calculate minimum number of samples in which each OTU appears.
bactarch_otu_PA <- as.data.frame(1*((rare_core_16s>0)==1))
bactarch_otu_PA$rowsum <- rowSums(bactarch_otu_PA)
min(bactarch_otu_PA$rowsum)
nrow(bactarch_otu_PA[bactarch_otu_PA$rowsum >= 12,])

#subset rarefied fungal phyllosphere table for only samples also in the bactarch data.
rare_core_fungal<- rare_core_fungal[,colnames(rare_core_fungal)%in%colnames(rare_core_16s)] 
#calculate number of samples in which each OTU appears.
fungal_otu_PA <- as.data.frame(1*((rare_core_fungal>0)==1))
fungal_otu_PA$rowsum <- rowSums(fungal_otu_PA)
min(fungal_otu_PA$rowsum)
nrow(fungal_otu_PA[fungal_otu_PA$rowsum >= 12,])



# Order the samples in the fungal phyllosphere core otu table
rare_core_fungal <- rare_core_fungal[,order(colnames(rare_core_fungal))]
# Order the 16s phyllosphere core otu table the same way
rare_core_16s <- rare_core_16s[,order(colnames(rare_core_16s))]
# check that samples are in the same order in both tables
colnames(rare_core_fungal) == colnames(rare_core_16s)
#bind the otu tables
rare_core_all <- rbind(rare_core_fungal, rare_core_16s)




#replace zeros with an empty space (for input to MENA analysis)
#zeroes meaning not detectable.
rare_core_all[rare_core_all == 0] <- ''

dim(rare_core_fungal)
dim(rare_core_16s)




#save as .txt file.
#write.table(rare_core_all, file="rare_core_all.txt", sep = '\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
#manually add a tab as the first character of that file. then conduct
#Molecular Ecological Network Analysis as described in manuscript.


#######
## Calculate various results from MENA analysis
## Create data needed for Cytoscape to produce Figure 5 (network)

###Note to open MENA's output Network file in Cytoscape, 
#you must change the file extension to .SIF rather than .TXT
#For cytoscape, you also need the 'nodes_withtax.csv' created below.

#######

#upload initial network nodes from MENA and add metadata.
initialnodes <- read.table("MENA_Node_attrib_file.txt", sep='\t',header=TRUE)
initialnodes <- subset(initialnodes, select=c(Name, node.degree))

#get list of fungal OTUs, and include Fungi as kingdom name
fungal_tax <- subset(keep_otus, select=c(Kingdom, otu))
#get list of 16s OTUs, and include Kingdom name (all are bacteria, no archaea)
tax_16s <- read.table("tax_16s.txt", header=TRUE, sep='\t')
bactarch_tax <- subset(tax_16s, select=c(Kingdom, otu))




nodes_withtax <- merge(initialnodes,bactarch_tax, by.x='Name', by.y='otu',all.x=TRUE)
nodes_withtax <- merge(nodes_withtax, fungal_tax, by.x='Name', by.y='otu',all.x=TRUE)
#combine the two columns
nodes_withtax$Kingdom.x <-as.character(nodes_withtax$Kingdom.x)
nodes_withtax$Kingdom.x <- replace_na(nodes_withtax$Kingdom.x, rep('Fungi'))

#clean up the table
nodes_withtax <- subset(nodes_withtax, select=-c(Kingdom.y))
nodes_withtax$Kingdom.x <- gsub(pattern = "k:", replacement = "", x = nodes_withtax$Kingdom.x )

colnames(nodes_withtax)[colnames(nodes_withtax) == 'Kingdom.x'] <- 'Kingdom'

#Write table to text file: will be used by Cytoscape to prepare Figure 5.
#write.csv(nodes_withtax, file="nodes_withtax.csv", quote=FALSE, sep=',', row.names=FALSE)


###find taxonomy of the most highly connected fungal OTUs in the microbial network
#first get the fungal OTUs (124/227 OTUs) and order by number of connections
fung_otus<-arrange(subset(nodes_withtax, nodes_withtax$Kingdom=='Fungi'), desc(node.degree))
#get taxonomy of the most highly-connected fungal OTUs.
subset(keep_otus, keep_otus$otu==fung_otus[1,1])
subset(keep_otus, keep_otus$otu==fung_otus[2,1])
subset(keep_otus, keep_otus$otu==fung_otus[3,1])
subset(keep_otus, keep_otus$otu==fung_otus[4,1])


###find taxonomy of the most highly connected bacterial OTU in the microbial network
#first get bacterial OTUs (23/227 OTUs) and order by number of connections
bac_otus<-arrange(subset(nodes_withtax, nodes_withtax$Kingdom=='Bacteria'), desc(node.degree))
#get taxonomy of the most highly-connected bacterial OTUs.
subset(tax_16s, as.character(tax_16s$otu)==bac_otus[1,1])

subset(tax_16s, as.character(tax_16s$otu)==bac_otus[2,1])

subset(tax_16s, as.character(tax_16s$otu)==bac_otus[3,1])


#Read in the MENA edge file. Will inform on number and type of connections between taxa.
mena_int<- read.csv('MENA_Edge_attrib_file.txt', sep='\t', header=F)

#Get list of network nodes (OTUs)
coreOTUs <- nodes_withtax$Name

#parse apart the Mena Int file into a dataframe.
interTAXA <- mena_int %>% separate(V1, c('start','connect', 'end','value'), sep='([\\(\\)\\=])', remove=T) %>%
  
  mutate(start=str_trim(start, side = "both"),
         end=str_trim(end, side = "both"),
         start.memb=if_else(start %in% coreOTUs, 1,0),
         end.memb=if_else(end %in% coreOTUs, 1, 0),
         core.int=start.memb+end.memb,
         
         domain.start= if_else(start %in% tax_16s$otu, '16S', 'ITS'),
         domain.end= if_else(end %in% tax_16s$otu, '16S', 'ITS'),
         fungus=if_else(start %in% fungal_tax$otu, 1,0),
         fungus=if_else(end %in% fungal_tax$otu,1,fungus),
         bacteria= if_else(start %in% tax_16s$otu, 1, 0),
         bacteria= if_else(end %in% tax_16s$otu, 1, bacteria),
         bact_fung=if_else(bacteria+fungus==2, 1, 0))

#number of fungal OTUs in the network
nrow(nodes_withtax[nodes_withtax$Kingdom =="Fungi",])
#number of bacterial OTUs in the network
nrow(nodes_withtax[nodes_withtax$Kingdom =="Bacteria",])

#total number of bacteria-fungi interactions
sum(interTAXA$bact_fung)
#or
nrow(interTAXA[interTAXA$bact_fung == "1",])
#positive bacteria-fungi interactions
nrow(interTAXA[interTAXA$bact_fung == "1" & interTAXA$connect=="pp",])
#negative bacteria-fungi interactions
nrow(interTAXA[interTAXA$bact_fung == "1" & interTAXA$connect=="np",])

#total number of fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0",])
#positive fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0" & interTAXA$connect=="pp",])
#negative fungi-fungi interactions
nrow(interTAXA[interTAXA$fungus == "1" & interTAXA$bacteria=="0" & interTAXA$connect=="np",])


#total number of bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1",])
#positive bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1" & interTAXA$connect=="pp",])
#negative bacteria-bacteria interactions
nrow(interTAXA[interTAXA$fungus == "0" & interTAXA$bacteria=="1" & interTAXA$connect=="np",])




##get list of bacterial OTUs that were associated with fungi
intbactOTUs <-interTAXA[interTAXA$bact_fung == "1",]
##subset bacterial taxonomy for those OTUs
intbactOTUtax <- tax_16s[tax_16s$otu %in% intbactOTUs$end,]
#inspect the above table to determine common bacterial taxa


##get list of fungal OTUs that were associated with bacteria
intfungiOTUs <-interTAXA[interTAXA$bact_fung == "1",]
##subset bacterial taxonomy for those OTUs
intfungiOTUtax <- keep_otus[keep_otus$otu %in% intfungiOTUs$start,]
#inspect the above table to determine common fungal taxa



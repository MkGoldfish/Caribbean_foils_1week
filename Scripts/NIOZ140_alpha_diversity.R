##%######################################################%##
#                                                          #
####              Exploratory Analysis of               ####
####            16S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####           Project: foils statia NIOZ 140           ####
####            Adapted by Maaike April 2022            ####
###%#                                                  #%###
####               Corrected Primers                    ####
#                                                          #
##%######################################################%##
##%########################################################%##
##  version.string R version 4.2.1 (2022-06-23 ucrt)
##  nickname       Funny-Looking Kid   


###%####################################################%###
#                                                          #
####                      Alpha Diversity              ####
#                                                          #
##%######################################################%##


###_______________________________________________________________________________________#%###
####                 Working Directory                 ####
####_______________________________________________________________________________________#%###
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")


####_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
####_______________________________________________________________________________________#%###
library(devtools)
library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(compositions)
library(gghalves)
library(ggrepel)
library(microbiome)
library("remotes")
library("ggh4x")
library(ggpubr)
library("plotrix")
library("FactoMineR")
library("factoextra")
library(patchwork)
library(usedist)
library("ggsci")
library(rvg)
library(officer)
library("cowplot")

####_______________________________________________________________________________________#%###
####                  Import Data &                    ####
####_______________________________________________________________________________________#%###
tax <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/asv_table.txt', row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- read.delim('../../../Data_corrected_barcodes/metadata/mapping_file_details.txt', row.names = 1, na.strings = c(""))

map <- map %>% mutate(timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6" ,
  )
)

map <- sample_data(map)

#Import earlier created tidy data tibble with pre-processed absolute abundance values
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c("")) 
head(t2)

#Make sure al T3 and negative controls are goooone
t2$revComp <- NULL
t2 <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc") 
t2 <- t2 %>% mutate( timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6" ,
)
)

set.seed(42)
####_______________________________________________________________________________________#%###
####           Create & check physeq object            ####
####_______________________________________________________________________________________#%###
physeq_object = merge_phyloseq(otu, tax, map) 

summarize_phyloseq(physeq_object)
ntaxa(physeq_object)  ##13169
nsamples(physeq_object) ##63?
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #it says 0 sequences here

physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90") #eliminate the sample with 0seq
max(sample_sums(physeq_object))
min(sample_sums(physeq_object)) #534

physeq_object_MinSamp90 = physeq_object

#####subset T3 & merge ####         ### No longer in map file
sub1<- subset_samples(physeq_object, timepoint %in% c("Day 1", "Day 6"))
sub2 <- subset_samples(physeq_object, surface!="negative_c") 
physeq_object <- merge_phyloseq(sub1, sub2)

summarize_phyloseq(physeq_object)
nsamples(physeq_object)
sample_names(physeq_object)

physeq_object_Min90T3 = physeq_object

############################################################################################%###
###_______________________________________________________________________________________#%###
####  Colors!!####
####_______________________________________________________________________________________#%###

####_Colors!!______________________________________________________________________________####
# Different colorpallettes to choose from
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

HCB.col <- c("#EE7733")
PDB.col <- c("#0077BB")

############################################################################################%###
####_______________________________________________________________________________________#%###
####           Alpha Diversity          ####
####_______________________________________________________________________________________#%###

####__Obtain alpha diversity indices with phyloseq package__________________________######
####________________________________________________________________________________###%###

# Phyloseq, Microbiome and Vegan give same values for different indices
# Both Microbiome and Phyloseq are based on Vegan? Vegan has less options
# Phyloseq and Microbiome can work with phyloseq objects so seem easiest
# Microbiome::alpha has many measures we don't need and is slow

####___Alpha diversity with phyloseq________________________________________________________________####
# Since there is many 0's in this data, not all indices work, so we can not just do index=all
# Pick 4 of the most used measures
Alpha.Foils <- phyloseq::estimate_richness(physeq_object, measures=c("Observed", "Simpson", "Shannon", "Chao1") )
head(Alpha.Foils)

write_csv(Alpha.Foils , file = "../alpha_div/alpha_div_no_prune.csv")
Alpha.Foils  <- Alpha.Foils  %>% cbind(data.frame(physeq_object@sam_data)) 
Alpha.Foils <- filter(Alpha.Foils, timepoint != "T3")
Alpha.Foils$InputFileName <- NULL
Alpha.Foils$BarcodeSequence <- NULL
Alpha.Foils$BarcodeSequence_1 <- NULL
Alpha.Foils$LinkerPrimerSequence <- NULL
Alpha.Foils$ReversePrimerSequence <- NULL
Alpha.Foils$revComp <- NULL
Alpha.Foils <- Alpha.Foils %>% rownames_to_column(var = "Sample") 
head(Alpha.Foils )

Alpha.Foils["category"] <- NA

for (i in 1:nrow(Alpha.Foils)){
  
  if (Alpha.Foils$Material[i] %in% c("Nylon","PET")) {
    Alpha.Foils$category[i] <- "Hetero"
  }
  
  else 
    
    if (Alpha.Foils$Material[i] %in% c("PE","PP","PS")) {
      Alpha.Foils$category[i] <- "Carbon"
      
    }
  
}

head(Alpha.Foils)

treat_time <- str_c(Alpha.Foils$treatment, "_", Alpha.Foils$timepoint)
Alpha.Foils<- Alpha.Foils %>%
  add_column(treat_time)

time_mat <- str_c(Alpha.Foils$timepoint, "_", Alpha.Foils$Material)
Alpha.Foils <- Alpha.Foils %>%
  add_column(time_mat)

cat_time_treat <- str_c(Alpha.Foils$category, "_", Alpha.Foils$timepoint, "_", Alpha.Foils$treatment)
Alpha.Foils <- Alpha.Foils %>%
  add_column(cat_time_treat)

#Two options for plotting
#1, select the different columns to obtain the different measures in different columns
obs.ft <- Alpha.Foils %>%  select(Sample, Observed, timepoint:category)
chao1 <- Alpha.Foils %>%  select(Sample, Chao1, timepoint:category)
Shannon <- Alpha.Foils %>%  select(Sample, Shannon, timepoint:category)
Simpson <- Alpha.Foils %>%  select(Sample, Simpson, timepoint:category)

#2, transform the complete df in long format
Alpha.foils.long <- Alpha.Foils %>% select(!se.chao1) %>%  pivot_longer(cols = c("Observed", "Chao1", "Shannon", "Simpson"), names_to = "Diversity_Index", 
                                                  values_to = "Value")


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=T) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  

  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
 

# SUmmarize per measure
obs.ft.c <-  summarySE(obs.ft, measurevar="Observed", groupvars=c("description","Material","timepoint","treatment","category"))
chao1.c <-  summarySE(chao1, measurevar="Chao1", groupvars=c("description","Material","timepoint","treatment","category"))
Shannon.c <- summarySE(Shannon, measurevar="Shannon", groupvars=c("description","Material","timepoint","treatment","category"))
Simpson.c <- summarySE(Simpson, measurevar="Simpson", groupvars=c("description","Material","timepoint","treatment","category"))
# And bint into one df
Alphafoils.c <- bind_rows(obs.ft.c, chao1.c, Shannon.c, Simpson.c)

#Summarize long format
Alpha.foils.long.c <-  summarySE(Alpha.foils.long, measurevar="Value", groupvars=c("description","Material","timepoint","treatment", "category", "Diversity_Index"))
Richness.long.c <- Alpha.foils.long.c %>% filter(Diversity_Index %in% c("Observed", "Chao1"))
Diversity.long.c <- Alpha.foils.long.c %>% filter(Diversity_Index %in% c("Shannon", "Simpson"))

#### Statistics, are Alpha-divs significantly different ####

# select the different columns to obtain the different measures in different columns
obs.ft.v <- Alpha.Foils %>%  select(Observed, timepoint:category)
chao1.v <- Alpha.Foils %>%  select(Chao1, timepoint:category)
Shannon.v <- Alpha.Foils %>%  select(Shannon, timepoint:category)
Simpson.v <- Alpha.Foils %>%  select(Simpson, timepoint:category)
alpha_div_data  <- Alpha.Foils %>%  select(!Observed:Simpson)
head(alpha_div_data)

kt.chao.time <- kruskal.test(chao1$Chao1 ~ timepoint, data = alpha_div_data )
kt.chao.timeUV <- kruskal.test(chao1$Chao1 ~ treat_time, data = alpha_div_data )
kt.chao.cattimetreat <- kruskal.test(chao1$Chao1 ~ cat_time_treat, data = alpha_div_data )
kt.chao.timemat <- kruskal.test(chao1$Chao1 ~ time_mat, data = alpha_div_data )

pt.chao.timeUV <- pairwise.wilcox.test(chao1$Chao1, 
                     g = alpha_div_data$treat_time,
                     p.adjust.method = 'BH')

pt.chao.time <- pairwise.wilcox.test(chao1$Chao1, 
                                       g = alpha_div_data$timepoint,
                                       p.adjust.method = 'BH')


pt.chao.timemat<- pairwise.wilcox.test(chao1$Chao1, 
                                       g = alpha_div_data$time_mat,
                                       p.adjust.method = 'BH')

kt.simp.time <- kruskal.test(Simpson$Simpson ~ timepoint, data = alpha_div_data )
kt.simp.timeUV <- kruskal.test(Simpson$Simpson ~ treat_time, data = alpha_div_data )
kt.simp.cattimetreat <- kruskal.test(Simpson$Simpson ~ cat_time_treat, data = alpha_div_data )
kt.simp.time_mat <- kruskal.test(Simpson$Simpson ~ time_mat, data = alpha_div_data )

pt.simp.timeUV <- pairwise.wilcox.test(Simpson$Simpson, 
                                       g = alpha_div_data$treat_time,
                                       p.adjust.method = 'BH')

pt.simp.time <- pairwise.wilcox.test(Simpson$Simpson, 
                                     g = alpha_div_data$timepoint,
                                     p.adjust.method = 'BH')


####_______________________________________________________________________________________#%###
####          Observed richness                     ####
####_______________________________________________________________________________________#%###
Observed.Foils<- read_pptx()
Observed.Foils <- read_pptx("../Reports/Observed_Foils_202211.pptx")

### alpha div plot 
P.Richness <- ggplot(Richness.long.c,         #Pick data to plot
                   aes(x=interaction(Material, category), y = Value, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se, width=.1)) +
  geom_point(size = 3) +
  facet_nested_wrap(vars(Diversity_Index, timepoint, treatment),
                    nrow = 1,  axes = "margins",
                    remove_labels = "x", shrink = T) +
  guides( x = "axis_nested") + 
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 13, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 13), 
        legend.text=element_text(size = 13),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
        ggh4x.axis.nesttext.x = element_blank()) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  labs( title = "Richness", fill = "Polymer") +
  xlab("Polymer") + 
  ylab("")

P.Richness

### alpha div plot 
P.Diversity <- ggplot(Diversity.long.c,         #Pick data to plot
                     aes(x = Material, y = Value, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se, width=.1)) +
  geom_point(size = 3) +
  facet_nested_wrap(vars(Diversity_Index, timepoint, treatment),
                    nrow = 1,  axes = "margins", scale = "free_y",
                    remove_labels = "x", shrink = T) +
  guides( x = "axis_nested") + 
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  labs( title = "Diversity") 


P.Diversity

## Make plot with GGPLOT
P3 <- ggplot(obs.ft.c,         #Pick data to plot
                     aes(x=interaction(Material, category), y = Observed, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se, width=.3)) +
 geom_point(size = 4) +
  facet_nested( ~ timepoint + treatment) +
  guides( x = "axis_nested") +
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
        ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5)) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  ylab("Observed features") +
  xlab("Polymers") +
  labs( title = "", color = "Polymers", fill = "Polymers")

P3

editable_graph <- dml(ggobj = P3)
Observed.Foils <- add_slide(Observed.Foils) 
Observed.Foils <- ph_with(x = Observed.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Observed.Foils, target = "../Reports/Observed_Foils_202211.pptx")

####_______________________________________________________________________________________#%###
####           Chao1 Index    ####
####_______________________________________________________________________________________#%###
Chao1.Foils<- read_pptx()
Chao1.Foils <- read_pptx("../Reports/Chao1_index_Foils.pptx")

## Make plot with GGPLOT
Chao1 <- ggplot(chao1.c,         #Pick data to plot
                     aes(x=interaction(Material, category), y = Chao1, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Chao1-se, ymax=Chao1+se, width=.3)) +
  geom_point(size = 4) +
  facet_nested( ~ timepoint + treatment) +
  guides( x = "axis_nested") +
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 15, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 15), 
        legend.text=element_text(size = 15),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size= 20),
        strip.text.x = element_text(size = 17),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
        ggh4x.axis.nesttext.x = element_blank()) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  xlab("") +
  labs( title = "", color = "Polymers", fill = "Polymers")

Chao1

editable_graph <- dml(ggobj = Chao1.Plot)
Chao1.Foils <- add_slide(Chao1.Foils) 
Chao1.Foils <- ph_with(x = Chao1.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Chao1.Foils, target = "../Reports/Chao1_index_Foils.pptx")

####_______________________________________________________________________________________#%###
####           Simpson   ####
####_______________________________________________________________________________________#%###
Simpson.Foils <- read_pptx()
Simpson.Foils<- read_pptx("../Reports/Simpson_evenness_Foils.pptx")

Simpson <- ggplot(Simpson.c,         #Pick data to plot
                 aes(x=interaction(Material, category), y = Simpson, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Simpson-se, ymax=Simpson+se, width=.3)) +
  geom_point(size = 4) +
  facet_nested( ~ timepoint + treatment) +
  guides( x = "axis_nested") +
  theme_pubclean()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 15, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 15), 
        legend.text=element_text(size = 15),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size= 20),
        strip.text.x = element_text(size = 17),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
        ggh4x.axis.nesttext.x = element_blank()) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  xlab("") +
  ylab("Gini-Simpson Index") +
  labs( title = "", color = "Polymers", fill = "Polymers")

  
Simpson

editable_graph <- dml(ggobj = plotP1)
Simpson.Foils <- add_slide(Simpson.Foils ) 
Simpson.Foils <- ph_with(x = Simpson.Foils , editable_graph,location = ph_location_type(type = "body") )
print(Simpson.Foils , target = "../Reports/Simpson_evenness_Foils.pptx")


####_______________________________________________________________________________________#%###
####           Shanon ####
####_______________________________________________________________________________________#$###
Shannon.Foils <- read_pptx()
Shannon.Foils <- read_pptx("../Reports/Shannon_Eukaryotes.pptx")

Shannon <- ggplot(Shannon.c,         #Pick data to plot
                 aes(x=interaction(Material, category), y = Shannon, fill = Material, color = Material)) + #Pick factors to use
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se, width=.3)) +
  geom_point(size = 4) +
  facet_nested( ~ timepoint + treatment) +
  guides( x = "axis_nested") +
  theme_pubclean()+
  theme(legend.position = "bottom",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.text.x=element_text(size = 15, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 15), 
        legend.text=element_text(size = 15),
        legend.title = element_text(size=18),
        axis.title.y = element_text(size= 20),
        strip.text.x = element_text(size = 17),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank(),
        ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
        ggh4x.axis.nesttext.x = element_blank()) +
  scale_colour_manual(values = Pal.plast) +
  scale_fill_manual(values = Pal.plast) +
  xlab("") +
  ylab("Shannon index") +
  labs( title = "", color = "Polymers", fill = "Polymers")


Shannon

editable_graph <- dml(ggobj = plotP2)
Shannon.Foils <- add_slide(Shannon.Foils) 
Shannon.Foils <- ph_with(x = Shannon.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Shannon.Foils, target = "../Reports/Shannon_Eukaryotes.pptx")

legend.a <- get_legend(Shannon+
                       theme(legend.direction = "horizontal",
                             legend.title.align = 0.5))



plot_grid(Chao1 + theme(legend.position ="none"),
          Simpson + theme(legend.position ="none"),
          Shannon, 
          ncol = 2,
          nrow = 2,
          align = 'v',
          axis = "l",
          rel_heights = c(1,1.1),
          rel_widths = c(1,1))


####_______________________________________________________________________________________#%###
####           Summary plots for SupInfo to show averages ####
####_______________________________________________________________________________________#%###
Chao1.Foils<- read_pptx()
Chao1.Foils <- read_pptx("../Reports/Chao1_index_Foils.pptx")

## Make plot with GGPLOT
Chao1 <- ggplot(Alpha.Foils,          #Pick data to plot
                 aes(x = treatment, y = Chao1,  color = treatment)) + #Pick factors to use
  geom_boxplot(stat = "boxplot", outlier.colour =  NULL, linewidth = 1) +
  geom_point( position = "jitter", size = 3, alpha = 0.8) +
  theme_pubclean()+
  theme(axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
        axis.text.y=element_text(size= 12), 
        legend.text=element_text(size = 12),
        legend.title = element_text(size=15, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 13),
        plot.title = element_blank(),
        panel.border = element_rect(color = "grey90", fill = NA),
        panel.grid.major.y = element_line(color = "grey90", linetype = 3),
        panel.grid.major.x = element_blank()) +
  facet_grid( ~ timepoint) +
  guides(alpha = "none") +
  xlab("Treatment") +
  scale_colour_manual(values = pal.uv) +
  scale_fill_manual(values = pal.uv) 

Simpson

legend.a <- get_legend(Shannon+
                         theme(legend.direction = "vertical",
                               legend.title.align = 0.5))



plot_grid(Chao1 + theme(legend.position ="none"),
          Simpson + theme(legend.position ="none"),
          Shannon + theme(legend.position ="none"),
          legend.a,
          labels = c("A", "B", "C"),
          ncol = 4,
          nrow = 1,
          align = 'v',
          axis = "l",
          rel_heights = c(1,1,1,1),
          rel_widths = c(1,1,1,0.2))


##############################################################################################
###%#______________________________________________________________________________________#%###
####                Rarefaction                                                            ####
###%#______________________________________________________________________________________#%###
ASV.ABS.count <-  t2 %>%  select(Sample, Description, OTU, Abundance, timepoint, Material, treatment) %>%
  distinct() 
head(ASV.ABS.count)

#Adjust the number of significant digets to get from vegan in rarefactions
old <- options(pillar.sigfig = 8)
options(old)

# Transform count table to get samples per row and ASV per column
# Not per se necessary, but this means we can use standard MARGIN in transfromation and ordination
ASV.ABS_df1<- ASV.ABS.count %>% select(Sample, OTU, Abundance)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>% replace(is.na(.), 0)  %>% 
  arrange(Sample) %>% column_to_rownames(var = "Sample") %>% as.data.frame()
str(ASV.ABS_df1)
rownames(ASV.ABS_df1)
colnames(ASV.ABS_df1)

# Select the metadata we need for plotting
meta_df1 <- ASV.ABS.count %>% select(Sample, OTU, Abundance, timepoint, treatment, Material)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>% replace(is.na(.), 0)  %>% 
  arrange(Sample) %>% column_to_rownames(var = "Sample") %>% as.data.frame() %>% select(timepoint, treatment, Material)
rownames(meta_df1)
colnames(meta_df1)

# # Find the lowest amount of sequences in all samples with this
# # Can be used for rarefaction depth
# min_seqs_1 <- ASV.ABS.count%>%
#   group_by(Sample) %>%
#   summarize(n_seqs = sum(Abundance))  %>% 
#   summarize(min = min(n_seqs)) %>%
#   pull(min)

## The min amount of reads is 844, but that's very low. I will pick the 2nd one of 3360
# # A tibble: 5 × 2
# Sample      n_seqs
# <chr>        <int>
# 1 NIOZ140.108    844
# 2 NIOZ140.91    3360
# 3 NIOZ140.76    3814
# 4 NIOZ140.92    3826
# 5 NIOZ140.87    3936
# 6 NIOZ140.77    4010
# 7 NIOZ140.89    4782
# 8 NIOZ140.65    5106
# 9 NIOZ140.83    6960
# 10 NIOZ140.78   7974

# ASV.ABS_df2 <- ASV.ABS.count %>% select(Description, OTU, Abundance)  %>%
#   group_by(Description, OTU) %>% summarize(count = as.integer(mean(Abundance))) %>% 
#   mutate(across(c(OTU),factor)) %>%
#   pivot_wider(names_from = OTU, values_from = count, values_fill = 0)  %>%
#   column_to_rownames(var = "Description") %>% as.data.frame()
# str(ASV.ABS_df2)
# 
# 
# min_seqs_2 <- ASV.ABS.count%>%
#   group_by(Sample, Description) %>%
#   summarize(n_seqs = sum(Abundance)) %>% ungroup() %>% 
#   group_by(Description) %>% 
#   summarize(av_seqs = mean(n_seqs)) %>% 
#   summarize(min = min(av_seqs)) %>%
#   pull(min) %>% as.integer()
# 
#Remove the first column for rarefaction 
ASV.ABS.df1 <- ASV.ABS_df1[,-1]
# #another way to get a tibble with rarefied data
# samples <- vegan::rarefy(ASV.ABS.df1, 3360) %>%
#   as_tibble(rownames = "Sample") %>%
#   select(Sample, samp.val = value)

# #another way to get a tibble with rarefied data
# ASV.ABS.df2 <- ASV.ABS_df2[,-1]
# reps<- rarefy(ASV.ABS_df2, min_seqs_2) %>% 
#   as_tibble(rownames = "Description") %>% 
#   select(Description, samp.val = value)

# Create the rarefaction curves for both individual samples, and the replicates combined
samples.rarecurve <- vegan::rarecurve(ASV.ABS.df1, step = 250) 
# reps.rarecurve <- vegan::rarecurve(ASV.ABS.df2, step = 100)
#check structure of the object
str(samples.rarecurve)

samp.rare.data <- map_dfr(samples.rarecurve[-99], bind_rows) %>% 
  bind_cols(meta_df1[-c(99),], Sample = rownames(ASV.ABS_df1[-c(99),])) %>% 
  pivot_longer(-c("Sample", "treatment", "timepoint", "Material")) %>% 
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>% 
  select(-name)

head(samp.rare.data)
colnames(samp.rare.data)

# We add an extra column with samplenames to be used for the label
# Group by Sample and make sure there is only a label for the last value of the sample
samp.rare.data <- samp.rare.data %>% group_by(Sample) %>% mutate(label = if_else(value == max(value), as.character(Material), NA_character_))


rareplot.1 <- ggplot(samp.rare.data, aes(x = n_seqs, y = value, group = Sample)) +
  geom_line(aes(color = treatment), linewidth = 1.5) +
  # xlim(0, 50000) +
  facet_grid(.~ timepoint) +
  geom_label_repel(aes(label = label, color = treatment), nudge_x = 1, nudge_y = -0.5, na.rm = T) +
  scale_colour_manual(values = pal.uv) +
  theme_classic()
rareplot.1

rep.rare.data <- map_dfr(reps.rarecurve, bind_rows) %>% 
  bind_cols(Description = rownames(ASV.ABS_df2), .) %>% 
  pivot_longer(-Description) %>% 
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>% 
  select(-name)

head(rep.rare.data)

rareplot.2 <- ggplot(rep.rare.data, aes(x = n_seqs, y = value, color = Description)) +
  geom_line( linewidth = 1.5) +
  geom_vline(xintercept = min_seqs_2, linewidth = 1.5) +
  xlim(0, 50000) +
  theme_classic()


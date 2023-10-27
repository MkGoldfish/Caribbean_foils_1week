##%######################################################%##
#                                                          #
####             Diversity Analysis of                   ####
####            16S  Amplicon Sequencing Data           ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####           Project: OpenBio_ Prokaryotes            ####
####            Adapted by Maaike April 2022            ####
####              (this is the right one)               ####
####               Corrected Primers                    ####
#                                                          #
##%######################################################%##
############################################################
##  version.string R version 4.2.1 (2022-06-23 ucrt)
##  nickname       Funny-Looking Kid   

###%#_______________________________________________________________________________________#%###
####                 Working Directory                                                    ####
###%#_______________________________________________________________________________________#%###
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")

###%#_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
###%#_______________________________________________________________________________________#%###
library(devtools)
library(vegan)
library(compositions)
library(ggforce)
library(scales)
library("remotes")
library("GGally")
library(cowplot)
library(tidyverse)

###%#_______________________________________________________________________________________#%###
####                  Import Data &                    
###%#_______________________________________________________________________________________#%###
#Set seed for randomized processes
set.seed(42)
# set number of permutations for SIMPER
perm <- 9999

#Countdata
#Import earlier created tidy data tibble with pre-processed absolute abundance values
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c("")) 
head(t2)

#Make sure al T3 and negative controls are goooone
t2 <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc") 
#Subset t1 and t6 for tests within timepoints
t1 <- t2 %>% filter( timepoint %in% c("T1"))
t6 <- t2 %>% filter( timepoint %in% c("T6"))

#ASV count table, select relevant info
ASV.count <-  t2 %>%  select(Sample, OTU, Sample_rel_abund) %>%
  distinct() %>% arrange(Sample)
head(ASV.count)

# Transform count table to get samples per row and ASV per column
ASV.count.t <- ASV.count %>% select(Sample, OTU, Sample_rel_abund)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Sample_rel_abund) %>% replace(is.na(.), 0)  %>% 
  column_to_rownames(var = "Sample") %>% as.data.frame()
rownames(ASV.count.t)

#GEnus count tables
## Maybe filter genera with a minimum RA
Genus.count <- t2 %>% select(Sample, Genus, Genus_rel_abund_Sample) %>% filter(Genus != "unassigned") %>%  distinct() %>% arrange(Sample) 

Genus.count.t <- Genus.count %>% select(Sample, Genus, Genus_rel_abund_Sample)  %>%  mutate(across(c(Genus),factor)) %>%
  pivot_wider(names_from = Genus, values_from = Genus_rel_abund_Sample) %>% replace(is.na(.), 0) %>%   arrange(Sample) %>%
  column_to_rownames(var = "Sample") %>% as.data.frame() %>% sqrt()
rownames(Genus.count.t)
colnames(Genus.count.t)
head(Genus.count.t)

####Metadata####
#Remove timepoint T3 and negative controls, in case they are still there
#select relevant columns
#Make sure metadata is sorted by samplename in the same way as countdata and make samplenames rownames
meta <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc")  %>% 
  select(Sample, timepoint, replicate, Material, treatment, Description) %>%  unique() %>%  arrange(Sample) %>% column_to_rownames(var="Sample")
rownames(meta) 

#Omit columns we don't need for analysis
meta$LinkerPrimerSequence <- NULL
meta$ReversePrimerSequence <- NULL
meta$InputFileName <- NULL
meta$BarcodeSequence <- NULL
meta$BarcodeSequence_1 <- NULL
meta$revComp <- NULL

## add category of plastic for nesting in new column to dataframe
meta["category"] <- NA

for (i in 1:nrow(meta)){
  
  if (meta$Material[i] %in% c("Nylon","PET")) {
    meta$category[i] <- "Hetero"
  }
  
  else 
    
    if (meta$Material[i] %in% c("PE","PP","PS")) {
      meta$category[i] <- "Carbon"
      
    }
  
}

#add extra colums we do need
treat_mat <- str_c(meta$treatment, "_", meta$Material)
meta <- meta %>%
  add_column(treat_mat)

treat_time <- str_c(meta$treatment, "_", meta$timepoint)
meta <- meta %>%
  add_column(treat_time)

time_mat <- str_c(meta$timepoint, "_", meta$Material)
meta <- meta %>%
  add_column(time_mat)

cat_time <- str_c(meta$category, "_", meta$timepoint)
meta <- meta %>%
  add_column(cat_time)

cat_treat <- str_c(meta$category, "_", meta$treatment)
meta <- meta %>%
  add_column(cat_treat)

cat_time_treat <- str_c(meta$category, "_", meta$timepoint, "_", meta$treatment)
meta <- meta %>%
  add_column(cat_time_treat)

#Check column and rownames
colnames(meta)
rownames(meta)
metadata <- meta

metadata.t1 <- metadata %>% filter( timepoint %in% c("T1"))
metadata.t6 <- metadata %>% filter( timepoint %in% c("T6"))

head(metadata)

###%#_____________________________________________________________________________________#%###
####           Beta diversity: SIMilarity PErcentage _ Genus/Order level                       ####
###%#_____________________________________________________________________________________#%###
## DIfferent Genera between Timepoints
Simper<- simper(
 Genus.count.t, #Community matrix, samples are rows
   group = metadata$time_mat,#grouping factor
  permutations = perm) # Set number of permutations

summary(Simper, ordered = T)

#What's the structure of the Simper output?
str(Simper)

# # there are 3 lists, one for each pairwise comparison
# # order of OTU names in each of the output lists
# all.equal(colnames(Order.count.t), Simper.habitat$Benthic_Eulittoral$species)
# all.equal(colnames(Order.count.t), Simper.habitat$Benthic_Pelagic$species)
# all.equal(colnames(Order.count.t), Simper.habitat$Eulittoral_Pelagic$species)

#extract overall dissimilarity between habitats
Simper.overall <- lapply( # apply to each element in simper output
  Simper, # simper output
  function(x) { data.frame(
    # create data frame consisting of
    overall = x$overall ) # overall differences between habitats
  }
)

#add timepoint to name of lists in list of lists in case you test timepoint
# to know what is what :) 
names(Simper) <- paste("time_mat_",names(Simper.overall))

# use bind_rows to get a df instead of a list
overall.dissimilarity.time_mat <- bind_rows(Simper.overall,
                                           .id = paste('test'))

overall.dissimilarity.time_mat

# creating simper output tables
# remove $overall (only 1 element, not of the same length as other vectors in list)
# We don't want to print that value into the results files for further procsseing, that's why we remove it 
simper.df <- lapply(
  Simper,
  function(x) {
    do.call(data.frame,
            x[-3])
  }
)

str(simper.df)
summary(simper.df, ordered = T)


# write simper output to file
# one file for each comparison
# we use the created df for this, to make further processing of the files easier 
# only look at OTUs contributing at least to 70% of the differences between groups?? 
for(i in 1:length(simper.df)) {
  write.table(
  simper.df[[i]],
    paste("../SIMPER_202303/SIMPER_genus_time_mat_", # <---- Change filenames here to the test you did
          names(simper.df)[i], 
          ".txt", 
          sep = ""),
    sep = "\t",
    quote = F,
    dec = "."
  )
}

listenames <- lapply(simper.df, function(x) x$listnames)

# Print the extracted listenames
print(listenames)


#This piece of code is to store ungrouped results. 
#Due to the structure of the ungrouped results, missing p value, we can not use the loop.

simper.ungrouped <- lapply( # apply to each element in simper output
  Simper, # simper output
  function(x) {
    data.frame( # create data frame consisting of
      species = x$species, # OTU names
      average= x$average, 
      sd = x$sd,
      ratio = x$ratio,
      ord = x$ord,
      cusum = x$cusum
      # average contribution
    )
  }
)

str(simper.ungrouped)

ungrouped <- bind_rows(simper.ungrouped)

write.table(ungrouped, "../SIMPER_202303/SIMPER_Foils_Orders_Ungrouped.txt",
            na = "NA", dec = '.', quote = F, sep = "\t")

###%#-----------------------------------#%###

###%#----------------------------------#%###
##### Overall dissimilarities    ####
Overall.dissimilarities.genus <- bind_rows(overall.dissimilarity.ungrouped,
                                            overall.dissimilarity.time,
                                            overall.dissimilarity.category,
                                            overall.dissimilarity.treat_time,
                                            overall.dissimilarity.cat_time,
                                            overall.dissimilarity.cat_time_treat,
                                            overall.dissimilarity.polymers,
                                            overall.dissimilarity.time_mat)


write.csv(Overall.dissimilarities.orders, "../SIMPER_202303/Overall_dissimiliraties_Simper_genus.csv",
           na = "NA", dec = ".", quote = F)

write.table(Overall.dissimilarities.genus, "../SIMPER_202303/Overall_dissimiliraties_Simper_genus.txt",
            na = "NA", dec = '.', quote = F, sep = "\t")

###%#_____________________________________________________________________#%###
# convert lists into dataframes
# and store for later use in plotting
# als alternative for the automatic file generation above
#EXAMPLE
T1.T6 <- bind_rows(timepoint.simper.df, .id = "Test") %>% # bind rows to extract df's from list
  remove_rownames() %>% 
  select(Test, species, average, sd, ord, cusum, p) %>%
  filter(!species %in% c("NA")) %>%  # select relevant columns
  slice_min(order_by = cusum, n =10) #seect the first 10 genera, that have the most influence

write.table(T1.T6, "../SIMPER/T1.T6.top10.txt",
            na = "NA", dec = '.', quote = F, sep = "\t")

###%#______________________________________#%###
#### Make Heatmap plots Genera             #####
###%#______________________________________#%###
###### Import datafiles ####
Ungrouped <- read.delim("../SIMPER_202303/Simper_Foils_genus_Ungrouped.txt", header = T, row.names = NULL, na.strings = c(" "))
T1.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_timepoints_timepoints_alltime_total.txt", header = T, row.names = NULL, na.strings = c(" "))
T6.treat <- read.delim("../SIMPER_202303/SIMPER_genus_treat_time_treat_time_UV_T6_noUV_T6.txt", header = T, row.names = NULL, na.strings = c(" "))
UV.T1.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_treat_time_treat_time_UV_T6_UV_T1.txt", header = T, row.names = NULL, na.strings = c(" "))
noUV.T1.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_treat_time_treat_time_noUV_T6_no UV_T1.txt", header = T, row.names = NULL, na.strings = c(" "))

Carbon.Hetero <- read.delim("../SIMPER_202303/SIMPER_genus_category_category_Carbon_Hetero.txt", header = T, row.names = NULL, na.strings = c(" "))
Carbon.Hetero.T6<- read.delim("../SIMPER_202303/SIMPER_genus_cat_time_cat_time_Carbon_T6_Hetero_T6.txt", header = T, row.names = NULL, na.strings = c(" "))
Carbon.T1.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_cat_time_cat_time_Carbon_T6_Carbon_T1.txt", header = T, row.names = NULL, na.strings = c(" "))
Hetero.T1.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_cat_time_cat_time_Hetero_T6_Hetero_T1.txt", header = T, row.names = NULL, na.strings = c(" "))
Carbon.UV.noUV.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_cat_time_treat_cat_time_treat_Carbon_T6_UV_Carbon_T6_noUV.txt", header = T, row.names = NULL, na.strings = c(" "))
Hetero.UV.noUV.T6<- read.delim("../SIMPER_202303/SIMPER_genus_cat_time_treat_cat_time_treat_Hetero_T6_noUV_Hetero_T6_UV.txt", header = T, row.names = NULL, na.strings = c(" "))

PE.Nylon.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_time_mat_T6_Nylon_T6_PE.txt", header = T, row.names = NULL, na.strings = c(" "))
PE.PET.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_time_mat_T6_PE_T6_PET.txt", header = T, row.names = NULL, na.strings = c(" "))
PS.Nylon.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_time_mat_T6_PS_T6_Nylon.txt", header = T, row.names = NULL, na.strings = c(" "))
PS.PE.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_time_mat_T6_PS_T6_PE.txt", header = T, row.names = NULL, na.strings = c(" "))
PS.PP.T6 <- read.delim("../SIMPER_202303/SIMPER_genus_time_mat_T6_PS_T6_PP.txt", header = T, row.names = NULL, na.strings = c(" "))

Nylon.PE <- read.delim("../SIMPER_202303/SIMPER_genus_polymers_Nylon_PE.txt", header = T, row.names = NULL, na.strings = c(" "))
PE.PET <- read.delim("../SIMPER_202303/SIMPER_genus_polymers_PE_PET.txt", header = T, row.names = NULL, na.strings = c(" "))
PP.PET <- read.delim("../SIMPER_202303/SIMPER_genus_polymers_PP_PET.txt", header = T, row.names = NULL, na.strings = c(" "))
PS.PE <- read.delim("../SIMPER_202303/SIMPER_genus_polymers_PS_PE.txt", header = T, row.names = NULL, na.strings = c(" "))
PS.PP <- read.delim("../SIMPER_202303/SIMPER_genus_polymers_PS_PP.txt", header = T, row.names = NULL, na.strings = c(" "))

hcb <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
pdb <- as.character(read_lines("../data/PlasticDB_Prokaryotic_genera.txt"))

#Bund into dataframe
simper.df <- bind_rows("Day1 vs Day6" = T1.T6, 
                       "Day6 UV vs NoUV" = T6.treat,
                       "UV Day1 vs Day6" =  UV.T1.T6,
                       "noUV Day1 vs Day6" = noUV.T1.T6,
                       
                       "C-C vs H-A backbone" =  Carbon.Hetero, 
                       "Day 6 C-C vs H-A backbone" = Carbon.Hetero.T6, 
                       "C-C backbone Day 1 vs Day6" = Carbon.T1.T6,
                       "H-A backbone Day 1 vs Day6" =  Hetero.T1.T6, 
                       "C-C backbone Day 6 UV vs noUV" = Carbon.UV.noUV.T6,
                       "H-A backbone Day 6 UV vs noUV" = Hetero.UV.noUV.T6,
                       
                       "Day 6 PE vs Nylon" = PE.Nylon.T6, 
                       "Day 6 PE vs PET"=  PE.PET.T6, 
                       "Day 6 PS vs Nylon" = PS.Nylon.T6 , 
                       "Day 6 PS vs PE" = PS.PE.T6, 
                       "Day 6 PS vs PP" = PS.PP.T6,  
                       
                       "PE vs Nylon" = Nylon.PE, 
                       "PE vs PET" =  PE.PET , 
                       "PP vs PET" = PP.PET,
                       "PS vs PE" = PS.PE, 
                       "PS vs PP"= PS.PP,  
                       .id = "test") 


#Filter the df on the first XX species based on the cumsum values. 
simper.cumsum <- simper.df %>%  select(test, species, average, sd, ord, cusum) %>% # select relevant columns
  filter(!species %in% c("NA")) %>%  #remove NA genera 
  group_by(test) %>% slice_min(order_by = cusum, n = 10)  %>% ungroup() %>% 
  select(test,species, average, cusum) 

write.csv(simper.cumsum, "../SIMPER_202303/Cumulative_contributions_genera_top20.csv",
          na = "NA", dec = '.', quote = F, sep = "\t")

# Determine individual contribution per genus to cusum
# since cusum is cumulative, we have to get out indivudal contributios
simper.diffsum <- simper.cumsum %>% 
 arrange(test, cusum) %>% 
  mutate(diffsum = cusum - dplyr::lag(cusum, default = dplyr::first(cusum))) %>% ungroup()

head(simper.diffsum)
unique(simper.diffsum$test)

# With the previous piece of code, we got 0 for the first genus per test
# because we took the differce between current and previous
# this way, we get a column with the value for the first genus included. 
simper.diffsum <- simper.diffsum %>% mutate(contribution = ifelse(diffsum <= 0, cusum, diffsum)) 

#expand the dataframe for plotting
simper.diffsum.xp <- simper.diffsum %>%  select(test, species, average, cusum, diffsum, contribution) %>%
   complete(test, species)

#add HCB/PDB info to genus names
simper.diffsum.xp.symb <- simper.diffsum.xp %>% mutate(species = if_else(
  species %in% intersect(pdb,hcb), paste("#+", sep = "_", species), 
  if_else(species %in% pdb, paste("#", sep = "_", species),
          if_else(species %in% hcb, paste("+", sep = "_", species), species))))

viridis <- c("#fde725", "#b5de2b", "#6ece58", "#35b779", "#1f9e89", "#26828e", "#31688e", "#3e4989", "#482878")

Heatmap <- ggplot(simper.diffsum.xp.symb, # plot expanded df
                  aes(y=fct_reorder(species,cusum), x = test)) + #pick values to plot, in this case species vs test, and average dissimilarity as fill. 
  geom_tile(colour = "grey30", fill = 'white') + #color of tile borders
  geom_text(aes(label =  sprintf("%0.1f", round(100 * contribution, digits = 1))), size = 5) + #geom text, the amount of digits
  scale_fill_gradientn(name  = "Cumulative Contribution (%)",  colours = viridis, limits = c(0,0.301), na.value ="#bbbbbb", labels = scales::percent) + #set-up gradient, in this case magma from viridis package
  theme_classic() +
  scale_x_discrete(limits = c( "Day1 vs Day6", "UV Day1 vs Day6", "noUV Day1 vs Day6", "C-C backbone Day 1 vs Day6", "H-A backbone Day 1 vs Day6",   
                               "C-C vs H-A backbone" , "PE vs Nylon", "PE vs PET", "PP vs PET", "PS vs PE", "PS vs PP", 
                                   "Day6 UV vs NoUV", "C-C backbone Day 6 UV vs noUV", "H-A backbone Day 6 UV vs noUV", "Day 6 C-C vs H-A backbone", 
                               "Day 6 PE vs Nylon", "Day 6 PE vs PET", "Day 6 PS vs Nylon", 
                               "Day 6 PS vs PE",  "Day 6 PS vs PP" )) +
  theme(
    axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 12, face = "italic", color = "black"), 
    legend.text=element_text(size = 11),
    legend.title = element_text(size=11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey", fill = NA),
    #strip.background = element_rect( color = "#FFFFFF"),
    legend.position = "bottom"
     )+
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1.5, title.position = "bottom"), size = "none") +
  xlab(label = "SIMPER test") +
  ylab(label = "Genera contributing to dissimilarity")

Heatmap

## Plot overall dissimilarities and contributions in a table
Simper.summ <- read.delim( "../SIMPER_202303/Overall_dissimilarities_cusum_genus_summary.txt", header = T, row.names = NULL) %>%  na.omit()
Simper.summ.p <- Simper.summ %>% pivot_longer(cols = !Test,
                                              names_to = "variable",
                                              values_to = "value")

Table <- ggplot(Simper.summ.p, # plot expanded df
                aes(y=variable, x = Test, fill = value )) + #pick values to plot, in this case species vs test, and average dissimilarity as fill. 
  geom_tile(colour = "grey50", fill = "white") + #color of tile borders
  geom_text(aes(label = sprintf("%0.1f", round(100 * value, digits = 1))), size = 5, fontface = 'bold') + #geom text, the amount of digits
  theme_classic() +
  scale_x_discrete(limits = c( "T1.T6", "UV.T1.T6" ,"noUV.T1.T6", "Carbon.T1.T6", "Hetero.T1.T6",  
                                "Carbon.Hetero", "Nylon.PE", "PE.PET", "PP.PET", "PS.PE", "PS.PP", 
                               "T6.UV.noUV", "Carbon.UV.noUV.T6", "Hetero.UV.noUV.T6", "Carbon.Hetero.T6", 
                               "PE.Nylon.T6", "PE.PET.T6", "PS.Nylon.T6", 
                               "PS.PE.T6",  "PS.PP.T6" )) + 
   scale_y_discrete(breaks = c( "Cummulative.Difference", "Avg.dissimilarity"),
                   limits = c( "Cummulative.Difference", "Avg.dissimilarity"),
                   labels = c( "Cumulative difference\n caused by 10 genera", "Average Dissimilarity" ))+
  theme(
    axis.text.x=element_blank(), 
    axis.text.y=element_text(size = 12, face = "bold", color = "black"), 
    legend.text=element_text(size = 10),
    legend.title = element_text(size=11),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey", fill = NA),
    #strip.background = element_rect( color = "#FFFFFF")
    legend.position = "top" )+
  guides(size = "none") +
  xlab(label = "") +
  ylab(label = "")
 
Table

legend <- get_legend(Heatmap +
                       theme(legend.direction = "horizontal",
                             legend.position = "bottom",
                             legend.justification = "right",
                             legend.margin = margin(0,20,0,0),
                             legend.title.align = 0.5))

plot_grid(Table, 
          Heatmap + theme(legend.position ="none"),
          ncol = 1,
          align = 'v',
          axis = "l",
          rel_heights = c(0.18,1))

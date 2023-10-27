###Differential abundance analyses on genus level

setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/essai_microbimeR")

library(mia)
library(devtools)
library(phyloseq)
library(tidyverse)
library(vegan)
library("TreeSummarizedExperiment")
library(ggpubr)
library("plotrix")
library("FactoMineR")
library("factoextra")
library(usedist)
library(patchwork)
library(ALDEx2)
library(ggh4x)
library("ggpubr")
library(viridis)
library(cowplot)

#####data_import######
genera <- read_lines("./Analysis/genera.txt")
tax_gen <- read.csv2('./Analysis/tax_table_genus_agglom.csv', na.strings = "NA")
rownames(tax_gen) <- genera
#Generate tax table. 
tax_gen <- tax_table(as.matrix(tax_gen))
mat_gen <- as.matrix(read.csv2("./Analysis/genus_matrix.csv"))
rownames(mat_gen) <- genera
mat_gen<- otu_table(mat_gen, taxa_are_rows = T)
map <- sample_data(read.delim('../../Data_corrected_barcodes/metadata/mapping_file_details.txt', row.names = 1, na.strings = c("")))
map <- sample_data(map)
physeq_gen <- merge_phyloseq(tax_gen,mat_gen, map)

### Adding a column to devide plastic types based on backbone structure  
map["category"] <- NA

for (i in 1:nrow(map)){
  if (map$Material[i] %in% c("Nylon","PET")) {
    map$category[i] <- "heteratomes"
  }
  else
    if (map$Material[i] %in% c("PE","PP","PS")) {
      map$category[i] <- "C_backbone"
    }
}

### add columns to mapfile to make combinations of factors for analysis
map["cat_uv"] <- str_c(map$category,"_",map$treatment)
map["time_UV"] <- str_c(map$timepoint, "_", map$treatment)
map["treat_time"] <- str_c(map$treatment, "_", map$timepoint)
map["time_mat"] <- str_c(map$timepoint, "_", map$Material)
map["mat_time"] <- str_c(map$Material, "_", map$timepoint)
map["BarcodeSequence"] <- NULL
map["BarcodeSequence_1"] <- NULL
map["LinkerPrimerSequence"] <- NULL
map["ReversePrimerSequence"] <- NULL
map["InputFileName"] <- NULL

colnames(map)

map <- sample_data(map)

physeq_gen <- merge_phyloseq(tax_gen,mat_gen, map) 
physeq_object <-physeq_gen

####basic_info#####
physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90") #eliminate the sample with 0seq
max(sample_sums(physeq_object))
min(sample_sums(physeq_object)) #534


summarize_phyloseq(physeq_object) # 1,605,618 reads
ntaxa(physeq_object)  ##13169
nsamples(physeq_object) ##63?
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #580

   
physeq_gen <-  subset_samples(physeq_object, !timepoint %in% c("T3"))     


####################TSE######################
tse <- makeTreeSummarizedExperimentFromPhyloseq(physeq_gen) 
rowData(tse)$Genus %>% table()

count(as.data.frame(colData(tse)), category) %>% kable()
set.seed(42)

########## Timepoint #################
x <- aldex.clr(
  reads = assay(tse),
  conds = colData(tse)$category, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)
x_tt <- aldex.ttest(  x, paired.test = FALSE, verbose = FALSE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)


par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05)

aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05)

#Check results
out <- aldex_out  %>%
  filter(wi.eBH <= 0.05)   # here we chose the wilcoxon output rather than tt

#Filter results and add ASV taxanomy to know which ASV is who
out <- aldex_out %>% rownames_to_column(var= "ASV") %>%
   filter(wi.eBH <= 0.05)   # here we chose the wilcoxon output rather than tt
  left_join(asv_taxa, by = "ASV") 

# Store results
write_csv2(out, "./Analysis/aldex_genus_clr_category_output.csv")

############subsetting######
physeq_object <-physeq_gen
####____t6___####
sub_t6 <- subset_samples(physeq_object, timepoint %in% c("T6"))

sub_t6_uv <- subset_samples(sub_t6, treatment == "UV")
sub_t6_no.uv <- subset_samples(sub_t6, treatment == "no UV")

sub_t6_pe <- subset_samples(sub_t6, Material == "PE")
sub_t6_pp <- subset_samples(sub_t6, Material == "PP")
sub_t6_ps <- subset_samples(sub_t6, Material == "PS")
sub_t6_pet <- subset_samples(sub_t6, Material == "PET")
sub_t6_nylon <- subset_samples(sub_t6, Material == "Nylon")

sub_t6_cbackbone <- subset_samples(sub_t6, category =="C_backbone")
sub_t6_hatoms <- subset_samples(sub_t6, category =="heteratomes")

####____t1___#####
sub_t1 <- subset_samples(physeq_object, timepoint %in% c("T1"))

sub_t1_uv <- subset_samples(sub_t1, treatment == "UV")
sub_t1_no.uv <- subset_samples(sub_t1, treatment == "no UV")

sub_t1_pe <- subset_samples(sub_t1, Material == "PE")
sub_t1_pp <- subset_samples(sub_t1, Material == "PP")
sub_t1_ps <- subset_samples(sub_t1, Material == "PS")
sub_t1_pet <- subset_samples(sub_t1, Material == "PET")
sub_t1_nylon <- subset_samples(sub_t1, Material == "Nylon")

sub_t1_cbackbone <- subset_samples(sub_t1, category =="C_backbone")
sub_t1_hatoms <- subset_samples(sub_t1, category =="heteratomes")

###____category____####
sub_cbackbone <- subset_samples(physeq_object, category =="C_backbone")
sub_cbackbone_UV <- subset_samples(sub_cbackbone, treatment == "UV")
sub_cbackbone_noUV <- subset_samples(sub_cbackbone, treatment == "no UV")

sub_hatoms <- subset_samples(physeq_object, category =="heteratomes")
sub_hatoms_UV <- subset_samples(sub_hatoms, treatment == "UV")
sub_hatoms_noUV <- subset_samples(sub_hatoms, treatment == "no UV")

####____polymer_____########
sub_pe <- subset_samples(physeq_object, Material == "PE")
sub_pp <- subset_samples(physeq_object, Material == "PP")
sub_ps <- subset_samples(physeq_object, Material == "PS")
sub_pet <- subset_samples(physeq_object, Material == "PET")
sub_nylon <- subset_samples(physeq_object, Material == "Nylon")

####____polymer_pairwise_comparions_____####
sub_pe.ps <- subset_samples(physeq_object, Material %in% c("PE", "PS"))
sub_pp.ps <- subset_samples(physeq_object, Material %in% c("PP", "PS"))
sub_pe.nyl <- subset_samples(physeq_object, Material %in% c("PE", "Nylon"))
sub_pe.pet <- subset_samples(physeq_object, Material %in% c("PE", "PET"))
sub_pp.pet <- subset_samples(physeq_object, Material %in% c("PP", "PET"))

sub_t6_pe.ps <- subset_samples(sub_t6, Material %in% c("PE", "PS"))
sub_t6_pp.ps <- subset_samples(sub_t6, Material %in% c("PP", "PS"))
sub_t6_pe.nyl <- subset_samples(sub_t6, Material %in% c("PE", "Nylon"))
sub_t6_pe.pet <- subset_samples(sub_t6, Material %in% c("PE", "PET"))
sub_t6_pp.pet <- subset_samples(sub_t6, Material %in% c("PET", "PP"))

####____subsets UV vs non-UV_____######
sub_UV <- subset_samples(physeq_gen, treatment %in% c("UV"))
sub_UV_pe <- subset_samples(sub_UV, Material == "PE")
sub_UV_pp <- subset_samples(sub_UV, Material == "PP")
sub_UV_ps <- subset_samples(sub_UV, Material == "PS")
sub_UV_pet <- subset_samples(sub_UV, Material == "PET")
sub_UV_nylon <- subset_samples(sub_UV, Material == "Nylon")

sub_noUV <- subset_samples(physeq_gen, treatment %in% c("no UV"))
sub_noUV_pe <- subset_samples(sub_noUV, Material == "PE")
sub_noUV_pp <- subset_samples(sub_noUV, Material == "PP")
sub_noUV_ps <- subset_samples(sub_noUV, Material == "PS")
sub_noUV_pet <- subset_samples(sub_noUV, Material == "PET")
sub_noUV_nylon <- subset_samples(sub_noUV, Material == "Nylon")

###%#_________________________#%##
###### TSE subsets ############
tse <- makeTreeSummarizedExperimentFromPhyloseq(sub_pe.ps) 
rowData(tse)$Genus %>% table()

count(as.data.frame(colData(tse)), Material) %>% kable()

x <- aldex.clr(
  reads = assay(tse),
  conds = colData(tse)$Material,
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 150, 
  denom = "all",
  verbose = FALSE
)

x_tt <- aldex.ttest(  x, paired.test = FALSE, verbose = FALSE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05)

aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05)

out <- aldex_out %>% rownames_to_column(var= "Genus") %>%
  filter(wi.eBH <= 0.05)  # %>% # here we chose the wilcoxon output rather than tt
  # left_join(asv_taxa, by = "Genus") 
out

# Store the data
# Change filename according to analysis
write_csv2(out, "./Analysis/aldex_genus_clr_cbakcbone_UV_time_output.csv")

###%#___________________________________________________________#%###
#####  Differential abundant Genera plots                       #####
###%#___________________________________________________________#%###
## Read data
Genus.noUV.time <- read.csv2("./Analysis/aldex_genus_clr_noUV_time_output.csv")
Genus.UV.time <- read.csv2("./Analysis/aldex_genus_clr_UV_time_output.csv")
Genus.hatoms.time <- read.csv2("./Analysis/aldex_genus_clr_hatoms_time_output.csv")
Genus.cback.time <- read.csv2("./Analysis/aldex_genus_clr_cbackbone_time_output.csv")
Genus.cback.UV.time <- read.csv2("./Analysis/aldex_genus_clr_cbakcbone_UV_time_output.csv")
Genus.time<- read.csv2("./Analysis/aldex_genus_clr_timepoint_output.csv", sep = ',', dec = '.', header = T )

#Read HCB and PDB data for plotting
HCB <- as.character(read_lines("./data/Hydrocarbon_degraders_sorted_22_08.txt"))
PDB <- as.character(read_lines("./data/PlasticDB_Prokaryotic_genera.txt"))
plast.HCB <- read_lines("./data/PlasticDB_HCB_genera.txt")
Genus <- read.csv2("./../amplicon-analysis_plots_alfa/Analysis/genus.csv")


#Set categories for facet plotting
HCB_only <- setdiff(HCB, PDB)
PDB_only <- setdiff(PDB, HCB)
plast.hcb.intersect <- intersect(PDB, HCB)


#Combine all data together in one df for plotting
Diff.Ab.Genus.T1.vs.T6 <- bind_rows("Day1 vs Day6" = Genus.time, 
                           "noUV polymers Day 1 vs Day6" = Genus.noUV.time,
                           "UV polymers Day 1 vs Day6" = Genus.UV.time,
                           "H-A backbone Day 1 vs Day6" = Genus.hatoms.time,
                           "C-C backbone Day 1 vs Day6" = Genus.cback.time,
                           "C-C backbone UV Day 1 vs Day6" = Genus.cback.UV.time,
                           .id = "test")

Diff.Ab.Genus <- Diff.Ab.Genus.T1.vs.T6$Genus %>% unique()

write_lines(Diff.Ab.Genus, "Analysis/Diff.Ab.Genera.txt")

# ####_______________________________#%##
# #### Differentially Abundant HCBs   
# Diff.Ab.HCB <- Diff.Ab.Genus.T1.vs.T6$Genus
# Diff.Ab.HCB <- Diff.Ab.HCB[Diff.Ab.HCB%in%plast.HCB] %>% unique(Diff.Ab.HCB) 
# length(Diff.Ab.HCB)
# 
# Oil.gens <- read_csv("./../amplicon-analysis_plots_alfa/Analysis/oil_degraders_genera.csv") 
# Diff.Ab.oil <-  Genus %>% filter(Genus%in%Diff.Ab.HCB) 
# Diff.Ab.oil$mat_time <- str_c(Diff.Ab.oil$Material, "_", Diff.Ab.oil$timepoint)
# 
# #Store results 
# write_csv(Genus_oil, "./Analysis/Diff_ab_oil_degraders_genera.csv")

##%#___________________________________________________________#%###
#####  Differential abundant top10 genera                        #####
###%#___________________________________________________________#%###
## read top genera info (from bubbleplot)
Genus.top10 <- read_lines("../amplicon-analysis_plots_alfa/Analysis/Genus_top5_intersect_1%_20230313.txt")

Diff.Ab.Genus.top10 <- Diff.Ab.Genus.T1.vs.T6 %>%  filter(Genus%in%Genus.top10)

#bind different tests together into one dataframe for plotting
#Add symbols to indicate HCB and PDB in the plot
Diff.Ab.Genus.symb.top10 <- Diff.Ab.Genus.top10 %>% mutate(Genus = if_else(
  Genus %in% intersect(pdb,hcb), paste(Genus, sep = "   ", "#+"), 
  if_else(Genus %in% pdb, paste(Genus, sep = "  ", "#"),
          if_else(Genus %in% hcb, paste(Genus, sep = "  ", "+"), Genus))))

# what are in the ends the results of all time tests? 
time_tests <- Diff.Ab.Genus.symb.top10%>% dplyr::select(test, Genus, we.eBH, wi.eBH, effect)

head(time_tests)
write.csv(time_tests, "Analysis/Overview_Aldex2_genera.csv",
          na = "NA", dec = '.', quote = F)

time_tests.xp <- time_tests %>% complete(test, Genus) 

####___________________________#%###
####__Heatmaps__________________####
#### Differentially abundant genera in the top-10 intersection ####
time_tests.top10 <- Diff.Ab.Genus.symb.top10 %>% dplyr::select(test, Genus, we.eBH, wi.eBH, effect) %>% complete(test, Genus)
head(time_tests.top10)
unique(time_tests.top10$test)

Diff.Ab.Genus.rest <- Diff.Ab.Genus.T1.vs.T6 %>%  filter(!Genus%in%Genus.top10)
Rest <- unique(Diff.Ab.Genus.rest$Genus)

Diff.Ab.Genus.symb.rest <- Diff.Ab.Genus.rest %>% mutate(Genus = if_else(
  Genus %in% intersect(pdb,hcb), paste(Genus, sep = "  ", "#+"), 
  if_else(Genus %in% pdb, paste(Genus, sep = "  ", "#"),
          if_else(Genus %in% hcb, paste(Genus, sep = "  ", "+"), Genus))))

time_tests.rest <- Diff.Ab.Genus.symb.rest %>% dplyr::select(test, Genus, we.eBH, wi.eBH, effect) %>% complete(test, Genus)
head(time_tests.rest)

# pivot data longer to be able to plot q-values
q_vals = time_tests.top10 %>% pivot_longer(cols = we.eBH:wi.eBH,
                                           names_to = "variable",
                                           values_to = "value")

time_tests.top10$Genus <- factor(time_tests.top10$Genus, levels=rev(sort(unique(time_tests.top10$Genus))))

Heatmap.top10 <- ggplot(time_tests.top10, aes(y=Genus, x = test, fill = effect )) +
  geom_tile(color = "grey50")+
  geom_text(aes(label = ifelse(is.na(effect), "", sprintf("%0.2f", round(effect, digits = 2))) , size = 4)) +
  scale_fill_gradientn(name  = "Effect size", limits = c(-3, 3), colours = c( '#00767B', '#238F9D', '#42A7C6', '#60BCE9', '#9DCCEF',  
                                                                                           '#DEE6E7', '#ECEADA', '#F9D576', '#FFB954', '#FD9A44', 
                                                                                           '#F57634', '#E94C1F', '#D11807'), na.value = "#bbbbbb") +
 scale_x_discrete(limits = c("Day1 vs Day6", "UV polymers Day 1 vs Day6", 
                             "noUV polymers Day 1 vs Day6", "C-C backbone Day 1 vs Day6", 
                             "C-C backbone UV Day 1 vs Day6", "H-A backbone Day 1 vs Day6")) +
 theme_classic() +
  theme(
    axis.text.x=element_blank(), 
    axis.text.y=element_text(size= 13, face = "italic"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 15),
    panel.border = element_rect(color = "grey50", fill = NA),
    plot.margin = unit(c(0,0,0,0), "cm"),
    #strip.background = element_rect( color = "#FFFFFF")
    legend.position = "none")+
  guides(size = "none") +
  labs(title = "",
       subtitle = "",
       x = "", y= "", 
       fill = "Effect size")

Heatmap.top10

#### Rest of the differentially abundant genera  ####
time_tests.rest$Genus <- factor(time_tests.rest$Genus, levels=rev(sort(unique(time_tests.rest$Genus))))

Heatmap.rest <- ggplot(time_tests.rest, aes(y=Genus, x = test, fill = effect )) +
  geom_tile(color = "grey50")+
  geom_text(aes(label = ifelse(is.na(effect), "", sprintf("%0.2f", round(effect, digits = 2))) , size = 4)) +
  scale_fill_gradientn(name  = "Effect size", limits = c(-3, 3), colours = c( '#00767B', '#238F9D', '#42A7C6', '#60BCE9', '#9DCCEF',  
                                                                                           '#DEE6E7', '#ECEADA', '#F9D576', '#FFB954', '#FD9A44', 
                                                                                           '#F57634', '#E94C1F', '#D11807'), na.value = "#bbbbbb") +
                                                                                             theme_classic() +
  scale_x_discrete(limits = c("Day1 vs Day6", "UV polymers Day 1 vs Day6", 
                              "noUV polymers Day 1 vs Day6", "C-C backbone Day 1 vs Day6", 
                              "C-C backbone UV Day 1 vs Day6", "H-A backbone Day 1 vs Day6")) +
  theme(
    axis.text.x=element_text(size = 15, angle = 45, hjust = 1),  
    axis.text.y=element_text(size= 13, face = "italic"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 15),
    panel.border = element_rect(color = "grey50", fill = NA),
    #strip.background = element_rect( color = "#FFFFFF")
    legend.position = "bottom",
    plot.margin = unit(c(0,0,0,0), "cm"))+
  labs(title = "",
       subtitle = "",
       x = "", y= "") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "top", ticks.colour = "black"), size = "none")+
  xlab(label = "Pairwise test")

Heatmap.rest

# Combine the two heatmaps as two panels. 
legend <- get_legend(Heatmap.rest +
                       theme(legend.direction = "horizontal",
                             legend.position = "bottom"))
                              

plot_grid(Heatmap.top10 + theme(legend.position ="none"),
          Heatmap.rest + theme(legend.position ="none"),
          legend,
          ncol = 1,
          align = 'v',
          axis = "lrtb",
          labels = c("A", "B", ""),
          rel_heights = c(0.3,1,0.15))


## Facetted heatmap of both data subsets ----------------------------------------
#Combine all data together in one df for plotting
Diff.Ab.Genus.T1.vs.T6 <- bind_rows("Day 1 vs Day 6" = Genus.time, 
                                    "noUV-treated polymers: Day 1 vs 6" = Genus.noUV.time,
                                   "UV-treated polymers: Day 1 vs 6" = Genus.UV.time,
                                    "Hetero atoms: Day 1 vs 6" = Genus.hatoms.time,
                                   "C-backbone: Day 1 vs 6" = Genus.cback.time,
                                    "C-backbone UV-treated: Day 1 vs 6" = Genus.cback.UV.time,
                                    .id = "test")

time_tests <- Diff.Ab.Genus.T1.vs.T6 %>% dplyr::select(test, Genus, we.eBH, wi.eBH, effect) %>% complete(test, Genus) 

time_tests.subs <-time_tests %>% mutate(Subset = case_when(
  Genus %in% Genus.top10 ~ "Top 5 intersection",
  Genus %in% Rest ~ "Other genera",
))

time_tests.subs.symb <- time_tests.subs %>% mutate(Genus = if_else(
  Genus %in% intersect(pdb,hcb), paste(Genus, sep = "   ", "#+"), 
  if_else(Genus %in% pdb, paste(Genus, sep = "  ", "#"),
          if_else(Genus %in% hcb, paste(Genus, sep = "  ", "+"), Genus))))



time_tests.subs.symb$Genus <- factor(time_tests.subs.symb $Genus, levels=rev(sort(unique(time_tests.subs.symb $Genus))))

Heatmap.all <- ggplot(time_tests.subs.symb , aes(y=Genus, x = test, fill = effect )) +
  geom_tile(color = "grey50")+
  geom_text(aes(label = ifelse(is.na(effect), "", sprintf("%0.2f", round(effect, digits = 2))) , size = 4)) +
  scale_fill_gradientn(name  = "Effect size", limits = c(-3, 3), colours = c( '#00767B', '#238F9D', '#42A7C6', '#60BCE9', '#9DCCEF',  
                                                                              '#DEE6E7', '#ECEADA', '#F9D576', '#FFB954', '#FD9A44', 
                                                                              '#F57634', '#E94C1F', '#D11807'), na.value = "#bbbbbb") +
  theme_classic() +
  scale_x_discrete(limits = c("Day 1 vs Day 6", "UV-treated polymers: Day 1 vs 6",
                              "noUV-treated polymers: Day 1 vs 6",  "C-backbone: Day 1 vs 6", 
                              "C-backbone UV-treated: Day 1 vs 6","Hetero atoms: Day 1 vs 6" )) +
  facet_nested(fct_relevel(Subset, "Top 5 intersection", "Other genera") ~ ., drop = F,
                scales = "free_y", space = "free_y",
                axes = 'margins', as.table = T, 
                nest_line = element_line(),
                strip = strip_nested(background_y =  elem_list_rect(color = "grey30",
                                                                     fill = "white",  linewidth = 1),
                                      text_y = elem_list_text(size = 13, 
                                                              color =  "grey30", by_layer_y = F))) + 
  guides(y="axis_nested") +
  theme(
    axis.text.x=element_text(size = 12, angle = 45, hjust = 1),  
    axis.text.y=element_text(size= 12, face = "italic"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=15),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 15),
    panel.border = element_rect(color = "grey50", fill = NA),
    #strip.background = element_rect( color = "#FFFFFF")
    legend.position = "bottom",
   )+
  labs(title = "",
       subtitle = "",
       x = "", y= "") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "top", ticks.colour = "black"), size = "none")+
  xlab(label = "Pairwise test")

Heatmap.all
#### Heatmap qvalues ####
# pivot data longer to be able to plot q-values
q_vals = time_tests.xp %>% pivot_longer(cols = we.eBH:wi.eBH,
                                        names_to = "variable",
                                        values_to = "value") 

#Plot heatmap
Heatmap <- ggplot(q_vals, aes(y=Genus, x = test, fill = value )) +
  geom_tile(color = "grey50")+
  geom_text(aes(label = round(value, digits = 3)), color = "black", size = 4) +
  scale_fill_gradientn(name  = "q-values", limits = c(0,0.2),  colours = rev(viridis(20, option = "viridis")), na.value ="#bbbbbb") +
  facet_grid(.~variable ) +
  theme_classic() +
  
  theme(
    axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 12, face = "italic"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    #strip.text.x = element_text(size = 15),
    #strip.placement = "outside", 
    #strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 20),
    panel.border = element_rect(color = "grey50", fill = NA),
    #strip.background = element_rect( color = "#FFFFFF")
    legend.position = "top" )+
  labs(title = "q-values per genus Aldex2 tests") +
  xlab(label = "Pairwise test")

Heatmap 
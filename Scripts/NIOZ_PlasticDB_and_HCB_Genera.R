##%#_______________________________________________________________________________________#%###
####                 Working Directory                                                      ####       
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")

###%#_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
###%#_______________________________________________________________________________________#%###
library(devtools)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(stringr)
library(ggh4x)
library("ggthemes")

set.seed(42)
###%#_______________________________________________________________________________________#%###
####                  Import Data                                                            ####
###%#_______________________________________________________________________________________#%###
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c(""))

# Remove irrelevant columns
t2$LinkerPrimerSequence <- NULL
t2$ReversePrimerSequence <- NULL
t2$InputFileName <- NULL
t2$BarcodeSequence <- NULL
t2$BarcodeSequence_1 <- NULL
t2$revComp <- NULL

#Add relevant columns for plotting and facetting
t2 <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc") 
time_UV <- str_c(t2$timepoint, "_", t2$treatment)
t2 <- t2 %>%
  add_column(time_UV, .before ="Kingdom")
treat_time <- str_c(t2$treatment, "_", t2$timepoint)
t2 <- t2 %>%
  add_column(treat_time, .before ="Kingdom")
time_mat <- str_c(t2$timepoint, "_", t2$Material)
t2 <- t2 %>%
  add_column(time_mat, .before ="Kingdom")
t2 <- t2 %>% mutate( timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6"
)
)

###%#_______________________________________________________________________________________#%###
####_Colors!!______________________________________________________________________________####
# Different colorpallettes to choose from
pal.cb.tol.4 <- c('#77AADD', '#EE8866', '#AAAA00', '#99DDFF', '#9970AB', "#332288" ,
                           '#44BB99', '#5AAE61', '#222255', '#225555', '#225522', "#362C28", "#F4978E",
                           '#666633', '#663333', '#555555', '#EC7014', '#0077BB',  "#7A9B76", "#143642",
                           '#EE3377', '#33BBEE', '#CC3311', '#009988', '#FFB954', '#993404', "#5171A5", '#F7F056',  
                           '#1B7837', '#A01813', '#762A83', '#E94C1F',  '#b09e99' )
                           
pal.cb.mk.1        <- c( '#009175', '#FF66FD', '#FFB935', '#DE0D2E', 
                                  '#00EBC1', '#8E06CD', '#FF8735', '#FFE239',
                                  '#005FCC', '#FFACC6', '#005A01',  
                                  '#00C2F9', '#810D49', '#00F407',
                                  '#005745', '#ED0DFD', '#FF4235',
                                  '#00AF8E', '#6B069F', '#009503',
                                  '#0079FA', '#FF78AD', '#00D302',
                                  '#00E5F8', '#AB0D61', '#AFFF2A',
                                  '#00735C', '#FFA3FC', '#86081C',
                                  '#00CBA7', '#B40AFC', '#B20725',  
                                  '#00489E', '#FF2E95', '#00B408', 
                                  '#009FFA', '#D80D7B', '#007702', "#999999"  )
                                  
                        

Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

#___Genus Relative Abundance________________________________________________________________####
#select or upload data
# Genus <- read.csv("../Analysis/genus.csv")
Genus <-t2  %>%  dplyr::select(Sample, timepoint, Material, treatment, surface, 
                               Description, time_UV, treat_time, time_mat,
                               Genus, Genus_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("Day 1", "Day 6"))

###%#____________________________________#%###
#### Genera in Plastic DB                 ####
###%#____________________________________#%###
PDB <- read.delim("../data/PlasticDB_species_taxonomy.txt", header = T, na.strings = c(""))
HCB <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
plast.deg <- as.character(read_lines("../data/PlasticDB_Prokaryotic_genera.txt"))
plast.HCB <- read_lines("../data/PlasticDB_HCB_genera.txt")

#Set categories for facet plotting
HCB_only <- setdiff(HCB, plast.deg)
PDB_only <- setdiff(plast.deg, HCB)
plast.hcb.intersect <- intersect(plast.deg, HCB)

#Which genera do we find per cateory and how many?
Genus_plastic <- Genus %>%  filter(Genus%in%PDB$Genus)
Genus_oil <-  Genus %>%  filter(Genus%in%HCB)
Genus_plastoil <- Genus %>%  filter(Genus%in%plast.HCB) 
Genus_HCBonly <-  Genus %>%  filter(Genus%in%HCB_only)

length(unique(Genus_oil$Genus))
length(unique(Genus_plastic$Genus))
length(unique(Genus_plastoil$Genus))
length(unique(Genus_HCBonly$Genus))

#Define genera found in our dataset
found.genera <- Genus$Genus %>% unique() %>%  sort()
str(found.genera)

gen.in.pdb <- PDB %>% filter(Genus %in% found.genera)
pdb.gens <- PDB$Genus %>% unique()
length(pdb.gens)
length(unique(gen.in.pdb$Genus))

plastics <- PDB$Plastic %>% unique() %>% sort()
pdb.gen.filt <-gen.in.pdb %>% mutate(Polymer = if_else(str_detect(Plastic, "PET"), "PET",
                                                       (if_else(str_detect(Plastic, "PS"), "PS", 
                                                                (if_else(str_detect(Plastic, "PP"), "PP", 
                                                                         (if_else(Plastic %in% c("PEA", "PEC", "PEF", "PEG", "PES","PESU" ), "Other", 
                                                                            (if_else(str_detect(Plastic, "PE"), "PE", 
                                                                              (if_else(str_detect(Plastic, "Nylon"), "Nylon", 
                                                                            "Other"))))))))))))
PE <- pdb.gen.filt %>% filter(Polymer == "PE") %>% pull(Genus) %>% unique()
PP <- pdb.gen.filt %>% filter(Polymer == "PP") %>% pull(Genus) %>% unique()
PS <- pdb.gen.filt %>% filter(Polymer == "PS") %>% pull(Genus) %>% unique()
PET <- pdb.gen.filt %>% filter(Polymer == "PET") %>% pull(Genus) %>% unique()
Nylon <- pdb.gen.filt %>% filter(Polymer == "Nylon") %>% pull(Genus) %>% unique()

# Put genera in alphabetcial order
pdb.gen.filt$Genus <- factor(pdb.gen.filt$Genus, levels=rev(sort(unique(pdb.gen.filt$Genus))))

# set colors and shape
shapes = c(15,19,8,17,18,4)

#### Plot genera from our dataset present in plasticDB, and the polymers upon which they are found in plasticDB ####
plot.pdb <- ggplot(pdb.gen.filt, aes(x = Polymer, y = Genus)) +
  geom_point(aes(shape = Polymer), size = 3) +
   scale_shape_manual(values = shapes, limits = c("PE", "PP", "PS", "PET", "Nylon", "Other")) +
  scale_x_discrete(limits = c("PE", "PP", "PS", "PET", "Nylon", "Other")) +
  theme_pubclean() +
  theme(
    axis.text.x=element_text(size = 11, angle = 60, hjust = 0.9), 
    axis.text.y=element_text(size= 11, face = "italic"), 
    legend.text=element_text(size = 11),
    legend.title = element_text(size=12),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(),
    panel.border = element_rect(color = "grey", fill = NA),
    legend.key = element_rect(fill = NA), 
    legend.position = "right") +
  xlab("Plastics in PlasticDB")+ 
  ylab("Genera in dataset found in PlasticDB") + 
  labs(title = "")

plot.pdb


#### Bubbleplot of PlaticDB and HCB genera ####
# All genera that are both PDB and HCB, filter RA>0.25%
HCB_genera <- Genus_plastoil %>% dplyr::select( Description, Genus, Genus_rep_rel_abund, 
                                                timepoint, Material, treatment) %>% filter(Genus_rep_rel_abund > 0.0025)

#Add categories for plottingt to dataset
#HCB, PDB. Intersect category included
HCB_genera <- HCB_genera %>% mutate(Degrading = case_when(
  Genus %in% HCB_only ~ "HCB",
  Genus %in% PDB_only ~ "PlasticDB genus",
  Genus %in% plast.hcb.intersect ~ "HCB in PlasticDB",
) )

# Add backbone category for plotting
HCB_genera["category"] <-NA
for (i in 1:nrow(HCB_genera)){
  
  if (HCB_genera$Material[i] %in% c("Nylon","PET")) {
    HCB_genera$category[i] <- "Hetero"
  }
  else {HCB_genera$category[i] <- "Carbon"
  }
}

# Add the info on which polymer category the genera are found in plasticDB
HCB_genera_plasts <- HCB_genera %>% mutate(Genus.pol = 
                                             if_else(Genus %in% PE, paste(Genus, "I", sep = " "), Genus)) 
HCB_genera_plasts <- HCB_genera_plasts %>% mutate(Genus.pol =                                                 
                                                    if_else(Genus %in% PP, paste(Genus.pol, "II", sep = " "),
                                                            if_else(Genus %in% PS, paste(Genus.pol, "III", sep = " "),
                                                                    if_else(Genus %in% PET, paste(Genus.pol, "IV", sep = " "),
                                                                            if_else(Genus %in% Nylon, paste(Genus.pol, "V", sep = " "),
                                                                                    Genus.pol)))))
# #Put genera in alphabetical order
# HCB_bubble$Genus <- factor(HCB_bubble$Genus, levels=rev(sort(unique(HCB_bubble$Genus))))

# Bubbleplot with top n HCB crossection of all samples
HCB_bubble <- ggplot(HCB_genera_plasts, aes(x= interaction(Material, category), y= fct_reorder(Genus.pol,Genus_rep_rel_abund) )) +
  geom_point(aes(size=Genus_rep_rel_abund, fill = factor(Material)), alpha = 1, 
             shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,8))+
  guides(y="axis_nested", x = "axis_nested", fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("Polymer") +  
  facet_nested( Degrading ~ timepoint + treatment, drop = T,
                scales = "free_y", space = "free_y",
                axes = 'margins',switch = "y", 
                nest_line = element_line(),
                strip = strip_nested(
                  background_y =  elem_list_rect(color = c("grey10", "grey50", "grey30"),
                                                 fill = rep_len("white", 3),  linewidth = 1),
                  text_y = elem_list_text(size = rep_len(14,3), angle = rep_len(90,3), 
                                          color = c("grey10", "grey50", "grey30"), by_layer_y = F),
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  theme_bw()+
  theme(
    axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 12, face = 'italic'), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=14),
    legend.position = "right",
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.placement = "outside", 
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(size = 20, hjust = 0),
    panel.border = element_rect(color = "grey"),
      ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(size = 12, angle = 0, color = c("black", "grey"))) +
  labs(title = "", fill = "Polymer", size = "Relative \nAbundance",
       subtitle = "") 

HCB_bubble

### Barplot PDB-only genera ####
Genus.plast <- Genus_plastic

Genus.plast["category"] <- NA

for (i in 1:nrow(Genus.plast)){
  if (Genus.plast$Material[i] %in% c("Nylon","PET")) {
    Genus.plast$category[i] <- "Hetero"
  }
  
  else 
    if (Genus.plast$Material[i] %in% c("PE","PP","PS")) {
      Genus.plast$category[i] <- "Carbon"
    }
}


Plast.filtab <- Genus.plast %>% mutate(Genus = ifelse(Genus_rep_rel_abund<0.0025, "others<0.25%", Genus)) 
length(unique(Plast.filtab$Genus))

PDB.bar <- ggplot(Plast.filtab, aes(x=interaction(Material, category), y= Genus_rep_rel_abund, 
                                      fill=forcats::fct_reorder(Genus,Genus_rep_rel_abund )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(pal.cb.tol.4)) + 
  guides( x = "axis_nested",  fill = guide_legend(ncol = 1)) +
  facet_nested( ~ timepoint + treatment, drop = T, 
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
                )) + 
  
  theme_classic2() + 
  theme(
    axis.text.x=element_text( size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text( size= 12), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Relative Abundance")+xlab("Material") + 
  labs(fill = "Genus", title = "Genera found in PlasticDB") 

PDB.bar

# Create Barplot with relative abundance data of all HCB
ggplot(Genus_HCBonly, aes(x=Material, y= Genus_rep_rel_abund, fill=forcats::fct_reorder(Genus, Genus_rep_rel_abund)))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = pal.cb.mk.1) + 
  theme_classic2() +
  theme(
    legend.text=element_text(size = 10),
    legend.title = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15)) +
  facet_grid (treatment~timepoint ) + 
  guides(fill = guide_legend(ncol = 1, reverse = T)) +
  ylab("Relative Abundance")+ xlab("Material") + labs(fill = "HCB")


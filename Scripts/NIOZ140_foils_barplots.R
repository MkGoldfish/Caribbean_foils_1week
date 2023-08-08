
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
library(ggh4x)
library(rvg)
library(officer)
library("ggthemes")
library("wesanderson")
library("IslamicArt")
library("Polychrome")
library("colorblindr")



###%#_______________________________________________________________________________________#%###
####                  Import Data                                                            ####
###%#_______________________________________________________________________________________#%###
map <- sample_data(read.delim('../v460_data/mapping_file_details_corrections290322_no_t3.txt', row.names = 1, na.strings = c("")))
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv")

# Remove unnecessary columns
t2$LinkerPrimerSequence <- NULL
t2$ReversePrimerSequence <- NULL
t2$InputFileName <- NULL
t2$BarcodeSequence <- NULL
t2$BarcodeSequence_1 <- NULL
t2$revComp <- NULL

# Concatonate information to new columns, to create new categories for analysis. 
time_UV <- str_c(t2$timepoint, "_", t2$treatment)
t2 <- t2 %>% 
  add_column(time_UV, .before ="Kingdom")

treat_time <- str_c(t2$treatment, "_", t2$timepoint)
t2 <- t2 %>% 
  add_column(treat_time, .before ="Kingdom")

time_mat <- str_c(t2$timepoint, "_", t2$Material)
t2 <- t2 %>% 
  add_column(time_mat, .before ="Kingdom")

set.seed(42)
                                                                   
###%#_______________________________________________________________________________________#%###
####_Colors!!______________________________________________________________________________####
# Different colorpallettes to choose from
colors_M1 <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
                        "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
                        "#570000", "#FFF5B2")

colors_M3 <- c(   "#339989", "#faf3dd", "#04724d", "#98B9AB",
                           "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                           "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", 
                           "#662C91")
                           
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

##%#_______________________________________________________________________________________________#%###
##__Creating Barplots____________________________________________________________________________####

####___Kingdom Absolute Abundance________________________________________________________________####
#Kingdom_abs <- read.csv("../Analysis/Kingdom.csv")
Kingdom_abs <- t2  %>%  select(Abundance,timepoint, Material, treatment, surface, description, time_UV, treat_time, 
                               time_mat, Kingdom, Kingdom_rep_rel_abund, Kingdom_st_dev_abund_samples)%>% 
  #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))
  distinct() %>% filter( timepoint %in% c("T1", "T6") )

ggplot(Kingdom_abs, aes(x=Material , y= Abundance, fill=Kingdom))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = colors_M3)+
  theme_classic2()+  
  facet_grid (treatment~timepoint) +
  ylab("Abundance")+xlab("Samplepoint") + labs(fill = "Polymer")


####___Kingdom Relative Abundance________________________________________________________________####
#Kingdom <- read.csv("../Analysis/Kingdom.csv")
Kingdom <- t2  %>%  select(timepoint, Material, treatment, surface, description, time_UV, treat_time, time_mat, Kingdom, Kingdom_rep_rel_abund, Kingdom_st_dev_abund_samples)%>% 
  distinct() 
Kingdom <- Kingdom %>% filter ( timepoint %in% c("T1", "T6")) #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))

ggplot(Kingdom, aes(x=Material, y= Kingdom_rep_rel_abund, fill=Kingdom))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = colors_M3)+
  theme_classic2()+  
  facet_grid (timepoint~treatment) +
  ylab("Relative Abundance")+xlab("Material") + labs(fill = "Kingdom")


####___Phylum Relative Abundance________________________________________________________________####
#creating table phyla
# Phyla <- read.csv("../Analysis/phyla.csv")
Phyla <- t2 %>% select(treatment, timepoint, Material, Description, Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
  distinct() 

# Mutate table so taxa with an RA below certain treshold are lumped in one category
# to make plots easier to read. 
Phyla_filtAb <- Phyla %>% mutate(Phylum = ifelse(Phyla_rep_rel_abund<0.01, "others<1.0%", Phylum)) %>%
  filter ( timepoint %in% c("T1", "T6")) 

Phyla_filtAb["category"] <- NA

for (i in 1:nrow(Phyla_filtAb)){
  
  if (Phyla_filtAb$Material[i] %in% c("Nylon","PET")) {
    Phyla_filtAb$category[i] <- "Hetero"
  }
  
  else 
    
    if (Phyla_filtAb$Material[i] %in% c("PE","PP","PS")) {
      Phyla_filtAb$category[i] <- "Carbon"
      
    }
  
}


Phyla_bar <- ggplot(Phyla_filtAb, aes(x=interaction(Material, category), y= Phyla_rep_rel_abund, 
                         fill=forcats::fct_reorder(Phylum,Phyla_rep_rel_abund )))+
  geom_bar(stat="identity", position="stack")+ 
  scale_fill_manual(values = rev(colors_M1)) + 
  guides( x = "axis_nested",  guide_legend(ncol = 1)) +
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
   ylab("Relative Abundance")+xlab("Material") + labs(fill = "Phylum")

Phyla_bar

# select top6 phyla to plot
Prot <- Phyla %>% filter(Phylum %in% c("Proteobacteria", "Bacteroidetes", "Firmicutes", "Actinobacteria", "Cyanobacteria", "Verrucomicrobiota"))

p2 <- ggplot(Prot, aes(x=description, y= Phyla_rep_rel_abund, fill=Phylum)) +
  geom_col(stat="identity", position="dodge")+ scale_fill_manual(values = colors_by_Maaike) +
  theme_pubclean()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  labs(title = "OTU", subtitle = NULL, x =NULL , y = NULL ) + ylim(0,1)+  ##Labels
  geom_errorbar(aes(ymin=Phyla_rep_rel_abund, ymax=Phyla_rep_rel_abund+st_dev_Phylum_abund), width=.2,
                position=position_dodge(.9))
p2






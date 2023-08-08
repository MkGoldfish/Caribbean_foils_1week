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
library(rvg)
library(officer)
library("ggthemes")
library("stringr")
library(cowplot)
library (emojifont)

set.seed(42)
###%#_______________________________________________________________________________________#%###
####                  Import Data                                                            ####
###%#_______________________________________________________________________________________#%###
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c(""))

t2$LinkerPrimerSequence <- NULL
t2$ReversePrimerSequence <- NULL
t2$InputFileName <- NULL
t2$BarcodeSequence <- NULL
t2$BarcodeSequence_1 <- NULL
t2$revComp <- NULL

t2 <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc") 
t2 <- t2 %>% mutate( timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6"
)
)

time_UV <- str_c(t2$timepoint, "_", t2$treatment)
t2 <- t2 %>%
  add_column(time_UV, .before ="Kingdom")

treat_time <- str_c(t2$treatment, "_", t2$timepoint)
t2 <- t2 %>%
  add_column(treat_time, .before ="Kingdom")

time_mat <- str_c(t2$timepoint, "_", t2$Material)
t2 <- t2 %>%
  add_column(time_mat, .before ="Kingdom")


###%#_______________________________________________________________________________________#%###
####_Colors!!______________________________________________________________________________####
# Different colorpallettes to choose from
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")
pal.phylum <- c("#F2AF29","#ecc8af", "#c18c5d", "#495867","#436436", "#004e64")

###%#_______________________________________________________________________________________#%###
####                    Bubbleplots                                                         ####
###%#_______________________________________________________________________________________#%###

#___Genus Relative Abundance________________________________________________________________####
#select or upload data
# Genus <- read.csv("../Analysis/genus.csv")
Genus <- t2 %>% dplyr::select(timepoint, Material, treatment, surface, Description, Phylum, 
                               Order, Genus, Genus_rep_rel_abund)%>% 
  distinct() %>% filter ( timepoint %in% c("Day 1", "Day 6")) 

#Somehow Prevotella_7 was a result from SIlva in our dataset, so we renamed this annotation
Genus <- Genus %>%  mutate(Genus = ifelse(Genus == "Prevotella_7", "Prevotella", Genus))

#Cout how many unique genera we have. 
Genus  %>% select(Genus) %>%  unique() %>% count()

# select top-5 genera per sample Genera for plotting
Genus_top <- Genus %>% dplyr::select(Description, Genus, Genus_rep_rel_abund) %>% 
  filter ( !Genus %in% c("NA", "unassigned")) %>%  #First remove unassigned genera
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rep_rel_abund, n = 5) %>% ungroup()

#Check how much/which genera we find
unique(Genus_top$Genus)
# and with an RA above 1%
Genus_top %>% filter(Genus_rep_rel_abund > 0.01) %>% select(Genus) %>%  unique()

#Filter genera with RA>1%, to avoid bubbles with 0 value
top_genus= Genus %>% filter(Genus%in%unique((c(Genus_top$Genus)))) %>%  
  filter(Genus_rep_rel_abund > 0.01) 
top_genus %>% select(Genus) %>% unique() 

#Store top-5 to compare to differential abundance analysis
top.5 <- unique(top_genus$Genus)
write_lines(top.5, '../Analysis/Genus_top5_intersect_1%_20230313.txt')

# add category as column for plotting this category
top_genus["category"] <- NA

for (i in 1:nrow(top_genus)){
  
  if (top_genus$Material[i] %in% c("Nylon","PET")) {
    top_genus$category[i] <- "Hetero"
  }
  
  else 
    
    if (top_genus$Material[i] %in% c("PE","PP","PS")) {
      top_genus$category[i] <- "Carbon"
      
    }
  
}

head(top_genus)
####____Genus_per_Order_with_nested_axis_and_facets_______________________________________####
# hcb <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
# pdb <- as.character(read_lines("../data/PlasticDB_Prokaryotic_genera.txt"))
# 
# pdbgen <- Genus %>% filter(Genus %in% pdb)
# pdbgen$Genus  %>% unique()
# hcbgen <- top_genus %>% filter(Genus %in% hcb) 
# hcbgen$Genus %>% unique()
# intgen <- top_genus %>% filter(Genus %in% intersect(pdb,hcb)) 
# intgen$Genus %>% unique()

top_genus$Genus <- factor(top_genus$Genus, levels=rev(sort(unique(top_genus$Genus))))
top_genus <- top_genus %>% arrange(Order)

Genus_bubble <- ggplot(top_genus,aes(x=interaction(Material, category),y= Genus)) +
  geom_point(aes(size=Genus_rep_rel_abund, fill = factor(Material)), shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides( x = "axis_nested",  fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("") +  
  facet_nested(Phylum + Order ~ timepoint + treatment, drop = T, 
               scales = "free_y", space = "free_y",
               axes = 'margins',
               nest_line = element_line(),
               strip = strip_nested(
                 background_y =  elem_list_rect(color = "grey25",
                                                fill = "white", linewidth = 1),
                 text_y = elem_list_text(size = 12, angle = 0, 
                                         color = "grey25", by_layer_y = F),
                 background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0))
               )) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text( size = 14, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 13, face = "italic", color = "black"), 
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
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "", subtitle = "",
       fill = "Polymer", size = "Relative Abundance") 

Genus_bubble

####____Orders_per_phylum_nested_facets/axes_______________________________________####
Orders <-t2  %>%  dplyr::select(timepoint, Material, treatment, surface, 
                               Description, time_UV, treat_time, time_mat, Phylum,  
                               Order, Order_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("Day 1", "Day 6"))

head(Orders)

#Select top-5 orders per sample for plotting, only use Orders with RA>0.05%
Orders_top <- Orders %>% dplyr::select(Description, Order, Order_rep_rel_abund) %>% 
  mutate(across(c(Description),factor))%>% distinct() %>% filter(Order_rep_rel_abund > 0.005) %>% 
  group_by(Description) %>% slice_max(order_by = Order_rep_rel_abund, n = 5)
head(Orders_top)

#Store info for further use
top.5 <-unique(Orders_top$Order)
write_lines(top.6, '../Analysis/Orders_top5_intersect.txt')

sort(top.5)

#Remove unassigned orders, and only keep Orders with RA>0.1%
top_orders = Orders %>% filter(Order%in%unique(c(Orders_top$Order)))  %>% 
  filter( Order!="unassigned") %>% 
  filter(Order_rep_rel_abund > 0.01) 
unique(top_orders$Order)

#Add backbone category for plotting
top_orders["category"] <- NA

for (i in 1:nrow(top_orders)){
  if (top_orders$Material[i] %in% c("Nylon","PET")) {
    top_orders$category[i] <- "Hetero"
  }
  else 
    if (top_orders$Material[i] %in% c("PE","PP","PS")) {
      top_orders$category[i] <- "Carbon"
    }
}

#Put Orders in alphabetical order
top_orders$Order <- factor(top_orders$Order, levels=rev(sort(unique(top_orders$Order))))

Order_bubble <- ggplot(top_orders,aes(x=interaction(Material, category),y= Order )) +
  geom_point(aes(size=Order_rep_rel_abund, fill = factor(Material)), alpha = 1, shape = "circle filled", stroke = 1, colour = "black") +
  scale_fill_manual(values = Pal.plast) +
  scale_size(range = c(2,9.5))+
  guides(y="axis_nested", x= "axis_nested", fill = guide_legend(override.aes = list(size = 10))) +
  ylab("") +
  xlab("Polymer") +  
  facet_nested(Phylum ~ timepoint + treatment, drop = T, 
                scales = "free_y", space = "free_y",
                axes = 'margins',
                nest_line = element_line(),
                strip = strip_nested(
                  background_y =  elem_list_rect(color = c("#F2AF29","#ecc8af", "#c18c5d", "#495867","#436436", "#004e64"
                                                            ),
                                                 fill = rep_len("white", 6), linewidth = rep_len(1,6)),
                  text_y = elem_list_text(size = rep_len(13,6), angle = rep_len(0,6), face = rep_len("bold", 6),
                                          color = c("#F2AF29","#ecc8af", "#c18c5d", "#495867","#436436", "#004e64"
                                          )),
                  background_x = (element_rect(fill = "grey90", color = "grey90", linetype = 0)
                ))) + 
  theme_minimal()+
  theme(
    axis.text.x=element_text(size = 12, angle = 60, hjust = 1), 
    axis.text.y=element_text(size= 12, face = "italic"), 
    legend.text=element_text(size = 12),
    legend.title = element_text(size=13),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    strip.text.x = element_text(size = 13),
    plot.title = element_text(size = 10),
    panel.border = element_rect(color = "grey90", fill = NA),
    ggh4x.axis.nestline.x = element_line(linetype = c(6,1), linewidth = 1, color = c("black", "darkgrey")),
    ggh4x.axis.nesttext.x = element_text(angle = 0, color = c("black", "darkgrey"), hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linetype = 3),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom") +
  labs(title = "Intersection of the top-5 orders per sample with relative abundance >1%",
       fill = "Polymer", size = "Relative Abundance")

Order_bubble




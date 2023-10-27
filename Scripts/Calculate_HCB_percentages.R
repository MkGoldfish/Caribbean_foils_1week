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


t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv")

Genus <-t2  %>%  dplyr::select( Description, timepoint, Material, treatment, 
                               Genus, Genus_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("T1", "T6"))

Sub_t1 <- Genus %>% filter ( timepoint %in% c("T1"))
Sub_t6 <- Genus %>% filter ( timepoint %in% c("T6"))

HCB <- as.character(read_lines("../data/Hydrocarbon_degraders_sorted_22_08.txt"))
plast.deg <- as.character(read_lines("../data/PlasticDB_Prokaryotic_genera.txt"))
plast.HCB <- read_lines("../data/PlasticDB_HCB_genera.txt")
HCB_only <- setdiff(HCB, plast.deg)

Genus_plastic <- Genus %>%  filter(Genus%in%plast.deg)
Genus_oil <-  Genus %>%  filter(Genus%in%HCB)
Genus_plastoil <- Genus %>%  filter(Genus%in%plast.HCB) 
Genus_HCBonly <-  Genus %>%  filter(Genus%in%HCB_only)

length(unique(Genus_oil$Genus))
length(unique(Genus_plastic$Genus))
length(unique(Genus_plastoil$Genus))
length(unique(Genus_HCBonly$Genus))

#unique(Genus_HCBonly$Genus) %>% write_lines("../Analysis/HCB_genera_on_foils.txt")

HCB <- Genus_plastic %>% filter(Genus_rep_rel_abund > 0.0025)
length(unique(HCB$Genus))

HCBs<- Genus_plastoil %>% filter(Genus_rep_rel_abund > 0.0025) #%>% 
  #filter( timepoint %in% c("T6")) %>% filter ( treatment %in% c("no UV"))
plast <- HCBs %>%  filter(Genus%in%plast.deg)
length(unique(HCBs$Genus))

Gens_plastoil <- as.character(Genus_HCBonly$Genus) 
Gens_low_ab <- HCB$Genus %>% as.character()
setdiff(Genus_plastic$Genus, HCB$Genus) %>% write_lines("../Analysis/Low_abundance_PDB_genera.txt")


HCB_top <- Genus_plastic %>% dplyr::select(Description, 
                                       Genus, Genus_rep_rel_abund ) %>% filter(Genus_rep_rel_abund > 0.005) %>% 
  mutate(across(c(Description),factor))%>% distinct() %>% 
  group_by(Description) %>% slice_max(order_by = Genus_rep_rel_abund, n = 10) %>% arrange(desc(Genus_rep_rel_abund))

unique(HCB_top$Genus)

unique(HCBs$Genus)
length(HCBs$Genus)
length(unique(HCBs$Genus))
head(HCBs)

mean_HCB <- sum(mean(HCB$Genus_rep_rel_abund))*100
sd_HCB <- sd(HCB$Genus_rep_rel_abund)*100
mean_HCB
sd_HCB

mean_HCB <- Genus_HCBonly %>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(mean = mean(sum)) *100
sd_HCB <- Genus_HCBonly %>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(sd = sd(sum))*100

mean_plastic <- Genus_plastic %>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(mean = mean(sum)) *100
sd_plastic <- Genus_plastic %>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(sd = sd(sum))*100

mean_plastoil <- Genus_plastoil%>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(mean = mean(sum)) *100
sd_plastoil <- Genus_plastoil %>%  group_by(Description) %>% summarise(sum = sum(Genus_rep_rel_abund)) %>% summarise(sd = sd(sum))*100

Gens <- Genus %>%  filter(Genus == "Alcanivorax" ) %>% filter(Genus_rep_rel_abund > 0)
length(unique(Gens$Description))

Gen <- Gens %>%  filter(timepoint == "T6" ) 
length(unique(Gen$Description))
mean(Gen$Genus_rep_rel_abund)*100
sd(Gen$Genus_rep_rel_abund)*100


## Genera
Genera <-Genus   %>% 
  filter(Genus_rep_rel_abund > 0.0025)
length(unique(Genera$Genus))
 
Genera.avg<- Genus %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd <- Genus %>%  group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T1<- Genus %>%  filter(timepoint == 'T1')%>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T1 <- Genus %>%  filter(timepoint == 'T1')%>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T6<- Genera %>%  filter(timepoint == 'T6') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T6 <- Genera %>%  filter(timepoint == 'T6') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.UV<- Genera %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.UV <- Genera %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.noUV<- Genera %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.noUV <- Genera %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T1.UV<- Genus %>%  filter(timepoint == 'T1') %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T1.UV <- Genus %>%  filter(timepoint == 'T1') %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T1.noUV<- Genus %>%  filter(timepoint == 'T1') %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T1.noUV <- Genus %>%  filter(timepoint == 'T1') %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T6.UV<- Genus %>%  filter(timepoint == 'T6') %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T6.UV <- Genus %>%  filter(timepoint == 'T6') %>% filter(treatment  == 'UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))

Genera.avg.T6.noUV<- Genus %>%  filter(timepoint == 'T6') %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(avg = mean(Genus_rep_rel_abund))
Genera.sd.T6.noUV <- Genus %>%  filter(timepoint == 'T6') %>% filter(treatment  == 'no UV') %>% group_by(Genus) %>%  summarise(sd = sd(Genus_rep_rel_abund))



top12_Genera <- c("Enhydrobacter", "Cutibacterium",                                   
                  "Oleiphilus",  "Delftia",  "Bradyrhizobium" , "Tenacibaculum",          
                  "Staphylococcus" ,   "Erythrobacter",    "Stenotrophomonas" ,  "Aureicoccus" )                                                                   

sum(Genera.avg$avg)*100
Genera.avg.10 <- Genera.avg %>% filter(Genus %in% top12_Genera) 
sum(Genera.avg.10$avg)*100
Genera.sd.10 <- Genera.sd%>% filter(Genus %in% top12_Genera)
sqrt(sum((Genera.sd.10$sd)^2)) *100

avg.HCB <- Genera.avg %>% filter(Genus %in% HCB_only) 
sum(avg.HCB$avg)*100
sd.HCB <- Genera.sd %>% filter(Genus %in% HCB_only)
sqrt(sum((sd.HCB$sd)^2)) *100

avg.PDB <- Genera.avg %>% filter(Genus %in% plast.deg) 
sum(avg.PDB$avg)*100
sd.PDB <- Genera.sd %>% filter(Genus %in% plast.deg)
sqrt(sum((sd.PDB$sd)^2)) *100

avg.plastoil <- Genera.avg.T6 %>% filter(Genus %in% plast.HCB) 
sum(avg.plastoil$avg)*100
sd.plastoil <- Genera.sd.T6 %>% filter(Genus %in% plast.HCB) %>% filter (sd != " ") 
sqrt(sum((sd.plastoil$sd)^2))*100





## Phyla
Phyla <-t2  %>%  dplyr::select(Sample, Description,timepoint, treatment, 
                              Phylum, Phyla_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("T1", "T6"))

top6_phyla <- top_taxa(Phyla, 6)

top6_phyla <- c('Proteobacteria', 'Bacteroidota', 'Actinobacteriota', 'Firmicutes','Verrucomicrobiota', 'Cyanobacteria')

Phyls <- Phyla %>%  filter(Phylum %in% top6_phyla ) %>% filter(Phyla_rep_rel_abund > 0)  %>% group_by(Description)
length(unique(Phyls$Description))
unique(Phyls$Phylum)



Phyls <- Phyla %>%  filter(Phylum %in% top6_phyla ) %>% filter(Phyla_rep_rel_abund > 0)  
length(unique(Phyls$Description))
unique(Phyls$Phylum)

Proteobacteria <- Phyla %>%  filter(Phylum == 'Proteobacteria' ) %>% filter(Phyla_rep_rel_abund > 0)  
Bacteroidota <- Phyla %>%  filter(Phylum == 'Bacteroidota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Actinobacteriota <- Phyla %>%  filter(Phylum == 'Actinobacteriota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Firmicutes <- Phyla %>%  filter(Phylum == 'Firmicutes' ) %>% filter(Phyla_rep_rel_abund > 0)  
Verrucomicrobiota <- Phyla %>%  filter(Phylum == 'Verrucomicrobiota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Cyanobacteria <- Phyla %>%  filter(Phylum == 'Cyanobacteria' ) %>% filter(Phyla_rep_rel_abund > 0)  

sum(mean(Proteobacteria$Phyla_rep_rel_abund), mean(Bacteroidota$Phyla_rep_rel_abund), mean(Actinobacteriota$Phyla_rep_rel_abund), 
    mean(Firmicutes$Phyla_rep_rel_abund), mean(Verrucomicrobiota$Phyla_rep_rel_abund), mean(Cyanobacteria$Phyla_rep_rel_abund)) * 100
sqrt(sum(sd(Proteobacteria$Phyla_rep_rel_abund)^2, sd(Bacteroidota$Phyla_rep_rel_abund)^2, sd(Actinobacteriota$Phyla_rep_rel_abund)^2, 
         sd(Firmicutes$Phyla_rep_rel_abund)^2, sd(Verrucomicrobiota$Phyla_rep_rel_abund)^2, sd(Cyanobacteria$Phyla_rep_rel_abund)^2)) * 100

sqrt(sum(sd(Phyls$Phyla_rep_rel_abund)^2)) * 100

sqrt(sum(var(Proteobacteria$Phyla_rep_rel_abund), var(Bacteroidota$Phyla_rep_rel_abund), var(Actinobacteriota$Phyla_rep_rel_abund), 
         var(Firmicutes$Phyla_rep_rel_abund), var(Verrucomicrobiota$Phyla_rep_rel_abund), var(Cyanobacteria$Phyla_rep_rel_abund))) * 100

sqrt(sum(var(Phyls$Phyla_rep_rel_abund))) * 100

Phyls <- Phyla %>%   filter(Phyla_rep_rel_abund > 0.001)  
length(unique(Phyls$Phylum))
unique(Phyls$Phylum)

mean(Phyls$Phyla_rep_rel_abund) *100
sd(Phyls$Phyla_rep_rel_abund)*100

Phyls <- Phyla %>%  filter(Phylum %in% top6_phyla ) %>% filter(Phyla_rep_rel_abund > 0)  
length(unique(Phyls$Description))
unique(Phyls$Phylum)

Proteobacteria <- Phyla %>%  filter(Phylum == 'Proteobacteria' ) %>% filter(Phyla_rep_rel_abund > 0)  
Bacteroidota <- Phyla %>%  filter(Phylum == 'Bacteroidota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Actinobacteriota <- Phyla %>%  filter(Phylum == 'Actinobacteriota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Firmicutes <- Phyla %>%  filter(Phylum == 'Firmicutes' ) %>% filter(Phyla_rep_rel_abund > 0)  
Verrucomicrobiota <- Phyla %>%  filter(Phylum == 'Verrucomicrobiota' ) %>% filter(Phyla_rep_rel_abund > 0)  
Cyanobacteria <- Phyla %>%  filter(Phylum == 'Cyanobacteria' ) %>% filter(Phyla_rep_rel_abund > 0)  

sum(mean(Proteobacteria$Phyla_rep_rel_abund), mean(Bacteroidota$Phyla_rep_rel_abund), mean(Actinobacteriota$Phyla_rep_rel_abund), 
    mean(Firmicutes$Phyla_rep_rel_abund), mean(Verrucomicrobiota$Phyla_rep_rel_abund), mean(Cyanobacteria$Phyla_rep_rel_abund)) * 100
sqrt(sum(sd(Proteobacteria$Phyla_rep_rel_abund)^2, sd(Bacteroidota$Phyla_rep_rel_abund)^2, sd(Actinobacteriota$Phyla_rep_rel_abund)^2, 
         sd(Firmicutes$Phyla_rep_rel_abund)^2, sd(Verrucomicrobiota$Phyla_rep_rel_abund)^2, sd(Cyanobacteria$Phyla_rep_rel_abund)^2)) * 100

sqrt(sum(sd(Phyls$Phyla_rep_rel_abund)^2)) * 100

sqrt(sum(var(Proteobacteria$Phyla_rep_rel_abund), var(Bacteroidota$Phyla_rep_rel_abund), var(Actinobacteriota$Phyla_rep_rel_abund), 
         var(Firmicutes$Phyla_rep_rel_abund), var(Verrucomicrobiota$Phyla_rep_rel_abund), var(Cyanobacteria$Phyla_rep_rel_abund))) * 100

sqrt(sum(var(Phyls$Phyla_rep_rel_abund))) * 100

Phyls <- Phyla %>%  filter(Phylum == 'Bacteroidota')  %>% filter(Phyla_rep_rel_abund > 0.01)  %>%  filter(timepoint == 'T6')
length(unique(Phyls$Description))
unique(Phyls$Phylum)

mean(Phyls$Phyla_rep_rel_abund) *100
sd(Phyls$Phyla_rep_rel_abund)*100

## Orders
Orders <-t2  %>%  dplyr::select(Sample, Description,timepoint, treatment, 
                               Order, Order_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("T1", "T6"))  %>% filter(Order_rep_rel_abund > 0) 

Order <- Orders   %>% 
  filter(Order_rep_rel_abund > 0.01)
length(unique(Order$Order))

Orders.avg<- Orders %>% group_by(Order) %>%  summarise(avg = mean(Order_rep_rel_abund))
Orders.sd <- Orders %>%  group_by(Order) %>%  summarise(sd = sd(Order_rep_rel_abund))

Orders.avg.T1<- Orders %>%  filter(timepoint == 'T1')%>% group_by(Order) %>%  summarise(avg = mean(Order_rep_rel_abund))
Orders.sd.T1 <- Orders %>%  filter(timepoint == 'T1')%>% group_by(Order) %>%  summarise(sd = sd(Order_rep_rel_abund))

Orders.avg.T6<- Orders %>%  filter(timepoint == 'T6') %>% group_by(Order) %>%  summarise(avg = mean(Order_rep_rel_abund))
Orders.sd.T6 <- Orders %>%  filter(timepoint == 'T6') %>% group_by(Order) %>%  summarise(sd = sd(Order_rep_rel_abund))

Orders.avg.UV<- Orders %>%  filter(treatment  == 'UV')%>% group_by(Order) %>%  summarise(avg = mean(Order_rep_rel_abund))
Orders.sd.UV <- Orders %>%  filter(treatment  == 'UV')%>% group_by(Order) %>%  summarise(sd = sd(Order_rep_rel_abund))

Orders.avg.noUV<- Orders %>%  filter(treatment  == 'no UV')%>% group_by(Order) %>%  summarise(avg = mean(Order_rep_rel_abund))
Orders.sd.noUV <- Orders %>%  filter(treatment  == 'no UV')%>% group_by(Order) %>%  summarise(sd = sd(Order_rep_rel_abund))

top10_orders <- c("Pseudomonadales",     "Flavobacteriales",    "Propionibacteriales", "Rhodobacterales",     "Rhizobiales" ,        
                  "Chitinophagales" ,    "Caulobacterales"  ,   "Burkholderiales"  ,   "Cytophagales"   ,     "Arenicellales")

sum(Orders.avg$avg)*100
Orders.avg.10 <- Orders.avg %>% filter(Order %in% top10_orders) 
sum(Orders.avg.10$avg)*100
Orders.sd.10 <- Orders.sd%>% filter(Order %in% top10_orders)
sqrt(sum((Orders.sd.10$sd)^2)) *100

Orders.avg.T1.10 <- Orders.avg.T1 %>% filter(Order %in% top10_orders) %>% arrange(desc(avg))
Orders.sd.T1.10 <- Orders.sd.T1%>% filter(Order %in% top10_orders)

Orders.avg.T6.10 <- Orders.avg.T6 %>% filter(Order %in% top10_orders) %>% arrange(desc(avg))
Orders.sd.T6.10 <- Orders.sd.T6 %>% filter(Order %in% top10_orders)

Orders.avg.UV.10 <- Orders.avg.UV %>% filter(Order %in% top10_orders) 
Orders.sd.UV.10 <- Orders.sd.UV%>% filter(Order %in% top10_orders)

Orders.avg.noUV.10 <- Orders.avg.noUV%>% filter(Order %in% top10_orders) 
Orders.sd.noUV.10 <- Orders.sd.noUV %>% filter(Order %in% top10_orders)


Ords <- Orders %>%  filter(Order ==  "Arenicellales") %>%  filter(timepoint == 'T1') %>% filter(Order_rep_rel_abund > 0)  
length(unique(Ords$Description))
unique(Ords$Order)

mean(Ords$Order_rep_rel_abund) *100
sd(Ords$Order_rep_rel_abund)*100

# Kingdom
Kingdom <-t2  %>%  dplyr::select(Sample, Description,timepoint, treatment, 
                               Kingdom,Kingdom_rep_rel_abund)  %>% 
  distinct() %>%   filter ( timepoint %in% c("T1", "T6"))  %>% filter(Kingdom_rep_rel_abund > 0) 

archeae.t1 <- Kingdom %>%  filter(Kingdom == 'Archaea') %>%   filter(timepoint == 'T1')
archeae.t6 <- Kingdom %>%  filter(Kingdom == 'Archaea') %>%  filter(timepoint == 'T6')

length(unique(archeae.t6$Description))

mean(archeae.t1$Kingdom_rep_rel_abund) *100
sd(archeae.t1$Kingdom_rep_rel_abund)*100

mean(archeae.t6$Kingdom_rep_rel_abund) *100
sd(archeae.t6$Kingdom_rep_rel_abund)*100

archeae.t6.avg <- archeae.t6 %>% group_by(Description) %>%  summarise(avg = mean(Kingdom_rep_rel_abund)*100) 
archeae.t6.sd <- archeae.t6 %>% group_by(Description) %>%  summarise(sd = sd(Kingdom_rep_rel_abund)*100) 
archeae.t6.avg %>% top_n(1,avg)

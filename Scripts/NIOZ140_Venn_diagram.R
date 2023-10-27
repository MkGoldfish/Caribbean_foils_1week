##%#_______________________________________________________________________________________#%###
####                 Working Directory                                                      ####       
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")


###%#_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
###%#_______________________________________________________________________________________#%###
library(devtools)
library("ggplot2")
library("ggVennDiagram")
library("stringr")
library("tidyverse")
library("ggvenn")
library(rvg)
library(officer)
library("MicEco")
library(nVennR)
library(eulerr)
library(cowplot)

set.seed(42)
###%#_______________________________________________________________________________________#%###
####                  Import Data                                                            ####
###%#_______________________________________________________________________________________#%###
tax <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/asv_table.txt', row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- read.delim('../../../Data_corrected_barcodes/metadata/mapping_file_details.txt', row.names = 1, na.strings = c("")) %>% filter(timepoint %in% c("T1", "T6"))


t1 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c("")) 
t2 <- t1 %>% filter( timepoint %in% c("T1", "T6")) %>% select(OTU:Description, Order, Genus)

time_UV <- str_c(t2$timepoint, "_", t2$treatment)
t2 <- t2 %>%
  add_column(time_UV)

treat_time <- str_c(t2$treatment, "_", t2$timepoint)
t2 <- t2 %>%
  add_column(treat_time)

time_mat <- str_c(t2$timepoint, "_", t2$Material)
t2 <- t2 %>%
  add_column(time_mat)

t2 <-t2 %>% filter(Abundance > 1)

## add category of plastic for nesting in new column to dataframe
t2["category"] <- NA

for (i in 1:nrow(t2)){
  
  if (t2$Material[i] %in% c("Nylon","PET")) {
    t2$category[i] <- "heteratomes"
  }
  
  else 
    
    if (t2$Material[i] %in% c("PE","PP","PS")) {
      t2$category[i] <- "C_backbone"
    }
}

####Creating Phyloseq object for use with MicEco ####
map$LinkerPrimerSequence <- NULL
map$ReversePrimerSequence <- NULL
map$InputFileName <- NULL
map$BarcodeSequence <- NULL
map$BarcodeSequence_1 <- NULL

time_UV <- str_c(map$timepoint, "_", map$treatment)
map <- map %>%
  add_column(time_UV)

treat_time <- str_c(map$treatment, "_", map$timepoint)
map <- map %>%
  add_column(treat_time)

time_mat <- str_c(map$timepoint, "_", map$Material)
map <- map %>%
  add_column(time_mat)

## add category of plastic for nesting in new column to dataframe
map["category"] <- NA

for (i in 1:nrow(t2)){
  
  if (map$Material[i] %in% c("Nylon","PET")) {
    map$category[i] <- "heteratomes"
  }
  
  else 
    
    if (map$Material[i] %in% c("PE","PP","PS")) {
      map$category[i] <- "C_backbone"
      
    }
  
}

map <- map %>% mutate(timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6" ,
)
)

map <- sample_data(map)
physeq_object = merge_phyloseq(otu, tax, map)
physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90")
sub1<- subset_samples(physeq_object, timepoint %in% c("Day 1", "Day 6"))
sub2 <- subset_samples(physeq_object, surface!="negative_c") 
physeq_object <- merge_phyloseq(sub1, sub2)

### Create subset dataframes for plotting ####
##Pull out only the ASV IDs to character vector
colnames(t2)

T1 <- t2 %>% filter(timepoint =="T1")
T6 <- t2 %>% filter(timepoint =="T6")
UV <-t2 %>% filter(treatment =="UV")
noUV <-t2 %>% filter(treatment =="no UV")
PE <- t2 %>% filter(Material =="PE")
PP <- t2 %>% filter(Material =="PP")
PS <- t2 %>% filter(Material =="PS")
PET <- t2 %>% filter(Material =="PET")
Nylon <- t2 %>% filter(Material =="Nylon")

sub_t1 <- t2 %>% filter(timepoint =="T1") %>% pull(OTU) %>% unique()
sub_t6 <- t2 %>% filter(timepoint =="T6") %>% pull(OTU) %>% unique()
sub_hetero <-t2 %>% filter(category =="heteratomes") %>% pull(OTU)
sub_cback <-t2 %>% filter(category =="C_backbone") %>% pull(OTU)
sub_UV <-t2 %>% filter(treatment =="UV") %>% pull(OTU)
sub_noUV <-t2 %>% filter(treatment =="no UV") %>% pull(OTU)

sub_T1.UV <- T1 %>% filter(treatment =="UV") %>% pull(OTU)
sub_T1.noUV <-  T1 %>% filter(treatment =="no UV") %>% pull(OTU)
sub_T6.UV <- T6 %>% filter(treatment =="UV") %>% pull(OTU) 
sub_T6.noUV <- T6 %>% filter(treatment =="no UV") %>% pull(OTU) 

sub_het.UV <-UV %>% filter(category =="heteratomes") %>% pull(OTU)
sub_het.noUV <-noUV %>% filter(category =="heteratomes")  %>% pull(OTU)
sub_c.UV <-UV %>% filter(category =="C_backbone") %>%  pull(OTU)
sub_c.noUV <-noUV %>% filter(category =="C_backbone") %>% pull(OTU)

sub_t1.het <- T1 %>% filter(category =="heteratomes") %>% pull(OTU)
sub_t1.c <- T1 %>% filter(category =="C_backbone") %>% pull(OTU)
sub_t6.het <- T6 %>% filter(category =="heteratomes") %>% pull(OTU)
sub_t6.c <- T6 %>% filter(category =="C_backbone") %>% pull(OTU)

sub_PE <- t2 %>% filter(Material =="PE") %>% pull(OTU)
sub_PP <- t2 %>% filter(Material =="PP") %>% pull(OTU)
sub_PS <- t2 %>% filter(Material =="PS")%>% pull(OTU)
sub_PET <- t2 %>% filter(Material =="PET")%>% pull(OTU)
sub_Nylon <- t2 %>% filter(Material =="Nylon")%>% pull(OTU)

T1_PE <- T1 %>% filter(Material =="PE") %>% pull(OTU) %>% unique() %>% length()
T1_PP <-T1 %>% filter(Material =="PP") %>% pull(OTU) %>% unique() %>% length()
T1_PS <-T1 %>% filter(Material =="PS")%>% pull(OTU) %>% unique() %>% length()
T1_PET <- T1 %>% filter(Material =="PET")%>% pull(OTU) %>% unique() %>% length()
T1_Nylon <- T1 %>% filter(Material =="Nylon")%>% pull(OTU) %>% unique() %>% length()

T6_PE <- T6 %>% filter(Material =="PE") %>% pull(OTU) %>% unique() %>% length()
T6_PP<- T6 %>% filter(Material =="PP") %>% pull(OTU) %>% unique() %>% length()
T6_PS<- T6 %>% filter(Material =="PS")%>% pull(OTU) %>% unique() %>% length()
T6_PET<- T6 %>% filter(Material =="PET")%>% pull(OTU) %>% unique() %>% length()
T6_Nylon <- T6 %>% filter(Material =="Nylon")%>% pull(OTU) %>% unique() %>% length()

UV_PE <- PE %>% filter(treatment =="UV") %>% pull(OTU)
UV_PP <-PP %>% filter(treatment =="UV")%>% pull(OTU)
UV_PS <-PS %>%filter(treatment =="UV")%>% pull(OTU)
UV_PET <- PET %>% filter(treatment =="UV")%>% pull(OTU)
UV_Nylon <- Nylon %>% filter(treatment =="UV")%>% pull(OTU)

noUV_PE <- PE %>% filter(treatment =="no UV") %>% pull(OTU)
noUV_PP<- PP %>% filter(treatment =="no UV") %>% pull(OTU)
noUV_PS<- PS %>% filter(treatment =="no UV")%>% pull(OTU)
noUV_PET<- PET %>% filter(treatment =="no UV")%>% pull(OTU)
noUV_Nylon <- Nylon %>% filter(treatment =="no UV")%>% pull(OTU)

### create lists for plotting the venns ####
Time <- list("T1"= sub_t1, "T6" = sub_t6)
Category<- list("Heteroatomes"= sub_hetero, "C_backbone" = sub_cback)
Treatment <- list("UV" = sub_UV, "no UV" = sub_noUV)

Time_treat <- list("T1.UV" = sub_T1.UV, "T1.no UV" = sub_T1.noUV, "T6.UV" =sub_T6.UV,  "T6.noUV" = sub_T6.noUV)
T1_treat <- list("T1.UV" = sub_T1.UV, "T1.no UV" = sub_T1.noUV)
T6_treat <- list("T6.UV" =sub_T6.UV,  "T6.noUV" = sub_T6.noUV)
UV_time <- list("T1.UV" = sub_T1.UV, "T6.UV" = sub_T1.noUV)
noUV_time <- list("T1.noUV" =sub_T6.UV,  "T6.noUV" = sub_T6.noUV)

Time_cat <- list("T1.Het" = sub_t1.het, "T1.C" = sub_t1.c, "T6.Het" =sub_t6.het,  "T6.C" = sub_t6.c) 
Time_C <- list( "T1.C" = sub_t1.c,  "T6.C" = sub_t6.c) 
Time_het <- list("T1.Het" = sub_t1.het,  "T6.Het" =sub_t6.het) 

Cat_treat <- list("Het.UV" = sub_het.UV, "Het.no UV" = sub_het.noUV, "Cback.UV" =sub_c.UV,  "Cback.noUV" = sub_c.noUV )
Het_treat <- list("Het.UV" = sub_het.UV, "Het.no UV" = sub_het.noUV )
C_treat <- list("Cback.UV" =sub_c.UV,  "Cback.noUV" = sub_c.noUV )

Materials <- list("PE" =sub_PE,  "PP" = sub_PP, "PS" = sub_PS, "PET" = sub_PET, "Nylon" = sub_Nylon )
T1_Materials <- list("PE" =T1_PE,  "PP" = T1_PP,  "PS" = T1_PS, "PET" = T1_PET,"Nylon" = T1_Nylon )
T6_Materials <- list("PE" =T6_PE,  "PP" = T6_PP,  "PS" = T6_PS, "PET" = T6_PET, "Nylon" = T6_Nylon )

PE_time <- list("PE T1" = T1_PE, "PE T6" = T6_PE)
PP_time<- list("PP T1" = T1_PP, "PP T6" = T6_PP)
PS_time<- list("PS T1" = T1_PS, "PS T6" = T6_PS)
PET_time<- list("PET T1" = T1_PET, "PET T6" = T6_PET)
Nylon_time<- list("Nylon T1" = T1_Nylon, "Nylon T6" = T6_Nylon)

PE_treatment <- list("PE UV" = UV_PE, "PE no UV" = noUV_PE)
PP_treatment<- list("PP UV" = UV_PP, "PP no UV" = noUV_PP)
PS_treatment<- list("PS UV" = UV_PS, "PS no UV" = noUV_PS)
PET_treatment<- list("PET UV" = UV_PET, "PET no UV" = noUV_PET)
Nylon_treatment<- list("Nylon UV" = UV_Nylon, "Nylon no UV" = noUV_Nylon)

### Plot the venns and store in .pptx ####
Venndiagrams.Foils<- read_pptx()
Venndiagrams.Foils <- read_pptx("../Reports/Venndiagrams.Foils.abund+1_202213.pptx")

Venn <- ggvenn(Time) + labs( title = "Venn Foils ASVs")

Venn

editable_graph <- dml(ggobj = Venn)
Venndiagrams.Foils <- add_slide(Venndiagrams.Foils) 
Venndiagrams.Foils<- ph_with(x = Venndiagrams.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Venndiagrams.Foils, target = "../Reports/Venndiagrams.Foils.abund+1_202213.pptx")

### Use ggVennDiagram for polymers, since this has the option for more than 4 groups, ggvenn does not ####

Venn <- ggVennDiagram(T1_Materials[4:5]) + labs( title = "Venn Foils T1 ASV")

Venn

editable_graph <- dml(ggobj = Venn)
Venndiagrams.Foils <- add_slide(Venndiagrams.Foils) 
Venndiagrams.Foils<- ph_with(x = Venndiagrams.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Venndiagrams.Foils, target = "../Reports/Venndiagrams.Foils.abund+1_202213.pptx")



### Lists and Venns on Genus Level ####
t2 <- t2 %>% filter(!Genus =="NA")
colnames(t2)
T1 <- t2 %>% filter(timepoint =="T1")
T6 <- t2 %>% filter(timepoint =="T6")
UV <-t2 %>% filter(treatment =="UV")
noUV <-t2 %>% filter(treatment =="no UV")
PE <- t2 %>% filter(Material =="PE")
PP <- t2 %>% filter(Material =="PP")
PS <- t2 %>% filter(Material =="PS")
PET <- t2 %>% filter(Material =="PET")
Nylon <- t2 %>% filter(Material =="Nylon")

sub_t1 <- t2 %>% filter(timepoint =="T1") %>% pull(Genus) %>% unique()
sub_t6 <- t2 %>% filter(timepoint =="T6") %>% pull(Genus)  %>% unique()
sub_hetero <-t2 %>% filter(category =="heteratomes") %>% pull(Genus)  %>% unique()
sub_cback <-t2 %>% filter(category =="C_backbone") %>% pull(Genus)  %>% unique()
sub_UV <-t2 %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_noUV <-t2 %>% filter(treatment =="no UV") %>% pull(Genus)  %>% unique()

sub_T1.UV <- T1 %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_T1.noUV <-  T1 %>% filter(treatment =="no UV") %>% pull(Genus)  %>% unique()
sub_T6.UV <- T6 %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_T6.noUV <- T6 %>% filter(treatment =="no UV") %>% pull(Genus)  %>% unique()

sub_het.UV <-UV %>% filter(category =="heteratomes") %>% pull(Genus)  %>% unique()
sub_het.noUV <-noUV %>% filter(category =="heteratomes")  %>% pull(Genus)  %>% unique()
sub_c.UV <-UV %>% filter(category =="C_backbone") %>%  pull(Genus)  %>% unique()
sub_c.noUV <-noUV %>% filter(category =="C_backbone") %>% pull(Genus) %>% unique()

sub_t1.het <- T1 %>% filter(category =="heteratomes") %>% pull(Genus)  %>% unique()
sub_t1.c <- T1 %>% filter(category =="C_backbone") %>% pull(Genus)  %>% unique()
sub_t6.het <- T6 %>% filter(category =="heteratomes") %>% pull(Genus)  %>% unique()
sub_t6.c <- T6 %>% filter(category =="C_backbone") %>% pull(Genus)  %>% unique()

sub_t1.het.UV <- T1 %>% filter(category =="heteratomes") %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_t1.het.noUV <- T1 %>% filter(category =="heteratomes") %>% filter(treatment =="no UV")%>% pull(Genus)  %>% unique()
sub_t1.c.UV <- T1 %>% filter(category =="C_backbone") %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_t1.c.noUV <- T1 %>% filter(category =="C_backbone")%>% filter(treatment =="no UV") %>% pull(Genus)  %>% unique()

sub_t6.het.UV <- T6 %>% filter(category =="heteratomes") %>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_t6.het.noUV <- T6 %>% filter(category =="heteratomes") %>% filter(treatment =="no UV")%>% pull(Genus)  %>% unique()
sub_t6.c.UV <- T6 %>% filter(category =="C_backbone")%>% filter(treatment =="UV") %>% pull(Genus)  %>% unique()
sub_t6.c.noUV <- T6 %>% filter(category =="C_backbone")%>% filter(treatment =="no UV") %>% pull(Genus)  %>% unique()

sub_PE <- t2 %>% filter(Material =="PE") %>% pull(Genus)  %>% unique()
sub_PP <- t2 %>% filter(Material =="PP") %>% pull(Genus)  %>% unique()
sub_PS <- t2 %>% filter(Material =="PS")%>% pull(Genus)  %>% unique()
sub_PET <- t2 %>% filter(Material =="PET")%>% pull(Genus)  %>% unique()
sub_Nylon <- t2 %>% filter(Material =="Nylon")%>% pull(Genus) %>% unique()

T1_PE <- T1 %>% filter(Material =="PE") %>% pull(Genus) %>% unique()
T1_PP <-T1 %>% filter(Material =="PP") %>% pull(Genus) %>% unique()
T1_PS <-T1 %>% filter(Material =="PS")%>% pull(Genus) %>% unique()
T1_PET <- T1 %>% filter(Material =="PET")%>% pull(Genus) %>% unique()
T1_Nylon <- T1 %>% filter(Material =="Nylon")%>% pull(Genus) %>% unique()

T6_PE <- T6 %>% filter(Material =="PE") %>% pull(Genus) %>% unique()
T6_PP<- T6 %>% filter(Material =="PP") %>% pull(Genus) %>% unique()
T6_PS<- T6 %>% filter(Material =="PS")%>% pull(Genus) %>% unique()
T6_PET<- T6 %>% filter(Material =="PET")%>% pull(Genus) %>% unique()
T6_Nylon <- T6 %>% filter(Material =="Nylon")%>% pull(Genus) %>% unique()

UV_PE <- PE %>% filter(treatment =="UV") %>% pull(Genus) %>% unique()
UV_PP <-PP %>% filter(treatment =="UV")%>% pull(Genus) %>% unique()
UV_PS <-PS %>%filter(treatment =="UV")%>% pull(Genus) %>% unique()
UV_PET <- PET %>% filter(treatment =="UV")%>% pull(Genus) %>% unique()
UV_Nylon <- Nylon %>% filter(treatment =="UV")%>% pull(Genus) %>% unique()

noUV_PE <- PE %>% filter(treatment =="no UV") %>% pull(Genus) %>% unique()
noUV_PP<- PP %>% filter(treatment =="no UV") %>% pull(Genus) %>% unique()
noUV_PS<- PS %>% filter(treatment =="no UV")%>% pull(Genus) %>% unique()
noUV_PET<- PET %>% filter(treatment =="no UV")%>% pull(Genus) %>% unique()
noUV_Nylon <- Nylon %>% filter(treatment =="no UV")%>% pull(Genus) %>% unique()

##### create lists for plotting the venns ####
Time <- list("T1"= sub_t1, "T6" = sub_t6)
Category<- list("Heteroatomes"= sub_hetero, "C_backbone" = sub_cback)
Treatment <- list("UV" = sub_UV, "no UV" = sub_noUV)

Time_treat <- list("T1.UV" = sub_T1.UV, "T1.no UV" = sub_T1.noUV, "T6.UV" =sub_T6.UV,  "T6.noUV" = sub_T6.noUV)
T1_treat <- list("T1.UV" = sub_T1.UV, "T1.no UV" = sub_T1.noUV)
T6_treat <- list("T6.UV" =sub_T6.UV,  "T6.noUV" = sub_T6.noUV)
UV_time <- list("T1.UV" = sub_T1.UV, "T6.UV" = sub_T6.UV)
noUV_time <- list("T1.noUV" =sub_T1.noUV,  "T6.noUV" = sub_T6.noUV)

Time_cat <- list("T1.Het" = sub_t1.het, "T1.C" = sub_t1.c, "T6.Het" =sub_t6.het,  "T6.C" = sub_t6.c) 
T1_cat <- list("T1.Het" = sub_t1.het, "T1.C" = sub_t1.c) 
T6_cat <- list("T6.Het" =sub_t6.het,  "T6.C" = sub_t6.c) 
Time_C <- list( "T1.C" = sub_t1.c,  "T6.C" = sub_t6.c) 
Time_het <- list("T1.Het" = sub_t1.het,  "T6.Het" =sub_t6.het) 

T1_cat_treat <-list("t1.het.UV" = sub_t1.het.UV, "t1.het.noUV" = sub_t1.het.noUV, "t1.c.UV" = sub_t1.c.UV, "t1.c.noUV" = sub_t1.c.noUV)
T6_cat_treat <-list("t6.het.UV" = sub_t6.het.UV, "t6.het.noUV" = sub_t6.het.noUV, "t6.c.UV" = sub_t6.c.UV, "t6.c.noUV" = sub_t6.c.noUV)
T1_C_treat <-list( "t1.c.UV" = sub_t1.c.UV, "t1.c.noUV" = sub_t1.c.noUV) 
T6_C_treat <-list("t6.c.UV" = sub_t6.c.UV, "t6.c.noUV" = sub_t6.c.noUV)
T1_H_treat<-list("t1.het.UV" = sub_t1.het.UV, "t1.het.noUV" = sub_t1.het.noUV)
T6_H_treat <-list("t6.het.UV" = sub_t6.het.UV, "t6.het.noUV" = sub_t6.het.noUV)

Cat_treat <- list("Het.UV" = sub_het.UV, "Het.no UV" = sub_het.noUV, "Cback.UV" =sub_c.UV,  "Cback.noUV" = sub_c.noUV )
Het_treat <- list("Het.UV" = sub_het.UV, "Het.no UV" = sub_het.noUV )
C_treat <- list("Cback.UV" =sub_c.UV,  "Cback.noUV" = sub_c.noUV )

Materials <- list("PE" =sub_PE,  "PP" = sub_PP, "PS" = sub_PS, "PET" = sub_PET, "Nylon" = sub_Nylon )
T1_Materials <- list("PE" =T1_PE,  "PP" = T1_PP,  "PS" = T1_PS, "PET" = T1_PET,"Nylon" = T1_Nylon )
T6_Materials <- list("PE" =T6_PE,  "PP" = T6_PP,  "PS" = T6_PS, "PET" = T6_PET, "Nylon" = T6_Nylon )

PE_time <- list("PE T1" = T1_PE, "PE T6" = T6_PE)
PP_time<- list("PP T1" = T1_PP, "PP T6" = T6_PP)
PS_time<- list("PS T1" = T1_PS, "PS T6" = T6_PS)
PET_time<- list("PET T1" = T1_PET, "PET T6" = T6_PET)
Nylon_time<- list("Nylon T1" = T1_Nylon, "Nylon T6" = T6_Nylon)

PE_treatment <- list("PE UV" = UV_PE, "PE no UV" = noUV_PE)
PP_treatment<- list("PP UV" = UV_PP, "PP no UV" = noUV_PP)
PS_treatment<- list("PS UV" = UV_PS, "PS no UV" = noUV_PS)
PET_treatment<- list("PET UV" = UV_PET, "PET no UV" = noUV_PET)
Nylon_treatment<- list("Nylon UV" = UV_Nylon, "Nylon no UV" = noUV_Nylon)

##### Plot the venns and store in .pptx ####
Venndiagrams.Foils<- read_pptx()
Venndiagrams.Foils <- read_pptx("../Reports/Venndiagrams.Foils.Genera.abund+1.unique_202213.pptx")

Venn <- ggvenn(T6_H_treat) + labs( title = "Venn Foils Genuslevel")

Venn

editable_graph <- dml(ggobj = Venn)
Venndiagrams.Foils <- add_slide(Venndiagrams.Foils) 
Venndiagrams.Foils<- ph_with(x = Venndiagrams.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Venndiagrams.Foils, target = "../Reports/Venndiagrams.Foils.Genera.abund+1.unique_202213.pptx")

##Use ggVennDiagram for polymers, since this has the option for more than 4 groups, ggvenn does not
Venn <- ggVennDiagram(Materials[4:5]) + labs( title = "Venn Foils Genuslevel")

Venn

editable_graph <- dml(ggobj = Venn)
Venndiagrams.Foils <- add_slide(Venndiagrams.Foils) 
Venndiagrams.Foils<- ph_with(x = Venndiagrams.Foils, editable_graph,location = ph_location_type(type = "body") )
print(Venndiagrams.Foils, target = "../Reports/Venndiagrams.Foils.Genera.abund+1.unique_202213.pptx")

#### We can use the package eulerr to create Venns ####
venn.time <- plot(euler(Time), fraction = 0, 
                  weight = F, plot = T, relative = F, 
                  fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                  edges = list(col = "white", lwd = 2), 
                  labels = list(labels =  c("day 1", "day 6"),fontsize = 15, col = "white"),
                  quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                  legend = F)

venn.time

venn.treat <- plot(euler(Treatment), fraction = 0, 
                       weight = F, relative = F, 
                       fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                       edges = list(col = "white", lwd = 2), 
                       labels = list(fontsize = 15, col = "white"),                       
                       quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                       legend = F)
venn.treat

venn.cat <-  plot(euler(Category), fraction = 0, 
                     weight = F, relative = F, 
                     fills = list(fill = c("black", "grey90")),
                     edges = list(col = "white", lwd = 2), 
                  labels = list( labels =  c("Carbon\n backbone", "Heteroatom\n backbone"), fontsize = 13, col = c("white", "black", "white")),
                     quantities = list(type=c ("percent", "counts"), fontsize = 11, col = c("white", "black", "white"), fontface = "bold"),
                     legend = F)

venn.cat

venn.t6.treat <- plot(euler(T6_treat), fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                          edges = list(col = "white", lwd = 2), 
                          labels = list(labels =  c("day 6\n UV", "day 6\nno UV"), fontsize = 15, col = "white"),                         
                          quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                          legend = F)

venn.t6.treat

venn.t1.treat <- plot(euler(T1_treat), fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                          edges = list(col = "white", lwd = 2), 
                          labels = list(labels =  c("day 1\n UV", "day 1\nno UV"), fontsize = 15, col = "white"),
                          quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                          legend = F)

venn.t1.treat

venn.UV.time <- plot(euler(UV_time), fraction = 0, 
                         weight = F, relative = F, 
                         fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                         edges = list(col = "white", lwd = 2), 
                         labels = list(labels =  c("UV\n day 1", "UV\n day 6"), fontsize = 15, col = "white"),
                         quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                         legend = F)

venn.UV.time

venn.noUV.time <- plot(euler(noUV_time), fraction = 0, 
                           weight = F, relative = F, 
                           fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                           edges = list(col = "white", lwd = 2), 
                           labels = list(labels =  c("noUV\nday1", "noUV\nday6"), fontsize = 15, col = "white"),
                           quantities = list(type=c ("percent", "counts"), fontsize = 12, col = "white", fontface = "bold"),
                           legend = F)

venn.noUV.time

venn.t6.cat <- plot(euler(T6_cat), fraction = 0, 
                    weight = F, relative = F, 
                    fills = list(fill = c("black", "grey90")),
                    edges = list(col = "white", lwd = 2), 
                    labels = list( labels =  c("day 6 C\n backbone", "day 6 H\n backbone"), fontsize = 13, col = c("white", "black", "white")),
                    quantities = list(type=c("percent", "counts"), fontsize = 11, col = c("white", "black", "white"), fontface = "bold"),
                    legend = F)

venn.t6.cat

venn.t1.cat <- plot(euler(T1_cat),  fraction = 0, 
                    weight = F, relative = F, 
                    fills = list(fill = c("black", "grey90")),
                    edges = list(col = "white", lwd = 2), 
                    labels = list( labels =  c("day 1 C\n backbone", "day 1 H\n backbone"), fontsize = 13, col = c("white", "black", "white")),
                    quantities = list(type=c("percent", "counts"), fontsize = 11, col = c("white", "black", "white"), fontface = "bold"),
                    legend = F)

venn.t1.cat

venn.cback.time <- plot(euler(Time_C), fraction = 0, 
                        weight = F, relative = F, 
                        fills = list(fill = c("black", "grey90")),
                        edges = list(col = "white", lwd = 2), 
                        labels = list( labels =  c("Carbon\n backbone\n day 1", "Carbon\n backbone\n day 6"), fontsize = 13, col = c("white", "black", "white")),
                        quantities = list(type=c("percent", "counts"),  col = c("white", "black", "white"), fontface = "bold"),
                        legend = F)

venn.cback.time

venn.het.time <- plot(euler(Time_het),  fraction = 0, 
                      weight = F, relative = F, 
                      fills = list(fill = c("black", "grey90")),
                      edges = list(col = "white", lwd = 2), 
                      labels = list( labels =  c("H\n backbone\n day 1", "H backbone\n day 6"), fontsize = 13, col = c("white", "black", "white")),
                      quantities = list(type=c("percent", "counts"), fontsize = 11, col = c("white", "black", "white"), fontface = "bold"),
                      legend = F)

venn.het.time

venn.t6.pols <- plot(euler(T6_Materials), fraction = 0, 
                         weight = F, relative = F, 
                         fills = list(fill = c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )),
                         edges = list(col = "white", lwd = 2), 
                         labels = list( fontsize = 13, col = "white"),
                         quantities = list(type=c ("percent", "counts"), fontsize = 11, col = "white", fontface = "bold"),
                         legend = F)

venn.t6.pols

venn.t1.pols <- plot(euler(T1_Materials),  fraction = 0, 
                         weight = F, relative = F, 
                         fills = list(fill = c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )),
                         edges = list(col = "white", lwd = 2), 
                         labels = list( fontsize = 13, col = "white"),
                        shape = "ellipse", 
                         quantities = list(type=c ("percent", "counts"), fontsize = 11, col = "white", fontface = "bold"),
                         legend = F)

venn.t1.pols

plot_grid(venn.time,
          venn.UV.time,
          venn.noUV.time,
          venn.treat,
          venn.t1.treat,
          venn.t6.treat, 
          venn.cat,
          venn.cback.time,
          venn.het.time,
          ncol = 3,
          align = 'v',
          axis = "tbrl",
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L"),
          label_size = 15,
          rel_heights = c(1,1,1,1.8),
          rel_widths = c(1,1,1))


#### Use MicEco to create Euler diagram from physeq object ####
colnames(sample_data(physeq_object))
genus_physeq <- phyloseq::tax_glom(physeq_object, taxrank="Genus", NArm = TRUE)
ntaxa(genus_physeq)

#Make data subset
ps_day1 <- subset_samples(genus_physeq, timepoint %in% c("Day 1"))
ps_day6 <- subset_samples(genus_physeq, timepoint %in% c("Day 6"))
ps_UV <- subset_samples(genus_physeq, treatment %in% c("UV"))
ps_noUV <- subset_samples(genus_physeq, treatment %in% c("no UV"))
ps_het <- subset_samples(genus_physeq, category %in% c("heteratomes"))
ps_cbac <- subset_samples(genus_physeq, category %in% c("C_backbone"))


Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

HCB.col <- c("#EE7733")
PDB.col <- c("#0077BB")

venn.time <- ps_euler(physeq_object, "timepoint", fraction = 0, 
                      weight = F, plot = T, relative = F, 
                      fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                      edges = list(col = "white", lwd = 3), 
                      labels = list(fontsize = 18, col = "white"),
                      shape = "ellipse", 
                      quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                      legend = F)

venn.time

venn.treat <- ps_euler(genus_physeq, "treatment", fraction = 0, 
                      weight = F, relative = F, 
                      fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                      edges = list(col = "white", lwd = 3), 
                      labels = list(fontsize = 18, col = "white"),
                      shape = "ellipse", 
                      quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                      legend = F)
venn.treat

venn.cat <- ps_euler(genus_physeq, "category", fraction = 0, 
                       weight = F, relative = F, 
                       fills = list(fill = c("black", "grey90")),
                       edges = list(col = "white", lwd = 3), 
                       labels = list( labels =  c("Carbon\n backbone", "Heteroatom\n backbone"), fontsize = 18, col = c("white", "black", "white")),
                       shape = "ellipse", 
                       quantities = list(type=c ("percent", "counts"), fontsize = 15, col = c("white", "black", "white"), fontface = "bold"),
                       legend = F)

venn.cat

venn.t6.treat <- ps_euler(ps_day6, "time_UV", fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                          edges = list(col = "white", lwd = 3), 
                          labels = list(labels =  c("T6 noUV", "T6 UV"), fontsize = 18, col = "white"),
                          shape = "ellipse", 
                          quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                          legend = F)

venn.t6.treat

venn.t1.treat <- ps_euler(ps_day1, "time_UV", fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#999933", "#CC6677", "slateblue1")),
                          edges = list(col = "white", lwd = 3), 
                          labels = list(labels =  c("T1 noUV", "T1 UV"), fontsize = 18, col = "white"),
                          shape = "ellipse", 
                          quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                          legend = F)

venn.t1.treat

venn.UV.time <- ps_euler(ps_UV, "time_UV", fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                          edges = list(col = "white", lwd = 3), 
                          labels = list(labels =  c("T1 UV", "T6 UV"), fontsize = 18, col = "white"),
                          shape = "ellipse", 
                          quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                          legend = F)

venn.UV.time

venn.noUV.time <- ps_euler(ps_noUV, "time_UV", fraction = 0, 
                         weight = F, relative = F, 
                         fills = list(fill = c("#44AA99", "#882255", "#ee8866" )),
                         edges = list(col = "white", lwd = 3), 
                         labels = list(labels =  c("T1 noUV", "T6 noUV"), fontsize = 18, col = "white"),
                         shape = "ellipse", 
                         quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                         legend = F)

venn.noUV.time

venn.t6.pols <- ps_euler(ps_day6, "Material", fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )),
                          edges = list(col = "white", lwd = 3), 
                          labels = list( fontsize = 18, col = "white"),
                          shape = "ellipse", 
                          quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                          legend = F)

venn.t6.pols

venn.t1.pols <- ps_euler(ps_day1, "Material", fraction = 0, 
                          weight = F, relative = F, 
                          fills = list(fill = c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )),
                          edges = list(col = "white", lwd = 3), 
                          labels = list( fontsize = 18, col = "white"),
                          shape = "ellipse", 
                          quantities = list(type=c ("percent", "counts"), fontsize = 15, col = "white", fontface = "bold"),
                          legend = F)

venn.t1.pols



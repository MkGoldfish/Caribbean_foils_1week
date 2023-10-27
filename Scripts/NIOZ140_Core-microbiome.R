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
####                          #%###
####               Corrected Primers                    ####
#                                                          #
##%######################################################%##
#########################################################%###
##  version.string R version 4.2.1 (2022-06-23 ucrt)
##  nickname       Funny-Looking Kid   


##%######################################################%##
#                                                          #
####           Core community analysis            ####
#                                                          #
##%######################################################%##
BiocManager::install("MicrobiotaProcess")

####%#_______________________________________________________________________________________#%###
###                 Working Directory                ####
####%#_______________________________________________________________________________________#%###
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")

###%#_______________________________________________________________________________________#%###
####                 Load libraries       ####                                               
###%#_______________________________________________________________________________________#%###
library(devtools)
library(phyloseq)
library(microbiome)
library(grid)
library(tidyverse)
library("ggh4x")
library(ggpubr)
library(viridis)
library(ggvenn)
library("ggVennDiagram")

###%#_______________________________________________________________________________________#%###
####                  Import Data &                    ####
###%#_______________________________________________________________________________________#%###
tax <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- as.matrix(read.delim('../../../Data_corrected_barcodes/asv/asv_table.txt', row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('../../../Data_corrected_barcodes/metadata/mapping_file_details.txt', row.names = 1, na.strings = c("")))

set.seed(42)

## add new columns to map
map["cat_uv"] <- str_c(map$category,"_",map$treatment)
map["time_UV"] <- str_c(map$timepoint, "_", map$treatment)
map["treat_time"] <- str_c(map$treatment, "_", map$timepoint)
map["time_mat"] <- str_c(map$timepoint, "_", map$Material)

map$LinkerPrimerSequence <- NULL
map$ReversePrimerSequence <- NULL
map$InputFileName <- NULL
map$BarcodeSequence <- NULL
map$BarcodeSequence_1 <- NULL

colnames(map)

map["category"] <- NA
for (i in 1:nrow(map)){

  if (map$Material[i] %in% c("Nylon","PET")) {
    map$category[i] <- "Hetero"
  }
  else
    if (map$Material[i] %in% c("PE","PP","PS")) {
      map$category[i] <- "Carbon"
    }
}

map <- sample_data(map)

###%#_______________________________________________________________________________________#%###
####           Create & check physeq object          ####  
###%#_______________________________________________________________________________________#%###
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
sub1 <- subset_samples(physeq_object, timepoint != "T3")
sub2 <- subset_samples(physeq_object, replicate!="NA") 
physeq_object <- merge_phyloseq(sub1, sub2)
summarize_phyloseq(physeq_object)
nsamples(physeq_object)
sample_names(physeq_object)

physeq_object_Min90T3 = physeq_object

get_taxa_unique(physeq_object, "Kingdom") # unassigned in Kingdom
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c(" ", "Unassigned", "NA")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Kingdom") # all good now

get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
length(get_taxa_unique(physeq_object,"Phylum"))  #36
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA", " ")) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum")) #35

length(get_taxa_unique(physeq_object,"Order"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c(" Chloroplast", "Chloroplast", "chloroplast", " chloroplast"))
length(get_taxa_unique(physeq_object,"Order"))

length(get_taxa_unique(physeq_object,"Family"))
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria", " Mitochondria"))
length(get_taxa_unique(physeq_object,"Family"))
get_taxa_unique(physeq_object,"Family")


####Pruning to get rid of singletons#####
physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
summarize_phyloseq(physeq_object)
ntaxa(physeq_object) ##10115
nsamples(physeq_object)  #63
rank_names(physeq_object)
max(taxa_sums(physeq_object)) 
min(taxa_sums(physeq_object)) #2
max(sample_sums(physeq_object))
min(sample_sums(physeq_object)) #422

####loops to redefine weird taxonomy to a single common character string "unassigned"#####
taxo <- as.data.frame(physeq_object@tax_table)


for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA")))) {taxo[i,y] <- "unassigned" }
  }
} 

taxo <- tax_table(as.matrix(taxo)) #re-define as tax table object

physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map) # merge updated taxonomy
ntaxa(physeq_object)
nsamples(physeq_object) #59
min(taxa_sums(physeq_object)) #0
max(taxa_sums(physeq_object))
min(sample_sums(physeq_object)) #850
max(sample_sums(physeq_object))

physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE)#no singletons
summarize_phyloseq(physeq_object)

physeq_object_TaxCorrd = physeq_object




###%#_______________________________________________________________________________________#%###
####  Colors!!   ####
####_______________________________________________________________________________________#%###

# Different colorpallettes to choose from

colors_by_Maaike <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
                      "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
                      "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497",
                      "#8D6346", "#FF4B3E", "#149911", "#472d30", "#ce796b",
                      "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
                      "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
                      "#F2AF29", "#606c38", "#3d405b", "#94d2bd", "#772e25",
                      "#344e41", "#0047E0", "#6c584c", "#5f0f40", "#D7F171", 
                      "#c89f9c", "#339989", "#faf3dd", "#04724d", "#98B9AB",
                      "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                      "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", 
                      "#662C91")

colors_M1 <- c("#004e64", "#ecc8af", "#F2AF29", "#436436", "#00a5cf", 
               "#c18c5d", "#5f0f40", "#DC602E", "#495867", "#A29F15", 
               "#570000", "#FFF5B2", "#20221B", "#9fffcb", "#c08497", 
               "#8D6346", "#FF4B3E", "#149911", "#472d30")


colors_M2 <- c( "#ce796b", 
                "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
                "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
                "#F2AF29", "#606c38", "#3d405b", "#94d2bd", "#772e25",
                "#344e41", "#0047E0", "#6c584c", "#5f0f40", "#D7F171", "#c89f9c" )


colors_M3 <- c(   "#339989", "#faf3dd", "#04724d", "#98B9AB",
                  "#b09e99", "#AD343E", "#F2AF29", "#362C28", "#5171A5",
                  "#F7FE72", "#F4978E", "#7A9B76", "#8A7E72", "#143642", 
                  "#662C91")

pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

HCB.col <- c("#EE7733")
PDB.col <- c("#0077BB")




###%#_______________________________________________________________________________________#%###
####            Prevalence plots      ####
####_______________________________________________________________________________________#%###
Prevalence <- read_pptx()
Prevalence <- read_pptx("../Reports/Prevalence_OpenBio_Prokaryotes.pptx")

# Aggregate taxa on taxonomic level of choice for analysis on taxonomic level
Physeq_genus <- aggregate_taxa(physeq_object, "Genus")
Physeq_genus <- subset_taxa(Physeq_genus, !Order%in% c("unassigned"))
Physeq_order <-aggregate_taxa(physeq_object, "Order")
Physeq_phylum <-aggregate_taxa(physeq_object, "Phylum") 
taxa_names(Physeq_order)


#Prevalence on all samples on different taxonomic levels
# Prevalence is measurement, which describes in how many samples certain microbes were detected
# detection gives a minimum RA to filter with
# Prevalence is the fraction of samples 
#select prevalence  with detection argument
Prevalence_plot <- microbiome::plot_taxa_prevalence(physeq_object, "Genus", detection = 0.1/100) + theme(legend.position = "none") + 
  ylab("Prevalence") + labs(title = "Prevalence Openbio Eukaryotic Phylum, prevalence 1/100") #prevalence
Prevalence_plot

Prevalence_plot <- microbiome::plot_taxa_prevalence(Physeq_genus, "Genus", detection = 0.1/100) + theme(legend.position = "none") + 
  ylab("Prevalence") + labs(title = "Prevalence Openbio Eukaryotic Phylum, prevalence 1/100") #prevalence
Prevalence_plot

prevalence_genus <- prevalence(Physeq_genus, detection = 0.01, sort = T)
length(prevalence_genus)  #854

editable_graph <- dml(ggobj = Prevalence_plot)
Prevalence <- add_slide(Prevalence) 
Prevalence <- ph_with(x = Prevalence, editable_graph,location = ph_location_type(type = "body") )
print(Prevalence, target = "../Reports/Prevalence_Eukaryotes.pptx")

####_______________________________________________________________________________________#%###
####           Core taxa  on GENUS level      ####
####_______________________________________________________________________________________#%###
# change abundance into RA in physeq object
# RA is enough, since this function does not work with statistics
# it selects taxa based on minimum RA and prevalence
ps.rel <- microbiome::transform(physeq_object, "compositional")
ps.rel.genus <- microbiome::transform(Physeq_genus, "compositional")

# 
# det <- c(0, 0.1, 0.5, 2, 5, 20)/100
# prevalences <- seq(.05, 1, .05)
# 
# 
# core.line <- plot_core(ps.rel, prevalences = prevalences, 
#                        detections = det, plot.type = "lineplot") + 
#   xlab("Relative Abundance (%)") + 
#   theme_bw()
# 
# core.line

# determine prevalence of core genera
# this function gives the % of samples in which the Genus is detected
# with given detection limit (in this cae RA>0.05)
core.gen.prev <- prevalence(ps.rel.genus, detection = 0.01, sort = T) 
head(core.gen.prev) 
str(core.gen.prev)
#Convert named numeric vector into dataframe to write as a file. 
top.gen.prev.df <- data.frame(core.gen.prev) %>% top_n(20)
write.table(top.gen.prev.df, '../core_community/top.10_core_genera_prevalence_overall_0224.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")


# determine core taxa
# detection at least 1.0%, prevalence at least 10 of all the samples samples
# gives the names of the core genera under given conditions
core.genera.taxa <- core_members(ps.rel.genus, detection = 0.01, prevalence = 1/3)
length(core.genera.taxa) #9 core taxa

# gives the proportion of core genera for each sample
core.abund <- microbiome::core_abundance(ps.rel.genus, detection = 0.01, 
                                         prevalence = 2/3)

# filter complete phyloseq object to contain only prevalent taxa 
core.gens.pseq <- core(ps.rel.genus, detection = 0.01, prevalence = 2/3)

# #determine sample sums in filtered physeq
core.abundance <- core_abundance(ps.rel.genus, detection = 0.01, prevalence = 2/3)
core.gen.samp.abund <- data.frame(core.abundance)

write.table(core.gen.samp.abund, '../core_community/top.10_core_genera_prev0.67_samplepopulation_freq_0228.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")

#you can aggregate the physeq, to collapse rare genera in "Other"category
#core.gens2 <- aggregate_rare(core.gens, "Genus", detection = 0.01, prevalence = 3/59)

#Select taxonomy and relative abundance of core genera and store
gens.tax <- tax_table(core.gens.pseq) %>% as.data.frame()
gens.otu <- otu_table(core.gens.pseq) %>% as.data.frame()
top10.gens.abund.df <-cbind(gens.tax,gens.otu) %>% top_n(10)

write.table(top10.gens.abund.df, '../core_community/top.10_core_genera_RelAb_overall_0223.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")

#Set sequence of prevalences 
prevalences <- seq(0.01, 1, .05)
#set sequence of minimal detection level
detections <- round(10^seq(log10(1e-5), log10(.2), length = 20), 3)

cp.1 <- plot_core(ps.rel.genus,
                  plot.type = "heatmap",
                  colours = hcl.colors(5,palette = "temps"),
                  prevalences = prevalences, 
                  detections = detections, min.prevalence = 0.5) +
  xlab("Detection Threshold (Relative Abundance (%))")
cp.1 <- cp.1 + theme_bw() + ylab("Genera")
cp.1

#some long genera names, change as follows 
taxa_names(ps.rel.genus) <- gsub("Bacteria_Proteobacteria_Alphaproteobacteria_", "" , taxa_names(ps.rel.genus))
taxa_names(ps.rel.genus) <- gsub("Bacteria_Bacteroidota_Bacteroidia_"," ", taxa_names(ps.rel.genus))

cp.1 <- plot_core(ps.rel.genus,
                  plot.type = "heatmap",
                  colours = hcl.colors(10,palette = "temps"),
                  prevalences = prevalences, 
                  detections = detections, min.prevalence = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90),
        axis.text.y = element_text(face = "italic")) +
  labs(title = "Core community genera overall", subtitle = "min.prevalence 0.5") +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  ylab("Genera")
cp.1

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.01, 1, .05)

core.line <- plot_core(ps.rel.genus, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

core.line

###%# Subsetting to detrmine core for different timepoints and polymers  #%###

#T6 ####
sub_t6 <- subset_samples(ps.rel.genus , timepoint == "T6")
sub_t6_cbackbone <- subset_samples(sub_t6, category =="Carbon")
sub_t6_hatoms <- subset_samples(sub_t6, category =="Hetero")

#T1 ####
sub_t1 <- subset_samples(ps.rel.genus , timepoint == "T1")
sub_t1_cbackbone <- subset_samples(sub_t1, category =="Carbon")
sub_t1_hatoms <- subset_samples(sub_t1, category =="Hetero")

###____category____####
sub_cbackbone <- subset_samples(ps.rel.genus, category =="Carbon")
sub_hatoms <- subset_samples(ps.rel.genus, category =="Hetero")

####____polymer_____########
sub_pe <- subset_samples(ps.rel.genus , Material == "PE")
sub_pp <- subset_samples(ps.rel.genus , Material == "PP")
sub_ps <- subset_samples(ps.rel.genus , Material == "PS")
sub_pet <- subset_samples(ps.rel.genus , Material == "PET")
sub_nylon <- subset_samples(ps.rel.genus , Material == "Nylon")

####______________________________________#%###
####  Test core genera                     ####

# determine prevalence of core genera
# this function gives the % of samples in which the Genus is detected
# with given detection limit (in this cae RA>0.05)
core.gen.prev <- prevalence(sub_t6_hatoms, detection = 0.01, sort = T) 
head(core.gen.prev) 
str(core.gen.prev)
#Convert named numeric vector into dataframe to write as a file. 
top.gen.prev.df <- data.frame(core.gen.prev) %>% top_n(10)
write.table(top.gen.prev.df, '../core_community/top.10_core_genera__prevalence_hatoms_T6_0224.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")


# determine core taxa
# detection at least 1.0%, prevalence at least 10 of all the samples samples
# gives the names of the core genera under given conditions
core.genera.taxa <- core_members(sub_t6_hatoms, detection = 0.01, prevalence = 1/3)
length(core.genera.taxa) #11 core taxa

# filter complete phyloseq object to contain only prevalent taxa 
core.gens.pseq <- core(sub_t6_hatoms, detection = 0.01, prevalence = 1/3)

# #determine sample sums in filtered physeq
core.abundance <- core_abundance(sub_pe, detection = 0.01, prevalence = 1/3)
sample.data <- (sample_data(sub_pe))[,1]
core.gen.samp.abund <- cbind(sample.data,data.frame(core.abundance))


write.table(core.gen.samp.abund, '../core_community/core_genera_prev0.67_samplepopulation_PE_freq_0228.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")

#you can aggregate the physeq, to collapse rare genera in "Other"category
#core.gens2 <- aggregate_rare(core.gens, "Genus", detection = 0.01, prevalence = 3/59)

#Select taxonomy and relative abundance of core genera and store
gens.tax <- tax_table(core.gens.pseq) %>% as.data.frame()
gens.otu <- otu_table(core.gens.pseq) %>% as.data.frame()
top10.gens.abund.df <-cbind(gens.tax,gens.otu) 

write.table(top10.gens.abund.df, '../core_community/core_genera_RelAb_hatoms_T6_0224.txt',
            na = "NA", dec = '.', quote = F, sep = "\t")

#Set sequence of prevalences 
prevalences <- seq(0.01, 1, .05)
#set sequence of minimal detection level
detections <- round(10^seq(log10(1e-5), log10(.2), length = 20), 3)

cp.1 <- plot_core(sub_t6_hatoms,
                  plot.type = "heatmap",
                  colours = hcl.colors(5,palette = "temps"),
                  prevalences = prevalences, 
                  detections = detections, min.prevalence = 2/3) +
  xlab("Detection Threshold (Relative Abundance (%))")
cp.1 <- cp.1 + theme_bw() + ylab("Genera")
cp.1

#some long genera names, change as follows 
taxa_names(sub_t6_hatoms) <- gsub("Bacteria_Proteobacteria_Alphaproteobacteria_", "" , taxa_names(sub_t6_hatoms))
taxa_names(sub_t6_hatoms) <- gsub("Bacteria_Bacteroidota_Bacteroidia_"," ", taxa_names(sub_t6_hatoms))
taxa_names(sub_t6_hatoms) <- gsub("Bacteria_Proteobacteria_*proteobacteria_"," ", taxa_names(sub_t6_hatoms))
taxa_names(sub_t6_hatoms) <- gsub("Bacteria_Proteobacteria_Gammaproteobacteria_"," ", taxa_names(sub_t6_hatoms))

cp.1 <- plot_core(sub_t6_hatoms,
                  plot.type = "heatmap",
                  colours = hcl.colors(10,palette = "temps"),
                  prevalences = prevalences, 
                  detections = detections, min.prevalence = 2/3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90),
        axis.text.y = element_text(face = "italic")) +
  labs(title = "Core community genera heteroatome backbone T6", subtitle = "min.prevalence 0.667") +
  xlab("Detection Threshold (Relative Abundance (%))") + 
  ylab("Genera")
cp.1

####_______________________________________________________________________________________#%###
####         Shared cores by Venn                                         ####
####_______________________________________________________________________________________#%###
# convert to RA
ps.rel.genus <- microbiome::transform(Physeq_genus, "compositional")
# simple way to count number of samples in each group
table(meta(ps.rel.genus)$Material, useNA = "always")

#Make a list of timepoints
material <- unique(as.character(meta(ps.rel.genus)$Material))
print(material)

# loop through physeq, and one by one conbine identified core taxa
list_core <- c() # an empty object to store information

for (n in material){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.rel.genus, Material == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.01, # 0.001 in atleast 90% samples 
                         prevalence = 2/3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)
core <- table(stack(list_core))
write.table(core, "../core_community/Overlap_core_genera_polymers.txt", 
            na = "NA", dec = '.', quote = F, sep = "\t")


ggvenn(list_core,
       fill_color = Pal.plast,
       stroke_size = 0.5) + 
  labs( title = "Overlap Core genera Category")

plot(venn(list_core),
     fills = Pal.plast)


####_______________________________________________________________________________________#%###
####         Calculate average abundance of taxa levels                                   ####
####_______________________________________________________________________________________#%###
top.16.genus <- top_taxa(ps.rel.genus, 16)
top6.phyl <- top_taxa(Physeq_phylum, 6)
top10.ord <- top_taxa(Physeq_order, 10)
ps.rel.order <- microbiome::transform(Physeq_order, "compositional")

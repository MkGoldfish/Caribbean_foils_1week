# Script for curation of the amplicon data for NIOZ140
# Create phyloseq object, filter/clean data, turn into tiday format. 

##%#_______________________________________________________________________________________#%###
####                 Working Directory                                                      ####       
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")

#load needed libraries
library(devtools)
library(BiocManager)
library(phyloseq)
library(tidyverse)
library(microbiome)

#####data_import######
tax <- as.matrix(read.delim('../v460_data/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- as.matrix(read.delim('../v460_data/asv_table.txt', row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('../v460_data/mapping_file_details_corrections290322_no_t3.txt', row.names = 1, na.strings = c("")))

physeq_object = merge_phyloseq(otu, tax, map)
summarize_phyloseq(physeq_object) # Total number of reads = 1,892,002

sub1 <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))
sub2 <- subset_samples(physeq_object, surface!="negative_c") 
physeq_object <- merge_phyloseq(sub1, sub2) 

####basic_info#####
summarize_phyloseq(physeq_object)
ntaxa(physeq_object)
nsamples(physeq_object)  
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #
physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90") 
max(sample_sums(physeq_object))


get_taxa_unique(physeq_object, "Kingdom") # unassigned in Domains
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c("", "Unassigned")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Kingdom") # all good now
get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA"," NA" )) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c("Chloroplast")) 
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria"))

physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
physeq_object <- filter_taxa(physeq_object, function(x) sum(x) > 1, TRUE)#no singletons
#loops to redefine weird taxonomy to a single common character string "unassigned" 
taxo <- as.data.frame(physeq_object@tax_table)
  taxo <-taxo %>%    mutate_at(c(1:ncol(taxo)), ~replace_na(.,"unassigned")) %>%
  mutate_at(c(1:ncol(taxo)), ~str_replace_all(., c("uncultured","Uncultured",
                                                   "metagenome","Metagenome","unknown","Unknown","NA"), "unassigned"))

taxo <- tax_table(as.matrix(taxo))
physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map)

#getting the phyloseq object as tidy tibble 
tidy_psmelt <- function(physeq) {
  ### INSERT Initial variable and rank name checking and modding from `psmelt`
  # Get the OTU table with taxa as rows
  rankNames = rank_names(physeq, FALSE)
  sampleVars = sample_variables(physeq, FALSE) 
  otutab <- otu_table(physeq)
  if (!taxa_are_rows(otutab)) {
    otutab <- t(otutab)
  }
  # Convert the otu table to a tibble in tidy form
  tb <- otutab %>% 
    as("matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    tidyr::gather("Sample", "Abundance", -OTU)
  # Add the sample data if it exists
  if (!is.null(sampleVars)) {
    sam <- sample_data(physeq) %>%
      as("data.frame") %>% 
      tibble::as_tibble(rownames = "Sample")
    tb <- tb %>%
      dplyr::left_join(sam, by = "Sample")
  }
  # Add the tax table if it exists
  if (!is.null(rankNames)) {
    tax <- tax_table(physeq) %>%
      as("matrix") %>%
      tibble::as_tibble(rownames = "OTU")
    tb <- tb %>%
      dplyr::left_join(tax, by = "OTU")
  }
  tb %>%
    arrange(desc(Abundance))
  # Optional conversion to a data frame doesn't affect the speed/memory usage
  # %>% as.data.frame
}



getwd()
source("tidy_tibble_maker.R")
tidy_physeq<- tidy_psmelt(physeq_object)
tidy_physeq<- tidy_physeq %>% filter(!Phylum == "unassigned") 
tidy_physeq$Species <- if_else((!tidy_physeq$Genus=="unassigned"&!tidy_physeq$Species=="unassigned"), 
                                   str_c(tidy_physeq$Genus," ",tidy_physeq$Species), str_c(tidy_physeq$Genus," sp."))
tidy_physeq <- tidy_psmelt(physeq_object)
tidy_physeq$LinkerPrimerSequence <- NULL
tidy_physeq$ReversePrimerSequence <- NULL
tidy_physeq$InputFileName <- NULL
tidy_physeq$BarcodeSequence <- NULL
tidy_physeq$BarcodeSequence_1 <- NULL
tidy_physeq$g_species <- str_c(tidy_physeq$Genus, " ", tidy_physeq$Species)


##############relative_abund_calc###########
t3 <- tidy_physeq  %>% group_by(Sample) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup() %>%
  group_by(Description) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
  ungroup() %>% 
  #Kingdom_section
  group_by(Sample, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  group_by(Description, Kingdom) %>% 
  mutate(Kingdom_st_dev_abund_samples = sd(Kingdom_rel_abund_Sample)) %>% # standard dev of Kingdom relative abundances between replicates of Description (ployner_timepoint.days._treatment)
  mutate(Kingdom_rep_rel_abund = sum(rep_rel_abund)) %>% #Kingdom relative abundance per samples of desc 
  ungroup() %>%
  #Phylum_section
  group_by(Sample, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Phylum) %>% 
  mutate(st_dev_Phylum_abund = sd(Phylum_rel_abund_Sample)) %>%
  mutate(Phyla_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Class_section
  group_by(Sample, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Class) %>% 
  mutate(st_dev_Class_abund = sd(Class_rel_abund_Sample)) %>%
  mutate(Class_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Order_section
  group_by(Sample, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Order) %>% 
  mutate(st_dev_Order_abund = sd(Order_rel_abund_Sample)) %>%
  mutate(Order_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Family_section
  group_by(Sample, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Family) %>% 
  mutate(st_dev_Family_abund = sd(Family_rel_abund_Sample)) %>%
  mutate(Family_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Genus_section
  group_by(Sample, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Genus) %>% 
  mutate(st_dev_Genus_abund = sd(Genus_rel_abund_Sample)) %>%
  mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Species_section
  group_by(Sample, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(Description, Species) %>% 
  mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup()  

write_csv(t3,"../Analysis/tidy_dataset_NIOZ140_nov22.csv")

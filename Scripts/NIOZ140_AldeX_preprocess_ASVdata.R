####_______________________________________________________________________________________#%###
####                 Working Directory                 ####
####_______________________________________________________________________________________#%###
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/essai_microbimeR")

library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(rmarkdown)
library("microbiome")
library("TreeSummarizedExperiment")
library(ggpubr)
#library("wesanderson")
#library("plotrix")
library("FactoMineR")
library("factoextra")
library(usedist)
#library("heatmaply")
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(ALDEx2)
library(knitr)



#####        Data_import                 ######
tax <- as.matrix(read.delim('../../Data_corrected_barcodes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt', row.names = 1, na.strings = c(" ")))
tax <- tax_table(tax)
otu <- as.matrix(read.delim('../../Data_corrected_barcodes/asv/asv_table.txt', row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('../../Data_corrected_barcodes/metadata/sampleList_mergedBarcodes_NIOZ140_corrected.txt', row.names = 1, na.strings = c("")))

set.seed(42)

####______Adding a column to mapfile divide plastic type into their structure____####
map["category"] <- NA #initialize empty column

for (i in 1:nrow(map)) {
  if (map$Material[i] %in% c("PE","PP","PS")) {
    map$category[i] <- "C_backbone"
  }
  else
    if (map$Material[i] %in% c("Nylon","PET")) {
      map$category[i] <- "heteratomes"  
    }
}

### add columns to mapfile to make combinations of factors for analysis
map["cat_uv"] <- str_c(map$category,"_",map$treatment)
map["time_UV"] <- str_c(map$timepoint, "_", map$treatment)
map["treat_time"] <- str_c(map$treatment, "_", map$timepoint)
map["time_mat"] <- str_c(map$timepoint, "_", map$Material)

#Create physeq object
map <- sample_data(map)
physeq_object = merge_phyloseq(otu, tax, map)    

####check basic_info phyloseq object#####
summarize_phyloseq(physeq_object) # Total number of reads = 1892002
ntaxa(physeq_object)  ##13169
nsamples(physeq_object) ##64
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #it says 0 sequences here

physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90") #eliminate the sample with 0seq
max(sample_sums(physeq_object))
min(sample_sums(physeq_object)) #534

#####subset T3 & merge ####         ### No longer in map file
physeq_object <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))

sub1 <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))
sub2 <- subset_samples(physeq_object, surface!="negative_c")
physeq_object <- merge_phyloseq(sub1, sub2)
sample_data(physeq_object)

summarize_phyloseq(physeq_object) # Total number of reads = 1,215,189
ntaxa(physeq_object)
####_______________________________________________________________________________________#%###
####      Check & correct taxonomy assignments         ####
####_______________________________________________________________________________________#%###
get_taxa_unique(physeq_object, "Kingdom") # unassigned in Kingdom
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c(" ", "Unassigned", "NA")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Kingdom") # all good now
summarize_phyloseq(physeq_object) # Total number of reads = 1,215,087

get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
length(get_taxa_unique(physeq_object,"Phylum"))  #36
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA", " ")) 
found_phyla <- get_taxa_unique(physeq_object, "Phylum") 
length(get_taxa_unique(physeq_object,"Phylum")) #35
summarize_phyloseq(physeq_object)# Total number of reads = 1,213,980

length(get_taxa_unique(physeq_object,"Order"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c(" Chloroplast", "Chloroplast", "chloroplast", " chloroplast"))
length(get_taxa_unique(physeq_object,"Order")) #184
found_orders <- get_taxa_unique(physeq_object, "Order")
summarize_phyloseq(physeq_object) # Total number of reads = 970,884

length(get_taxa_unique(physeq_object,"Family"))
physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria", " Mitochondria"))
length(get_taxa_unique(physeq_object,"Family")) #274
get_taxa_unique(physeq_object,"Family")
found_families <- get_taxa_unique(physeq_object, "Family") 
summarize_phyloseq(physeq_object) #Total number of reads =  928,418

length(get_taxa_unique(physeq_object,"Genus"))
physeq_object <- subset_taxa(physeq_object, !Genus%in% c("NA", ""))

length(get_taxa_unique(physeq_object,"Genus"))

found_families <- get_taxa_unique(physeq_object, "Order") 
summarize_phyloseq(physeq_object) #1,834,003 reads
length(get_taxa_unique(physeq_object, "Genus")) #628
found_genera <- get_taxa_unique(physeq_object, "Genus") 

####Pruning to get rid of singletons#####
physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons
physeq_object <-  subset_samples(physeq_object, !timepoint %in% c("T3"))
summarize_phyloseq(physeq_object) # Total number of reads = 891,917
ntaxa(physeq_object) ##9370
nsamples(physeq_object)  #63
rank_names(physeq_object)
max(taxa_sums(physeq_object)) 
min(taxa_sums(physeq_object)) #2
max(sample_sums(physeq_object))
min(sample_sums(physeq_object)) #422


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

tidy_physeq <- tidy_psmelt(physeq_object)
tidy_physeq$LinkerPrimerSequence <- NULL
tidy_physeq$ReversePrimerSequence <- NULL
tidy_physeq$InputFileName <- NULL
tidy_physeq$BarcodeSequence <- NULL
tidy_physeq$BarcodeSequence_1 <- NULL

#### Make genus tables
tidy_genus <- tidy_physeq  %>% group_by(Sample, Genus) %>% mutate(genus_sample_abundance = sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup()
tidy_genus$OTU <- NULL ## Getting rid of asv mentions to avoid repeitions
tidy_genus <- distinct(tidy_genus)
length(unique(tidy_genus$Genus))
tidy_genus %>% dplyr::select(Sample, Genus, genus_sample_abundance)# %>%  rownames_to_column(var = "Sample") 

genus_table <- tidy_genus %>% rownames_to_column() %>% dplyr::select(Sample, Genus, genus_sample_abundance)  

mat_genus = data.frame(matrix(0, nrow = length(unique(tidy_genus$Genus)),
                              ncol = length(unique(tidy_genus$Sample))),
                       row.names = unique(tidy_genus$Genus), stringsAsFactors = T)
# group_by(description) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
# ungroup() %>% 

#mat_genus = data.frame(matrix(0, nrow = length(unique(tidy_genus$Genus)), ncol = length(unique(tidy_genus$Sample))), row.names = unique(tidy_genus$Genus), stringsAsFactors = T)
# group_by(description) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
# ungroup() %>% 

colnames(mat_genus) <-  unique(tidy_genus$Sample)
for (s in 1:nrow(mat_genus)) { 
  for (h in 1:ncol(mat_genus)) {
    v <- filter(genus_table, Genus == rownames(mat_genus)[s] & Sample == colnames(mat_genus)[h]) %>% distinct()
    mat_genus[s,h] <- v$genus_sample_abundance 
    
  }
}

write_csv2(mat_genus, "./Analysis/genus_matrix.csv")

p2 <- physeq_object%>% subset_taxa(Genus %in% c(rownames(mat_genus)))
m <- as_tibble(p2@tax_table@.Data)
m$Species <- NULL
m <- distinct(m)
m <- data.frame(m)
rownames(m) <- m$Genus
genera <- m$Genus
write_lines(genera, "./Analysis/genera.txt")
write_csv2(m,"./Analysis/tax_table_genus_agglom.csv")
tax_genus <- tax_table(as.matrix(m))



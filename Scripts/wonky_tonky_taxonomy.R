# This function check taxonomy assignments in physeq-object, and removes, things we don't want: NA's, Chloroplast and mitochondria

wonky_tonky_taxonomy <- function(physeq_object){ 
  
  get_taxa_unique(physeq_object, "Kingdom") # unassigned in Kingdom
  physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & !Kingdom%in% c(" ", "Unassigned", "NA")) #let's eliminate those otus
  get_taxa_unique(physeq_object, "Kingdom") # all good now
  
  # get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
  phyla_before <- length(get_taxa_unique(physeq_object,"Phylum"))  
  physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA", " ")) 
  get_taxa_unique(physeq_object, "Phylum")
  phyla_after <- length(get_taxa_unique(physeq_object,"Phylum")) 
  present_phyla <- get_taxa_unique(physeq_object,"Phylum") %>% sort(decreasing = FALSE) %>% 
                  write.csv( "../Analysis/Found_phyla.csv" )
  
  order_before <- length(get_taxa_unique(physeq_object,"Order"))
  physeq_object <- subset_taxa(physeq_object, !Order%in% c(" Chloroplast", "Chloroplast", "chloroplast", " chloroplast"))
  order_after <- length(get_taxa_unique(physeq_object,"Order"))
  present_orders <- get_taxa_unique(physeq_object,"Order") %>% sort(decreasing = FALSE) %>%
                   write.csv( "../Analysis/Found_orders.csv" )
  
  
  family_before <- length(get_taxa_unique(physeq_object,"Family"))
  physeq_object <- subset_taxa(physeq_object, !Family%in% c("Mitochondria", " Mitochondria"))
  family_after <- length(get_taxa_unique(physeq_object,"Family"))
  present_fams <- get_taxa_unique(physeq_object,"Family") %>% sort(decreasing = FALSE) %>% 
                  write.csv( "../Analysis/Found_families.csv" )
  
  
  present_gena <- get_taxa_unique(physeq_object,"Genus") %>% sort(decreasing = FALSE) %>% 
                  write.csv( "../Analysis/Found_genera.csv" )
  genera_amount <- length(get_taxa_unique(physeq_object,"Genus"))
  
  
  present_species <- get_taxa_unique(physeq_object,"Species") %>% sort(decreasing = FALSE) %>% 
                    write.csv( "../Analysis/Found_species.csv" )
  species_amount<- length(get_taxa_unique(physeq_object,"Species"))
  
 
  
  
  taxa_sums_before_after <- list (sprintf("Phyla amount before is %s", phyla_before),
                                  sprintf("Phyla amount after is %s", phyla_after),
                                  sprintf("Order amount before is %s", order_before),
                                  sprintf("Order amount after is %s", order_after),
                                  sprintf("Family amount before is %s", family_before),
                                  sprintf("Family amount after is %s", family_after),
                                  sprintf("Genera amount is %s", genera_amount),
                                  sprintf("Species amount is %s", species_amount)
                            
  )

  return(taxa_sums_before_after)
  
  
  }

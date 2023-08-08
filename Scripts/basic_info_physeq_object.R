#Function to summarize and check physeq-objects

basic_info_physeq_object <- function(physeq_object) { 
  
  total_taxa <- ntaxa(physeq_object) 
  total_samples <- nsamples(physeq_object) 
  # sample_names <- sample_names(physeq_object)  
  # taxa <-taxa_names(physeq_object) 
  ranks <- rank_names(physeq_object)
  sample_sums <- sample_sums(physeq_object)
  taxa_sums <- taxa_sums(physeq_object)
  min_readnumber <- min(sample_sums(physeq_object))
  max_readnumber <- max(sample_sums(physeq_object))
  min_taxa <- min(taxa_sums(physeq_object)) 
  max_taxa <- max(taxa_sums(physeq_object))
  
  
  interesting_info <- list (sprintf("Total taxa is %s", total_taxa),
  sprintf("Total samples is %s", total_samples),
  sprintf("Rank name is %s", ranks),
  sprintf("Samplesum %s", sample_sums),
  sprintf("Taxasum %s", taxa_sums),
  sprintf("Lowest readnumber is %s", min_readnumber),
  sprintf("Highest readnumber is %s", max_readnumber),
  sprintf("Lowest taxa sum is %s", min_taxa),
  sprintf("Highest taxa sum is %s", max_taxa)
  ) 
  
  return(interesting_info)
 
  
}



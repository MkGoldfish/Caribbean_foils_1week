# This function is to transform the phyloseq object into tidy tibble 


tidy_tibble_maker <- function(physeq) {
  
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

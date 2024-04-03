make_table_2 <- function() {
  
  xlsx <- "data/MetaMLST_CF.xlsx"
  
  hi <- readxl::read_xlsx(xlsx, sheet = 1, range = "A1:C6") %>%
    mutate(species = "Haemophilus influenzae")
  
  pa <- readxl::read_xlsx(xlsx, sheet = 2, range = "A1:C2") %>%
    mutate(species = "Pseudomonas aeruginosa")
  
  sa <- readxl::read_xlsx(xlsx, sheet = 3, range = "A1:C15") %>%
    mutate(species = "Staphylococcus aureus")
  
  sm <- readxl::read_xlsx(xlsx, sheet = 4, range = "A1:C8") %>%
    mutate(species = "Stenotrophomonas maltophilia")
  
  sp <- readxl::read_xlsx(xlsx, sheet = 5, range = "A1:C2") %>%
    mutate(species = "Streptococcus pyogenes")
  
  metadata <- get_metadata()
  
  cols <- c("Species", "Patient ID", "Sample ID", "Sequence Type (ST)")
  
  rbind(hi, pa, sa, sm, sp) %>%
    dplyr::rename(sample_id = 3) %>%
    mutate(sample_id = gsub("_.*", "", sample_id)) %>%
    inner_join(metadata, by = "sample_id") %>%
    select(species, subject_id, alias, ST) %>%
    setNames(cols) %>%
    arrange(`Patient ID`)
  
}

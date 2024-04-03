run_maaslin <- function() {
  
  dir.create("maaslin", showWarnings = FALSE)
  
  metadata <- get_metadata() %>%
    arrange(sample_id) %>%
    column_to_rownames("sample_id")
  
  mpa <- read_tsv("data/metaphlan.tsv") %>%
    arrange(sample_id) %>%
    column_to_rownames("sample_id")
  
  Maaslin2::Maaslin2(input_data = mpa,
                     input_metadata = metadata,
                     output = "maaslin",
                     fixed_effects = c("fev1", "antibiotics", "sex"),
                     reference = "fev1,Mi/N",
                     random_effects = "subject_id",
                     correction = "holm",
                     plot_heatmap = FALSE,
                     plot_scatter = FALSE)
  
}

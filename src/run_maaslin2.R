library(tidyverse)
library(Maaslin2)

setwd("~/Desktop/CF_SHOTGUN/2022/")

metadata <- read_tsv("data/maaslin2_metadata.tsv", show_col_types = FALSE) %>%
  arrange(sample_id) %>%
  column_to_rownames("sample_id")

mpa <- read_tsv("data/metaphlan.2022.tsv", show_col_types = FALSE) %>%
  arrange(sample_id) %>%
  column_to_rownames("sample_id")

run_maaslin <- function(random_effect = c(NULL, "subject_id")) {
  
  if (is.null(random_effect)) {
    out_dir <- "maaslin/maaslin2_random_effect_false"
  } else {
    out_dir <- "maaslin/maaslin2_random_effect_true"
  }
  
  maaslin2_models <- Maaslin2::Maaslin2(
    input_data = mpa,
    input_metadata = metadata,
    output = out_dir,
    fixed_effects = c("fev1", "antibiotics", "sex"),
    reference = "fev1,Mi/N",
    random_effects = random_effect,
    plot_heatmap = FALSE,
    plot_scatter = FALSE)
  
}

run_maaslin(random_effect = NULL)
run_maaslin(random_effect = "subject_id")

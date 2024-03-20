library(readxl)
library(tidyverse)

setwd("~/Desktop/CF_SHOTGUN/")

out_dir <- "2022/"
suffix <- ".2022.tsv"

write_table <- function(data, name, out_dir = "2022/", suffix = ".2022.tsv") {
  write.table(data, paste0(out_dir, name, suffix), sep = "\t", quote = FALSE, row.names = FALSE)
}

######################
# MetaPhlAn + HUMAnN #
######################

mpa <- read_tsv("MetaPhlAn2/MetaPhlAn2_CF.txt", show_col_types = FALSE) %>%
  filter(ID != "#SampleID" & grepl("s__", ID) & !(grepl("t__", ID))) %>%
  mutate(ID = gsub(".*\\|", "", ID)) %>%
  setNames(gsub("_S.*", "", names(.))) %>%
  pivot_longer(!ID, names_to = "sample_id", values_to = "abundance") %>%
  pivot_wider(names_from = "ID", values_from = "abundance")

write_table(data = mpa, name = "metaphlan")

humann <- read_excel("CF_humann2_fast/CF_HUMAnN2_fast.xlsx") %>%
  setNames(gsub("_S*", "", names(.)))

write_table(humann, "humann")

#####################
# tidy the metadata #
#####################

sample_sheet <- read_csv("SampleSheet_JR.csv", skip = 18, col_select = c(1, 2), show_col_types = FALSE)

raw <- openxlsx::read.xlsx("Metadata_CF.xlsx", fillMergedCells = TRUE, colNames = FALSE, sheet = 1)

header <- as.data.frame(t(raw[c(1, 2), ])) %>%
  rename(V1 = 1, V2 = 2) %>%
  mutate(name = case_when(
    V2 == V1 ~ V1,
    V2 != V1 ~ paste0(V1, ": ", V2)
  )) %>%
  .$name

metadata <- raw[-c(1, 2),] %>%
  setNames(header) %>%
  as_tibble() %>%
  mutate_at(vars(contains(c("Resistance genes mechanisms", "Antibiotics in prior 2 months"))), funs(recode(., "Y" = 1, `0` = 0))) %>%
  mutate(Mucoidy = na_if(Mucoidy, 0))

numeric_columns <- as.numeric(suppressWarnings(
  metadata %>%
    mutate_if(is.character, as.numeric) %>%
    setNames(c(1:ncol(.))) %>%
    select(-1, 11) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    colnames(.)
))

tidied <- metadata %>%
  rename(sample_id = 1) %>%
  mutate(sample_id = gsub("\\.", "_", as.character(round(as.numeric(gsub("\\.$", "", sample_id)), 3)))) %>%
  mutate(across(.cols=numeric_columns, .fns=as.numeric)) %>%
  left_join(., sample_sheet, by = c("sample_id" = "Sample_ID")) %>%
  select(1, ncol(.), 2:(ncol(.) - 1))

write_table(tidied, "metadata")

############
# PanPhlAn #
############

# gene family Prokka annotations

centroid_headers <- list.files("PanPhlAn_out", full.names = T, pattern = ".headers")

panphlan_annotations <- c()

for (i in centroid_headers) {
  
  species <- gsub(".*panphlan_", "", gsub("_centroids.ffn.headers", "", i))
  
  tmp <- data.frame(species = species, header = readLines(i)) %>%
    as_tibble() %>%
    mutate(prokka_id = substr(gsub(".*\\:", "", header), 1, 12)) %>%
    mutate(prokka_annotation = substring(gsub(".*\\:", "", header), 14)) %>%
    mutate(gene = gsub("\\:.*", "", gsub(".*\\:g", "g", header))) %>%
    mutate(reference_genome = gsub("\\:.*", "", gsub(".*\\:G", "G", header))) %>%
    select(-header)
  
  panphlan_annotations <- rbind(panphlan_annotations, tmp)
  
}

write_table(panphlan_annotations, "panphlan_annotations")

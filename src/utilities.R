get_metadata <- function(full = FALSE) {
  
  metadata <- read_tsv("data/metadata.tsv")
  
  antibiotics <- metadata %>%
    select(2, contains("Antibiotics in prior 2 months")) %>%
    rename(sample_id = 1) %>%
    pivot_longer(!sample_id) %>%
    group_by(sample_id) %>%
    summarise(antibiotics = sum(value)) %>%
    ungroup()
  
  cols <- c("sample_id", "alias", "sex", "fev1", "patient")
  
  metadata %>%
    select(2, 1, 4, 6) %>%
    setNames(cols) %>%
    inner_join(antibiotics, by = "sample_id") %>%
    mutate(fev1 = recode(fev1, "Mi" = "Mi/N", "N" = "Mi/N")) %>%
    mutate(subject_id = str_pad(gsub("_.*", "", alias), 2, "left", "0")) %>%
    mutate(tmp = as.numeric(gsub(".*_", "", alias))) %>%
    arrange(tmp) %>%
    group_by(subject_id) %>%
    mutate(visit = row_number()) %>%
    ungroup()
  
}

get_top_n_spp <- function(n, mpa) {
  mpa %>%
    group_by(species) %>%
    summarise(abundance = mean(abundance)) %>%
    top_n(n, abundance) %>%
    pull(species)
}

pos_neg <- function(pathogen, mpa, alpha) {
  
  tmp <- unlist(str_split(gsub(".*__", "", pathogen), "_"))
  genus <- paste0(str_sub(tmp[1], 1, 1), ". ")
  bug <- paste0(genus, tmp[2])
  
  mpa %>%
    filter(species == pathogen) %>%
    rename(pathogen = 2) %>%
    mutate(status = case_when(abundance != 0 ~ "Positive",
                              abundance == 0 ~ "Negative")) %>%
    inner_join(alpha, by = "sample_id") %>%
    mutate(species = bug) %>%
    select(species, status, shannon)
  
}

get_dat <- function() {
  
  dat <- read_csv("data/resistome/resistome_culture_seq.csv") %>%
    mutate(Class = recode(Class,
                          "betalactams" = "Beta-lactams",
                          "Cationic antimicrobial peptides" = "Cationic AMPs")) %>%
    pivot_longer(cols = !c("Method", "Class"), names_to = "sample_id") %>%
    mutate(present = case_when(value > 0 ~ 1, value == 0 ~ 0)) %>%
    select(1, 2, 3, 5)
  
  abx <- dat %>%
    filter(Method == "Culture") %>%
    pull(Class) %>%
    unique()
  
  categories <- get_categories()
  
  dat %>%
    filter(Class %in% abx) %>%
    pivot_wider(names_from = Method, values_from = present, values_fill = 0) %>%
    mutate(value = case_when(
      Culture == 0 & Sequencing == 0 ~ categories[1],
      Culture == 1 & Sequencing == 1 ~ categories[2],
      Culture == 0 & Sequencing == 1 ~ categories[3],
      Culture == 1 & Sequencing == 0 ~ categories[4]
    )) %>%
    mutate(subject_id = str_pad(gsub("_.*", "", sample_id), 2, "left", "0"))
  
}

get_hit <- function(x) {
  
  x %>%
    mutate(hit = paste0(Class, "|", sample_id)) %>%
    select(hit) %>%
    .$hit
  
}

get_set <- function(method, x) {
  
  both_neg <- x %>%
    filter(value == "Both -") %>%
    get_hit()
  
  pos <- x %>%
    select(Class, sample_id, all_of(method)) %>%
    dplyr::rename("Method" = 3) %>%
    filter(Method > 0) %>%
    get_hit()
  
  output <- c(both_neg, pos)
  
  return(output)
  
}

get_categories <- function() {
  
  neg <- "\U2212"
  pos <- "+"
  
  c(paste("Both", neg),
    paste("Both", pos),
    paste("Sequencing", pos),
    paste("Culture", pos))
  
}

get_pal <- function(ref_gen = FALSE) {
  
  subject_ids <- get_metadata() %>%
    pull(subject_id) %>%
    unique() %>%
    sort()
  
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"),
           RColorBrewer::brewer.pal(3, "Set2"))[1:11]
  
  if (ref_gen) {
    
    pal <- c(pal, "ivory")
    
    names(pal) <- c(subject_ids, "Ref. genome")
    
  } else {
    
    names(pal) <- c(subject_ids)
    
  }
  
  return(pal)
  
}

run_pca <- function(inp_dir, sp = c("HI", "PA", "SA", "SM")) { 
  
  abbreviations <- c("HI", "PA", "SA", "SM")
  species_names <- c("***Haemophilus influenzae***",
                     "***Pseudomonas aeruginosa***",
                     "***Staphylococcus aureus***",
                     "***Stenotrophomonas maltophilia***")
  sp_df <- tibble(abbreviation = abbreviations, species_name = species_names)
  
  sp_name <- sp_df %>% filter(abbreviation == sp) %>% pull(species_name)
  
  filename <- sprintf("%s/gene_presence_absence_%s.csv", inp_dir, sp)
  
  dat <- read_delim(filename) %>%
    rename(gene = 1) %>%
    select(gene, contains("bmtagger")) %>%
    setNames(gsub("_.*", "", names(.))) %>%
    column_to_rownames("gene") %>%
    t()
  
  distances <- vegan::vegdist(dat, method = "bray")
  
  pca <- vegan::wcmdscale(distances, k = 2) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    rename(PC1 = 2, PC2 = 3) %>%
    mutate(species = sp_name)
  
  return(pca)
  
}

make_tree <- function(tree, metadata, pal) {
  
  species <- tree %>%
    basename() %>%
    str_replace("s__", "") %>%
    str_replace("_", " ") %>%
    str_replace("\\..*", "")
  
  title <- sprintf("***%s***", species)
  
  tre <- ggtree::read.tree(tree)
  
  ann <- tibble(tip_label = tre$tip.label) %>%
    mutate(sample_id = tip_label %>%
             str_replace("'", "") %>%
             str_replace("_.*", "")) %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(subject_id = ifelse(is.na(subject_id), "Ref. genome", subject_id))
  
  p <- ggtree::ggtree(tre, layout = "circular") %<+%
    ann
  
  p +
    geom_tippoint(aes(fill = subject_id), pch = 21, colour = "black", size = 3.5) +
    theme(plot.title = ggtext::element_markdown(),
          legend.title = element_text(face = "bold")) +
    ggtitle(title) +
    labs(fill = "Patient ID") +
    scale_fill_manual(values = pal)
  
}

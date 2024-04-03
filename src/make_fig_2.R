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

make_fig_2a <- function(metaphlan, metadata) {
  
  mpa <- metaphlan %>%
    pivot_longer(!sample_id, names_to = "species", values_to = "abundance")
  
  top_25_spp <- get_top_n_spp(25, mpa)
  
  heatmap_inp <- mpa %>%
    filter(species %in% top_25_spp) %>%
    pivot_wider(names_from = "sample_id", values_from = "abundance", values_fill = 0) %>%
    column_to_rownames("species") %>%
    as.matrix()
  
  ann_col <- metadata %>%
    select(sample_id, subject_id) %>%
    column_to_rownames("sample_id") %>%
    rename(Patient = 1)
  
  pal <- get_pal()
  
  annotation_colours <- list("Patient" = pal)
  
  ggplotify::as.ggplot(
    pheatmap::pheatmap(
      log(heatmap_inp + 1),
      color = viridis::viridis(option = "turbo", n = 100),
      border_color = "black",
      annotation_col = ann_col,
      annotation_colors = annotation_colours
    )
  )
  
}

make_fig_2b <- function(metaphlan, metadata) {
  
  alpha <- metaphlan %>%
    column_to_rownames("sample_id") %>%
    vegan::diversity(index = "shannon") %>%
    data.frame() %>%
    rownames_to_column("sample_id") %>%
    rename(shannon = 2)
  
  mpa <- metaphlan %>%
    pivot_longer(!sample_id, names_to = "species", values_to = "abundance") %>%
    inner_join(metadata, by = "sample_id")
  
  pathogens <- c(
    "s__Stenotrophomonas_maltophilia", 
    "s__Staphylococcus_aureus",
    "s__Prevotella_melaninogenica",
    "s__Pseudomonas_aeruginosa"
  )
  
  p_inp <- bind_rows(lapply(pathogens, pos_neg, mpa, alpha))
  
  max_y <- max(p_inp$shannon) * 1.1
  
  ann <- p_inp %>%
    group_by(species) %>%
    rstatix::wilcox_test(formula = shannon ~ status, data = .) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p, method = "fdr")) %>%
    filter(fdr <= 0.05) %>%
    select(species, fdr) %>%
    mutate(label = paste0("q=", round(fdr, 3))) %>%
    mutate(shannon = max_y, status = NA)
  
  p <- ggplot(p_inp, aes(x = species, y = shannon, fill = status)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.25), pch = 21, colour = "black") +
    geom_boxplot(position = position_dodge(), alpha = 0.25, outlier.shape = NA) +
    geom_text(data = ann, aes(x = species, y = shannon, label = label)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    theme(axis.text.x = element_text(face = "italic")) +
    labs(x = "Species", y = "Shannon index", fill = "Status") +
    scale_fill_manual(values = wesanderson::wes_palette("Royal1"))
  
  cowplot::plot_grid(plot.new(), p, plot.new(),
                     nrow = 1,
                     rel_widths = c(1/6, 2/3, 1/6))
  
}

make_fig_2 <- function() {
  
  metaphlan <- read_tsv("data/metaphlan.tsv")
  
  metadata <- get_metadata()
  
  fig_2a <- make_fig_2a(metaphlan, metadata)
  fig_2b <- make_fig_2b(metaphlan, metadata)
  
  fig_2 <- cowplot::plot_grid(fig_2a, fig_2b,
                              nrow = 2,
                              rel_heights = c(1, 0.5),
                              scale = 0.9)
  
  return(fig_2)
  
}
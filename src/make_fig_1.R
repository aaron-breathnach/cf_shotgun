make_fig_1a <- function(metaphlan, metadata, subject_ids) {
  
  alpha <- metaphlan %>%
    column_to_rownames("sample_id") %>%
    vegan::diversity(index = "shannon") %>%
    data.frame() %>%
    rownames_to_column("sample_id") %>%
    dplyr::rename(shannon = 2)
  
  p_inp <- inner_join(alpha, metadata, by = "sample_id") %>%
    filter(subject_id %in% subject_ids)
  
  ggplot(p_inp, aes(x = visit, y = shannon)) +
    geom_line(aes(group = 1)) +
    facet_grid(~subject_id, space = "free", scales = "free") +
    theme_bw(base_size = 12.5) +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    labs(x = "Visit", y = "Shannon index") +
    scale_x_continuous(breaks = 1:10)
  
}

make_fig_1b <- function(metaphlan, metadata, subject_ids) {
  
  mpa <- metaphlan %>%
    pivot_longer(!sample_id, names_to = "species", values_to = "abundance")
  
  top_10_spp <- get_top_n_spp(10, mpa)
  
  p_inp <- inner_join(mpa, metadata, by = "sample_id") %>%
    filter(subject_id %in% subject_ids & species %in% top_10_spp)
  
  ggplot(p_inp, aes(x = visit, y = abundance)) +
    geom_area(aes(fill = species), colour = "black")  +
    facet_grid(~subject_id, space = "free", scales = "free") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = 1:10) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_brewer(palette = "Spectral") +
    theme_bw(base_size = 12.5) +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    labs(x = "Visit", y = "Relative abundance (%)", fill = "Species")
  
}

make_fig_1 <- function() {
  
  metadata <- get_metadata()
  
  subject_ids <- metadata %>%
    group_by(subject_id) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 3) %>%
    pull(subject_id)
  
  metaphlan <- read_tsv("data/metaphlan.tsv")
  
  fig_1a <- make_fig_1a(metaphlan, metadata, subject_ids)
  fig_1b <- make_fig_1b(metaphlan, metadata, subject_ids)
  
  fig_1 <- cowplot::plot_grid(fig_1a, fig_1b, nrow = 2, align = "v", axis = "lr")
  
  return(fig_1)
  
}

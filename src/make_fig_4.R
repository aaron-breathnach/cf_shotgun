make_fig_4a <- function(metadata, pal) {

  ids <- c("HI", "PA", "SA", "SM")
  
  inp_dir <- "data/panphlan"
  
  p_inp <- purrr::map(ids, function(x) run_pca(inp_dir, x)) %>%
    bind_rows() %>%
    inner_join(metadata, by = "sample_id")
  
  p <- ggplot(p_inp, aes(PC1, PC2)) + 
    facet_wrap(~species, scales = "free") +
    geom_point(aes(fill = subject_id), pch = 21, colour = "black", size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text.x = ggtext::element_markdown()) +
    ggtitle("PCA plots based on PanPhlAn outputs") + 
    labs(x = "PC1", y = "PC2", fill = "Patient ID")
  
  cowplot::plot_grid(plot.new(), p, plot.new(), nrow = 1, rel_widths = c(1, 2, 1))
  
}

make_fig_4b <- function(metadata, pal) {
  
  trees <- list.files("data/strainphlan", full.names = TRUE)
  
  plot_list <- purrr::map(trees, function(tree) make_tree(tree, metadata, pal))
  
  p <- cowplot::plot_grid(plotlist = plot_list, scale = 0.95, nrow = 1)
  
  return(p)
  
}

make_fig_4 <- function() {
  
  metadat <- get_metadata()
  palette <- get_pal(ref_gen = TRUE)
  
  fig_4a <- make_fig_4a(metadat, palette)
  fig_4b <- make_fig_4b(metadat, palette)
  
  fig_4 <- cowplot::plot_grid(fig_4a, fig_4b,
                              nrow = 2,
                              scale = 0.95,
                              labels = c("A", "B"))
  
  return(fig_4)
  
}

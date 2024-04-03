make_fig_3a <- function() {
  
  pal <- get_pal()
  
  dat <- read_csv("data/resistome/resistome_class_cpm.csv") %>%
    pivot_longer(cols = !(c(Pseudomonas, Sample, Individual)),
                 names_to = "ARG",
                 values_to = "CPM") %>%
    filter(CPM > 0) %>%
    mutate(Individual = as.character(str_pad(Individual, 2, pad = "0"))) %>%
    mutate(ARG = recode(ARG,
                        "Multi-drug resistance" = "MDR",
                        "Cationic antimicrobial peptides" = "Cationic AMPs"))
  
  ggplot(dat, aes(x = ARG, y = log1p(CPM), fill = Individual)) +
    geom_jitter(colour = "grey", width = 0.2, pch = 21) +
    geom_boxplot(colour = "black", fill = NA, outlier.shape = NA) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 12.5, face = "bold"),
      axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12.5),
      legend.title = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 12.5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 25) 
    ) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = pal) +
    labs(x = "", y = "log(copies per million)", fill = "Patient ID") +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  
}

make_fig_3b <- function() {
  
  dat <- get_dat()
  
  culture <- get_set("Culture", dat)
  seq <- get_set("Sequencing", dat)
  approaches <- list("Culture" = culture, "Sequencing" = seq)
  
  ggvenn::ggvenn(approaches, fill_alpha = 0.5, text_size = 5) +
    scale_fill_manual(values = c("red", "blue"))
  
}

make_fig_3c <- function() {
  
  dat <- get_dat()
  
  categories <- get_categories()
  
  pal <- c("lightgrey", "purple", "skyblue", "salmon")
  
  names(pal) <- categories
  
  ggplot(dat, aes(x = sample_id, y = Class, fill = value)) +
    geom_tile(colour = "black") +
    facet_grid(~subject_id, scales = "free", space = "free") +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 12.5, face = "bold"),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12.5),
      strip.text = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 12.5),
      legend.position = "top",
      plot.margin = margin(l = 25)
    ) +
    scale_x_discrete(
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0)),
      position = "right"
    ) +
    scale_fill_manual(values = pal) +
    labs(x = "Sample ID", y = "Class", fill = "")
  
}

make_fig_3 <- function() {
  
  fig_3a <- make_fig_3a()
  fig_3b <- make_fig_3b()
  fig_3c <- make_fig_3c()
  
  fig_3 <- cowplot::plot_grid(fig_3a, fig_3b,
                              nrow = 1,
                              rel_widths = c(1, 0.75),
                              labels = c("A", "B"))
  
  fig_3 <- cowplot::plot_grid(fig_3, fig_3c,
                              nrow = 2,
                              rel_heights = c(1, 1),
                              labels = c("", "C"))
  
  return(fig_3)
  
}

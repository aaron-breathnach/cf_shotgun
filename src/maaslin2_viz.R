#!/usr/bin/env Rscript

library(tidyverse)
library(scales)
library(cowplot)

setwd("~/Desktop/CF_shotgun/2022/")

palette <- c("Mi/N"="#FF0000", "Mo"="#00A08A", "Se"="#F2AD00", "M"="skyblue", "F"="hotpink", "≤3"="orangered", ">3"="midnightblue")

theme_aw <- function() {
  
  theme_classic() +
    theme(
      plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12.5, face = "bold"),
      axis.text = element_text(size = 12.5),
      strip.text = element_text(size = 12.5, face = "bold")
    )
  
}

boxplot_and_density_plot <- function(feat, var, data, n_per_group, annotation, palette, symbol, out_dir) {
  
  nudge <- (nlevels(as.factor(data$variable)) - 1) * 0.5
  
  plot_input <- data %>%
    filter(feature == feat)
  
  ann_n <- n_per_group %>%
    filter(feature == feat)
  
  ann_p <- annotation %>%
    filter(feature == feat)
  
  # boxplot
  p_boxplot <- ggplot(plot_input, aes(x = variable, y = abundance, colour = variable)) +
    geom_jitter(aes(colour = variable, fill = variable), pch = 21, alpha = 0.5, width = 0.25) +
    geom_boxplot(colour = "black", fill = NA, outlier.shape = NA) +
    geom_text(
      data = ann_n,
      aes(x = variable, y = abundance, colour = variable, label = paste0("|n>0|=", n)),
      size = 5
    ) +
    geom_text(
      data = ann_p, aes(x = variable, y = abundance, label = qval),
      size = 5,
      colour = "black",
      nudge_x = nudge,
      show.legend = F
    ) +
    theme_aw() +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none") +
    labs(title = feat, x = var, y = paste0("Abundance (", symbol, ")"))
  
  ggsave(plot = p_boxplot, file = paste0(out_dir, "/", var, "/boxplots/boxplot_", var, "_", feat, ".png"), dpi = 300, width = 5, height = 5, bg = "white")
  
  # density plot + bar
  
  delta <- as.numeric((
    plot_input %>%
      filter(abundance > 0) %>%
      mutate(delta = (max(sqrt(abundance)) + min(sqrt(abundance))) / 2))[1, 5])
  
  ann_p_density <- ann_p %>%
    mutate(abundance = exp(delta))
  
  p_density_in <- plot_input %>% filter(abundance > 0)
  
  p_density <- ggplot(p_density_in, aes(x = sqrt(abundance))) +
    geom_density(aes(colour = variable, fill = variable), alpha = 0.5) +
    geom_text(data = ann_p_density, aes(x = log(abundance), y = Inf, label = qval), vjust = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_aw() +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "top") +
    labs(x = paste0("sqrt(Abundance [", symbol, "])"), y = "Density", colour = var, fill = var)
  
  p_bar <- ggplot(ann_n, aes(x = variable, y = percent_non_zero, label = paste0(round(percent_non_zero), "%"))) +
    geom_bar(aes(colour = variable, fill = variable), stat = "identity", width = 0.75, alpha = 0.5) +
    geom_text(aes(colour = variable), vjust = -0.25) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_aw() +
    theme(legend.position = "none") +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    labs(x = var, y = "non-0 (%)")
  
  p_title <- ggdraw(get_title(p_boxplot))
  
  p_bd <- plot_grid(p_bar, p_density, axis="tblr", align="h", rel_widths=c(1, 3))
  p_title <- plot_grid(p_title, rel_widths = c(1, 1, 1))
  p_combined <- plot_grid(p_title, p_bd, nrow = 2, rel_heights = c(0.05, 1))
  
  ggsave(plot = p_combined, file=paste0(out_dir, "/", var, "/density_plots/density_plot_", var, "_", feat, ".png"), dpi = 300, width = 10, height = 5, bg = "white")
  
  output <- list("boxplot" = p_boxplot, "density_plot" = p_combined)
  
  return(output)
  
}

make_plots <- function(variable, df, metadata, sig_res, palette, symbol, out_dir) {
  
  dir.create(paste0(out_dir, "/", variable, "/boxplots"), recursive = TRUE)
  dir.create(paste0(out_dir, "/", variable, "/density_plots"), recursive = TRUE)
  
  q_values <- sig_res %>%
    filter(metadata == variable) %>%
    select(feature, qval) %>%
    group_by(feature) %>%
    filter(qval == min(qval)) %>%
    ungroup() %>%
    mutate(qval = paste0("q=", scientific(qval)))
  
  data <- df %>%
    pivot_longer(!SampleID, names_to = "feature", values_to = "abundance") %>%
    filter(feature %in% q_values$feature) %>%
    mutate(feature = factor(feature)) %>%
    inner_join(., metadata, by = c("SampleID" = "sample_id")) %>%
    select(SampleID, all_of(variable), feature, abundance) %>%
    rename(variable=2)
  
  zeros <- data %>%
    mutate(zero = case_when(
      abundance == 0 ~ "zero",
      abundance > 0 ~ "non-zero"
    )) %>%
    group_by(feature, variable, zero) %>%
    tally() %>%
    ungroup() %>%
    group_by(feature, variable) %>%
    mutate(total = sum(n)) %>%
    filter(zero == "non-zero") %>%
    mutate(percent_non_zero = 100 * (n / total)) %>%
    ungroup() %>%
    select(feature, variable, percent_non_zero) %>%
    pivot_wider(names_from = variable, values_from = percent_non_zero, values_fill = 0) %>%
    pivot_longer(!feature, names_to = "variable", values_to = "percent_non_zero")
  
  n_per_group <- data %>%
    filter(abundance > 0) %>%
    group_by(feature, variable) %>%
    tally() %>%
    ungroup() %>%
    merge(., data, by = c("feature", "variable")) %>%
    group_by(feature) %>%
    mutate(temp = max(abundance)) %>%
    ungroup() %>%
    select(feature, variable, n, temp) %>%
    rename(abundance = temp) %>%
    mutate(abundance = 1.1 * abundance) %>%
    distinct() %>%
    pivot_wider(names_from = variable, values_from = n, values_fill = 0) %>%
    pivot_longer(cols = !(c(feature, abundance)), names_to = "variable", values_to = "n") %>%
    merge(., zeros, by=c("feature", "variable"))
  
  group <- levels(as.factor(data$variable))[1]
  
  my_pal <- palette[levels(as.factor(data$variable))]
  
  annotation <- data %>%
    group_by(feature) %>%
    filter(abundance == max(abundance)) %>%
    ungroup() %>%
    mutate(abundance = 1.25 * abundance) %>%
    merge(., q_values, by = "feature") %>%
    mutate(variable = group)
  
  features <- levels(as.factor(annotation$feature))
  
  plots <- lapply(features, boxplot_and_density_plot, variable, data, n_per_group, annotation, my_pal, symbol, out_dir)
  
}

extract_plots <- function(list, plot) {
  
  plots <- c()
  
  n_variables <- length(list)
  
  for (i in 1:n_variables) {
    
    n <- length(list[[i]])
    
    for (j in 1:n) {
      
      p <- list[[i]][[j]][plot]
      
      plots <- c(plots, p)
      
    }
  }
  
  return(plots)
  
}

make_pdf <- function(plots, filename, h, w) {
  pdf(filename, height = h, width = w, family = "sans")
  print(plots)
  dev.off()
}

make_maaslin_plots <- function(maaslin_in, maaslin_out, metadata, out_dir, symbol = c("%", "CPM"), palette, q = 0.25) {
  
  if (!(dir.exists(out_dir))) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  maaslin2_input <- read_tsv(maaslin_in, show_col_types = F) %>%
    rename(SampleID = 1)
  
  # metadata
  metadata <- read_tsv(metadata, show_col_types = F)
  
  # the output from MaAsLin2
  maaslin2_out <- read_tsv(paste0("./", maaslin_out, "/significant_results.tsv"), show_col_types = F) %>%
    filter(qval <= q) %>%
    mutate(coef_sign = case_when( coef > 0 ~ 1 , coef < 0 ~ -1)) %>%
    mutate(value = factor(value)) %>%
    mutate(value.var = -log(qval) * coef_sign) %>%
    mutate(feature = factor(feature)) %>%
    mutate(sig = case_when(
      qval <= 0.001 ~ "***", 
      qval > 0.001 & qval <= 0.01 ~ "**",
      qval > 0.01 & qval <= 0.05 ~ "*",
      qval > 0.05 & qval <= 0.25 ~ "."
    ))
  
  # significance level
  significance <- maaslin2_out %>%
    group_by(metadata, feature) %>%
    filter(qval == min(qval))

  # summary heatmap
  p_tile <- ggplot(data = maaslin2_out, aes(x = feature, y = metadata, label = sig)) +
    geom_tile(aes(fill = value.var), colour = "black") +
    coord_equal() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "darksalmon") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12.5, face = "bold"),
      legend.title = element_text(size = 12.5, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "top",
      strip.text.y = element_text(angle = 0, size = 12.5, face = "bold"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank()
    ) +
    labs(x = "feature", y = "effect", fill = "sign(coef) * -log(qval)") +
    geom_text(data = significance, aes(x = feature, y = metadata, label = sig), size = 5, fontface = "bold")

  width <- nlevels(as.factor(maaslin2_out$feature)) / 1.5
  height <- nlevels(as.factor(maaslin2_out$metadata)) * 2
  
  ggsave(plot = p_tile, file = paste0(out_dir, "/maaslin2_heatmap.png"), dpi = 300, height = height, width = width)
  
  # diagnostic plots
  variables <- levels(as.factor(maaslin2_out$metadata))
  
  plots <- lapply(variables, make_plots, maaslin2_input, metadata, maaslin2_out, palette, symbol, out_dir)
  
  boxplots <- extract_plots(plots, "boxplot")
  make_pdf(boxplots, paste0(out_dir, ".boxplots.pdf"), h = 5, w = 5)
  
  density_plots <- extract_plots(plots, "density_plot")
  make_pdf(density_plots, paste0(out_dir, ".density_plots.pdf"), h = 6, w = 10)

}

list.dirs()

make_maaslin_plots(
  maaslin_in = "data/metaphlan.2022.tsv",
  maaslin_out = "maaslin/maaslin2_random_effect_true/",
  metadata = "data/maaslin2_metadata.tsv",
  out_dir = "maaslin/maaslin2_f",
  palette = palette,
  symbol = "%",
  q = 0.25
)

smaaslin_in = "metaphlan.2022.tsv"
maaslin_out = "maaslin2_random_effect_false"
metadata = "maaslin2_metadata.tsv"
out_dir = "maaslin2_figures_re_false"
palette = palette
symbol = "%"

variable = variables[1]
df = maaslin2_input
metdata = metadata
sig_res = maaslin2_out
palette = palette
symbol = symbol

marrangeGrob(grobs = plotz, nrow = 1, ncol = 1)

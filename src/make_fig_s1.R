make_fig_s1 <- function() {
  
  p_inp <- read_delim("data/reads_per_sample.tsv") %>%
    select(1, 5:6) %>%
    pivot_longer(!sample_id, names_to = "stage", values_to = "reads") %>%
    mutate(stage = factor(stage, levels = c("Raw", "Processed")))
  
  ggplot(p_inp, aes(x = stage, y = reads)) +
    geom_violin(aes(fill = stage), trim = FALSE, scale = "width", width = 0.75) +
    geom_boxplot(width = 0.1,  color = "gray20", fill = "white") +
    theme_bw(base_size = 12.5) +
    theme(legend.title = element_text(face="bold"),
          legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = wesanderson::wes_palette("Royal1")) +
    labs(y = "No. Reads", fill = "Stage") +
    scale_y_sqrt(label = scales::comma, breaks = c(1000000, 5000000, 1000000, 20000000))
  
}
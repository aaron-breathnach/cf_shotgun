library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggvenn)

data <- read_csv("~/Desktop/CF_SHOTGUN/resistome_out/resistome_culture_seq.csv") %>%
  mutate(Class = recode(
  Class,
  "betalactams" = "Beta-lactams"  
))

levels(as.factor(data$Class))

df <- data %>%
  pivot_longer(cols = !c("Method", "Class"), names_to = "sample_id", values_to = "value") %>%
  mutate(present = case_when(
    value > 0 ~ 1,
    value == 0 ~ 0
  )) %>%
  select(-value)

abx <- df %>%
  filter(Method == "Culture") %>%
  select(Class) %>%
  distinct() %>%
  .$Class

hmap_in <- df %>%
  filter(Class %in% abx) %>%
  mutate(Class = recode(Class, "Cationic antimicrobial peptides" = "Cationic AMPs")) %>%
  pivot_wider(names_from = Method, values_from = present, values_fill = 0) %>%
  mutate(value = case_when(
    Culture == 0 & Sequencing == 0 ~ "Both -",
    Culture == 1 & Sequencing == 1 ~ "Both +",
    Culture == 0 & Sequencing == 1 ~ "Sequencing +",
    Culture == 1 & Sequencing == 0 ~ "Culture +"
  )) %>%
  mutate(subject_id = str_pad(gsub("_.*", "", sample_id), 2, "left", "0"))

agreement <- hmap_in %>%
  group_by(Class) %>%
  mutate(agree = case_when(
    grepl("Both", value) ~ "yes",
    !grepl("Both", value) ~ "no"
  )) %>%
  group_by(agree) %>%
  tally() %>%
  ungroup() %>%
  mutate(percent = 100 * n / sum(n))

title <- paste0(round(agreement[2, 3], 2), "% agreement")

# hc_in <- hmap_in %>%
#   mutate(value = case_when(
#     value == "Sequencing +" ~ 4,
#     value == "Both +" ~ 3,
#     value == "Culture +" ~ 2,
#     value == "Both -" ~ 1
#   )) %>%
#   select(sample_id, Class, value) %>%
#   pivot_wider(names_from = Class, values_from = value, values_fill = 0) %>%
#   column_to_rownames("sample_id")

# get_order <- function(x) {
#   dendro <- as.dendrogram(hclust(dist(x)))
#   order <- labels(dendro)
#   return(order)
# }
# 
# x_order <- get_order(hc_in)
# y_order <- get_order(t(hc_in))
# 
# hmap_in$sample_id <- factor(hmap_in$sample_id, levels = x_order)
# hmap_in$Class <- factor(hmap_in$Class, levels = y_order)

amr_hmap <- ggplot(hmap_in, aes(x = sample_id, y = Class, fill = value)) +
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
  scale_fill_manual(
    values = c(
      "Both -" = "lightgrey",
      "Both +" = "purple",
      "Sequencing +" = "skyblue",
      "Culture +" = "salmon"
    )
  ) +
  # ggtitle(title) +
  labs(x = "Sample ID", y = "Class", fill = "")
amr_hmap

##############################
## ResitomeAnalyzer boxplot ##
##############################

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

amr_boxplot <- read_csv("~/Desktop/CF_SHOTGUN/resistome_out/resistome_class_cpm.csv") %>%
  pivot_longer(cols = !(c(Pseudomonas, Sample, Individual)), names_to = "ARG", values_to = "CPM") %>%
  filter(CPM > 0) %>%
  mutate(Individual = as.character(str_pad(Individual, 2, pad = "0"))) %>%
  mutate(ARG = recode(ARG,
                        "Multi-drug resistance" = "MDR",
                        "Cationic antimicrobial peptides" = "Cationic AMPs")) %>%
  ggplot(., aes(x = ARG, y = CPM, fill = Individual)) +
  geom_jitter(colour = "grey", width = 0.2, pch = 21) +
  geom_boxplot(colour = "black", fill = NA, outlier.shape = NA) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12.5, face = "bold"),
    axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12.5),
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 12.5)
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "Copies per million (CPM)", fill = "Patient ID")

##################
## venn diagram ##
##################

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
    rename("Method" = 3) %>%
    filter(Method > 0) %>%
    get_hit()
  
  output <- c(both_neg, pos)
  
  return(output)
  
}

culture <- get_set("Culture", hmap_in)
seq <- get_set("Sequencing", hmap_in)
approaches <- list("Culture" = culture, "Sequencing" = seq)

venn_diagram <- ggvenn::ggvenn(approaches, fill_alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"))

#######################
## combine the plots ##
#######################

# a <- plot_grid(plot.new(), amr_boxplot, plot.new(), rel_widths = c(1/6, 2/3, 1/6), nrow = 1, labels = c("", "A", ""))
# 
# b_c <- plot_grid(amr_hmap, venn_diagram, nrow = 1, rel_widths = c(1, 0.5), labels = c("B", "C"))
# 
# a_b_c <- plot_grid(a, b_c, nrow = 2)

a_b <- plot_grid(amr_boxplot, venn_diagram, nrow = 1, rel_widths = c(1, 0.5), labels = c("A", "B"))

a_b_c <- plot_grid(a_b, amr_hmap, nrow = 2, rel_heights = c(1, 1), labels = c("", "C"))

ggsave(a_b_c, file = "~/Desktop/CF_SHOTGUN/2022/plots/amr_fig.png", height = 9, width = 13.5, dpi = 300, bg = "white")

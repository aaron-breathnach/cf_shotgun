library(tidyverse)
library(vegan)
library(RColorBrewer)

setwd("~/Desktop/CF_SHOTGUN/2022/")

theme_aw <- function() {
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
}

make_subsets <- function(combo, data) {
  
  data_subset <- data %>% 
    filter(grepl(combo, variable)) %>%
    mutate(facet = gsub("\\|", " v ", combo))
  
  return(data_subset)
  
}

make_mds_plots <- function(mds, var, variable_name, palette, fixed = TRUE) {
  
  data <- mds %>%
    select(all_of(var), MDS1, MDS2) %>%
    rename(variable = 1)
  
  var_to_incl <- data %>%
    group_by(variable) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 5) %>%
    .$variable
  
  data <- data %>% filter(variable %in% var_to_incl)
  
  if (length(var_to_incl) > 2) {
    
    combos <- apply(combn(var_to_incl, 2), 2, paste, collapse = '|')
    
    data <- bind_rows(lapply(combos, make_subsets, data))
    
    tmp <- ggplot(data, aes(x = MDS1, y = MDS2, colour = variable, fill = variable)) +
      facet_grid(~ facet)
    
  } else {
    
    tmp <- ggplot(data, aes(x = MDS1, y = MDS2, colour = variable, fill = variable))
    
  }
  
  p <- tmp +
    stat_ellipse(geom = "polygon", alpha = 0.25) +
    geom_point(pch = 21, colour = "black") +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(colour = variable_name, fill = variable_name)
  
  if (fixed) {
    p <- ggplotify::as.ggplot(egg::set_panel_size(p, width = unit(5, "cm"), height = unit(5, "cm")))
  }
  
  return(p)
  
}

pos_neg <- function(pathogen, mpa, alpha) {
  
  tmp <- unlist(str_split(gsub(".*__", "", pathogen), "_"))
  genus <- paste0(str_sub(tmp[1], 1, 1), ". ")
  bug <- paste0(genus, tmp[2])
  print(bug)
  
  info <- mpa %>%
    filter(species == pathogen) %>%
    rename(pathogen = 2) %>%
    mutate(status = case_when(
      abundance != 0 ~ "Positive",
      abundance == 0 ~ "Negative"
    )) %>%
    inner_join(., alpha, by = "sample_id") %>%
    mutate(species = bug) %>%
    select(species, status, shannon)
  
  return(info)
  
}

make_alpha_box <- function(var, p_vals, input, palette) {
  
  tmp <- suppressWarnings(any(unlist(lapply(unlist(strsplit(var, "")), function(x) { !is.na(as.numeric(x)) }))))
  
  if (tmp) {
    x_lab <- str_to_upper(var)
  } else {
    x_lab <- str_to_title(var)
  }
  
  df <- input %>%
    select(sample_id, subject_id, all_of(var), shannon) %>%
    rename(variable = 3)
  
  p_value <- p_vals %>% filter(covariate == var) %>% pull(p)
  
  pal <- palette[levels(as.factor(df$variable))]
  
  x <- df$variable %>%
    unique() %>%
    sort() %>%
    nth(1)
  
  p <- ggplot(df, aes(x = variable, y = shannon, fill = variable)) +
    geom_jitter(width = 0.25, pch = 21, colour = "black") +
    geom_boxplot(alpha = 0.25, outlier.shape = NA) +
    geom_text(data = tibble(variable = x,
                            shannon = Inf,
                            label = paste0("p=", round(p_value, 3))),
              aes(x = variable, y = shannon, label = label),
              vjust = 1) +
    theme_aw() +
    theme(legend.position = "none") +
    labs(x = x_lab, y = "Shannon index") +
    scale_fill_manual(values = palette)
  
  n <- nlevels(as.factor(df$variable))
  
  width <- (n / 3) * 9
  
  p_fixed <- ggplotify::as.ggplot(egg::set_panel_size(p, height = unit(6, "cm"), width = unit(width, "cm")))
  
  return(p_fixed)
  
}

############
# metadata #
############

metadata <- read_tsv("data/metadata.2022.tsv", show_col_types = FALSE)

antibiotics <- metadata %>%
  select(2, contains("Antibiotics in prior 2 months")) %>%
  rename(sample_id = 1) %>%
  pivot_longer(!sample_id) %>%
  group_by(sample_id) %>%
  summarise(antibiotics = sum(value)) %>%
  ungroup()

meta <- metadata %>%
  select(2, 1, 4, 6) %>%
  rename(alias = 2) %>%
  mutate(patient = str_pad(gsub("_.*", "", alias), 2, "left", "0")) %>%
  setNames(c("sample_id", "alias", "sex", "fev1", "subject_id", "patient")) %>%
  inner_join(., antibiotics, by = "sample_id") %>%
  mutate(fev1 = recode(fev1, "Mi" = "Mi/N", "N" = "Mi/N")) %>%
  mutate(time = as.numeric(gsub(".*_", "", alias))) %>%
  arrange(time) %>%
  group_by(subject_id) %>%
  mutate(time = row_number()) %>%
  arrange(subject_id)

write.table(meta, "data/maaslin2_metadata.tsv", sep = "\t", quote = F, row.names = F)

#############
# MetaPhlAn #
#############

mpa_raw <- read_tsv("data/metaphlan.2022.tsv", show_col_types = FALSE)

mpa_veg <- mpa_raw %>%
  column_to_rownames("sample_id")

##

mpa <- mpa_raw %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "abundance")

get_top_n_spp <- function(n, mpa) {
  mpa %>%
    group_by(species) %>%
    summarise(abundance = mean(abundance)) %>%
    top_n(n, abundance) %>%
    pull(species)
}

top_25_spp <- get_top_n_spp(25, mpa)

top_10_spp <- get_top_n_spp(10, mpa)

heatmap_inp <- mpa %>%
  filter(species %in% top_25_spp) %>%
  pivot_wider(names_from = "sample_id", values_from = "abundance", values_fill = 0) %>%
  column_to_rownames("species") %>%
  as.matrix()

ann_col <- meta %>%
  select(sample_id, subject_id) %>%
  rename(Patient = 2) %>%
  column_to_rownames("sample_id")

patients <- meta %>%
  pull(subject_id) %>%
  unique()

patient_pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))[1:length(patients)]
names(patient_pal) <- sort(patients)

annotation_colours <- list("Patient" = patient_pal)

pheatmap::pheatmap(log(heatmap_inp + 1),
                   color = viridis::viridis(option = "turbo", n = 100),
                   border_color = "black",
                   annotation_col = ann_col,
                   annotation_colors = annotation_colours,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "plots/top_25_spp.png")

###################
# alpha diversity #
###################

alpha <- data.frame(vegan::diversity(mpa_veg, index = "shannon")) %>%
  rownames_to_column("sample_id") %>%
  rename(shannon = 2)

mpa <- mpa_raw %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "abundance")

mpa_meta <- inner_join(mpa, meta, by = "sample_id")

# alpha diversity: pathogens

pathogens <- c(
  "s__Stenotrophomonas_maltophilia", 
  "s__Staphylococcus_aureus",
  "s__Prevotella_melaninogenica",
  "s__Pseudomonas_aeruginosa")

input <- bind_rows(lapply(pathogens, pos_neg, mpa, alpha))

max_y <- max(input$shannon) * 1.1

ann <- input %>%
  group_by(species) %>%
  rstatix::wilcox_test(formula = shannon ~ status, data = .) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "fdr")) %>%
  filter(fdr <= 0.05) %>%
  select(species, fdr) %>%
  mutate(label = paste0("q=", round(fdr, 3))) %>%
  mutate(shannon = max_y, status = NA)

alpha_box <- ggplot(input, aes(x = species, y = shannon, fill = status)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25), pch = 21, colour = "black") +
  geom_boxplot(position = position_dodge(), alpha = 0.25, outlier.shape = NA) +
  geom_text(data = ann, aes(x = species, y = shannon, label = label)) +
  theme_aw() +
  theme(axis.text.x = element_text(face = "italic")) +
  labs(x = "Species", y = "Shannon index", fill = "Status") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1"))

ggsave(alpha_box, file = "plots/alpha_div.png", height = 4, width = 6)

# alpha diversity: sex, fev1, antibiotics

alpha_meta <- alpha %>%
  inner_join(., meta, by = "sample_id") %>%
  mutate(antibiotics = ifelse(antibiotics <= 3, paste0("\U2264", 3), ">3"))

p_vals <- lmerTest::lmer(shannon ~ antibiotics + fev1 + sex + (1 | subject_id), alpha_meta) %>%
  anova() %>%
  as.data.frame() %>%
  rownames_to_column("covariate") %>%
  select(1, ncol(.)) %>%
  rename(p = 2) %>%
  mutate(p = round(p.adjust(p, "BH"), 3))

vars <- c("sex", "fev1", "antibiotics")

palette <- c("Mi/N"="#FF0000", "Mo"="#00A08A", "Se"="#F2AD00", "M"="skyblue", "F"="hotpink", "≤3"="orangered", ">3"="midnightblue")

plots <- lapply(vars, make_alpha_box, p_vals, alpha_meta, palette)

alpha_combo <- cowplot::plot_grid(plotlist = plots, rel_widths = c(1, (1 + (1 / 3)), 1), nrow = 1)

ggsave(alpha_combo, file = "plots/alpha_div_metadata.png", width = 10, height = 4, dpi = 300)

########
# pcoa #
########

vd <- vegan::vegdist(mpa_veg, method = "bray")

df_md <- meta %>%
  column_to_rownames("sample_id")

perm <- permute::how(nperm = 999)
permute::setBlocks(perm) <- with(df_md, subject_id)

set.seed(321)
adonis_out <- vegan::adonis2(vd ~ fev1 + antibiotics + sex,
                             data = df_md,
                             by = "margin",
                             permutations = perm) %>%
  as.data.frame() %>%
  drop_na() %>%
  rownames_to_column("variable") %>%
  select(variable, 4, 6) %>%
  rename(r2 = 2, p = 3) %>%
  # mutate(p = p.adjust(p, "BH")) %>%
  mutate(signif = case_when(
    p < 0.050 & p > 0.010 ~ "*",
    p < 0.010 & p > 0.001 ~ "**",
    p == 0.001 ~ "***"
  )) %>%
  mutate(variable = recode(variable, "fev1" = "FEV1", "antibiotics" = "Antibiotics", "sex" = "Sex"))
adonis_out

adonis_bar <- ggplot(adonis_out, aes(x = r2, y = reorder(variable, r2), label = paste0("p=", p))) +
  geom_bar(stat = "identity", colour = "black", fill = "steelblue", width = 0.5) +
  geom_text(hjust = -0.125) +
  scale_x_continuous(expand = expansion(mult = c(0, 1/3))) +
  theme_aw() +
  theme(axis.title.x = ggtext::element_markdown()) +
  labs(x = "_R_<sup>2</sup>", y = "Variable")

adonis_bar_fixed <- ggplotify::as.ggplot(egg::set_panel_size(adonis_bar, width = unit(5, "cm"), height = unit(5, "cm")))
ggsave(adonis_bar_fixed, file = "permanova.png", width = 7.5, height = 5)

#############
# MDS plots #
#############

patient_pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

set.seed(123)
mds <- vegan::metaMDS(vd)$points %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  inner_join(., meta, by = "sample_id")

## plots coloured by subject_id

mds_patient_id <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = subject_id)) +
  geom_point(colour = "black", pch = 21, size = 3) +
  scale_fill_manual(values = patient_pal) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "Patient ID")

mds_abx <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = antibiotics)) +
  geom_point(colour = "black", pch = 21, size = 3) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "# antibiotics")

mds_sex <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = sex)) +
  geom_point(colour = "black", pch = 21, size = 3) +
  scale_fill_manual(values = c("F" = "pink", "M" = "lightblue")) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "Sex")
mds_sex

mds_fev <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = fev1)) +
  geom_point(colour = "black", pch = 21, size = 3) +
  scale_fill_manual(values = palette) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(fill = "FEV1")
mds_fev

fig_s2 <- cowplot::plot_grid(bc_dists_plot, adonis_bar, mds_patient_id,
                   mds_fev, mds_abx, mds_sex, ncol = 3,
                   align = "v", axis = "lr")

ggsave("plots/fig_s2.png", fig_s2, width = 15, height = 7.5)

## calculate inter- and intra-individual Bray-Curtis distances between samples

sub_meta <- meta %>%
  select(sample_id, subject_id)

bc_dists <- vd %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("sample_1") %>%
  pivot_longer(!sample_1, names_to = "sample_2", values_to = "distance") %>%
  inner_join(., sub_meta, by = c("sample_1" = "sample_id")) %>%
  rename(subject_1 = subject_id) %>%
  inner_join(., sub_meta, by = c("sample_2" = "sample_id")) %>%
  rename(subject_2 = subject_id) %>%
  filter(sample_1 != sample_2) %>%
  mutate(pair = case_when(
    sample_1 < sample_2 ~ paste0(sample_1, "-", sample_2),
    sample_2 < sample_1 ~ paste0(sample_2, "-", sample_1)
  ))

bc_dists_inter <- bc_dists %>%
  filter(subject_1 != subject_2) %>%
  group_by(subject_1) %>%
  summarise(distance = mean(distance)) %>%
  ungroup() %>%
  mutate(measure = "Inter-individual")

bc_dists_intra <- bc_dists %>%
  filter(subject_1 == subject_2) %>%
  group_by(subject_1) %>%
  summarise(distance = mean(distance)) %>%
  ungroup() %>%
  mutate(measure = "Intra-individual")

long_subs <- meta %>%
  group_by(subject_id) %>%
  tally() %>%
  filter(n > 3) %>%
  pull(subject_id)

meta_time <- meta %>%
  ungroup() %>%
  select(sample_id, time)

meta

bc_dist_init <- bc_dists %>%
  filter(subject_1 == subject_2 & subject_1 %in% long_subs) %>%
  inner_join(meta_time, by = c("sample_1" = "sample_id")) %>%
  filter(time == min(time)) %>%
  as.data.frame() %>%
  select(sample_1, sample_2, distance) %>%
  inner_join(meta_time, by = c("sample_2" = "sample_id"))

ggplot(bc_dist_init, aes(x = time, y = distance)) +
  geom_line(aes(colour = sample_1, group = sample_1)) +
  geom_point()

bcd <- rbind(bc_dists_inter, bc_dists_intra)

tmp <- inner_join(mpa, meta, by = "sample_id") %>%
  filter(subject_id %in% long_subs & species %in% top_10_spp)

area_chart <- ggplot(tmp, aes(x = time, y = abundance)) +
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

ggsave("plots/area_chart.png", area_chart, height = 4, width = 10, dpi = 300)

alpha_stability <- alpha_meta %>%
  filter(subject_id %in% long_subs)

alpha_stab <- ggplot(alpha_stability, aes(x = time, y = shannon)) +
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

temp_stab <- cowplot::plot_grid(alpha_stab, area_chart, nrow = 2, align = "v", axis = "lr")

ggsave("plots/temporal_stability.png", temp_stab, height = 8, width = 12, dpi = 300)

bc_dists_plot <- rbind(bc_dists_inter, bc_dists_intra) %>%
  ggplot(., aes(x = measure, y = distance, fill = subject_1)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.25, pch = 21, size = 3) +
  scale_fill_manual(values = patient_pal) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Measure", y = "Mean Bray-Curtis distance", fill = "Patient ID")

b_c <- cowplot::plot_grid(
  bc_dists_plot,
  adonis_bar,
  rel_widths = c(0.75, 1),
  scale = 0.95,
  nrow = 1,
  labels = c("B", "C")
  )

a <- cowplot::plot_grid(
  plot.new(),
  mds_patient_id,
  plot.new(),
  rel_widths = c(1/6, 2/3, 1/6),
  nrow = 1,
  labels = c("", "A", "")
  )

subject_id_plot <- cowplot::plot_grid(
  a,
  b_c,
  scale = 0.95,
  nrow = 2
)

ggsave(subject_id_plot, file = "~/Desktop/CF_SHOTGUN/2022/plots/individuals.png",
       height = 7.5, width = 7.5, dpi = 300, bg = "white")

mds_fev <- make_mds_plots(mds, "fev1", "FEV1", wesanderson::wes_palette("Darjeeling1"))
mds_sex <- make_mds_plots(mds, "sex", "Sex", c("M"="skyblue", "F"="hotpink"))
mds_abx <- make_mds_plots(mds, "antibiotics", "Antibiotics", c("orangered", "midnightblue"))

plot_list <- list(mds_sex, mds_fev, mds_abx)

combo <- cowplot::plot_grid(
  plotlist = plot_list,
  rel_widths = c(1, 3, 1),
  nrow = 1
)

ggsave(combo, file = "beta_div.png", height = 5, width = 15, dpi = 300)

# adonis + mds_fev1
fev <- make_mds_plots(mds, "fev1", "FEV1", wesanderson::wes_palette("Darjeeling1"), fixed = F)

p_beta <- cowplot::plot_grid(adonis_bar, fev, align = "h", axis = "tb", rel_widths = c(2, 4))
ggsave(p_beta, file = "adonis_and_fev1_mds.png", height = 3, width = 9, dpi = 300, scale = 0.9)
  
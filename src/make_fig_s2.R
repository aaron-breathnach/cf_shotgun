make_fig_s2a <- function(vd, metadata) {
  
  bc_dists <- vd %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("sample_1") %>%
    pivot_longer(!sample_1, names_to = "sample_2", values_to = "distance") %>%
    inner_join(metadata[, c(1, 6)], by = c("sample_1" = "sample_id")) %>%
    rename(subject_1 = subject_id) %>%
    inner_join(metadata[, c(1, 6)], by = c("sample_2" = "sample_id")) %>%
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
  
  p_inp <- rbind(bc_dists_inter, bc_dists_intra) 
  
  ggplot(p_inp, aes(x = measure, y = distance, fill = subject_1)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    geom_jitter(width = 0.25, pch = 21, size = 3) +
    scale_fill_manual(values = get_pal()) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.position = "bottom") +
    labs(x = "Measure", y = "Mean Bray-Curtis distance", fill = "Patient ID")
  
}

make_fig_s2b <- function(vd, metadata) {
  
  df_md <- metadata %>%
    arrange(sample_id) %>%
    column_to_rownames("sample_id")
  
  perm <- permute::how(nperm = 999)
  permute::setBlocks(perm) <- with(df_md, subject_id)
  
  set.seed(123)
  p_inp <- vegan::adonis2(vd ~ fev1 + antibiotics + sex,
                          data = df_md,
                          by = "margin",
                          permutations = perm) %>%
    as.data.frame() %>%
    drop_na() %>%
    rownames_to_column("variable") %>%
    select(variable, 4, 6) %>%
    rename(r2 = 2, p = 3) %>%
    mutate(variable = recode(variable,
                             "fev1" = "FEV1",
                             "antibiotics" = "Antibiotics",
                             "sex" = "Sex"))
  
  ggplot(p_inp, aes(x = r2, y = reorder(variable, r2))) +
    geom_bar(stat = "identity", colour = "black", fill = "steelblue", width = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 1/3))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.position = "bottom") +
    theme(axis.title.x = ggtext::element_markdown()) +
    labs(x = "_R_<sup>2</sup>", y = "Variable")
  
}

make_fig_s2 <- function() {
  
  pal <- get_pal()
  
  metaphlan <- read_delim("data/metaphlan.tsv")
  
  vd <- metaphlan %>%
    arrange(sample_id) %>%
    column_to_rownames("sample_id") %>%
    vegan::vegdist(method = "bray")
  
  metadata <- get_metadata()
  
  fig_s2a <- make_fig_s2a(vd, metadata)
  
  fig_s2b <- make_fig_s2b(vd, metadata)
  
  set.seed(123)
  
  mds <- vegan::metaMDS(vd)$points %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(metadata, by = "sample_id")
  
  fig_s2c <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = subject_id)) +
    geom_point(colour = "black", pch = 21, size = 3) +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Patient ID")
  
  fig_s2d <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = fev1)) +
    geom_point(colour = "black", pch = 21, size = 3) +
    scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1")) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(fill = "FEV1")
  
  fig_s2e <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = antibiotics)) +
    geom_point(colour = "black", pch = 21, size = 3) +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(fill = "# antibiotics")
  
  fig_s2f <- ggplot(mds, aes(x = MDS1, y = MDS2, fill = sex)) +
    geom_point(colour = "black", pch = 21, size = 3) +
    scale_fill_manual(values = c("F" = "pink", "M" = "lightblue")) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Sex")
  
  plot_list <- list(fig_s2a, fig_s2b, fig_s2c, fig_s2d, fig_s2e, fig_s2f)
  
  patchwork::wrap_plots(plot_list)
  
}

library(tidyverse)

truncate_x <- function(x, len = 20) {
  
  x_trunc <- str_trunc(x, width = len)
  
  y <- gsub(".* ", "", x)
  
  if (str_ends(x_trunc, "\\.\\.\\.") & str_starts(y, "\\*")) {
    x_trunc <- x_trunc %>%
      str_replace("\\.\\.\\.", "") %>%
      str_c("*...")
  }
  
  return(x_trunc)
  
}

cluster_axis <- function(data) {
  
  if (nrow(data) > 1) {
    
    data %>%
      dist() %>%
      hclust() %>%
      as.dendrogram() %>%
      labels()
    
  } else {
    
    data %>%
      rownames()
    
  }
}

setwd("~/Desktop/CF_SHOTGUN/2022")

mpa <- read_delim("data/metaphlan.2022.tsv")

spp <- mpa %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  group_by(species) %>%
  summarise(relab = mean(relab)) %>%
  ungroup() %>%
  top_n(10, relab) %>%
  pull(species)

tidy_species_name <- function(species_name) {
  
  tmp <- species_name %>%
    str_replace("s__", "") %>%
    str_split("_") %>%
    unlist()
  
  sp_id <- tmp[2:length(tmp)]
  
  any_num <- any(suppressWarnings(sp_id %>%
                                    tail(-1) %>%
                                    purrr::map(function(x) !is.na(as.numeric(x))) %>%
                                    unlist()))
  
  any_unclassified <- any(lapply(sp_id, function(x) x == "unclassified") %>% unlist())
  
  if (any_num | any_unclassified) {
    
    genus <- paste0("*", tmp[1], "*")
    
    species <- paste(sp_id, collapse = " ")
    
  } else {
    
    genus <- paste0("*", str_sub(tmp[1], end = 1), ".*")
      
    species <- paste0("*", paste0(sp_id, collapse = "/"), "*")
    
  }
  
  paste(genus, species)
  
}

cor_inp <- mpa %>%
  select(1, all_of(spp)) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

colnames(cor_inp) <- lapply(colnames(cor_inp), tidy_species_name) %>%
  unlist()

corr <- Hmisc::rcorr(x = cor_inp, type = "spearman")

rval <- corr$r

x_ord <- rval %>%
  cluster_axis()

y_ord <- rval %>%
  t() %>%
  cluster_axis()

pval <- corr$P %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  pivot_longer(!x, names_to = "y", values_to = "p") %>%
  mutate(fdr = p.adjust(p, method = "BH")) %>%
  filter(fdr <= 0.05) %>%
  mutate(sig = case_when(
    fdr <= 0.001 ~ "***",
    fdr <= 0.010 & fdr > 0.001 ~ "**",
    fdr <= 0.050 & fdr > 0.010 ~ "*"
  ))

p_inp <- rval %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  pivot_longer(!x, names_to = "y", values_to = "r") %>%
  mutate(x = factor(x, levels = x_ord)) %>%
  mutate(y = factor(y, levels = y_ord))

x_labels <- lapply(levels(p_inp$x), truncate_x) %>%
  unlist()

p <- ggplot(p_inp, aes(x = x, y = y)) +
  geom_tile(aes(fill = r), colour = "black", linewidth = 0.25) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red") +
  geom_text(data = pval,
            aes(x = x, y = y, label = sig),
            size = 3.75,
            colour = "black") +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown(vjust = 0.75),
        legend.position = "top") +
  scale_x_discrete(expand = expansion(mult = c(0, 0)),
                   labels = x_labels) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  coord_equal() +
  labs(x = "", y = "", fill = "*r*")
p

ggsave(filename = "plots/correlation_heatmap.png",
       plot = p,
       width = 5,
       height = 5,
       dpi = 300)

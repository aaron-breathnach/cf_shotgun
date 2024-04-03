library(tidyverse)
library(ggtree)

R.utils::sourceDirectory("src")

fig_1 <- make_fig_1()
fig_2 <- make_fig_2()
fig_3 <- make_fig_3()
fig_4 <- make_fig_4()

fig_s1 <- make_fig_s1()
fig_s2 <- make_fig_s2()

table_2 <- make_table_2()

run_maaslin()

ggsave("figures/fig_1.png", fig_1, width = 15, height = 7.50, bg = "white")
ggsave("figures/fig_2.png", fig_2, width = 15, height = 12.5, bg = "white")
ggsave("figures/fig_3.png", fig_3, width = 15, height = 10.0, bg = "white")
ggsave("figures/fig_4.png", fig_4, width = 15, height = 10.0, bg = "white")

ggsave("figures/fig_s1.png", fig_s1, width = 6.25, height = 5.0)
ggsave("figures/fig_s2.png", fig_s2, width = 12.5, height = 7.5)

write_tsv(table_2, "data/table_2.tsv")

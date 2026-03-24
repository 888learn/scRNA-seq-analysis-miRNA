
rm(list = ls())
set.seed(123)
library(qs2); library(Seurat); library(cowplot); library(stringr)
library(ggunchull); library(ggrepel); library(tidyverse)
library(ggplot2); library(dplyr); library(ggbeeswarm); library(ggsignif)
setwd("9.T cell proliferation")
sce.all <- qs_read("../3.sce.singlet_int_celltype.qs2")

t_cell_types <- c("Cytotoxic CD8 T", "Naive T", "Treg CD4 T",
                  "Helper CD4 T", "Memory CD4 T", "NKT", "Cycling T")
sce.tcell <- subset(sce.all, subset = celltype %in% t_cell_types)

s.genes.mouse   <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes.mouse <- str_to_title(cc.genes.updated.2019$g2m.genes)
sce.tcell <- CellCycleScoring(sce.tcell,
                              s.features   = s.genes.mouse,
                              g2m.features = g2m.genes.mouse,
                              set.ident    = FALSE)

group_colors <- c("Untreated_Group1" = "#3C5488FF",
                  "NC_siR_Group1"    = "#E64B35FF",
                  "miR146a5p_Group1" = "#00A087FF")
group_levels <- c("Untreated_Group1", "NC_siR_Group1", "miR146a5p_Group1")
comparisons <- list(c("Untreated_Group1", "NC_siR_Group1"),
                    c("NC_siR_Group1",    "miR146a5p_Group1"),
                    c("Untreated_Group1", "miR146a5p_Group1"))
x_labels <- c("Untreated_Group1" = "Untreated",
              "NC_siR_Group1"    = "NC siRNA",
              "miR146a5p_Group1" = "miR-146a-5p")
plot_theme <- theme_bw(base_size = 12) +
  theme(axis.text.x        = element_text(angle = 35, hjust = 1),
        strip.text         = element_text(face = "bold"),
        strip.background   = element_rect(fill = "#F0F0F0"),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", hjust = 0.5))

s_raw <- sce.tcell@meta.data %>%
  filter(Phase == "S", orig.ident %in% group_levels) %>%
  mutate(orig.ident = factor(orig.ident, levels = group_levels))
s_df <- s_raw %>%
  group_by(orig.ident, celltype) %>%
  summarise(Mean = mean(S.Score), SEM = sd(S.Score)/sqrt(n()), .groups = "drop")

p_s <- ggplot(s_raw, aes(x = orig.ident, y = S.Score, fill = orig.ident)) +
  geom_bar(data = s_df, aes(y = Mean),
           stat = "identity", width = 0.65, alpha = 0.85) +
  geom_errorbar(data = s_df,
                aes(y = Mean, ymin = Mean - SEM, ymax = Mean + SEM),
                width = 0.2, linewidth = 0.7) +
  geom_jitter(aes(color = orig.ident), width = 0.15, size = 0.5, alpha = 0.3) +
  geom_signif(comparisons = comparisons, test = "wilcox.test",
              step_increase = 0.12, tip_length = 0.02, textsize = 3) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_x_discrete(labels = x_labels) +
  labs(y = "S.Score (Mean ± SEM)", x = "",
       title = "S Phase Score by T Cell Subtype") +
  plot_theme
p_s
ggsave(p_s, file = "Tcell_Sphase_score_bar.pdf", width = 16, height = 14)

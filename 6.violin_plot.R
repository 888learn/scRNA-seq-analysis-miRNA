library(ggpubr)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(qs2)
rm(list = ls())

sce.all <- qs_read("../../3.sce.singlet_int_celltype.qs2")
sce.all$orig.ident <- factor(sce.all$orig.ident,
                             levels = c("Untreated_Group1", "NC_siR_Group1", "miR146a5p_Group1"))

sce.all_list <- SplitObject(sce.all, split.by = "celltype")

Macrophage_cell_seob <- sce.all_list$`Mature macrophage`
Idents(Macrophage_cell_seob) <- Macrophage_cell_seob$orig.ident

Macrophage_cell_seob$orig.ident <- factor(Macrophage_cell_seob$orig.ident,
                                          levels = c("Untreated_Group1", "NC_siR_Group1", "miR146a5p_Group1"))

genes_use <- c("Itm2b", "Thbs1", "Ccl4", "Cd74")

for (gene in genes_use) {
  p_gene <- VlnPlot(
    object = Macrophage_cell_seob,
    features = gene,
    group.by = "orig.ident",
    pt.size = 0.1,
    cols = c("gray", "red", "green"))
  
  p_gene[[1]]$layers[c(2, 1)] <- p_gene[[1]]$layers[c(1, 2)]
  
  ggsave(filename = paste0(gene, "_violinplot.pdf"),
         plot = p_gene, width = 6, height = 6)
}

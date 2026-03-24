
library(scRNAtoolVis)
library(Seurat)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(qs2)
rm(list = ls())
sce.all <- qs_read("../../3.sce.singlet_int_celltype.qs2")
sce.all$orig.ident <- factor(
  sce.all$orig.ident,
  levels = c("Untreated_Group1", "NC_siR_Group1", "miR146a5p_Group1"))
Idents(sce.all) <- "celltype"
sce.all <- subset(
  sce.all,
  idents = setdiff(levels(sce.all),
                   c("LSEC", "Cycling T", "Cycling macrophage", "Neutrophil"))
)
celltypes <- levels(sce.all)
markers.list <- list()
for (ct in celltypes) {
  message("compare: ", ct)
  sce.sub <- subset(sce.all, idents = ct)
  Idents(sce.sub) <- "orig.ident"
  markers.ct <- FindMarkers(
    sce.sub,
    ident.1 = "miR146a5p_Group1",
    ident.2 = "NC_siR_Group1",
    only.pos = FALSE,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1)
  markers.ct$gene <- rownames(markers.ct)
  markers.ct$celltype <- ct
  markers.list[[ct]] <- markers.ct
}
markers.all <- do.call(rbind, markers.list)
write.csv(markers.all, "markers_by_celltype.csv", row.names = FALSE)
markers.all <- read.csv("markers_by_celltype.csv")
keep_cells <- c("Cytotoxic CD8 T", "Naive T", "Treg CD4 T", "Helper CD4 T",
                "Memory CD4 T", "B cell", "NK", "Monocytic macrophage",
                "Mature macrophage", "DC")
markers.all <- subset(markers.all, celltype %in% keep_cells)
pbmc.markers <- markers.all
pbmc.markers$celltype <- factor(pbmc.markers$celltype,
                               levels = c("Cytotoxic CD8 T", 'Naive T','Treg CD4 T','Helper CD4 T',
                                          'Memory CD4 T',"B cell", "NK",'Monocytic macrophage',
                                          'Mature macrophage',"DC"))
colour=c("#E64B35FF","#4DBBD5FF","#00A087FF", "#3C5488FF", "#F39B7FFF","#8491B4FF","#91D1C2FF",
         "#DC0000FF", "#7E6148FF","#B09C85FF","#999999", "#FFD700", "#8A2BE2")
pbmc.markers$cluster <- pbmc.markers$celltype
p <- jjVolcano( diffData = pbmc.markers,
  log2FC.cutoff = 0.5,
  col.type = "adjustP",
  topGeneN = 5, 
  tile.col = colour)  
p
ggsave("1.jjvolcano_plot.pdf", width = 14, height = 10)

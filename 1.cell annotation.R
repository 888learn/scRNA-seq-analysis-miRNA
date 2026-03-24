
rm(list = ls())
set.seed(123)
library(qs2);library(Seurat);library(cowplot);library(stringr)
library(ggunchull)
library(ggrepel)
library(tidyverse)
sce.all <- qs_read("sce.all_int_2.1singlet.qs2")

library(SingleR)
library(celldex)
ref <- ImmGenData()
expr <- GetAssayData(sce.all, layer = "data")
meta <- sce.all@meta.data
singleR_results <- SingleR(test = expr, ref = ref, labels = ref$label.main)

sce.all <- AddMetaData(sce.all, metadata = singleR_results$labels, 
                          col.name = "SingleR_labels")
table(sce.all@meta.data$SingleR_labels)

DimPlot(sce.all, reduction = "umap", group.by = "SingleR_labels", label = TRUE,
        pt.size = 0.5) + 
  ggtitle("umap of SingleR Annotated Cells") + NoLegend()
ggsave("SingleR_umap.pdf", width = 8, height = 6)


Idents(sce.all) <- sce.all$RNA_snn_res.2.1
library(COSG)
marker_cosg <- cosg(
  sce.all,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1,
  n_genes_user = 100)

save(marker_cosg, file = "COSG_marker_result.Rdata")
write.csv(marker_cosg$names, file = "COSG_marker_genes.csv")

top_3 <- unique(as.character(apply(marker_cosg$names, 2, head, 3)))
DotPlot(sce.all, features = top_3) + coord_flip()
ggsave("COSG_DotPlot_Top3_markers.pdf", width = 12, height = 20)


celltype=data.frame(ClusterID=0:43, celltype= NA) 

celltype[celltype$ClusterID %in% c(5,7,13,14,23,37),2]= 'Monocytic macrophage'
celltype[celltype$ClusterID %in% c(0, 9,19,26, 36,40,41),2]='B cell'
celltype[celltype$ClusterID %in% c(29 ),2]='Cytotoxic CD8 T'
celltype[celltype$ClusterID %in% c(12 ),2]='Helper CD4 T'
celltype[celltype$ClusterID %in% c(21 ),2]='Memory CD4 T'
celltype[celltype$ClusterID %in% c(1,10,22,32 ),2]='NK'
celltype[celltype$ClusterID %in% c(8),2]='NKT'
celltype[celltype$ClusterID %in% c(6,20,27,30, 34,38,43),2]='DC'
celltype[celltype$ClusterID %in% c(2,3,4,15,17,18,33 ),2]='Mature macrophage'
celltype[celltype$ClusterID %in% c(11 ),2]='Neutrophil'
celltype[celltype$ClusterID %in% c(16,24),2]='Naive T'
celltype[celltype$ClusterID %in% c(28),2]='Treg CD4 T'
celltype[celltype$ClusterID %in% c(25),2]='Cycling macrophage'
celltype[celltype$ClusterID %in% c(31),2]='KC'
celltype[celltype$ClusterID %in% c(35,42),2]='Cycling T'
celltype[celltype$ClusterID %in% c(39),2]='LSEC'

celltype
table(celltype$celltype)

sce.all@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.2.1 == celltype$ClusterID[i]),
                    'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)

qs_save(sce.all, "sce.singlet_int_celltype.qs2")


sce.all <- qs_read("sce.singlet_int_celltype.qs2")
remove_types <- c("LSEC", "Cycling T", "Cycling macrophage",'Neutrophil')
sce.all@meta.data$celltype <- as.character(sce.all@meta.data$celltype)
sce.all <- subset(sce.all, cells = rownames(sce.all@meta.data)[!(sce.all@meta.data$celltype %in% remove_types)])

sce.all$celltype <- factor(sce.all$celltype)
sce.all$celltype <- droplevels(sce.all$celltype)

Idents(sce.all) <- sce.all@meta.data$celltype
cell_colors <- c("Cytotoxic CD8 T"     = "#E64B35FF",
                 "Naive T"             = "#4DBBD5FF",
                 "Treg CD4 T"          = "#00A087FF",
                 "Helper CD4 T"        = "#3C5488FF",
                 "Memory CD4 T"        = "#F39B7FFF",
                 "B cell"              = "#8491B4FF",
                 "NK"                  = "#91D1C2FF",
                 "NKT"                 = "#DC0000FF",
                 "Monocytic macrophage"= "#7E6148FF",
                 "Mature macrophage"   = "#B09C85FF",
                 "KC"                  = "#999999",
                 "DC"                  = "#FFD700")

tsnecelltype <- DimPlot(sce.all, reduction = "tsne", group.by = "celltype", pt.size = 0.1,
                        label = T,label.box = T,repel = TRUE) + 
  theme(legend.position = "none") +
  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values = cell_colors)
tsnecelltype

umapcelltype <- DimPlot(sce.all, reduction = "umap", group.by = "celltype", pt.size = 0.1,
                        label = T,label.box = T,repel = TRUE ) +
  theme(legend.position = "none") +
  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values = cell_colors)
umapcelltype

ggsave(tsnecelltype, file = "tsne_celltype.pdf", width = 10, height = 8)
ggsave(umapcelltype, file = "umap_celltype.pdf", width = 10, height = 8)

table(sce.all@meta.data$celltype)

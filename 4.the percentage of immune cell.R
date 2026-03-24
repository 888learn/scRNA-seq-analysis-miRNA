
library(qs2)
library(reshape2)
library(Seurat)
rm(list = ls())
sce.all <- qs_read("sce.singlet_int_celltype.qs2")

sce.all <- subset(sce.all, subset = orig.ident %in% 
                             c("Untreated_Group1", "NC_siR_Group1","miR146a5p_Group1"))
sce.all$orig.ident <- factor(sce.all$orig.ident, levels = c("Untreated_Group1", "NC_siR_Group1", "miR146a5p_Group1"))

library(dplyr)

cell_type_by_sample <- sce.all@meta.data %>%
  group_by(orig.ident, celltype) %>% 
  summarise(CellNumber = n()) %>%
  ungroup()    

cell_type_by_sample <- cell_type_by_sample %>%
  filter(celltype != "Cycling macrophage", celltype != "Cycling T", celltype != "LSEC",
         celltype != "Neutrophil")

cell_type_by_sample <- cell_type_by_sample %>%
  group_by(orig.ident) %>%
  mutate(Percentage = round(CellNumber / sum(CellNumber) * 100, 2))

cell_type_by_sample
  
write.csv(cell_type_by_sample, file = "immunecell_type_by_sample.csv", row.names = FALSE)

cell_type_by_sample$cell_type <- factor(cell_type_by_sample$celltype, levels = c(
  "Cytotoxic CD8 T", 'Naive T','Treg CD4 T','Helper CD4 T',
  'Memory CD4 T',"B cell", "NK", 'NKT','Monocytic macrophage',
  'Mature macrophage','KC',"DC"))  

library(ggplot2)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)
extra_colors <- c("#999999", "#FFD700", "#8A2BE2")
my_colors <- c(npg_colors, extra_colors)
length(unique(cell_type_by_sample$cell_type))
length(my_colors)

p <- ggplot(cell_type_by_sample, aes(x = orig.ident, y = Percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  ylab("Cell Type Percentage (%)") + xlab("Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 20, hjust = 1)) +
  guides(fill = guide_legend(title = "Cell Type")) +
  scale_fill_manual(values = my_colors)
print(p)

ggsave("cell_type_percentage1_deletecycling.pdf", plot = p, width = 8, height = 12)

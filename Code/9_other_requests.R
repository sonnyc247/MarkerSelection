#### average expression ####

hodge_average <- AverageExpression(new_Seu_AIBS_obj_for_test)
hodge_average <- hodge_average$RNA
hodge_average <- hodge_average[new_CgG_results$gene,]
hodge_average <- tibble::rownames_to_column(hodge_average)
colnames(hodge_average) <- paste0(colnames(hodge_average), "_mean_expression")
colnames(hodge_average)[1] <- "gene"
new_CgG_results_pairdata <- merge(new_CgG_results_pairdata, hodge_average, by = "gene")

write.csv(new_CgG_results_pairdata, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_CgG_results_pairdata.csv", row.names = FALSE) #save/export results

#### Jordan marker filtering ####

Jordan_CTX <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Inputs/SUPPL_S1.csv", stringsAsFactors=FALSE)
Jordan_CTX_filtered <- Jordan_CTX[Jordan_CTX$avg_logFC > 1,]
length(unique(Jordan_CTX_filtered$gene))
duplicated_genes <- unique(Jordan_CTX_filtered[duplicated(Jordan_CTX_filtered$gene), "gene"])
Jordan_CTX_filtered <- Jordan_CTX_filtered[!(Jordan_CTX_filtered$gene %in% duplicated_genes),]

write.csv(Jordan_CTX_filtered, "Jordan_Mouse_Markers_Filtered.csv")

#### Annotating cell types and colours ####

Celltypes_and_colours <- new_Seu_AIBS_obj@meta.data
Celltypes_and_colours <- Celltypes_and_colours[,c("subclass_label", "subclass_color", "class_label", "class_color")]
Celltypes_and_colours <- unique(Celltypes_and_colours)
row.names(Celltypes_and_colours) <- NULL
Celltypes_and_colours$Our_label <- Celltypes_and_colours$subclass_label
Celltypes_and_colours[Celltypes_and_colours$class_label == "GABAergic","Our_label"] <- paste0("Inh_", Celltypes_and_colours[Celltypes_and_colours$class_label == "GABAergic","Our_label"])
Celltypes_and_colours[Celltypes_and_colours$class_label == "Glutamatergic","Our_label"] <- paste0("Exc_", Celltypes_and_colours[Celltypes_and_colours$class_label == "Glutamatergic","Our_label"])
Celltypes_and_colours <- Celltypes_and_colours[,c(5,1:4)]
names(Celltypes_and_colours)[2:5] <- paste0("AIBS_", names(Celltypes_and_colours)[2:5])

write.csv(Celltypes_and_colours, "Name_and_colour_scheme.csv")

#### PVALB cells plot - PVALB and SST genes in different brain regions ####

Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj.rds")
names(Seu_AIBS_obj@meta.data)
Idents(Seu_AIBS_obj) <- "subclass_label"
Seu_AIBS_obj <- subset(Seu_AIBS_obj, idents = "PVALB", invert = F) #take PVALB cells only
names(Seu_AIBS_obj@meta.data)
Idents(Seu_AIBS_obj) <- "region_label"

PVALB_Matrix <- as.data.frame(Seu_AIBS_obj@assays$RNA@data)
PVALB_Matrix <- PVALB_Matrix[c("PVALB", "SST"),]
PVALB_Matrix <- t(PVALB_Matrix)
PVALB_Matrix <- as.data.frame(PVALB_Matrix)
PVALB_Meta <- as.data.frame(Seu_AIBS_obj@meta.data)

PVALB_graph <- merge(PVALB_Matrix, PVALB_Meta, by = "row.names")
unique(PVALB_graph[,c("region_label", "region_color")])

library(ggplot2)

#scale_fill_manual(values=c("MTG" = "#5CCCCC",
#                           "V1C" = "#CC0099",
#                           "CgG" = "#CCA83D",
#                           "M1lm" = "#00FF40",
#                           "S1ul" = "#2E4999",
#                           "S1lm" = "#9326FF",
#                           "M1ul" = "#589917",
#                           "A1C"  = "#FF7373")) +

ggplot(PVALB_graph[PVALB_graph$region_label %in% c("MTG", "CgG", "V1C"),], aes(x=region_label, y=PVALB, fill=region_label)) + 
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.3, size = 0.6) +
  ylab('Normalized PVALB mRNA expression (CPM)') +
  xlab('Brain region from AIBS snRNAseq SMART-Seq v4 dataset') +
  scale_fill_manual(values=c("MTG" = "#5CCCCC",
                             "V1C" = "#CC0099",
                             "CgG" = "#CCA83D")) +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size=8), 
        axis.title.x = element_text(size = 10, margin = margin(t = 7, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 7, b = 0, l = 0)),
        axis.text.y = element_text(size = 8)) 


ggsave(dpi = 300, 
       units = "mm", 
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Misc",
       filename = "PVALB.jpg",
       device = "jpeg")

#### Boxplots for Dan's rosmap presentation ####

library(Seurat)
library(ggplot2)
library(cowplot)

### load objects

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")
table(Seu_plot_object$NeuN_Region)
Seu_plot_object <- subset(Seu_plot_object, subset = NeuN_Region == "MTG_Neuronal", invert = TRUE)

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_22MAY21.rds") #load cain seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_22MAY21.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

### get metadata

metadata_for_plot <- Seu_plot_object@meta.data
metadata_for_plot <- tibble::rownames_to_column(metadata_for_plot)

# for hodge
metadata_for_plot <- metadata_for_plot[, c("rowname", "nFeature_RNA", "subclass_label")] 
colnames(metadata_for_plot)[3] <- "subclass"
metadata_for_plot$Dataset <- "Hodge"

# for snDatasets
metadata_for_plot <- metadata_for_plot[, c("rowname", "nFeature_RNA", "predicted.id.CgG")] 
colnames(metadata_for_plot)[3] <- "subclass"
metadata_for_plot$Dataset <- "Zhou" #as appropriate

#metadata_for_plot_combined <- metadata_for_plot # first time
metadata_for_plot_combined <- rbind(metadata_for_plot_combined, metadata_for_plot)

table(metadata_for_plot_combined$Dataset)

# plot

metadata_for_plot_combined$Dataset_backup <- metadata_for_plot_combined$Dataset
metadata_for_plot_combined$Dataset <- factor(metadata_for_plot_combined$Dataset, levels=c("Hodge", "Mathys", "Cain", "Zhou"), ordered=TRUE)

all_plots <- ggplot(metadata_for_plot_combined, aes(Dataset, nFeature_RNA)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5) +
  theme_classic() +
  labs(y = "Number of detected genes")

SST_plots <- ggplot(metadata_for_plot_combined[metadata_for_plot_combined$subclass == "SST",], aes(Dataset, nFeature_RNA)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5) +
  theme_classic() +
  labs(y = "Number of detected genes")

plot_grid(all_plots,
          SST_plots,
          labels = c("All subclasses", "SST subclasss cells"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

SST_holder <- metadata_for_plot_combined[metadata_for_plot_combined$subclass == "SST",]
metadata_for_plot_combined$Plot <- "All subclasses"
SST_holder$Plot <- "SST cells"
metadata_for_plot_dup <- rbind(metadata_for_plot_combined, SST_holder)

metadata_for_plot_dup$tech <- "10X Genomics"
metadata_for_plot_dup[metadata_for_plot_dup$Dataset == "Hodge", "tech"] <- "SMART-seq"
metadata_for_plot_dup$tech <- factor(metadata_for_plot_dup$tech, levels=c("SMART-seq", "10X Genomics"), ordered=TRUE)

ggplot(metadata_for_plot_dup, aes(Dataset, nFeature_RNA, fill = tech)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5) +
  scale_fill_manual(values = c("cyan", "white")) +
  theme_classic() +
  theme(legend.margin=margin(0, 0, 0, -10),
        axis.title.y=element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x=element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  labs(y = "Number of detected genes", fill = "Sequencing \ntechnology") +
  facet_wrap(~Plot)

ggsave(width = 180,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Misc/",
       filename = "Genetech_Boxplots.pdf",
       device = "pdf")

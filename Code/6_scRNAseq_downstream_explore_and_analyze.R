#### All involved packages ####

library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(cowplot)
library(ggbeeswarm)
library(broom)

#
#### Confusion matrix cgg vs m1 ####

# qc thresh cut for IT cells

ncol(Seu_plot_object)
table(Seu_plot_object$predicted.id)
Seu_plot_object <- subset(Seu_plot_object, subset = predicted.id == "IT" & prediction.score.max < 0.8, invert = T) 
table(Seu_plot_object$predicted.id)

# get numbers for confusion matrix

df_for_CgGvM1_plot <- merge(predictions_Zhou_CgG, predictions_Zhou_M1, by = "row.names")
df_for_CgGvM1_plot <- merge(predictions_Mathys_CgG[!(predictions_Mathys_CgG$prediction.score.max < 0.8
                                                 & predictions_Mathys_CgG$predicted.id == "IT"),], 
                                                 predictions_Mathys_M1, 
                                                 by = "row.names",
                                                 all.x = T,
                                                 all.y = F)

confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$predicted.id.y, df_for_CgGvM1_plot$predicted.id.x)) 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

#heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

mathys_M1pct <- ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("CgG-mapped subclass") + 
  ylab("M1-mapped subclass") +
  labs(fill = "% of M1- \n mapped cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

plot_grid(mathys_M1pct,
          cain_M1pct,
          zhou_M1pct,
          labels = c("Mathys", "Cain", "Zhou"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

### export plot (most recently plotted)

ggsave(width = 180,
       height = 400,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "CgGvsM1_M1pct_IT8filter.pdf",
       device = "pdf")

#### Confusion matrix original vs m1 ####

# load data

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only

# get numbers for confusion matrix

metadata_for_plot <- Seu_map_object@meta.data

metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "Subcluster")] # for Mathys
identical(rownames(metadata_for_plot), rownames(predictions_Mathys_M1))
df_for_CgGvM1_plot <- merge(metadata_for_plot, predictions_Mathys_M1, by = "row.names")

table(metadata_for_plot$subtype) # for Cain
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "subtype")] 
identical(rownames(metadata_for_plot), rownames(predictions_Cain_M1))
df_for_CgGvM1_plot <- merge(metadata_for_plot, predictions_Cain_M1, by = "row.names")
df_for_CgGvM1_plot <- merge(predictions_Cain_M1[!(predictions_Cain_M1$prediction.score.max < 0.8
                                                & predictions_Cain_M1$predicted.id == "L2/3 IT"),], 
                            metadata_for_plot, 
                            by = "row.names",
                            all.x = T,
                            all.y = F)


confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$Subcluster, df_for_CgGvM1_plot$predicted.id)) # for mathys 
confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$subtype, df_for_CgGvM1_plot$predicted.id)) # for cain 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

#heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

cain_M1pct_08 <- ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("Mapped subclass") + 
  ylab("Pre-mapping cell group") +
  labs(fill = "% of pre- \n mapping \n cell group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

plot_grid(mathys_M1pct,
          cain_M1pct,
          cain_M1pct_08,
          labels = c("Mathys", "Cain", "Cain L2/3 IT Filter .8"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

### export plot (most recently plotted)

ggsave(width = 180,
       height = 550,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "OrigvsM1.pdf",
       device = "pdf")

#### Confusion matrix original vs allhodgeregions_expandedITwithL3/5 ####

# load data

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_08Jun21.rds") #load Mathys Seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_11JUN21.rds") #load Cain Seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for Cain object only

# get numbers for confusion matrix

metadata_for_plot <- Seu_map_object@meta.data

table(metadata_for_plot$Subcluster) # for Mathys
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "Subcluster", "predicted.id.AllHodge_ExpSubclas")] # for Mathys

table(metadata_for_plot$subtype) # for Cain
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "subtype", "predicted.id.AllHodge_ExpSubclas")] 

confusion_martix_hold <- as.matrix(table(metadata_for_plot$Subcluster, metadata_for_plot$predicted.id.AllHodge_ExpSubclas)) # for mathys 
confusion_martix_hold <- as.matrix(table(metadata_for_plot$subtype, metadata_for_plot$predicted.id.AllHodge_ExpSubclas)) # for cain 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("Mapped subclass") + 
  ylab("Pre-mapping cell group") +
  labs(fill = "% of pre- \n mapping \n cell group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

### export plot (most recently plotted)

ggsave(width = 180,
       height = 300,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "OrigvsAllHodgeExpIT_Cain.pdf",
       device = "pdf")

#
#### Confusion matrix original vs mtgncgg_top20pct ####

# load data

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_26Jul21.rds") #load Mathys Seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_27JUL21.rds") #load Cain Seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for Cain object only

# get numbers for confusion matrix

metadata_for_plot <- Seu_map_object@meta.data

table(metadata_for_plot$Subcluster) # for Mathys
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "Subcluster", "predicted.id.MTGnCgG_20Pct", "qc_passing.MTGnCgG_20Pct")] # for Mathys

table(metadata_for_plot$subtype) # for Cain
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "subtype", "predicted.id.MTGnCgG_20Pct", "qc_passing.MTGnCgG_20Pct")] 

metadata_for_plot <- metadata_for_plot[metadata_for_plot$qc_passing.MTGnCgG_20Pct == T, ] # as appropriate

confusion_martix_hold <- as.matrix(table(metadata_for_plot$Subcluster, metadata_for_plot$predicted.id.MTGnCgG_20Pct)) # for mathys 
confusion_martix_hold <- as.matrix(table(metadata_for_plot$subtype, metadata_for_plot$predicted.id.MTGnCgG_20Pct)) # for cain 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

Cain_filtered <- ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("Mapped subclass") + 
  ylab("Pre-map cell group by Cain et al.") +
  labs(fill = "% of pre- \n map cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

plot_grid(Mathys_unfiltered,
          Cain_unfiltered,
          labels = c("a)", "b)"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

### export plot (most recently plotted)

ggsave(width = 180,
       height = 360,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "MTGnCgG_vs_Orig_Unfiltered.pdf",
       device = "pdf")

#
#### lm of snCTPs and for rosmap factors (CgG Unfiltered; pilot/initial use case) ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
mathys_cell_type_prop_df$dataset <- "Mathys"
zhou_cell_type_prop_df$dataset <- "Zhou"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) 

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, zhou_cell_type_prop_df, cain_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'Sst', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check subclasses and datasets
cell_type_list <- as.character(unique(ad_snrnaseq_df$subclass))
table(ad_snrnaseq_df$subclass, ad_snrnaseq_df$dataset)

cell_type_list = cell_type_list[! cell_type_list %in% c('L5/6 IT Car3')] #remove subclasses that are not in all datasets; CgG map
cell_type_list = cell_type_list[! cell_type_list %in% c('L5 IT', 'L5/6 NP', 'L6 IT', 'L6 IT Car3', 'L6b', 'Sst Chodl', 'VLMC', 'L5 ET')] #remove subclasses that are not in all datasets; M1 map

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD'), subclass %in% cell_type_list) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, "subclass"]
cell_type_list = cell_type_list[! cell_type_list %in% c('L6b')] #remove subclasses that are pure 0 in 1 dataset; CgG map
remove(test)

# calculate beta coefficients

beta_coefs_non_meta_df = lapply(levels(ad_snrnaseq_df$dataset), function(curr_dataset){
  print(curr_dataset)
  
  dataset_df = lapply(cell_type_list, function(curr_cell_type){
    print(curr_cell_type)
    df = ad_snrnaseq_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    
    return(model_df) 
  }) %>% bind_rows()
  return(dataset_df)  
  #beta_out = 
}) %>% bind_rows()

# for troubleshooting the above function (basically convert the above into for loops):

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = ad_snrnaseq_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
  }
  
}

### Plot results

# sets factor levels to match order we want

beta_coefs_non_meta_df$dataset = factor(beta_coefs_non_meta_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
beta_coefs_non_meta_df = merge(beta_coefs_non_meta_df, subclass_meta, by.x = 'subclass', by.y = 'subclass') #cgg mapping
Seu_M1AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_M1AIBS_obj.rds")
subclass_meta_M1 <- unique(Seu_M1AIBS_obj@meta.data[c("subclass_label", "class_label", "class_color")])
remove(Seu_M1AIBS_obj)
colnames(subclass_meta_M1)[1:2] <- c("subclass", "class")
beta_coefs_non_meta_df = merge(beta_coefs_non_meta_df, subclass_meta_M1, by.x = 'subclass', by.y = 'subclass') #M1 mapping

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_wrap(~subclass) + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('')

### save results

beta_coefs_non_meta_df$mapping <- "M1"
beta_coefs_non_meta_df$filtering <- "Filtered"

beta_coefs_non_meta_df_CgG_Unfiltered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_CgG_Filtered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_M1_Unfiltered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_M1_Filtered <- beta_coefs_non_meta_df

### combine results

beta_coefs_non_meta_df_M1_Filtered$subclass <-toupper(beta_coefs_non_meta_df_M1_Filtered$subclass)
beta_coefs_non_meta_df_M1_Unfiltered$subclass <-toupper(beta_coefs_non_meta_df_M1_Unfiltered$subclass)

identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_M1_Unfiltered))
identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_M1_Filtered))
identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_CgG_Filtered))

beta_coefs_non_meta_df_combined = bind_rows(beta_coefs_non_meta_df_CgG_Unfiltered, 
                                            beta_coefs_non_meta_df_CgG_Filtered,
                                            beta_coefs_non_meta_df_M1_Unfiltered,
                                            beta_coefs_non_meta_df_M1_Filtered)
beta_coefs_non_meta_df_combined[beta_coefs_non_meta_df_combined$class == "Non-neuronal", "class"] <- "Non-Neuronal"
beta_coefs_non_meta_df_combined$filtering = factor(beta_coefs_non_meta_df_combined$filtering, levels = c('Unfiltered', 'Filtered'))

### combined graph

beta_coefs_non_meta_df_combined %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(subclass ~ mapping*filtering, scales = "free") + 
  #facet_wrap(~mapping*filtering, scales = "free") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('')


#### template for All (even if a group is not mapped in a given dataset) mapped identity proportions' lm analysis ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### lm 

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
zhou_cell_type_prop_df$dataset <- "Zhou"
mathys_cell_type_prop_df$dataset <- "Mathys"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) #checks before rbind all together

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, cain_cell_type_prop_df, zhou_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Cain', 'Zhou'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'Sst', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD')) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, c("subclass", "dataset")]
remove(test)

# for troubleshooting the above function (basically convert the above into for loops):

beta_coefs_non_meta_df_ALL <- model_df[0,] # first time (run the loop below, get model df as a template)

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  temp_df <- ad_snrnaseq_df[ad_snrnaseq_df$dataset == curr_dataset, ]
  cell_type_list <- unique(temp_df$subclass)
  
  #if (curr_dataset == "Zhou"){
  #cell_type_list <- cell_type_list[cell_type_list != "L6b"]
  #}  #CgG only
  
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = temp_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    model_df$mapping = "M1" # as appropriate
    model_df$filtering = "Filtered" # as appropriate
    beta_coefs_non_meta_df_ALL <- bind_rows(beta_coefs_non_meta_df_ALL, model_df)
  }
  
}

beta_coefs_non_meta_df_ALL$mapping <-"CgG" # after first run
beta_coefs_non_meta_df_ALL$filtering <-"Filtered" # after first run

beta_coefs_non_meta_df_backup <- beta_coefs_non_meta_df_ALL # as desired
table(beta_coefs_non_meta_df_ALL$mapping, beta_coefs_non_meta_df_ALL$filtering, exclude = "ifany")

### Plot results

# final adjustments to metadata/factors

beta_coefs_non_meta_df_ALL$dataset = factor(beta_coefs_non_meta_df_ALL$dataset, levels = c('Mathys', 'Cain', "Zhou"))
beta_coefs_non_meta_df_ALL$filtering = factor(beta_coefs_non_meta_df_ALL$filtering, levels = c('Unfiltered', 'Filtered'))
subclass_meta_total <- bind_rows(subclass_meta, subclass_meta_M1)
subclass_meta_total[subclass_meta_total$class == "Non-neuronal", "class"] <- "Non-Neuronal"
beta_coefs_non_meta_df_ALL = merge(beta_coefs_non_meta_df_ALL, unique(subclass_meta_total[,c("subclass", "class")]), by.x = 'subclass', by.y = 'subclass') #cgg mapping

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = subclass, y = estimate, fill = filtering)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("blue", "cyan")) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(dataset ~ mapping*filtering, scales = "free", space = "free") + 
  #facet_wrap(~dataset*mapping*filtering, scales = "free_x") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(subclass ~ mapping*filtering, scales = "free", space = "free") + 
  #facet_wrap(~dataset*mapping*filtering, scales = "free_x") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

#save

write.csv(beta_coefs_non_meta_df_combined, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/LOADAD_lm_results/beta_coefs_non_meta_df_4in1_common_subclasses.csv")

#
#### lm of all-hodge (brain regions) snCTPs and for rosmap factors ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/AllHodge_Unfiltered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/AllHodge_Unfiltered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/AllHodge_Unfiltered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
mathys_cell_type_prop_df$dataset <- "Mathys"
zhou_cell_type_prop_df$dataset <- "Zhou"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) 

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, zhou_cell_type_prop_df, cain_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'SST', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD')) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, c("subclass", "dataset")]
remove(test)

# for troubleshooting the above function (basically convert the above into for loops):

beta_coefs_non_meta_df_ALL <- beta_coefs_non_meta_df[0,c(2:7,1)] # first time

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  temp_df <- ad_snrnaseq_df[ad_snrnaseq_df$dataset == curr_dataset, ]
  cell_type_list <- unique(temp_df$subclass)
  
  #if (curr_dataset == "Zhou"){
  #cell_type_list <- cell_type_list[cell_type_list != "L6b"]
  #}  #CgG only
  
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = temp_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$subclass == curr_cell_type, ]
    if (mean(df$cell_type_proportion == 0)) {
      next
    }
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    model_df$mapping = "M1" # as appropriate
    model_df$filtering = "Filtered" # as appropriate
    beta_coefs_non_meta_df_ALL <- bind_rows(beta_coefs_non_meta_df_ALL, model_df)
  }
  
}

beta_coefs_non_meta_df_ALL$mapping <-"All_Hodge" # after first run
beta_coefs_non_meta_df_ALL$filtering <-"Unfiltered" # after first run

### Plot results

# final adjustments to metadata/factors

beta_coefs_non_meta_df_ALL$dataset = factor(beta_coefs_non_meta_df_ALL$dataset, levels = c('Mathys', 'Cain', "Zhou"))
Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds")
subclass_class_holder <- unique(Seu_ref_object[[c("subclass_label_expanded_L35IT", "class_label")]])
names(subclass_class_holder) <- c("subclass", "class")
beta_coefs_non_meta_df_ALL = merge(beta_coefs_non_meta_df_ALL, subclass_class_holder, by.x = 'subclass', by.y = 'subclass') 

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD') %>% 
  ggplot(aes(x = subclass, y = estimate, fill = class)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("red", "blue", "grey")) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(dataset ~ class, scales = "free", space = "free") + 
  #facet_wrap(~dataset*mapping*filtering, scales = "free_x") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, margin=margin(t=15)))
  
#save

write.csv(beta_coefs_non_meta_df_ALL, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/LOADAD_lm_results/beta_coefs_non_meta_df_allhodge_unfiltered.csv")







#### lm of mapped identity proportions' (MTGnCgG top 20 pct most var genes) with rosmap factors ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/MTGnCgG_20Pct/Filtered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/MTGnCgG_20Pct/Filtered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/MTGnCgG_20Pct/Filtered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### lm 

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
zhou_cell_type_prop_df$dataset <- "Zhou"
mathys_cell_type_prop_df$dataset <- "Mathys"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) #checks before rbind all together

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, cain_cell_type_prop_df, zhou_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Cain', 'Zhou'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data for SST
ad_snrnaseq_df %>% 
  filter(subclass == 'SST', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check for cases where all counts of a cell type = zero (mean would = exactly 0) for C and AD LOAD values
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD')) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, c("subclass", "dataset")]
remove(test)

# for troubleshooting the above function (basically convert the above into for loops):

beta_coefs_non_meta_df_ALL <- model_df[0,] # first time (run the loop below, get model df as a template)

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  temp_df <- ad_snrnaseq_df[ad_snrnaseq_df$dataset == curr_dataset, ]
  cell_type_list <- unique(temp_df$subclass)

  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = temp_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    model_df$mapping = "MTGnCgG_Top20PctMostVar" # as appropriate
    model_df$filtering = "Unfiltered" # as appropriate
    beta_coefs_non_meta_df_ALL <- bind_rows(beta_coefs_non_meta_df_ALL, model_df)
  }
  
}

beta_coefs_non_meta_df_backup <- beta_coefs_non_meta_df_ALL # as desired
table(beta_coefs_non_meta_df_ALL$dataset, beta_coefs_non_meta_df_ALL$filtering, exclude = "ifany")

### Plot results

# final adjustments to metadata/factors

beta_coefs_non_meta_df_ALL$dataset = factor(beta_coefs_non_meta_df_ALL$dataset, levels = c('Mathys', 'Cain', "Zhou"))
beta_coefs_non_meta_df_ALL$filtering = factor(beta_coefs_non_meta_df_ALL$filtering, levels = c('Unfiltered', 'Filtered'))

Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds")
subclass_meta <- unique(Seu_ref_object@meta.data[,c("subclass_label", "class_label")])
length(intersect(beta_coefs_non_meta_df_ALL$subclass, subclass_meta$subclass_label)) #check
remove(Seu_ref_object)
beta_coefs_non_meta_df_ALL = merge(beta_coefs_non_meta_df_ALL, subclass_meta, by.x = 'subclass', by.y = 'subclass_label') #cgg mapping
colnames(beta_coefs_non_meta_df_ALL)[10] <- "class"

#save

write.csv(beta_coefs_non_meta_df_ALL, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/LOADAD_lm_results/beta_coefs_non_meta_df_MTGnCgG_20Pct.csv")

### plots beta coefficients for LOAD across each dataset faceted by cell type

Unfiltered_lm <- beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD', filtering == "Unfiltered") %>% 
  ggplot(aes(x = subclass, y = estimate, fill = class)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("red", "blue", "grey")) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(dataset ~ class, scales = "free", space = "free") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

# combined plot

plot_grid(Filtered_lm,
          Unfiltered_lm,
          labels = c("a)", "b)"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

# export plot (most recently plotted)

ggsave(width = 360,
       height = 360,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Lm_betas/",
       filename = "MTGnCgG20Pct_LOADADbeta.pdf",
       device = "pdf")

#
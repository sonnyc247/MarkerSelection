#### packages ####

library(Seurat)
library(dplyr)
library(MAST)
library(metap)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(tidyr)
library(reshape2)

#### DE expression testing ####

### Taking single-cell/nucleus expression datasets, we test for DE gene expression for only the markers that we are trying to validate ###

### prep/set-up

Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_test_object <- subset(Seu_test_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_test_object <- subset(Seu_test_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

names(Seu_test_object@meta.data)
Idents(Seu_test_object) <- "predicted.id"
table(Idents(Seu_test_object))

Result_df_MTGandCgG_lfct2.0 <- new_MTGnCgG_lfct2_results[,2:12]

marker_list <- intersect(Result_df_MTGandCgG_lfct2.0$gene, row.names(Seu_test_object))

### do MAST test of markers for all cell types 

validation_mast_results <- FindAllMarkers(Seu_test_object, 
                                          slot = "data", 
                                          logfc.threshold = 0, 
                                          min.pct = 0, 
                                          only.pos = FALSE, 
                                          return.thresh = 1,
                                          features = marker_list,
                                          test.use = "MAST") #find markers

### process and save output

validation_mast_results$group_gene <- paste0(validation_mast_results$cluster, "_", validation_mast_results$gene)

validation_results_Mathys <- validation_mast_results
validation_results_Cain <- validation_mast_results
validation_results_Zhou <- validation_mast_results

### validation

Result_df_MTGandCgG_lfct2.0$In_Mathys <- FALSE
Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$gene %in% validation_results_Mathys$gene, "In_Mathys"] <- TRUE

Result_df_MTGandCgG_lfct2.0$In_Cain <- FALSE
Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$gene %in% validation_results_Cain$gene, "In_Cain"] <- TRUE

Result_df_MTGandCgG_lfct2.0$In_Zhou <- FALSE
Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$gene %in% validation_results_Zhou$gene, "In_Zhou"] <- TRUE

Result_df_MTGandCgG_lfct2.0$group_gene <- paste0(Result_df_MTGandCgG_lfct2.0$subclass, "_", Result_df_MTGandCgG_lfct2.0$gene)

intersect(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Mathys == TRUE, "group_gene"],
          validation_results_Mathys$group_gene) %>% length()
setdiff(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Mathys == TRUE, "group_gene"],
        validation_results_Mathys$group_gene)
setdiff(Result_df_MTGandCgG_lfct2.0$subclass, validation_results_Mathys$cluster)

intersect(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Cain == TRUE, "group_gene"],
          validation_results_Cain$group_gene) %>% length()
sort(setdiff(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Cain == TRUE, "group_gene"],
        validation_results_Cain$group_gene))
setdiff(Result_df_MTGandCgG_lfct2.0$subclass, validation_results_Cain$cluster)

intersect(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Zhou == TRUE, "group_gene"],
          validation_results_Zhou$group_gene) %>% length()
outlier_markers <- setdiff(Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$In_Zhou == TRUE, "group_gene"],
        validation_results_Zhou$group_gene)
setdiff(Result_df_MTGandCgG_lfct2.0$subclass, validation_results_Zhou$cluster)

retest <- Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$group_gene %in% outlier_markers, ]
retest <- retest[retest$subclass != "L4 IT",]
retest <- retest$gene

evaluation_df <- validation_results_Cain[validation_results_Cain$group_gene %in% Result_df_MTGandCgG_lfct2.0$group_gene,]
evaluation_df <- validation_results_Mathys[validation_results_Mathys$group_gene %in% Result_df_MTGandCgG_lfct2.0$group_gene,]
evaluation_df <- validation_results_Zhou[validation_results_Zhou$group_gene %in% Result_df_MTGandCgG_lfct2.0$group_gene,]

### retest zhou cells (looking for why some markers don't show up - maybe due to things like min.cells.feature, for example)

validation_mast_results <- FindAllMarkers(Seu_test_object, 
                                          slot = "data", 
                                          logfc.threshold = -Inf, 
                                          min.pct = -Inf, 
                                          only.pos = FALSE, 
                                          return.thresh = 100,
                                          features = retest,
                                          min.cells.feature = 0,
                                          min.cells.group = 0,
                                          test.use = "MAST") #find markers

cluster.averages <- AverageExpression(Seu_test_object)
cluster.averages <- cluster.averages$RNA
cluster.averages <- cluster.averages[retest,]

#### meta-p value combination of results from DE expression testing ####

validation_mast_results <- validation_results_Mathys
validation_mast_results <- validation_results_Cain
validation_mast_results <- validation_results_Zhou


validation_mast_results <- validation_mast_results[validation_mast_results$group_gene %in% Result_df_MTGandCgG_lfct2.0$group_gene,]

table(validation_mast_results$cluster)

metap_holder <- as.data.frame(table(validation_mast_results$cluster))
metap_holder$metap <- 0

for (cellgroup in unique(validation_mast_results$cluster)){
  #print(cellgroup)
  imp_val <- min(validation_mast_results[(validation_mast_results$p_val > 0) &
                                         (validation_mast_results$cluster == cellgroup), "p_val"])
  validation_mast_results[(validation_mast_results$p_val == 0) & 
                          (validation_mast_results$cluster == cellgroup), "p_val"] <- imp_val
  test_result <- sumlog(validation_mast_results[validation_mast_results$cluster == cellgroup,"p_val"])
  metap_holder[metap_holder$Var1 == cellgroup, "metap"] <- test_result$p
  }


metap_holder$output <- paste0(metap_holder$metap, " (", metap_holder$Freq, ")")

Mathys_metap <- metap_holder
Cain_metap <- metap_holder
Zhou_metap <- metap_holder

metap_result <- merge(Mathys_metap[,c("Var1", "output")], Cain_metap[,c("Var1", "output")], by = "Var1", all.x = T, all.y = T)
metap_result <- merge(metap_result, Zhou_metap[,c("Var1", "output")], by = "Var1", all.x = T, all.y = T)
names(metap_result) <- c("Subclass", "Mathys_meta_p", "Cain_meta_p", "Zhou_meta_p")

metap_holder <- as.data.frame(table(Result_df_MTGandCgG_lfct2.0$subclass))
metap_holder <- merge(metap_holder, metap_result, by.x = "Var1", by.y = "Subclass", all.x = T, all.y = T)
names(metap_holder) <- c("Subclass", "Number_of_markers", "Mathys_meta_p(marker_n)", "Cain_meta_p(marker_n)", "Zhou_meta_p(marker_n)")

write.csv(metap_holder, "subclass_MTGandCgG_lfct2.0_marker_validation_metap.csv")

#source('https://raw.githubusercontent.com/yaowuliu/ACAT/master/R/ACAT.R')

#for (cellgroup in unique(test$subclass)){
#  print(cellgroup)
#  print(ACAT(test[test$subclass == cellgroup,"p_val"]))
#}


#### confusion matrix/heatmap + QC based on prediction scores ####

# getting data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# seeing columns to graph, as appropriate for dataset

unique(Seu_plot_object$subtype) # for cain
unique(Seu_plot_object$predicted.id)
unique(Seu_plot_object$Subcluster) # for mathys

# qc thresh cut for IT cells

ncol(Seu_plot_object)
table(Seu_plot_object$predicted.id)
Seu_plot_object <- subset(Seu_plot_object, subset = predicted.id == "IT" & prediction.score.max < 0.8, invert = T) 
table(Seu_plot_object$predicted.id)

# get numbers for confusion matrix

confusion_martix_hold <- as.matrix(table(Seu_plot_object$subtype, Seu_plot_object$predicted.id)) # for cain
confusion_martix_hold <- as.matrix(table(Seu_plot_object$Subcluster, Seu_plot_object$predicted.id)) # for mathys
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

#heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

# get order information

subclass_meta_info <- readr::read_csv("~/collabgit/AD_snRNAseq/data/subclass_meta_info.csv")

# reorder confusion matrix

col_order <- intersect(subclass_meta_info$subclass, colnames(confusion_martix_hold))
confusion_martix_hold <- confusion_martix_hold[,col_order]

# get clustered rownames
clustered_holder <- pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = T, cluster_cols = F)
clustered_holder <- rownames(confusion_martix_hold[clustered_holder$tree_row[["order"]],])

# reorder cain
clustered_holder <- clustered_holder[c(16:22,6:14,27,34:35,23:24,28:33,36,25:26,15,38,2:5,37,1)]
confusion_martix_hold <- confusion_martix_hold[clustered_holder,]

# reorder mathys
clustered_holder <- clustered_holder[c(4:12,3,2,1,37:41,33:36,25:27,42:44,28:30,15:16,20:22,31:32,13,23,14,17,24,19,18)]
confusion_martix_hold <- confusion_martix_hold[clustered_holder,]

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### save confusion matrix
confusion_martix_hold_Cain <- confusion_martix_hold
confusion_martix_hold_Cain_filt <- confusion_martix_hold
confusion_martix_hold_Cain_ordered <- confusion_martix_hold
confusion_martix_hold_Mathys <- confusion_martix_hold
confusion_martix_hold_Mathys_filt <- confusion_martix_hold
confusion_martix_hold_Mathys_ordered <- confusion_martix_hold

### plot in ggplot

confusion_matrix_hold <- confusion_martix_hold_Cain_ordered
melted_conf_mtx <- melt(t(confusion_martix_hold_Cain_ordered[nrow(confusion_martix_hold_Cain_ordered):1,]))
melted_conf_mtx <- melt(t(confusion_martix_hold_Mathys_ordered[nrow(confusion_martix_hold_Mathys_ordered):1,]))
melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

ggplot(melted_conf_mtx[melted_conf_mtx$Var2 != "None.NA",], aes(Var1,Var2, fill=value)) + 
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
       height = 210,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "Cain_mapping_conf_mtx_h210.pdf",
       device = "pdf")

#confusion_martix_hold <- confusion_martix_hold_Cain
#confusion_martix_hold <- confusion_martix_hold_Cain_filt
#confusion_martix_hold <- confusion_martix_hold_Mathys
#confusion_martix_hold <- confusion_martix_hold_Mathys_filt


#confusion_martix_hold_man <- confusion_martix_hold[,c(2,8,9,1,3,15,4,12,16,5,6,7,13,10,11,14)]
#pheatmap::pheatmap(confusion_martix_hold_man_pct, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows=T, cluster_cols=F)
#heatmap.2(log(confusion_martix_hold+1))

#### troubleshoot IT ####

Idents(Seu_plot_object) <- "predicted.id"
Seu_plot_object <- subset(Seu_plot_object, idents = c("IT", "Endothelial"))

IT_Troubleshooter <- Seu_plot_object@meta.data

IT_Troubleshooter <- IT_Troubleshooter[,-(5:9)]
IT_Troubleshooter <- IT_Troubleshooter[,-31]
IT_Troubleshooter <- IT_Troubleshooter[,-(29:30)]
IT_Troubleshooter$prediction.score.2nd <- 0 
IT_Troubleshooter$predicted.id.next <- "Unknown"

for (cell in (row.names(IT_Troubleshooter))){
  rowtorank <- IT_Troubleshooter[cell,9:27]
  rankedrows <- as.data.frame(rank(-rowtorank, ties.method = "min"))
  rankedrows <- tibble::rownames_to_column(rankedrows)
  second_type <- rankedrows[rankedrows[,2] == 2,1]
  IT_Troubleshooter[cell,"predicted.id.next"] <-  paste(unlist(second_type), collapse='_')
  IT_Troubleshooter[cell,"prediction.score.2nd"] <- rowtorank[,second_type[1]]
}

IT_Troubleshooter$top2_prediction_diff <- IT_Troubleshooter$prediction.score.max - IT_Troubleshooter$prediction.score.2nd
IT_Troubleshooter$prediction_consistency <- "Not_consistent"
table(IT_Troubleshooter[,c("subtype", "predicted.id")])
IT_Troubleshooter[(IT_Troubleshooter$predicted.id == "Endothelial") &
                  (IT_Troubleshooter$subtype %in% c("Endo.1", "Endo.2", "Endo.4", "Endo.3")),"prediction_consistency"] <- "Consistent"
IT_Troubleshooter[(IT_Troubleshooter$predicted.id == "IT") &
                    (IT_Troubleshooter$subtype %in% c("Exc.Exc.RORB_L5_IT_2", "Exc.Exc.L6.IT.THEMIS", "Exc.Exc.RORB_L5_IT_1")),"prediction_consistency"] <- "Consistent"
IT_Troubleshooter[(IT_Troubleshooter$predicted.id == "IT") &
                    (IT_Troubleshooter$subtype %in% c("Exc.Exc.L5", "Exc.Exc.L3", "Exc.Exc.FEZEF2.L5.ET")),"prediction_consistency"] <- "Maybe"
IT_Troubleshooter[(IT_Troubleshooter$predicted.id == "IT") &
                    (IT_Troubleshooter$subtype %in% c("None.NA")),"prediction_consistency"] <- "NA"

ggplot(IT_Troubleshooter, aes(x = prediction_consistency, y = prediction.score.max)) + 
  geom_violin(scale = "width") +
  geom_jitter(size = 0.2, alpha = 0.15) +
  #geom_boxplot(outlier.size = 2.5, outlier.alpha = 0.1) +
  facet_wrap(~predicted.id * cogdx, scales = "free") +
  scale_color_manual(values=c("red", "blue")) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(IT_Troubleshooter, aes(x = prediction_consistency, y = prediction.score.2nd)) + 
  #geom_violin(scale = "width") +
  #geom_jitter(size = 0.2, alpha = 0.15) +
  geom_boxplot(outlier.size = 2.5, outlier.alpha = 0.1) +
  facet_wrap(~predicted.id, scales = "free") +
  scale_color_manual(values=c("red", "blue")) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(IT_Troubleshooter, aes(x = prediction_consistency, y = top2_prediction_diff)) + 
  #geom_violin(scale = "width") +
  #geom_jitter(size = 0.2, alpha = 0.15) +
  geom_boxplot(outlier.size = 2.5, outlier.alpha = 0.1) +
  facet_wrap(~predicted.id, scales = "free") +
  scale_color_manual(values=c("red", "blue")) +
  theme_classic() +
  theme(legend.position = "none")

### heatmap of second best alignment ###

IT_inconsistent <- IT_Troubleshooter[IT_Troubleshooter$predicted.id == "IT",]
IT_inconsistent <- IT_inconsistent[IT_inconsistent$prediction_consistency == "Not_consistent",]
IT_conf_mtx <- as.matrix(table(IT_inconsistent$subtype, IT_inconsistent$predicted.id.next))
IT_conf_mtx <- IT_conf_mtx[,-18]
IT_conf_mtx <- (IT_conf_mtx/rowSums(IT_conf_mtx))*100
pheatmap::pheatmap(IT_conf_mtx, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

#test <- colRanks(as.matrix(IT_Troubleshooter[,9:27]),)
#IT_Troubleshooter_2ndmax <- apply(IT_Troubleshooter[,9:27], 1, FUN = function(x) which(x == sort(x, decreasing = TRUE)[2]))


#### investigation into prediction scores ####

# getting data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

metadata_for_predscore_breakdown <- Seu_plot_object@meta.data

Projid_info <- readRDS("~/collabgit/AD_snRNAseq/data/ROSmaster.rds") #load ROS metadata
metadata_for_predscore_breakdown <- merge(metadata_for_predscore_breakdown, Projid_info, by = "projid", all.x = T, all.y = F)
remove(Projid_info)

metadata_for_predscore_breakdown %<>% mutate(LOAD = case_when((braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                                              (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                                              TRUE ~ 'OTHER')) 
remove(Seu_plot_object)
metadata_for_predscore_breakdown$over.9 <- metadata_for_predscore_breakdown[, c("prediction.score.max")] > 0.9
metadata_for_predscore_breakdown$over.8 <- metadata_for_predscore_breakdown[, c("prediction.score.max")] > 0.8

#metadata_for_predscore_breakdown$cell_name <- paste0(metadata_for_predscore_breakdown$orig.ident, "_", 
#                                                     metadata_for_predscore_breakdown$nCount_RNA, "_",
#                                                     metadata_for_predscore_breakdown$nFeature_RNA) #for zhou only

#metadata_for_predscore_breakdown$cell_name <- metadata_for_predscore_breakdown$TAG #for mathys only

metadata_for_predscore_breakdown <- metadata_for_predscore_breakdown[,c("cell_name", "predicted.id", "prediction.score.max", "over.9", "over.8", "projid", "LOAD")]

metad_for_pb_Cain <- metadata_for_predscore_breakdown
metad_for_pb_Zhou <- metadata_for_predscore_breakdown
metad_for_pb_Mathys <- metadata_for_predscore_breakdown

metadata_for_predscore_breakdown <- metad_for_pb_Cain
metadata_for_predscore_breakdown <- metad_for_pb_Zhou
metadata_for_predscore_breakdown <- metad_for_pb_Mathys

### table

predbreakdown <- as.data.frame(table(metadata_for_predscore_breakdown[,c("over.9", "predicted.id", "projid")], useNA = "ifany"))
predbreakdown <- merge(predbreakdown, unique(metadata_for_predscore_breakdown[,c("projid", "LOAD")]), all.x = TRUE, by = "projid")
predbreakdown <- spread(predbreakdown, value = "Freq", key = "over.9")
names(predbreakdown)[4:5] <- c("Under_thresh", "Over_thresh")
predbreakdown$Total <- predbreakdown$Under_thresh + predbreakdown$Over_thresh
predbreakdown$Under_thresh_pct <- (predbreakdown$Under_thresh/predbreakdown$Total)*100
predbreakdown$Over_thresh_pct <- (predbreakdown$Over_thresh/predbreakdown$Total)*100
predbreakdown <- predbreakdown[predbreakdown$LOAD %in% c("AD", "C"),] #as desired

### plot

predforplot <- as.data.frame(table(metadata_for_predscore_breakdown[,c("over.9", "predicted.id", "projid")], useNA = "ifany"))
predforplot <- merge(predforplot, unique(metadata_for_predscore_breakdown[,c("projid", "LOAD")]), all.x = TRUE, by = "projid")
predforplot <- predforplot[predforplot$LOAD %in% c("AD", "C"),] #as desired

ggplot() +
  geom_bar(data=predforplot, aes(y = Freq, x = predicted.id, fill = over.9), stat="identity",
           position='fill') +
  theme_classic() + 
  facet_wrap(~LOAD*projid, scales = "fixed") +
  scale_fill_manual(values = c("grey", "red")) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Cell pass percentage") + xlab("Patient ID") + labs(fill = "Prediction score pass/fail breakdown") +
  theme(axis.text.x = element_text(angle = 90))

ggplot() +
  geom_bar(data=predforplot, aes(y = Freq, x = predicted.id, fill = over.9), stat="identity",
           position='fill') +
  theme_classic() + 
  facet_wrap(~LOAD, scales = "fixed", ncol = 1) +
  scale_fill_manual(values = c("grey", "red")) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Cell pass percentage") + xlab("Patient ID") + labs(fill = "Prediction score pass/fail breakdown") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(predbreakdown, aes(x = LOAD, y = Over_thresh_pct)) + 
  #geom_jitter(size = 0.2, alpha = 0.15) +
  geom_boxplot(outlier.size = 2.5, outlier.alpha = 0.1) +
  facet_wrap(~predicted.id, scales = "fixed") +
  scale_color_manual(values=c("red", "blue")) +
  theme_classic() +
  theme(legend.position = "none")

#### "Crappy" Marker proportions lfct 2 vs 2.5 in mathys ####

# look at the % of good markers/bad markers (defined by max marker expression in each group) in the lfct 2.0 and 2.5 marker lists

Seu_mathys_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds")
names(Seu_mathys_obj@meta.data)
Idents(Seu_mathys_obj)

Idents(Seu_mathys_obj) <- "predicted.id"

Mathys_subclass_cpm_avg <- AverageExpression(Seu_mathys_obj, slot = "data")
Mathys_subclass_cpm_avg <- Mathys_subclass_cpm_avg$RNA

Mathys_2.0_cpm_avg <- Mathys_subclass_cpm_avg[intersect(new_MTGnCgG_lfct2_results$gene, rownames(Mathys_subclass_cpm_avg)),]
Mathys_2.5_cpm_avg <- Mathys_subclass_cpm_avg[intersect(Result_df_MTGandCgG_final$gene, rownames(Mathys_subclass_cpm_avg)),]

Mathys_2.0_cpm_avg$Max_Express_Subclass <- colnames(Mathys_2.0_cpm_avg)[apply(Mathys_2.0_cpm_avg,1,function(x) which(x==max(x)))]
Mathys_2.0_cpm_avg <- tibble::rownames_to_column(Mathys_2.0_cpm_avg)
Mathys_2.0_cpm_avg <- merge(Mathys_2.0_cpm_avg, new_MTGnCgG_lfct2_results[,c("gene", "subclass")], by.x = "rowname", by.y = "gene")
Mathys_2.0_cpm_avg$Good_marker <- Mathys_2.0_cpm_avg$Max_Express_Subclass == Mathys_2.0_cpm_avg$subclass
Bad_marker_pct_Mathys2.0 <- as.data.frame(table(Mathys_2.0_cpm_avg[,c("subclass", "Good_marker")]))
Bad_marker_pct_Mathys2.0 <- tidyr::spread(Bad_marker_pct_Mathys2.0, Good_marker, Freq)
Bad_marker_pct_Mathys2.0$Total <- Bad_marker_pct_Mathys2.0[,"FALSE"] + Bad_marker_pct_Mathys2.0[,"TRUE"]
Bad_marker_pct_Mathys2.0$Pct_good_2 <- (Bad_marker_pct_Mathys2.0[,"TRUE"] / Bad_marker_pct_Mathys2.0[,"Total"])*100

Mathys_2.5_cpm_avg$Max_Express_Subclass <- colnames(Mathys_2.5_cpm_avg)[apply(Mathys_2.5_cpm_avg,1,function(x) which(x==max(x)))]
Mathys_2.5_cpm_avg <- tibble::rownames_to_column(Mathys_2.5_cpm_avg)
Mathys_2.5_cpm_avg <- merge(Mathys_2.5_cpm_avg, Result_df_MTGandCgG_final[,c("gene", "subclass")], by.x = "rowname", by.y = "gene")
Mathys_2.5_cpm_avg$Good_marker <- Mathys_2.5_cpm_avg$Max_Express_Subclass == Mathys_2.5_cpm_avg$subclass
Bad_marker_pct_Mathys2.5 <- as.data.frame(table(Mathys_2.5_cpm_avg[,c("subclass", "Good_marker")]))
Bad_marker_pct_Mathys2.5 <- tidyr::spread(Bad_marker_pct_Mathys2.5, Good_marker, Freq)
Bad_marker_pct_Mathys2.5$Total <- Bad_marker_pct_Mathys2.5[,"FALSE"] + Bad_marker_pct_Mathys2.5[,"TRUE"]
Bad_marker_pct_Mathys2.5$Pct_good_2.5 <- (Bad_marker_pct_Mathys2.5[,"TRUE"] / Bad_marker_pct_Mathys2.5[,"Total"])*100

Bad_marker_pct_combined <- merge(Bad_marker_pct_Mathys2.0, Bad_marker_pct_Mathys2.5, by = "subclass")
Bad_marker_pct_combined$Pct2.5_minus_Pct2.0 <- Bad_marker_pct_combined$Pct_good_2.5 - Bad_marker_pct_combined$Pct_good_2
colnames(Bad_marker_pct_combined) <- c("Subclass", "n_Bad_markers_lfct2", "n_Good_markers_lfct2", "n_Total_markers_lfct2", "pct_Good_markers_lfct2",
                                       "n_Bad_markers_lfct2.5", "n_Good_markers_lfct2.5", "n_Total_markers_lfct2.5", "pct_Good_markers_lfct2.5", "Pct2.5_minus_Pct2.0")

write.csv(Bad_marker_pct_combined, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Validation/Mathys_lfct_good_marker_pct.csv") #save/export results

#### QC based on sex-specific genes (Old/unused) ####

hodge_cpm <- new_Seu_AIBS_obj_for_test@assays$RNA@data[c("XIST", "KDM5D", "RPS4Y1"),]
hodge_cpm <- as.matrix(hodge_cpm)
hodge_cpm <- t(hodge_cpm)
hodge_cpm <- as.data.frame(hodge_cpm)
hodge_cpm <- tibble::rownames_to_column(hodge_cpm)

hodge_meta <- new_Seu_AIBS_obj_for_test@meta.data
hodge_meta <- tibble::rownames_to_column(hodge_meta)

hodge_qcbysex <- merge(hodge_cpm, hodge_meta, by = "rowname")

hodge_qcbysex <- tidyr::gather(hodge_qcbysex, "Gene", "Expression_Cpm", c(2:4))

ggplot(hodge_qcbysex, aes(x = donor_sex_label, y = Expression_Cpm, color = donor_sex_label)) + 
  geom_jitter(size = 0.5) +
  facet_wrap(~Gene, scales = "free") +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() +
  theme(legend.position = "none")

for (gene in unique(hodge_qcbysex$Gene)) {
  
  for (sex in unique(hodge_qcbysex$donor_sex_label)){
    
    quant_holder <- quantile(hodge_qcbysex[((hodge_qcbysex$donor_sex_label == sex) &
                                              (hodge_qcbysex$Gene == gene)), "Expression_Cpm"], 
                             probs = c(0.9, 0.95))
    print(c(gene, sex, quant_holder))
    
  }
  
}

ggplot(hodge_qcbysex, aes(x = donor_sex_label, y = Expression_Cpm, color = donor_sex_label)) + 
  geom_jitter(size = 0.5) +
  facet_wrap(~Gene, scales = "free") +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() 

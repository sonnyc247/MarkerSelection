#### packages ####

library(Seurat)
library(dplyr)
library(MAST)
library(metap)
library(gplots)
library(RColorBrewer)

#### DE expression testing ####

### prep/set-up

Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_test_object <- subset(Seu_test_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

names(Seu_test_object@meta.data)
Idents(Seu_test_object) <- "predicted.id"
table(Idents(Seu_test_object))

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

#### meta-p value combination ####

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


#### confusion matrix/heatmap ####

# getting data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_test_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# seeing columns to graph, as appropriate for dataset

unique(Seu_plot_object$subtype)
unique(Seu_plot_object$predicted.id)

# get numbers for confusion matrix

confusion_martix_hold <- as.matrix(table(Seu_plot_object$subtype, Seu_plot_object$predicted.id))
confusion_martix_hold <- as.matrix(confusion_martix_hold)

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)
heatmap(log(confusion_martix_hold+1), col=brewer.pal(9 ,"Blues"), keep.dendro = F)

pheatmap::pheatmap(log(confusion_martix_hold+1), treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"))

confusion_martix_hold_man <- confusion_martix_hold[,c(1,13,2,10,9,4,12,16,8,5,6,7,11,15,14,3)]
pheatmap::pheatmap(log(confusion_martix_hold_man+1), treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows=T, cluster_cols=F)

heatmap.2(log(confusion_martix_hold+1))


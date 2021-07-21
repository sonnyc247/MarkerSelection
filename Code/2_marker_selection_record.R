#### get excit/inhib markers MTG ####

### filter data
Idents(new_Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd to make it easy for us to remove neurons from MTG/CgG later

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj, idents = "CgG_Neuronal", invert = TRUE) #remove "CgG_Neuronal" cells

table(new_Seu_AIBS_obj_for_test$class_label) #double check what classes we have

Idents(new_Seu_AIBS_obj_for_test) <- "class_label" #assign "class_label" as the key grouping variable now

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_MTG_class <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers

new_AIBS_markers_roc_MTG_class <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates

dup_list <- unique(new_AIBS_markers_mast_MTG_class[duplicated(new_AIBS_markers_mast_MTG_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTG_class <- new_AIBS_markers_mast_MTG_class[!(new_AIBS_markers_mast_MTG_class$gene %in% dup_list),] #remove duplicated marker genes

dup_list <- unique(new_AIBS_markers_roc_MTG_class[duplicated(new_AIBS_markers_roc_MTG_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTG_class <- new_AIBS_markers_roc_MTG_class[!(new_AIBS_markers_roc_MTG_class$gene %in% dup_list),] #remove duplicated marker genes

remove(dup_list) #clean temporary object

### finalize df

length(intersect(new_AIBS_markers_mast_MTG_class$gene, new_AIBS_markers_roc_MTG_class$gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_MTG_class[,c(8,7,5,6,2,1,3)], new_AIBS_markers_mast_MTG_class[,c(1,5,7,6,3,4,2)], by = "gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
unique_genes <- setdiff(new_AIBS_markers_mast_MTG_class$gene, new_AIBS_markers_roc_MTG_class$gene) #genes unique to mast_MTG
Result_df <- Result_df[,1:9] #remove redundant columns

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:11)] <- c("ensembl_id", "class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

new_MTG_results_class <- Result_df #store the results in R
write.csv(new_MTG_results_class, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTG_results_class.csv") #save/export results



#### get excit/inhib markers CgG ####

### filter data
Idents(new_Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd to make it easy for us to remove neurons from MTG/CgG later

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj, idents = "MTG_Neuronal", invert = TRUE) #remove "CgG_Neuronal" cells

table(new_Seu_AIBS_obj_for_test$class_label) #double check what classes we have

Idents(new_Seu_AIBS_obj_for_test) <- "class_label" #assign "class_label" as the key grouping variable now

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_CgG_class <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers

new_AIBS_markers_roc_CgG_class <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates

dup_list <- unique(new_AIBS_markers_mast_CgG_class[duplicated(new_AIBS_markers_mast_CgG_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_CgG_class <- new_AIBS_markers_mast_CgG_class[!(new_AIBS_markers_mast_CgG_class$gene %in% dup_list),] #remove duplicated marker genes

dup_list <- unique(new_AIBS_markers_roc_CgG_class[duplicated(new_AIBS_markers_roc_CgG_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_CgG_class <- new_AIBS_markers_roc_CgG_class[!(new_AIBS_markers_roc_CgG_class$gene %in% dup_list),] #remove duplicated marker genes

remove(dup_list) #clean temporary object

### finalize df

length(intersect(new_AIBS_markers_mast_CgG_class$gene, new_AIBS_markers_roc_CgG_class$gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_CgG_class[,c(8,7,5,6,2,1,3)], new_AIBS_markers_mast_CgG_class[,c(1,5,7,6,3,4,2)], by = "gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
unique_genes <- setdiff(new_AIBS_markers_mast_CgG_class$gene, new_AIBS_markers_roc_CgG_class$gene) #genes unique to mast_CgG
Result_df <- Result_df[,1:9] #remove redundant columns

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:11)] <- c("ensembl_id", "class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

new_CgG_results_class <- Result_df #store the results in R
write.csv(new_CgG_results_class, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_CgG_results_class.csv") #save/export results

#### get Neu/Nonneu markers CgG ####

### filter data, modify grouping
Idents(new_Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd to make it easy for us to remove neurons from MTG/CgG later

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj, idents = "MTG_Neuronal", invert = TRUE) #remove "CgG_Neuronal" cells

table(new_Seu_AIBS_obj_for_test$class_label) #double check what classes we have

new_Seu_AIBS_obj_for_test$class_label[new_Seu_AIBS_obj_for_test$class_label=="Glutamatergic"] <- "Neuronal" #collapsing the excit. group
new_Seu_AIBS_obj_for_test$class_label[new_Seu_AIBS_obj_for_test$class_label=="GABAergic"] <- "Neuronal" #collapsing the inhib. group
levels(new_Seu_AIBS_obj_for_test) <- c("Neuronal", "Non-neuronal")

table(new_Seu_AIBS_obj_for_test$class_label, useNA = "always") #double check what classes we have now

Idents(new_Seu_AIBS_obj_for_test) <- "class_label" #assign "class_label" as the key grouping variable now

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_CgG_NeuNonN <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers

new_AIBS_markers_roc_CgG_NeuNonN <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates

dup_list <- unique(new_AIBS_markers_mast_CgG_NeuNonN[duplicated(new_AIBS_markers_mast_CgG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_CgG_NeuNonN <- new_AIBS_markers_mast_CgG_NeuNonN[!(new_AIBS_markers_mast_CgG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

dup_list <- unique(new_AIBS_markers_roc_CgG_NeuNonN[duplicated(new_AIBS_markers_roc_CgG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_CgG_NeuNonN <- new_AIBS_markers_roc_CgG_NeuNonN[!(new_AIBS_markers_roc_CgG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

remove(dup_list) #clean temporary object

### finalize df

length(intersect(new_AIBS_markers_mast_CgG_NeuNonN$gene, new_AIBS_markers_roc_CgG_NeuNonN$gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_CgG_NeuNonN[,c(8,7,5,6,2,1,3)], new_AIBS_markers_mast_CgG_NeuNonN[,c(1,5,7,6,3,4,2)], by = "gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
unique_genes <- setdiff(new_AIBS_markers_mast_CgG_NeuNonN$gene, new_AIBS_markers_roc_CgG_NeuNonN$gene) #genes unique to mast_CgG
Result_df <- Result_df[,1:9] #remove redundant columns

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:11)] <- c("ensembl_id", "class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

new_CgG_results_NeuNonN <- Result_df #store the results in R
write.csv(new_CgG_results_NeuNonN, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_CgG_results_NeuNonN.csv") #save/export results

#### get Neu/Nonneu markers MTG ####

### filter data, modify grouping
Idents(new_Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd to make it easy for us to remove neurons from MTG/CgG later

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj, idents = "CgG_Neuronal", invert = TRUE) #remove "CgG_Neuronal" cells

table(new_Seu_AIBS_obj_for_test$class_label) #double check what classes we have

new_Seu_AIBS_obj_for_test$class_label[new_Seu_AIBS_obj_for_test$class_label=="Glutamatergic"] <- "Neuronal" #collapsing the excit. group
new_Seu_AIBS_obj_for_test$class_label[new_Seu_AIBS_obj_for_test$class_label=="GABAergic"] <- "Neuronal" #collapsing the inhib. group
levels(new_Seu_AIBS_obj_for_test) <- c("Neuronal", "Non-neuronal")

table(new_Seu_AIBS_obj_for_test$class_label, useNA = "always") #double check what classes we have now

Idents(new_Seu_AIBS_obj_for_test) <- "class_label" #assign "class_label" as the key grouping variable now

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_MTG_NeuNonN <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers

new_AIBS_markers_roc_MTG_NeuNonN <- FindAllMarkers(new_Seu_AIBS_obj_for_test, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates

dup_list <- unique(new_AIBS_markers_mast_MTG_NeuNonN[duplicated(new_AIBS_markers_mast_MTG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTG_NeuNonN <- new_AIBS_markers_mast_MTG_NeuNonN[!(new_AIBS_markers_mast_MTG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

dup_list <- unique(new_AIBS_markers_roc_MTG_NeuNonN[duplicated(new_AIBS_markers_roc_MTG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTG_NeuNonN <- new_AIBS_markers_roc_MTG_NeuNonN[!(new_AIBS_markers_roc_MTG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

remove(dup_list) #clean temporary object

### finalize df

length(intersect(new_AIBS_markers_mast_MTG_NeuNonN$gene, new_AIBS_markers_roc_MTG_NeuNonN$gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_MTG_NeuNonN[,c(8,7,5,6,2,1,3)], new_AIBS_markers_mast_MTG_NeuNonN[,c(1,5,7,6,3,4,2)], by = "gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
unique_genes <- setdiff(new_AIBS_markers_mast_MTG_NeuNonN$gene, new_AIBS_markers_roc_MTG_NeuNonN$gene) #genes unique to mast_MTG
Result_df <- Result_df[,1:9] #remove redundant columns

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:11)] <- c("ensembl_id", "class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

new_MTG_results_NeuNonN <- Result_df #store the results in R
write.csv(new_MTG_results_NeuNonN, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTG_results_NeuNonN.csv") #save/export results


# comparing cell type overlaps

overlap_summary <- data.frame(Celltype=character(),
                              CgG_total_markers=double(), 
                              Overlapping_markers=double(),
                              MTG_total_markers=double(),
                              stringsAsFactors=FALSE) 

for (i in 1:length(unique(new_AIBS_markers_mast_MTG$cluster))) {
  
  cellgroup <- as.character(unique(new_AIBS_markers_mast_MTG$cluster)[i])
  temp_df <- new_AIBS_markers_mast_CgG[new_AIBS_markers_mast_CgG$cluster == cellgroup,] 
  temp_df2 <- new_AIBS_markers_mast_MTG[new_AIBS_markers_mast_MTG$cluster == cellgroup,]
  
  overlap_summary[i,"Celltype"] <- cellgroup
  overlap_summary[i,"CgG_total_markers"] <- nrow(temp_df)
  overlap_summary[i,"Overlapping_markers"] <- length(intersect(temp_df$gene, temp_df2$gene))
  overlap_summary[i,"MTG_total_markers"] <- nrow(temp_df2)
  
}

write.csv(overlap_summary, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/mast_overlap_summary.csv") #save/export results

#### get markers from CgG and MTG (but no other brain regions for non-neurons) ####

unique(new_Seu_AIBS_obj$NeuN_Region)
unique(new_Seu_AIBS_obj$outlier_call)

Idents(new_Seu_AIBS_obj) <- "subclass_label"
unique(Idents(new_Seu_AIBS_obj))
table(Idents(new_Seu_AIBS_obj))

new_AIBS_markers_mast_MTGandCgG <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_MTGandCgG <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### for mast
length(unique(new_AIBS_markers_mast_MTGandCgG$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_mast_MTGandCgG[duplicated(new_AIBS_markers_mast_MTGandCgG$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTGandCgG <- new_AIBS_markers_mast_MTGandCgG[!(new_AIBS_markers_mast_MTGandCgG$gene %in% dup_list),] #remove duplicated marker genes

### for roc
length(unique(new_AIBS_markers_roc_MTGandCgG$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_roc_MTGandCgG[duplicated(new_AIBS_markers_roc_MTGandCgG$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTGandCgG <- new_AIBS_markers_roc_MTGandCgG[!(new_AIBS_markers_roc_MTGandCgG$gene %in% dup_list),] #remove duplicated marker genes

length(intersect(new_AIBS_markers_mast_MTGandCgG$gene, new_AIBS_markers_roc_MTGandCgG$gene)) #see intersect of marker genes

### combine
new_AIBS_markers_mast_MTGandCgG$group_gene <- paste0(new_AIBS_markers_mast_MTGandCgG$cluster, "_", new_AIBS_markers_mast_MTGandCgG$gene)
new_AIBS_markers_roc_MTGandCgG$group_gene <- paste0(new_AIBS_markers_roc_MTGandCgG$cluster, "_", new_AIBS_markers_roc_MTGandCgG$gene)

Result_df_MTGandCgG <- merge(new_AIBS_markers_roc_MTGandCgG, new_AIBS_markers_mast_MTGandCgG, by = "group_gene", all.x = T, all.y = T) #combine the marker df

Result_df_MTGandCgG <- Result_df_MTGandCgG %>% mutate(gene = coalesce(gene.x, gene.y))
Result_df_MTGandCgG <- Result_df_MTGandCgG %>% mutate(cluster = coalesce(cluster.x, cluster.y))
Result_df_MTGandCgG <- Result_df_MTGandCgG %>% mutate(avg_logFC = coalesce(avg_logFC.x, avg_logFC.y))
Result_df_MTGandCgG <- Result_df_MTGandCgG %>% mutate(pct.1 = coalesce(pct.1.x, pct.1.y))
Result_df_MTGandCgG <- Result_df_MTGandCgG %>% mutate(pct.2 = coalesce(pct.2.x, pct.2.y))

### to add entrez and ensembl IDs to the output/result df

Result_df_MTGandCgG <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df_MTGandCgG, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno

Result_df_MTGandCgG$ensembl_gene_id[is.na(Result_df_MTGandCgG$ensembl_gene_id)] <- "NA"
Result_df_MTGandCgG$has_ensembl <- Result_df_MTGandCgG$ensembl_gene_id != "NA"

table(Result_df_MTGandCgG[, c("cluster", "has_ensembl")])

### clean and/or export

Result_df_MTGandCgG_final <- Result_df_MTGandCgG[, c(1:3,20,22:23,21,5,7,13,17,4,24)]

### Micaela's list overlap
humanMarkersCommon <- readRDS("~/git/MarkerSelection/Data/Inputs/humanMarkersCommon.rds")
humanMarkersCommon <- unique(unname(unlist(humanMarkersCommon)))
Result_df_MTGandCgG_final$Mic_Overlap <- Result_df_MTGandCgG_final$gene %in% humanMarkersCommon
###

colnames(Result_df_MTGandCgG_final)[c(4,8:11,14)] <- c("subclass", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj", "hMC_overlap") #rename some columns for clarity

write.csv(Result_df_MTGandCgG_final[,c(1:11,14)], "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_results.csv") #save/export results

write.csv(Result_df_MTGandCgG_final, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/") #save/export results

### Check for overlap between 2.5 cutoff markers and common gene from Micaela

humanMarkersCommon <- readRDS("~/git/MarkerSelection/Data/Inputs/humanMarkersCommon.rds")
humanMarkersCommon <- unique(unname(unlist(humanMarkersCommon)))

Result_df_MTGandCgG_final$Mic_Overlap <- Result_df_MTGandCgG_final$gene %in% humanMarkersCommon
Marker_overlap_holder <- data.frame(unclass(table(Result_df_MTGandCgG_final$subclass, Result_df_MTGandCgG_final$Mic_Overlap)))

names(Marker_overlap_holder) <- c("Not_in_list", "In_list")
Marker_overlap_holder$Total_n_markers <- Marker_overlap_holder$Not_in_list + Marker_overlap_holder$In_list

Marker_overlap_holder <- Marker_overlap_holder[,c(3,2,1)]

write.csv(Marker_overlap_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Validation/MTGnGgG_Lfct2.5_common_overlap.csv") #save/export results

#### get markers from CgG and MTG again - try lower lfct thresholds (2.0 used in the end) ####

unique(new_Seu_AIBS_obj$NeuN_Region)

Idents(new_Seu_AIBS_obj) <- "subclass_label"
unique(Idents(new_Seu_AIBS_obj))
table(Idents(new_Seu_AIBS_obj))

new_AIBS_markers_mast_MTGandCgG_lfct2.0 <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_MTGandCgG_lfct2.0 <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

new_AIBS_markers_mast_MTGandCgG_lfct1.5 <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 1.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_MTGandCgG_lfct1.5 <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 1.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc


### for mast
length(unique(new_AIBS_markers_mast_MTGandCgG_lfct1.5$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_mast_MTGandCgG_lfct1.5[duplicated(new_AIBS_markers_mast_MTGandCgG_lfct1.5$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTGandCgG_lfct1.5 <- new_AIBS_markers_mast_MTGandCgG_lfct1.5[!(new_AIBS_markers_mast_MTGandCgG_lfct1.5$gene %in% dup_list),] #remove duplicated marker genes

### for roc
length(unique(new_AIBS_markers_roc_MTGandCgG_lfct1.5$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_roc_MTGandCgG_lfct1.5[duplicated(new_AIBS_markers_roc_MTGandCgG_lfct1.5$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTGandCgG_lfct1.5 <- new_AIBS_markers_roc_MTGandCgG_lfct1.5[!(new_AIBS_markers_roc_MTGandCgG_lfct1.5$gene %in% dup_list),] #remove duplicated marker genes

length(intersect(new_AIBS_markers_mast_MTGandCgG_lfct1.5$gene, new_AIBS_markers_roc_MTGandCgG_lfct1.5$gene)) #see intersect of marker genes

### combine
new_AIBS_markers_mast_MTGandCgG_lfct1.5$group_gene <- paste0(new_AIBS_markers_mast_MTGandCgG_lfct1.5$cluster, "_", new_AIBS_markers_mast_MTGandCgG_lfct1.5$gene)
new_AIBS_markers_roc_MTGandCgG_lfct1.5$group_gene <- paste0(new_AIBS_markers_roc_MTGandCgG_lfct1.5$cluster, "_", new_AIBS_markers_roc_MTGandCgG_lfct1.5$gene)

Result_df_MTGandCgG_lfct1.5 <- merge(new_AIBS_markers_roc_MTGandCgG_lfct1.5, new_AIBS_markers_mast_MTGandCgG_lfct1.5, by = "group_gene", all.x = T, all.y = T) #combine the marker df

Result_df_MTGandCgG_lfct1.5 <- Result_df_MTGandCgG_lfct1.5 %>% mutate(gene = coalesce(gene.x, gene.y))
Result_df_MTGandCgG_lfct1.5 <- Result_df_MTGandCgG_lfct1.5 %>% mutate(cluster = coalesce(cluster.x, cluster.y))
Result_df_MTGandCgG_lfct1.5 <- Result_df_MTGandCgG_lfct1.5 %>% mutate(avg_logFC = coalesce(avg_logFC.x, avg_logFC.y))
Result_df_MTGandCgG_lfct1.5 <- Result_df_MTGandCgG_lfct1.5 %>% mutate(pct.1 = coalesce(pct.1.x, pct.1.y))
Result_df_MTGandCgG_lfct1.5 <- Result_df_MTGandCgG_lfct1.5 %>% mutate(pct.2 = coalesce(pct.2.x, pct.2.y))

### to add entrez and ensembl IDs to the output/result df

Result_df_MTGandCgG_lfct1.5 <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df_MTGandCgG_lfct1.5, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno

Result_df_MTGandCgG_lfct1.5$ensembl_gene_id[is.na(Result_df_MTGandCgG_lfct1.5$ensembl_gene_id)] <- "NA"
Result_df_MTGandCgG_lfct1.5$has_ensembl <- Result_df_MTGandCgG_lfct1.5$ensembl_gene_id != "NA"

table(Result_df_MTGandCgG_lfct2.0[, c("cluster", "has_ensembl")])
table(Result_df_MTGandCgG_lfct1.5[, c("cluster", "has_ensembl")])

### clean and/or export

Result_df_MTGandCgG_lfct2.0 <- Result_df_MTGandCgG_lfct2.0[, c(1:3,20,22:23,21,5,7,13,17,4,24)]
colnames(Result_df_MTGandCgG_lfct2.0)[c(4,8:11)] <- c("subclass", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

write.csv(Result_df_MTGandCgG_lfct2.0, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTGnCgG_lfct2_results.csv") #save/export results
write.csv(Result_df_MTGandCgG_lfct2.0[,1:11], "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTGnCgG_lfct2_results.csv") #save/export results

### compare

new_AIBS_markers_mast_MTGandCgG_screen <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 0.1, min.pct = 0.1, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers


new_comparator <- Result_df_MTGandCgG_lfct2.0
new_comparator <- new_comparator[new_comparator$has_ensembl == T,]
table(new_comparator$cluster)
intersect(new_CgG_results$comparator, new_comparator$group_gene) %>% length()

subclass_to_focus <- "PVALB"
intersect(new_CgG_results[new_CgG_results$subclass == subclass_to_focus, "gene"], new_comparator[new_comparator$cluster == subclass_to_focus, "gene"])

intersect(Result_df_MTGandCgG[Result_df_MTGandCgG$cluster == subclass_to_focus, "gene"], new_comparator[new_comparator$cluster == subclass_to_focus, "gene"])
union(Result_df_MTGandCgG[Result_df_MTGandCgG$cluster == subclass_to_focus, "gene"], new_comparator[new_comparator$cluster == subclass_to_focus, "gene"])

#### get excit/inhib markers from MTG + CgG (lfct 2.0) ####

### prep data
new_Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")

table(new_Seu_AIBS_obj$NeuN_Region) #double check that we have all the cells we want
table(new_Seu_AIBS_obj$class_label) #double check what classes we have

Idents(new_Seu_AIBS_obj) <- "class_label" #assign "class_label" as the key grouping variable now
table(Idents(new_Seu_AIBS_obj)) #see # of cells in each class

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_MTGandCgG_lfct2.0_class <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_MTGandCgG_lfct2.0_class <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc


### for mast
length(unique(new_AIBS_markers_mast_MTGandCgG_lfct2.0_class$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_mast_MTGandCgG_lfct2.0_class[duplicated(new_AIBS_markers_mast_MTGandCgG_lfct2.0_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTGandCgG_lfct2.0_class <- new_AIBS_markers_mast_MTGandCgG_lfct2.0_class[!(new_AIBS_markers_mast_MTGandCgG_lfct2.0_class$gene %in% dup_list),] #remove duplicated marker genes

### for roc
length(unique(new_AIBS_markers_roc_MTGandCgG_lfct2.0_class$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_roc_MTGandCgG_lfct2.0_class[duplicated(new_AIBS_markers_roc_MTGandCgG_lfct2.0_class$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTGandCgG_lfct2.0_class <- new_AIBS_markers_roc_MTGandCgG_lfct2.0_class[!(new_AIBS_markers_roc_MTGandCgG_lfct2.0_class$gene %in% dup_list),] #remove duplicated marker genes

length(intersect(new_AIBS_markers_mast_MTGandCgG_lfct2.0_class$gene, new_AIBS_markers_roc_MTGandCgG_lfct2.0_class$gene)) #see intersect of marker genes

### combine
Result_df_MTGandCgG_lfct2.0_class <- merge(new_AIBS_markers_roc_MTGandCgG_lfct2.0_class, new_AIBS_markers_mast_MTGandCgG_lfct2.0_class, by = "gene", all.x = T, all.y = T) #combine the marker df

identical(Result_df_MTGandCgG_lfct2.0_class[,c("avg_logFC.x", "pct.1.x", "pct.2.x", "cluster.x")],
          Result_df_MTGandCgG_lfct2.0_class[,c("avg_logFC.y", "pct.1.y", "pct.2.y", "cluster.y")])

identical(Result_df_MTGandCgG_lfct2.0_class$avg_logFC.x, Result_df_MTGandCgG_lfct2.0_class$avg_logFC.y)
identical(Result_df_MTGandCgG_lfct2.0_class$pct.1.x, Result_df_MTGandCgG_lfct2.0_class$pct.1.y)
identical(Result_df_MTGandCgG_lfct2.0_class$pct.2.x, Result_df_MTGandCgG_lfct2.0_class$pct.2.y)
identical(Result_df_MTGandCgG_lfct2.0_class$cluster.x, Result_df_MTGandCgG_lfct2.0_class$cluster.y)

### to add entrez and ensembl IDs to the output/result df

Result_df_MTGandCgG_lfct2.0_class <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df_MTGandCgG_lfct2.0_class, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno

Result_df_MTGandCgG_lfct2.0_class$ensembl_gene_id[is.na(Result_df_MTGandCgG_lfct2.0_class$ensembl_gene_id)] <- "NA"
Result_df_MTGandCgG_lfct2.0_class$has_ensembl <- Result_df_MTGandCgG_lfct2.0_class$ensembl_gene_id != "NA"

table(Result_df_MTGandCgG_lfct2.0_class[, c("cluster.x", "has_ensembl")])

### clean and/or export

Result_df_MTGandCgG_lfct2.0_class <- Result_df_MTGandCgG_lfct2.0_class[, c(1:3,10,8,9,7,4,6,11,15,17)]
colnames(Result_df_MTGandCgG_lfct2.0_class)[c(4:11)] <- c("class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

write.csv(Result_df_MTGandCgG_lfct2.0_class[,1:11], "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTGnCgG_lfct2_results_class.csv") #save/export results

#### get neu/nonneu markers from MTG + CgG (lfct 2.0) ####

### prep data
new_Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")

table(new_Seu_AIBS_obj$NeuN_Region) #double check that we have all the cells we want
table(new_Seu_AIBS_obj$NeuN) #double check what groups we have

Idents(new_Seu_AIBS_obj) <- "NeuN" #assign "NeuN" as the key grouping variable now
table(Idents(new_Seu_AIBS_obj)) #see # of cells in each group

### screen data

library("MAST") #as needed

### find markers

new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.0, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc


### for mast
length(unique(new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN[duplicated(new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN <- new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN[!(new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN$gene %in% dup_list),] #remove duplicated marker genes

### for roc
length(unique(new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN[duplicated(new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN <- new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN[!(new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN$gene %in% dup_list),] #remove duplicated marker genes

length(intersect(new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN$gene, new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN$gene)) #see intersect of marker genes

### combine
Result_df_MTGandCgG_lfct2.0_NeuN <- merge(new_AIBS_markers_roc_MTGandCgG_lfct2.0_NeuN, new_AIBS_markers_mast_MTGandCgG_lfct2.0_NeuN, by = "gene", all.x = T, all.y = T) #combine the marker df

identical(Result_df_MTGandCgG_lfct2.0_NeuN[,c("avg_logFC.x", "pct.1.x", "pct.2.x", "cluster.x")],
          Result_df_MTGandCgG_lfct2.0_NeuN[,c("avg_logFC.y", "pct.1.y", "pct.2.y", "cluster.y")])

identical(Result_df_MTGandCgG_lfct2.0_NeuN$avg_logFC.x, Result_df_MTGandCgG_lfct2.0_NeuN$avg_logFC.y)
identical(Result_df_MTGandCgG_lfct2.0_NeuN$pct.1.x, Result_df_MTGandCgG_lfct2.0_NeuN$pct.1.y)
identical(Result_df_MTGandCgG_lfct2.0_NeuN$pct.2.x, Result_df_MTGandCgG_lfct2.0_NeuN$pct.2.y)
identical(Result_df_MTGandCgG_lfct2.0_NeuN$cluster.x, Result_df_MTGandCgG_lfct2.0_NeuN$cluster.y)

### to add entrez and ensembl IDs to the output/result df

Result_df_MTGandCgG_lfct2.0_NeuN <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df_MTGandCgG_lfct2.0_NeuN, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno

Result_df_MTGandCgG_lfct2.0_NeuN$ensembl_gene_id[is.na(Result_df_MTGandCgG_lfct2.0_NeuN$ensembl_gene_id)] <- "NA"
Result_df_MTGandCgG_lfct2.0_NeuN$has_ensembl <- Result_df_MTGandCgG_lfct2.0_NeuN$ensembl_gene_id != "NA"

table(Result_df_MTGandCgG_lfct2.0_NeuN[, c("cluster.x", "has_ensembl")])

### clean and/or export

Result_df_MTGandCgG_lfct2.0_NeuN <- Result_df_MTGandCgG_lfct2.0_NeuN[, c(1:3,10,8,9,7,4,6,11,15,17)]
colnames(Result_df_MTGandCgG_lfct2.0_NeuN)[c(4:11)] <- c("nueron_or_not", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

write.csv(Result_df_MTGandCgG_lfct2.0_NeuN[,1:11], "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTGnCgG_lfct2_results_NeuN.csv") #save/export results

#### get markers from all of hodge ####

metadata <- new_Seu_AIBS_obj@meta.data
Idents(new_Seu_AIBS_obj) <- "subclass_label"

table(new_Seu_AIBS_obj$subclass_label) #double check what classes we have

Idents(new_Seu_AIBS_obj) #check we have subclass as active identity

### find markers

new_AIBS_markers_mast_ALL <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers

new_AIBS_markers_roc_ALL <- FindAllMarkers(new_Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates

dup_list <- unique(new_AIBS_markers_mast_MTG_NeuNonN[duplicated(new_AIBS_markers_mast_MTG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_MTG_NeuNonN <- new_AIBS_markers_mast_MTG_NeuNonN[!(new_AIBS_markers_mast_MTG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

dup_list <- unique(new_AIBS_markers_roc_MTG_NeuNonN[duplicated(new_AIBS_markers_roc_MTG_NeuNonN$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_MTG_NeuNonN <- new_AIBS_markers_roc_MTG_NeuNonN[!(new_AIBS_markers_roc_MTG_NeuNonN$gene %in% dup_list),] #remove duplicated marker genes

remove(dup_list) #clean temporary object

### finalize df

length(intersect(new_AIBS_markers_mast_MTG_NeuNonN$gene, new_AIBS_markers_roc_MTG_NeuNonN$gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_MTG_NeuNonN[,c(8,7,5,6,2,1,3)], new_AIBS_markers_mast_MTG_NeuNonN[,c(1,5,7,6,3,4,2)], by = "gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
unique_genes <- setdiff(new_AIBS_markers_mast_MTG_NeuNonN$gene, new_AIBS_markers_roc_MTG_NeuNonN$gene) #genes unique to mast_MTG
Result_df <- Result_df[,1:9] #remove redundant columns

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:11)] <- c("ensembl_id", "class", "pct.1", "pct.2", "avg_logFC", "roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

new_MTG_results_NeuNonN <- Result_df #store the results in R
write.csv(new_MTG_results_NeuNonN, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_MTG_results_NeuNonN.csv") #save/export results

#### get markers from all of hodge with new IT groups (and setting new IT Groups) ####

### initial setup

library(Seurat)
library(MAST)
library(dplyr)

Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds")
table(Seu_AIBS_obj$outlier_call, exclude = "ifany") #check for outliers
table(Seu_AIBS_obj$NeuN_Region, exclude = "ifany") #check our data composition

### make new IT groups

Seu_AIBS_obj$subclass_label_expanded <- Seu_AIBS_obj$subclass_label_expanded #make new variable, starting from subclasses

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 062",
                                                                                     "Neuron 063",
                                                                                     "Neuron 064",
                                                                                     "Neuron 065",
                                                                                     "Neuron 066",
                                                                                     "Neuron 067",
                                                                                     "Neuron 068",
                                                                                     "Neuron 069",
                                                                                     "Neuron 070",
                                                                                     "Neuron 071",
                                                                                     "Neuron 072")] <- "L2/3 IT" #set L2/3 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 073",
                                                                                     "Neuron 074",
                                                                                     "Neuron 075",
                                                                                     "Neuron 076")] <- "L6 IT" #set L6 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 079",
                                                                                     "Neuron 080",
                                                                                     "Neuron 081",
                                                                                     "Neuron 082",
                                                                                     "Neuron 083",
                                                                                     "Neuron 084",
                                                                                     "Neuron 085",
                                                                                     "Neuron 086",
                                                                                     "Neuron 087")] <- "L5 IT" #set L5 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 077",
                                                                                     "Neuron 078")] <- "Other IT" #set Leftover IT

Seu_AIBS_obj$subclass_label_expanded_L35IT[Seu_AIBS_obj$subclass_label_expanded %in% c("L5 IT",
                                                                                       "Other IT")] <- "L3/5 IT" #set Leftover IT

table(Seu_AIBS_obj$subclass_label_expanded_L35IT)
Idents(Seu_AIBS_obj) <- "subclass_label_expanded_L35IT"
table(Idents(Seu_AIBS_obj)) #double check what subclasses we have and that they're set as active identity

### subset as needed ###

Idents(Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd
Seu_AIBS_obj <- subset(Seu_AIBS_obj, subset = NeuN_Region %in% c("MTG_Neuronal", 
                                                                 "V1C_Neuronal", 
                                                                 "M1lm_Neuronal", 
                                                                 "S1ul_Neuronal", 
                                                                 "S1lm_Neuronal", 
                                                                 "M1ul_Neuronal", 
                                                                 "A1C_Neuronal"), invert = TRUE) #remove all non-CgG neurons

table(Seu_AIBS_obj$subclass_label_expanded) #double check what subclasses we have
Seu_AIBS_obj <- subset(Seu_AIBS_obj, subset = subclass_label_expanded == "L4 IT", invert = TRUE) #remove L4 IT for n = 1

### find markers

Idents(Seu_AIBS_obj) <- "subclass_label_expanded_L35IT" #assign proper labels
Idents(Seu_AIBS_obj) <- "class_label" #assign proper labels

new_AIBS_markers_mast_expIT_ALL <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 1.1, min.pct = .15, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_expIT_ALL <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 1.1, min.pct = .15, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### remove duplicates, if desired

#dup_list <- unique(new_AIBS_markers_mast_expIT_ALL[duplicated(new_AIBS_markers_mast_expIT_ALL$gene),"gene"]) #list of duplicated genes
#new_AIBS_markers_mast_expIT_ALL <- new_AIBS_markers_mast_expIT_ALL[!(new_AIBS_markers_mast_expIT_ALL$gene %in% dup_list),] #remove duplicated marker genes

#dup_list <- unique(new_AIBS_markers_roc_expIT_ALL[duplicated(new_AIBS_markers_roc_expIT_ALL$gene),"gene"]) #list of duplicated genes
#new_AIBS_markers_roc_expIT_ALL <- new_AIBS_markers_roc_expIT_ALL[!(new_AIBS_markers_roc_expIT_ALL$gene %in% dup_list),] #remove duplicated marker genes

#remove(dup_list) #clean temporary object

### finalize df

new_AIBS_markers_mast_expIT_ALL$group_gene <- paste0(new_AIBS_markers_mast_expIT_ALL$cluster, "_", new_AIBS_markers_mast_expIT_ALL$gene)
new_AIBS_markers_roc_expIT_ALL$group_gene <- paste0(new_AIBS_markers_roc_expIT_ALL$cluster, "_", new_AIBS_markers_roc_expIT_ALL$gene)

length(intersect(new_AIBS_markers_mast_expIT_ALL$group_gene, new_AIBS_markers_roc_expIT_ALL$group_gene)) #see intersect of marker genes
Result_df <- merge(new_AIBS_markers_roc_expIT_ALL, new_AIBS_markers_mast_expIT_ALL, by = "group_gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order

test <- Result_df[complete.cases(Result_df),]
identical(test$avg_log2FC.x, test$avg_log2FC.y)
identical(test$pct.1.x, test$pct.1.y)
identical(test$pct.2.x, test$pct.2.y)
identical(test$cluster.x, test$cluster.y)
identical(test$gene.x, test$gene.y) # as appropriate
remove(test)

Result_df <- Result_df %>% mutate(cluster = coalesce(cluster.x, cluster.y))
Result_df <- Result_df %>% mutate(avg_log2FC = coalesce(avg_log2FC.x, avg_log2FC.y))
Result_df <- Result_df %>% mutate(pct.1 = coalesce(pct.1.x, pct.1.y))
Result_df <- Result_df %>% mutate(pct.2 = coalesce(pct.2.x, pct.2.y))
Result_df <- Result_df %>% mutate(gene = coalesce(gene.x, gene.y))
Result_df <- Result_df[,c("gene", "cluster", "pct.1", "pct.2", "avg_log2FC", "avg_diff", "myAUC", "power", "p_val", "p_val_adj")]

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df)[c(3:4,8:12)] <- c("ensembl_id", "subclass", "roc_avg_diff","roc_myAUC", "roc_power", "MAST_p_val","MAST_p_val_adj") #rename some columns for clarity

write.csv(Result_df, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_class_ITexpand_WL35IT_lfct11_minpct15_dup.csv") #save/export results

#
#### get markers by pairwise comparisons between specific subclasses ####

### filter data
Idents(new_Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd to make it easy for us to remove neurons from MTG/CgG later

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj, idents = "MTG_Neuronal", invert = TRUE) #remove "MTG_Neuronal" cells

table(new_Seu_AIBS_obj_for_test$NeuN_Region) #double check what classes we have

Idents(new_Seu_AIBS_obj_for_test) <- "subclass_label" #assign "class_label" as the key grouping variable now

table(Idents(new_Seu_AIBS_obj_for_test)) #double check what classes we have

new_Seu_AIBS_obj_for_test <- subset(new_Seu_AIBS_obj_for_test, idents = "L4 IT", invert = TRUE) #remove "L4 IT" cell since there is an N of 1

### screen data

library("MAST") #as needed

### find markers

#make template
pairwise_results <- FindMarkers(new_Seu_AIBS_obj_for_test, 
                                slot = "data", 
                                ident.1 = "SST", 
                                ident.2 = "IT",
                                features = c("SST","CBLN2"),
                                test.use = "MAST")

pairwise_results$target_group <- "test"
pairwise_results$reference_group <- "test"
pairwise_results <- tibble::rownames_to_column(pairwise_results)


for (cellgroup in unique(Idents(new_Seu_AIBS_obj_for_test))){
  
  comparator_list <- as.character(unique(Idents(new_Seu_AIBS_obj_for_test)))
  comparator_list <- comparator_list[comparator_list!=cellgroup] 
  
  for (referencegroup in comparator_list) {
    
    holder <- FindMarkers(new_Seu_AIBS_obj_for_test, slot = "data", 
                          ident.1 = cellgroup, 
                          ident.2 = referencegroup,
                          features = new_CgG_results$gene,
                          logfc.threshold = 0, 
                          min.pct = 0, 
                          test.use = "MAST")
    holder$target_group <- cellgroup
    holder$reference_group <- referencegroup
    holder <- tibble::rownames_to_column(holder)
    pairwise_results <- rbind(pairwise_results, holder)
    
  }
  
}

### extract 2nd highest expressing group from df

pairwise_results <- pairwise_results[-(1:2),]
names(pairwise_results)[1] <- "gene"
pairwise_results$subclass_gene <- paste0(pairwise_results$target_group, "_", pairwise_results$gene)
new_CgG_results$subclass_gene <- paste0(new_CgG_results$subclass, "_", new_CgG_results$gene)
selected_pairwise_results <- pairwise_results[pairwise_results$subclass_gene %in% unique(new_CgG_results$subclass_gene),]

Result_df <- selected_pairwise_results[1,]
Result_df[,"subclass_gene"] <- "test"
Result_df_pct <- selected_pairwise_results[1,]
Result_df_pct[,"subclass_gene"] <- "test"

for (cellgroup in unique(selected_pairwise_results$target_group)){
  
  subclass_holder <- selected_pairwise_results[selected_pairwise_results$target_group == cellgroup,]
  
  for (gene_of_interest in unique(subclass_holder$gene)) {
    gene_holder <- selected_pairwise_results[selected_pairwise_results$gene == gene_of_interest,]
    holder <- gene_holder[which.min(gene_holder$avg_logFC),]
    Result_df <- rbind(Result_df, holder)
    holder <- gene_holder[which.max(gene_holder$pct.2),]
    Result_df_pct <- rbind(Result_df_pct, holder)
    
  }
  
}

### finalize + export results

Result_df <- Result_df[-1,]
Result_df <- Result_df[c(9,8,3,5,2)]
Result_df_pct <- Result_df_pct[-1,]
Result_df_pct <- Result_df_pct[c(9,8,3,5,2)]

colnames(Result_df)[2:5] <- paste0("2ByFC_",colnames(Result_df)[2:5])
colnames(Result_df_pct)[2:5] <- paste0("2ByPct_",colnames(Result_df_pct)[2:5])

new_CgG_results_pairdata <- merge(new_CgG_results, Result_df, by = "subclass_gene")
new_CgG_results_pairdata <- merge(new_CgG_results_pairdata, Result_df_pct, by = "subclass_gene")

new_CgG_results_pairdata <- new_CgG_results_pairdata[-1]

#### get markers from collapsing all instances of IT in Hodge (Old/unused) ####

unique(Idents(Seu_AIBS_obj))
Seu_AIBS_obj$subclass_label_collapsed <- Seu_AIBS_obj$subclass_label

Seu_AIBS_obj$subclass_label_collapsed[Seu_AIBS_obj$subclass_label_collapsed=="L4 IT"] <- "IT" #collapse to IT
Seu_AIBS_obj$subclass_label_collapsed[Seu_AIBS_obj$subclass_label_collapsed=="L5/6 IT Car3"] <- "IT" #collapse to IT

Idents(Seu_AIBS_obj) <- "subclass_label_collapsed"
unique(Idents(Seu_AIBS_obj))

table(Idents(Seu_AIBS_obj))
table(Seu_AIBS_obj$subclass_label)

new_AIBS_markers_mast_ITcol <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_ITcol <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### for mast
length(unique(new_AIBS_markers_mast_ITcol$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_mast_ITcol[duplicated(new_AIBS_markers_mast_ITcol$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_mast_ITcol <- new_AIBS_markers_mast_ITcol[!(new_AIBS_markers_mast_ITcol$gene %in% dup_list),] #remove duplicated marker genes

### for roc
length(unique(new_AIBS_markers_roc_ITcol$gene)) #check for unique marker genes
dup_list <- unique(new_AIBS_markers_roc_ITcol[duplicated(new_AIBS_markers_roc_ITcol$gene),"gene"]) #list of duplicated genes
new_AIBS_markers_roc_ITcol <- new_AIBS_markers_roc_ITcol[!(new_AIBS_markers_roc_ITcol$gene %in% dup_list),] #remove duplicated marker genes

length(intersect(new_AIBS_markers_mast_ITcol$gene, new_AIBS_markers_roc_ITcol$gene)) #see intersect of marker genes

### combine
new_AIBS_markers_mast_ITcol$group_gene <- paste0(new_AIBS_markers_mast_ITcol$cluster, "_", new_AIBS_markers_mast_ITcol$gene)
new_AIBS_markers_roc_ITcol$group_gene <- paste0(new_AIBS_markers_roc_ITcol$cluster, "_", new_AIBS_markers_roc_ITcol$gene)

Result_df_ITcol <- merge(new_AIBS_markers_roc_ITcol, new_AIBS_markers_mast_ITcol, by = "group_gene", all.x = T, all.y = T) #combine the marker df

Result_df_ITcol <- Result_df_ITcol %>% mutate(gene = coalesce(gene.x, gene.y))
Result_df_ITcol <- Result_df_ITcol %>% mutate(cluster = coalesce(cluster.x, cluster.y))
Result_df_ITcol <- Result_df_ITcol %>% mutate(avg_logFC = coalesce(avg_logFC.x, avg_logFC.y))
Result_df_ITcol <- Result_df_ITcol %>% mutate(pct.1 = coalesce(pct.1.x, pct.1.y))
Result_df_ITcol <- Result_df_ITcol %>% mutate(pct.2 = coalesce(pct.2.x, pct.2.y))

### to add entrez and ensembl IDs to the output/result df

Result_df_ITcol <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df_ITcol, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno

Result_df_ITcol$ensembl_gene_id[is.na(Result_df_ITcol$ensembl_gene_id)] <- "NA"
Result_df_ITcol$has_ensembl <- Result_df_ITcol$ensembl_gene_id != "NA"
#### get celltype markers from Mathys dataset (Old/unused - part of validation attempts) ####

# prep metadata

mathys_meta$matched_group <- "Unmatched" #add var for class conversion

for (cellgroup in best_mathys_subclasses_df$pre.cluster) {
  mathys_meta[mathys_meta$pre.cluster == cellgroup, "matched_group"] <- best_mathys_subclasses_df[best_mathys_subclasses_df$pre.cluster == cellgroup, "subclass"]
} #convert classes

table(mathys_meta[,c("pre.cluster", "matched_group")]) #see class conversion

rownames(mathys_meta) <- mathys_meta$TAG

###Seurat

Seu_mathys_obj <- CreateSeuratObject(counts = mathys_expr, meta.data = mathys_meta) #make seurat object

# checks

avg_expr_mat <- t(avg_expr_mat)
avg_expr_mat <- as.data.frame(avg_expr_mat)
colnames(avg_expr_mat) <- 1:21
avg_expr_mat <- avg_expr_mat[,-18]

Idents(Seu_mathys_obj) <- "pre.cluster"
test <- AverageExpression(Seu_mathys_obj, use.counts = TRUE)
test <- test$RNA
test <- test[matching_genes,]

identical(avg_expr_mat[order(names(avg_expr_mat))], test[order(names(test))])
remove(test)

test <- as.data.frame(Seu_mathys_obj@meta.data) #extracting metadata from Seurat object
test <- test[,4:ncol(test)] #removing extra columns made by Seurat
identical(test, mathys_meta) #identical test, should be TRUE

remove(test) #remove temporary test object

# screen proportion and pct

Seu_mathys_obj <- NormalizeData(Seu_mathys_obj, normalization.method = "LogNormalize", scale.factor = 1000000) #normalize data

Idents(Seu_mathys_obj) <- "matched_group" #set idents

mathys_markers_screen_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = .1, min.pct = .1, only.pos = TRUE, return.thresh = .05, test.use = "MAST")
mathys_markers_screen_mast[as.character(unique(Idents(Seu_mathys_obj))),] #see fold change and pct values, use for cut-off

# get markers

mathys_markers_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.2, min.pct = .25, only.pos = TRUE, test.use = "MAST") #find markers
mathys_markers_roc <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.2, min.pct = .25, only.pos = TRUE, test.use = "roc") #find markers using roc

#mathys_markers_stringent_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
#mathys_markers_stringent_roc <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

#check duplicates, overlaps, then consolidate

length(unique(mathys_markers_mast$gene)) #check for unique marker genes
length(unique(mathys_markers_roc$gene)) #check for unique marker genes

mathys_markers_mast$clust_gene <- paste0(mathys_markers_mast$cluster,"_",mathys_markers_mast$gene) #make cluster-gene combo for match check
mathys_markers_roc$clust_gene <- paste0(mathys_markers_roc$cluster,"_",mathys_markers_mast$gene)

length(intersect(unique(mathys_markers_mast$clust_gene), unique(mathys_markers_roc$clust_gene))) #check overlap between methods
identical(unique(mathys_markers_mast$clust_gene), unique(mathys_markers_roc$clust_gene)) #check overlap between methods

Result_df <- merge(mathys_markers_mast[,c(8,7,6,4,5,2,1,3)], mathys_markers_roc[,c(1,3,9)], by = "clust_gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
Result_df <- Result_df[,c(2,3,8,4,6,9,10,7,5)] #reorder

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df) <- colnames(new_CgG_results)

mathys_results <- Result_df #save results

write.csv(Result_df, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/mathys_results.csv") #export results

# get comparisons for DE genes with no thresholds
mathys_markers_mega_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 0, min.pct = 0, only.pos = TRUE, return.thresh = 1, test.use = "MAST")

# check and plot results
mathys_markers_mega_mast$clusterchar <- as.character(mathys_markers_mega_mast$cluster)
mathys_markers_mega_mast[mathys_markers_mega_mast$clusterchar == "L5.6.NP","clusterchar"] <- "L5/6 NP"
mathys_markers_mega_mast$markercheck <- paste0(mathys_markers_mega_mast$clusterchar, "_", mathys_markers_mega_mast$gene)
test <- new_CgG_results
test$markercheck <- paste0(test$subclass, "_", test$gene)
length(intersect(mathys_markers_mega_mast$markercheck, test$markercheck))

mathys_markers_mega_mast$volcano_group <- "Not significant"
mathys_markers_mega_mast[mathys_markers_mega_mast$p_val < 0.05, "volcano_group"] <- "Independently significant"
mathys_markers_mega_mast[mathys_markers_mega_mast$p_val_adj < 0.05, "volcano_group"] <- "Significant after correction"

for (cellgroup in unique(new_CgG_results$subclass)){
  temp_markers <- new_CgG_results[new_CgG_results$subclass == cellgroup, "gene"]
  mathys_markers_mega_mast[((mathys_markers_mega_mast$clusterchar == cellgroup) &
                              (mathys_markers_mega_mast$gene %in% temp_markers))  
                           , "volcano_group"] <- "Hodge marker"
} 

table(mathys_markers_mega_mast$volcano_group)

table(mathys_markers_mega_mast[mathys_markers_mega_mast$volcano_group == "Hodge marker", "cluster"])
table(new_CgG_results[new_CgG_results$subclass %in% unique(mathys_markers_mega_mast$clusterchar),"subclass"])

mathys_average <- AverageExpression(Seu_mathys_obj)
mathys_average <- mathys_average$RNA

mathys_markers_mega_mast[mathys_markers_mega_mast$volcano_group == "Hodge marker" & mathys_markers_mega_mast$cluster == "IT",]

intersect(mathys_markers_mega_mast[mathys_markers_mega_mast$cluster == "SST", "gene"],new_CgG_results[new_CgG_results$subclass == "SST", "gene"])
setdiff(new_CgG_results[new_CgG_results$subclass == "IT", "gene"], mathys_markers_mega_mast[mathys_markers_mega_mast$cluster == "IT", "gene"])

intersect(test[(test$gene %in% row.names(mathys_average)) &
                 (test$subclass %in% unique(mathys_markers_mega_mast$clusterchar))
               ,"markercheck"], mathys_markers_mega_mast$markercheck)

mathys_markers_mega_mast$p_val_adj_forplot <- log10(mathys_markers_mega_mast$p_val_adj)*(-1)
mathys_markers_mega_mast[mathys_markers_mega_mast$p_val_adj_forplot == Inf, "p_val_adj_forplot"] <- 400

library(ggplot2)

ggplot(mathys_markers_mega_mast[mathys_markers_mega_mast$volcano_group == "Hodge marker",], aes(x=avg_logFC, y= p_val_adj_forplot, color=volcano_group)) + 
  geom_point(size = 1) +
  facet_wrap(~cluster, scales = "fixed") +
  theme(legend.position = "none")


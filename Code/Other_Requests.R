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

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

#### get celltype markers from Mathys ####

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

mathys_markers_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 1.25, min.pct = .15, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
mathys_markers_roc <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 1.25, min.pct = .15, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

#mathys_markers_stringent_mast <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
#mathys_markers_stringent_roc <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

#check duplicates, overlaps, then consolidate

length(unique(mathys_markers_mast$gene)) #check for unique marker genes
length(unique(mathys_markers_stringent_mast$gene)) #check for unique marker genes

mathys_markers_mast$clust_gene <- paste0(mathys_markers_mast$cluster,"_",mathys_markers_mast$gene) #make cluster-gene combo for match check
mathys_markers_roc$clust_gene <- paste0(mathys_markers_roc$cluster,"_",mathys_markers_mast$gene)

length(intersect(unique(mathys_markers_mast$clust_gene), unique(mathys_markers_roc$clust_gene))) #check overlap between methods
identical(unique(mathys_markers_mast$clust_gene), unique(mathys_markers_roc$clust_gene)) #check overlap between methods

Result_df <- merge(mathys_markers_mast[,c(8,7,6,4,5,2,1,3)], mathys_markers_roc[,c(1,3,9)], by = "clust_gene", all.x = TRUE, all.y = TRUE) #combine the marker df; may want to double check the indices/order
Result_df <- Result_df[,c(2,3,8,4,6,9,10,7,5)] #reorder

Result_df <- merge(Gene_anno[,c("gene", "entrez_id", "ensembl_gene_id")], Result_df, by = "gene", all.y = TRUE) #add entrez and ensembl ids, keeping all results, even if they don't have a corresponding entry from Gene-Anno
colnames(Result_df) <- colnames(new_CgG_results)

write.csv(Result_df, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/mathys_results.csv") #save/export results

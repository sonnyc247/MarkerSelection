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

#### Seurat data integration and comparison between Mathys and Hodge ####

new_Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")
Seu_mathys_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds")

Idents(new_Seu_AIBS_obj) <- "subclass_label"
Idents(Seu_mathys_obj)

new_Seu_AIBS_obj <- FindVariableFeatures(new_Seu_AIBS_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Seu_mathys_obj <- FindVariableFeatures(Seu_mathys_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

tanchors <- FindTransferAnchors(reference = new_Seu_AIBS_obj, query = Seu_mathys_obj, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = new_Seu_AIBS_obj$subclass_label, dims = 1:30)
Seu_mathys_obj <- AddMetaData(Seu_mathys_obj, metadata = predictions)

test <- data.frame(unclass(table(Seu_mathys_obj$predicted.id, Seu_mathys_obj$matched_group)))

Idents(Seu_mathys_obj) <- "predicted.id"

mathys_markers_mega_mast_anchor <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 0, min.pct = 0, only.pos = TRUE, return.thresh = 1, test.use = "MAST")

predicted_id_match <- data.frame(unclass(table(Seu_mathys_obj$predicted.id, Seu_mathys_obj$matched_group)))

mathys_markers_selected_mast_anchor <- FindAllMarkers(Seu_mathys_obj, slot = "data", features=intersect(new_CgG_results$gene, row.names(mathys_average)), logfc.threshold = 0, min.pct = 0, return.thresh = 1, test.use = "MAST")

# checks and plots

mathys_markers_mega_mast_anchor$markercheck <- paste0(mathys_markers_mega_mast_anchor$cluster, "_", mathys_markers_mega_mast_anchor$gene)
test <- new_CgG_results
test$markercheck <- paste0(test$subclass, "_", test$gene)
length(intersect(mathys_markers_mega_mast_anchor$markercheck, test$markercheck))

mathys_markers_mega_mast_anchor$volcano_group <- "Not significant"
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val < 0.05, "volcano_group"] <- "Independently significant"
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val_adj < 0.05, "volcano_group"] <- "Significant after correction"

for (cellgroup in unique(new_CgG_results$subclass)){
  temp_markers <- new_CgG_results[new_CgG_results$subclass == cellgroup, "gene"]
  mathys_markers_mega_mast_anchor[((mathys_markers_mega_mast_anchor$cluster == cellgroup) &
                                   (mathys_markers_mega_mast_anchor$gene %in% temp_markers))  
                                    , "volcano_group"] <- "Hodge marker"
} 

table(mathys_markers_mega_mast_anchor$volcano_group)

table(mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$volcano_group == "Hodge marker", "cluster"])
table(new_CgG_results[new_CgG_results$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster),"subclass"])

mathys_average <- AverageExpression(Seu_mathys_obj)
mathys_average <- mathys_average$RNA

intersect(test[(test$gene %in% row.names(mathys_average)) &
               (test$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster))
               ,"markercheck"], mathys_markers_mega_mast_anchor$markercheck)

setdiff(test[(test$gene %in% row.names(mathys_average)) &
             (test$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster)),"markercheck"], 
             mathys_markers_mega_mast_anchor$markercheck)

mathys_markers_mega_mast_anchor$p_val_adj_forplot <- log10(mathys_markers_mega_mast_anchor$p_val_adj)*(-1)
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val_adj_forplot == Inf, "p_val_adj_forplot"] <- 400

library(ggplot2)

ggplot(mathys_markers_mega_mast_anchor, aes(x=avg_logFC, y= p_val_adj_forplot, color=volcano_group)) + 
  geom_point(size = .25) +
  facet_wrap(~cluster, scales = "fixed")

ggplot(mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$volcano_group == "Hodge marker",], aes(x=avg_logFC, y= p_val_adj_forplot, color=volcano_group)) + 
  geom_point(size = 1) +
  facet_wrap(~cluster, scales = "fixed") +
  theme(legend.position = "none")

Seu_mathys_obj$SST_or_Not <- Seu_mathys_obj$predicted.id
Seu_mathys_obj@meta.data[Seu_mathys_obj$SST_or_Not != "SST","SST_or_Not"] <- "Not_SST" #collapsing the non-sst groups
Idents(Seu_mathys_obj) <- "SST_or_Not"
VlnPlot(Seu_mathys_obj, features = "SST", slot = "data")
mathys_average_SST <- AverageExpression(Seu_mathys_obj)
mathys_average_SST <- mathys_average_SST$RNA

test <- t(mathys_average_SST["SST",])
test <- tibble::rownames_to_column(as.data.frame(test))
ggplot(data=test, aes(x=rowname, y=SST, fill=rowname)) + scale_fill_manual(values=c("white", "black")) +
  geom_bar(stat="identity", color = "black") + theme_bw() +theme(axis.title.x=element_blank()) + theme(legend.position = "none")

# gdp

library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis

gdp_plot <- t(as.data.frame(Seu_mathys_obj[["RNA"]]@data)) #get transposed lnCPM matrix
gdp_plot <- t(as.data.frame(new_Seu_AIBS_obj[["RNA"]]@data)) #get transposed lnCPM matrix

# format count/expression matrix for group_dot_plot smooth running
library(tibble)
gdp_plot <- rownames_to_column(as.data.frame(gdp_plot)) #get sample names as a column
colnames(gdp_plot)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_plot) <- gdp_plot$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(Seu_mathys_obj@meta.data) #create metadata copy for group_dot_plot
gdp_anno <- as.data.frame(new_Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot

colnames(gdp_anno)[4] <- "sample_name"
#gdp_anno$subclass_label <- gdp_anno$predicted.id
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"
gdp_anno[gdp_anno$subclass_label == "SST", "subclass_color"] <- "red"

gdp_anno$it2_query_subclass_id <- gdp_anno$it2_query_subclass_label #prepare "_id" needed for plotting
gdp_anno$it2_query_subclass_color <- gdp_anno$class_color #prepare "_color" needed for plotting; copy "class_color" for now
gdp_anno$it2_query_subclass_label <- factor(gdp_anno$it2_query_subclass_label, levels = c("Inh_SST", "Inh_VIP", "Inh_PVALB", "Inh_LAMP5", "Inh_PAX6", "Pyramidal", "Astro_FGFR3", "Oligo",  "OPC_MYT1", "Micro_C1QC", "Peri_MUSTN1", "Endo_CLDN5", "VLMC_CYP1B1")) #setting factor order

# set which genes to plot

gdp_markers <- new_CgG_results[new_CgG_results$subclass == "SST", "gene"] #get SST markers 
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

gdp_markers <- new_CgG_results[new_CgG_results$subclass == "SST", "gene"] #get SST markers
gdp_markers <- intersect(gdp_markers, row.names(mathys_average))
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

gdp_markers <- ACC_results[ACC_results$adapted_cluster_name == "PVALB", "gene"] #get PVALB markers 
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

# do the plot

gdp_anno_test <- merge(gdp_anno, Name_and_colour_scheme, by.x = "subclass_label", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
gdp_anno_test$old_labels <- gdp_anno_test$subclass_label
gdp_anno_test$subclass_label <- gdp_anno_test$Our_label
gdp_anno_test$subclass_label <- factor(gdp_anno_test$subclass_label, levels = c("Exc_IT",
                                                                              "Exc_L5 ET",
                                                                              "Exc_L5/6 IT Car3",
                                                                              "Exc_L5/6 NP",
                                                                              "Exc_L6 CT",  
                                                                              "Exc_L6b",
                                                                              "Inh_SST",
                                                                              "Inh_LAMP5",
                                                                              "Inh_PAX6",
                                                                              "Inh_PVALB",
                                                                              "Inh_VIP",
                                                                              "Astrocyte",
                                                                              "Endothelial",      
                                                                              "Microglia",
                                                                              "Oligodendrocyte",
                                                                              "OPC",
                                                                              "Pericyte",
                                                                              "VLMC"))

group_dot_plot(gdp_plot, 
               gdp_anno_test, 
               genes = gdp_markers, 
               grouping = "subclass", 
               log_scale = TRUE,
               font_size = 14,
               max_size = 29,
               rotate_counts = TRUE)

ggsave(dpi = 300, 
       limitsize = F,
       filename = "Hodge_SST_gdp.jpg",
       device = "jpeg")

ggsave(filename = "L5_6_NP_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)

### systemic gdp ###



#volcano

mathys_markers_selected_mast_anchor$volcano_group <- "Not significant"
mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$p_val < 0.05, "volcano_group"] <- "Significant"

for (cellgroup in unique(new_CgG_results$subclass)){
  temp_markers <- new_CgG_results[new_CgG_results$subclass == cellgroup, "gene"]
  mathys_markers_selected_mast_anchor[((mathys_markers_selected_mast_anchor$cluster == cellgroup) &
                                       (mathys_markers_selected_mast_anchor$gene %in% temp_markers))  
                                       , "volcano_group"] <- "Hodge marker"
} 

table(mathys_markers_selected_mast_anchor$volcano_group)

mathys_markers_selected_mast_anchor$p_val_forplot <- log10(mathys_markers_selected_mast_anchor$p_val)*(-1)
mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$p_val_forplot == Inf, "p_val_forplot"] <- 400

groups_of_interest <- c("SST", "VIP", "Oligodendrocyte", "Microglia")
groups_of_interest <- c("VLMC", "Pericyte", "Endothelial", "Microglia")

ggplot(mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster %in% groups_of_interest,], aes(x=avg_logFC, y= p_val_forplot, color=volcano_group)) + 
  geom_point(size = .75) +
  facet_wrap(~cluster, scales = "fixed") +
  scale_color_manual(values=c("red", "black", "orange")) +
  theme_bw() 

ggplot(mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster %in% groups_of_interest
                                           & mathys_markers_selected_mast_anchor$volcano_group == "Hodge marker",], aes(x=avg_logFC, y= p_val_forplot, color=volcano_group)) + 
  geom_point(size = .75) +
  facet_wrap(~cluster, scales = "fixed") +
  scale_color_manual(values=c("red", "black", "orange")) +
  theme_bw() +
  theme(legend.position = "none")

### log fc plot

mathys_markers_selected_mast_anchor$clust_gene <- paste0(mathys_markers_selected_mast_anchor$cluster,"_",mathys_markers_selected_mast_anchor$gene)
test2 <- new_CgG_results
test2$clust_gene <- paste0(test2$subclass,"_",test2$gene)
test <- merge(test2, mathys_markers_selected_mast_anchor, by = "clust_gene")

names(test)[c(8,14)] <- c("Hodge_logFC", "Mathys_logFC")

ggplot(test[test$cluster %in% groups_of_interest,], aes(x=Mathys_logFC, y= Hodge_logFC)) + 
  geom_point(size = .75) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~cluster, scales = "fixed") +
  theme_bw() 

test3 <- mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster == "SST",]

library(metap)

for (cellgroup in unique(test$subclass)){
  print(cellgroup)
  print(sumlog(test[test$subclass == cellgroup,"p_val"]))
}

source('https://raw.githubusercontent.com/yaowuliu/ACAT/master/R/ACAT.R')

for (cellgroup in unique(test$subclass)){
  print(cellgroup)
  print(ACAT(test[test$subclass == cellgroup,"p_val"]))
}

#### get markers by pairwise comparisons ####

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



#### QC based on sex-specific genes ####

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

#### systemic hodge gdp ####

library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis

gdp_plot <- t(as.data.frame(new_Seu_AIBS_obj[["RNA"]]@data)) #get transposed lnCPM matrix

# format count/expression matrix for group_dot_plot smooth running
library(tibble)
gdp_plot <- rownames_to_column(as.data.frame(gdp_plot)) #get sample names as a column
colnames(gdp_plot)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_plot) <- gdp_plot$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(new_Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot

#colnames(gdp_anno)[4] <- "sample_name"
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"

# set which genes to plot

for (cellgroup in names(humanMarkersCommon)[-5]) {
  gdp_anno[gdp_anno$subclass_label == cellgroup, "subclass_color"] <- "red"
  temp_markers <- unlist(humanMarkersCommon[cellgroup])
  gdp_markers <- new_CgG_results_pairdata[new_CgG_results_pairdata$gene %in% temp_markers, c("gene", "roc_power")] #get group markers 
  gdp_markers <- gdp_markers[order(gdp_markers$roc_power, decreasing = T),]
  
  group_dot_plot(gdp_plot, 
                 gdp_anno, 
                 genes = gdp_markers$gene, 
                 grouping = "subclass", 
                 log_scale = TRUE,
                 font_size = 10,
                 max_size = 20,
                 rotate_counts = TRUE)
  
  if (cellgroup == "L5/6 NP"){
  ggsave(filename = "L5_6_NP_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)
  } else if (cellgroup == "L5/6 IT Car3"){
    ggsave(filename = "L5_6_IT_Car3_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)  
  } else{
  ggsave(filename = paste0(cellgroup,"_gdp.pdf"), path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)  
  }
  
  gdp_anno[gdp_anno$subclass_label == cellgroup, "subclass_color"] <- "white"

  }

#### Other Misc ####

### average expression

hodge_average <- AverageExpression(new_Seu_AIBS_obj_for_test)
hodge_average <- hodge_average$RNA
hodge_average <- hodge_average[new_CgG_results$gene,]
hodge_average <- tibble::rownames_to_column(hodge_average)
colnames(hodge_average) <- paste0(colnames(hodge_average), "_mean_expression")
colnames(hodge_average)[1] <- "gene"
new_CgG_results_pairdata <- merge(new_CgG_results_pairdata, hodge_average, by = "gene")

write.csv(new_CgG_results_pairdata, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/new_CgG_results_pairdata.csv", row.names = FALSE) #save/export results

### Jordan marker filtering

Jordan_CTX <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Inputs/SUPPL_S1.csv", stringsAsFactors=FALSE)
Jordan_CTX_filtered <- Jordan_CTX[Jordan_CTX$avg_logFC > 1,]
length(unique(Jordan_CTX_filtered$gene))
duplicated_genes <- unique(Jordan_CTX_filtered[duplicated(Jordan_CTX_filtered$gene), "gene"])
Jordan_CTX_filtered <- Jordan_CTX_filtered[!(Jordan_CTX_filtered$gene %in% duplicated_genes),]

write.csv(Jordan_CTX_filtered, "Jordan_Mouse_Markers_Filtered.csv")

### Annotating cell types and colours

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

#### Pseudobulk quick attempt ####

### Add some metadata
test <- Seu_mathys_obj@meta.data

row.names(mathys_meta_df) <- mathys_meta_df$TAG
mathys_meta_df <- mathys_meta_df[,34:37]

Seu_mathys_obj <- AddMetaData(Seu_mathys_obj, mathys_meta_df)

### Pseudobulk analysis

celltype_of_interest <- "SST"
Idents(Seu_mathys_obj) <- "broad.cell.type"
celltype_of_interest <- "Ex"
temp_mathys_obj <- subset(Seu_mathys_obj, idents = celltype_of_interest)
Idents(temp_mathys_obj) <- "projid"
mathys_subject_count_mtx_holder <- AverageExpression(temp_mathys_obj)
mathys_subject_count_mtx_holder <- mathys_subject_count_mtx_holder$RNA

mathys_temp_metadata <- temp_mathys_obj@meta.data
mathys_temp_count_mtx <- temp_mathys_obj@assays$RNA@counts
mathys_temp_count_mtx <- as.matrix(mathys_temp_count_mtx)

for (subjectID in unique(mathys_temp_metadata$projid)) {
  cell_list <- as.character(mathys_temp_metadata[mathys_temp_metadata$projid == subjectID, "TAG"])
  count_matrix_forsum <- mathys_temp_count_mtx[,cell_list]
  gene_sums <- rowSums(count_matrix_forsum)
  mathys_subject_count_mtx_holder[,as.character(subjectID)] <- gene_sums
}

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS('/external/rprshnas01/netdata_kcni/dflab/data/rosmap/phenotype/ROSmaster.rds')
# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex)
# two samples are missing diagnoses - these are AD cases according to the mathys paper 
ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset
ros_meta_small = ros_meta_small %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))

mathys_temp_metadata <- merge(mathys_temp_metadata , ros_meta_small, by = 'projid')
mathys_temp_subject_metadata <- mathys_temp_metadata[,c(1,34:37)]
mathys_temp_subject_metadata <- unique(mathys_temp_subject_metadata)
row.names(mathys_temp_subject_metadata) <- mathys_temp_subject_metadata$projid

### Seurat method

temp_mathys_subject_obj <- CreateSeuratObject(counts = mathys_subject_count_mtx_holder, meta.data = mathys_temp_subject_metadata)
Idents(temp_mathys_subject_obj) <- "pathoAD"

library(DESeq2)
#Temp_DE_Results <- FindAllMarkers(temp_mathys_subject_obj, test.use = "DESeq2")
#names(Temp_DE_Results)[6] <- "pathoAD"

Temp_DE_Results <- FindMarkers(temp_mathys_subject_obj, ident.1 = "1", test.use = "DESeq2")
Test <- AverageExpression(temp_mathys_subject_obj, slot = "counts")
Test <- Test$RNA

write.csv(Temp_DE_Results, "Pseudobulk_SST_InitialRun1_ADvsControl.csv")

### direct use of deseq2

mathys_temp_subject_metadata$ncell <- 0

# get number of cells per subject

for (subjectID in unique(mathys_temp_metadata$projid)) {
  cell_list <- as.character(mathys_temp_metadata[mathys_temp_metadata$projid == subjectID, "TAG"])
  mathys_temp_subject_metadata[as.character(subjectID),"ncell"] <- length(cell_list)
}

# do deseq2

library(DESeq2)
#library(tidyverse)

mathys_temp_subject_metadata <- mathys_temp_subject_metadata %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))
str(mathys_temp_subject_metadata)

dds <- DESeqDataSetFromMatrix(countData = mathys_subject_count_mtx_holder, 
                              colData = mathys_temp_subject_metadata, 
                              design = ~ pathoAD + age_death + msex + ncell)

### first basic attempt

dds <- DESeq(dds)
resultsNames(dds)
Temp_DE_Results_noSeurat <- results(dds, name="pathoAD_1_vs_0")
Temp_DE_Results_noSeurat <- as.data.frame(Temp_DE_Results_noSeurat)

### second attempts from digging into seurat

dds <- DESeq2::estimateSizeFactors(object = dds)
dds <- DESeq2::estimateDispersions(object = dds, fitType = "local")
dds <- DESeq2::nbinomWaldTest(object = dds)
resultsNames(dds)
Temp_DE_Results_noSeurat <- DESeq2::results(object = dds,
                                            contrast = c("pathoAD", "1", "0"),
                                            alpha = 0.05)

#testresults <- lfcShrink(dds, coef="pathoAD_1_vs_0", type="apeglm")
#testresults <- results(dds, contrast=c("pathoAD","1","0"))

for (ADCondition in unique(mathys_temp_metadata$projid)) {
  cell_list <- as.character(mathys_temp_metadata[mathys_temp_metadata$projid == subjectID, "TAG"])
  count_matrix_forsum <- mathys_temp_count_mtx[,cell_list]
  gene_sums <- rowSums(count_matrix_forsum)
  mathys_subject_count_mtx_holder[,as.character(subjectID)] <- gene_sums
}

write.csv(Temp_DE_Results_noSeurat, "Pseudobulk_SST_InitialRun2_ADvsControl_DESEQ2withDesign.csv")

c("GOLT1B", "ATF6B", "DDRGK1", "TUBB2A", "BEX2", "ATPIF1", "RASGEF1B", "NGFRAP1", "LINGO1", "NTNG1")

table(mathys_temp_subject_metadata$pathoAD, mathys_temp_subject_metadata$msex)
t.test(mathys_temp_subject_metadata$age_death[mathys_temp_subject_metadata$pathoAD=="1"], 
       mathys_temp_subject_metadata$age_death[mathys_temp_subject_metadata$pathoAD=="0"])
t.test(mathys_temp_subject_metadata$ncell[mathys_temp_subject_metadata$pathoAD=="1"], 
       mathys_temp_subject_metadata$ncell[mathys_temp_subject_metadata$pathoAD=="0"])

#### collapse AIBS IT ####

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

#### CgG and MTG only ####

unique(new_Seu_AIBS_obj$NeuN_Region)

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

#### CgG and MTG only 2.0 lfct ####

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

#### PVALB plot ####

Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj.rds")
names(Seu_AIBS_obj@meta.data)
Idents(Seu_AIBS_obj) <- "subclass_label"
Seu_AIBS_obj <- subset(Seu_AIBS_obj, idents = "PVALB", invert = F) #remove "CgG_Neuronal" cells
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

#### get excit/inhib markers from MTG + CgG ####

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

#### get neu/nonneu markers from MTG + CgG ####


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

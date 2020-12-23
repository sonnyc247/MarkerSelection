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
gdp_anno$subclass_label <- gdp_anno$predicted.id
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

group_dot_plot(gdp_plot, 
               gdp_anno, 
               genes = gdp_markers, 
               grouping = "subclass", 
               log_scale = TRUE,
               font_size = 10,
               max_size = 20,
               rotate_counts = TRUE)

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


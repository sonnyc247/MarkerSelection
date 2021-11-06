#### Packages ####

library(Seurat)
library(tidyr)
library(dplyr)
library(Matrix)
library(magrittr)
library(tidyverse)

#### Hodge et al (copied from 1_Loading_and_selecting_data.Rmd and 2_Find_DE_genes_in_Seurat.Rmd (in Run_V2)) ####

### metadata
new_metadata <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/metadata.csv", stringsAsFactors=FALSE) #this file corresponds to https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv, we are using FALSE for strings as factors so that we can edit and order them later

names(new_metadata) #see columns that we have to work with, "outlier_call" and "subclass_label" are of particular interest to us here

table(new_metadata$outlier_call, useNA = "ifany") #see the outlier data; there are 1985 "outliers", likely corresponding to "tome_sample_meta_filtered$class_label == 'Exclude'" from Run_V1

new_metadata_filtered <- new_metadata[new_metadata$outlier_call == "False",] #filter out outliers, we do not want them in our analyses later on

table(new_metadata_filtered$subclass_label, useNA = "ifany") #for this demo/workbook, we will be working with subclass labels; this is to let us see what cell groups (subclasses) we have

table(new_metadata_filtered[,c("subclass_label", "region_label")], useNA = "ifany") #we can see from this that non-neuronal cell types have drastically decreased sample sizes when split to brain regions

table(new_metadata_filtered[,c("class_label", "region_label")], useNA = "ifany") #re-illustrate the point above

new_metadata_filtered$NeuN <- new_metadata_filtered$class_label #make a new variable for classifying cells as neuron/nonneuron
new_metadata_filtered[new_metadata_filtered$NeuN != "Non-neuronal", "NeuN"] <- "Neuronal" #change gabaergic/glutamatergic to Neuronal, store in "NeuN" variable
table(new_metadata_filtered[,c("NeuN", "region_label")], useNA = "ifany") #see our sample distribution now when looking at neurons and non-neurons

new_metadata_filtered$NeuN_Region <- paste0(new_metadata_filtered$region_label, "_", new_metadata_filtered$NeuN) #may be useful for filtering later on
table(new_metadata_filtered$NeuN_Region, useNA = "ifany") #final view of our sample distribution

### count_matrix
new_count_matrix <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv", stringsAsFactors=FALSE) #this file corresponds to https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv 

new_count_matrix[1:10,1:10] #preview of uploaded data
row.names(new_count_matrix) <- new_count_matrix$sample_name #set sample names as rownames - this is a formatting step important for later
new_count_matrix <- new_count_matrix[new_metadata_filtered$sample_name,] #filter for the cells we specified previously (removing outliers)
new_count_matrix <- new_count_matrix[,2:ncol(new_count_matrix)] #the df so far has samples as the first column, which is generally not desired for later on - the data needs to be in a numeric matrix. Thus, we set the row names to the first column and remove the first column; careful to run this only ONCE
new_count_matrix <- t(new_count_matrix) #count matrix will need to be transposed for later steps; count matrices should be samples as columns and genes as rows
new_count_matrix[1:10,1:7] #preview of uploaded data again

### Seurat
library(Seurat) #load the needed library

#check that the inputs are formatted correctly

new_metadata_filtered[1:10,1:5] #this is formatted correctly (in terms of rows and columns) but will need to set sample names as rownames
row.names(new_metadata_filtered) <- new_metadata_filtered$sample_name #setting row names to sample names
new_metadata_filtered[1:10,1:5] #checking the metadata df again
new_count_matrix[1:10,1:7] #checking count data again

Seu_AIBS_obj <- CreateSeuratObject(counts = t(new_count_matrix), meta.data = new_metadata_filtered) #this creates the Seurat object from our data

### comfirming count matrix
test <- as.data.frame(Seu_AIBS_obj@assays$RNA@counts) #extracting count matrix from Seurat object

#the follow need to be done to change the format of "test", not any data within it; this will allow us to use "identical" to compare "test" and the original count data in new_count_matrix
test <- as.data.frame(test) #change to df
temp_rownames <- row.names(test) #save rownames
test <- lapply(test, as.integer) #change all columns to integer
test <- as.data.frame(test) #change back to df
row.names(test) <- temp_rownames #reassign rownames

identical(test, new_count_matrix) #indentical check; will be TRUE

remove(temp_rownames) #remove temporary object

### confirming metadata

test <- as.data.frame(Seu_AIBS_obj@meta.data) #extracting metadata from Seurat object
test <- test[,4:ncol(test)] #removing extra columns made by Seurat
identical(test, new_metadata_filtered) #identical test, should be TRUE

remove(test) #remove temporary test object

### normalize data

Seu_AIBS_obj <- NormalizeData(Seu_AIBS_obj, normalization.method = "LogNormalize", scale.factor = 1000000)

Seu_AIBS_obj[["RNA"]]@data[1:10,1:7] #preview/see our normalized data
Seu_AIBS_obj[["RNA"]]@counts[1:10,1:7] #as opposed to our previous counts


#### Mathys et al (copied from finding markers in marker_finding_record.R and from mathys_cluster_annotation.R (in Run_V2)) ####

# define path to data
data_path = '/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/gene_expresssion_processed/'

# read datasets into environment
mathys_expr = readMM(file = paste0(data_path, 'filtered_count_matrix.mtx'))
mathys_genes = read.csv(file = paste0(data_path, 'filtered_gene_row_names.txt'), header = F)
mathys_meta = read.csv(file = paste0(data_path, 'filtered_column_metadata.txt'), sep = '\t')

# this assigns row and column names to mathys_expr sparse matrix
dimnames(mathys_expr) = list(mathys_genes %>% unlist %>% as.character(), mathys_meta %>% pull(TAG))

# pulls sonny's cell type specific markers annotated at the subclass level
sonny_markers = read.csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/new_CgG_results.csv'))

### define a simple correlation based classifier to map mathys pre.clusters onto AIBS /Hodge subclasses
matching_genes = intersect(mathys_genes %>% unlist, sonny_markers$gene) %>% make.names()

# get subclass level averaged data from sonny from github
allen_human_avg_expr = read.csv(url('https://github.com/sonnyc247/MarkerSelection/raw/master/Data/Outputs/new_CgG_avg_count.csv'))
rownames(allen_human_avg_expr) = allen_human_avg_expr$X
allen_human_avg_expr_small = allen_human_avg_expr[matching_genes, -1]

# define a smaller expression matrix consisting of just marker genes
mathys_expr_markers_only = mathys_expr[matching_genes, ] 

# this is a crappy bit of code that just averages counts per pre.cluster from mathys
# WARNING: it does NO normalization in it's current state
cluster_means = lapply(1:21, function(ind){
  print(ind)
  #subcluster_label = unique(mathys_meta$pre.cluster)[ind]
  keep_cells = mathys_meta %>% filter(pre.cluster == ind) %>% pull(TAG) 
  temp_expr_mat = mathys_expr_markers_only[, keep_cells]
  return(temp_expr_mat %>% rowMeans())
  #print(precluster_ind)
})
names(cluster_means) = 1:21
avg_expr_mat = cluster_means %>% bind_rows() %>% as.data.frame()
rownames(avg_expr_mat) = matching_genes
# avg_expr_mat is the cluster averaged version of the mathys dataset at the pre.cluster level

# correlate every mathys pre.cluster with the subclass averages from AIBS
corr_results = cor(t(avg_expr_mat), allen_human_avg_expr_small, method = 'spearman', use = "pairwise.complete.obs")

# for each mathys pre.cluster, find the AIBS subclass cluster with the highest correlation
best_mathys_subclasses = apply(corr_results, 1, function(vec){
  return(names(which.max(vec)))
}) %>% unlist() 

# this is the final mapping from mathys pre.clusters onto AIBS subclasses
best_mathys_subclasses_df = data.frame(pre.cluster = c(1:17,19:21), subclass = best_mathys_subclasses) #untested after edit by sonny

#write.csv(file = "~/Dropbox/mathys/best_mathys_subclasses_df.csv", best_mathys_subclasses_df)

remove(sonny_markers)
remove(matching_genes)
remove(allen_human_avg_expr)
remove(allen_human_avg_expr_small)
remove(avg_expr_mat)
remove(corr_results)
remove(best_mathys_subclasses)

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


#### Cain et al ####

### read in data
Cain_cell_metadata <- readr::read_csv("/external/rprshnas01/external_data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv")
Cain_gene_metadata <- readr::read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv")
Cain_COO_matrix <- readr::read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv.gz")

Cain_cell_metadata <- readr::read_csv("/external/rprshnas01/external_data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv.gz")
Cain_gene_metadata <- readr::read_csv("/external/rprshnas01/external_data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv")
Cain_COO_matrix <- readr::read_csv("/external/rprshnas01/external_data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv.gz")


identical(as.numeric(row.names(Cain_COO_matrix)), as.numeric(Cain_COO_matrix$X1)) #check that indexes/rownames are the same

### format/assemble count matrix
Cain_sparse_matrix <- Matrix::sparseMatrix(i = as.numeric(Cain_COO_matrix$i),
                                           j = as.numeric(Cain_COO_matrix$j),
                                           x = as.numeric(Cain_COO_matrix$x),
                                           dimnames = list(as.character(Cain_gene_metadata$x), 
                                                           as.character(Cain_cell_metadata$cell_name))) #make matrix

Cain_sparse_matrix[1:20,1:3] #see matrix sample

### format metadata
Cain_cell_metadata <- separate(Cain_cell_metadata, col = "specimenID", remove = F, into = c("U1", "U2", "U3", "U4", "Patho_AD_Maybe")) #try to extract metadata from specimen ID

unique(Cain_cell_metadata$U1) #not varying variable
Cain_cell_metadata <- Cain_cell_metadata[,c(1,2,4:9)] #remove it

unique(Cain_cell_metadata$Patho_AD_Maybe) #assume this is patho for now, collapse 1 and 0
Cain_cell_metadata$Patho_AD_assumed <- Cain_cell_metadata$Patho_AD_Maybe
Cain_cell_metadata[Cain_cell_metadata$Patho_AD_Maybe == "pAD0", "Patho_AD_assumed"] <- "Path0"
Cain_cell_metadata[Cain_cell_metadata$Patho_AD_Maybe == "pAD1", "Patho_AD_assumed"] <- "Path1"
unique(Cain_cell_metadata$Patho_AD_assumed)

### make seurat object
Cain_cell_metadata <- as.data.frame(Cain_cell_metadata)
row.names(Cain_cell_metadata) <- Cain_cell_metadata$cell_name
Seu_cain_obj <- CreateSeuratObject(counts = Cain_sparse_matrix, meta.data = Cain_cell_metadata)

Idents(Seu_cain_obj) <- "Patho_AD_assumed" #check that seurat object was created properly
Idents(Seu_cain_obj)

test <- Seu_cain_obj@meta.data #another check that seurat object was created properly
identical(Cain_cell_metadata, test[,4:12])
remove(test)

remove(Cain_COO_matrix)
remove(Cain_sparse_matrix)
remove(Cain_gene_metadata)
remove(Cain_cell_metadata)

Seu_cain_obj <- NormalizeData(Seu_cain_obj, normalization.method = "LogNormalize", scale.factor = 1000000) #preprocess/normalize seurat object
Seu_cain_obj[["RNA"]]@data[1:10,1:10] #see normalized data
Seu_cain_obj[["RNA"]]@counts[1:10,1:10] #see counts for comparison

### add some additional metadata

meta_holder <- as.data.frame(read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/phenotype/ROSMAP_clinical.csv"))
cain_rosmap_metadata <- as.data.frame(read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21670836/ROSMAP_clinical_metadata_24donors.csv"))
cain_rosmap_metadata <- merge(cain_rosmap_metadata[, c("projid", "Sample_ID")], meta_holder, by = "projid", all.x = T)

meta_holder <- Seu_cain_obj@meta.data[, c("specimenID", "group_assumed")]
meta_holder <- tibble::rownames_to_column(meta_holder)
meta_holder$convertednames <- stringr::str_extract(meta_holder$rowname, "[^_]+")
meta_holder <- merge(meta_holder[,c("convertednames", "rowname")], cain_rosmap_metadata, by.x = "convertednames", by.y = "Sample_ID", all.x = T, all.y = F)
row.names(meta_holder) <- meta_holder$rowname
meta_holder <- meta_holder[,-2]
Seu_cain_obj <- AddMetaData(Seu_cain_obj, metadata = meta_holder)

cain_subject_metadata <- Seu_cain_obj@meta.data
cain_subject_metadata <- cain_subject_metadata[, c("group_assumed", "individualID", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")]
cain_subject_metadata <- unique(cain_subject_metadata)

table(cain_subject_metadata$group_assumed, cain_subject_metadata$braaksc)
table(cain_subject_metadata$group_assumed, cain_subject_metadata$ceradsc)
table(cain_subject_metadata$group_assumed, cain_subject_metadata$cogdx)
table(cain_subject_metadata$group_assumed, cain_subject_metadata$dcfdx_lv)
#### Zhou et al ####

Zhou_sample_list <- list.dirs(path = "/external/rprshnas01/kcni/ychen/SCC-Bashwork/Temp/syn21670836/", full.names = FALSE, recursive = TRUE)
Zhou_sample_list <- Zhou_sample_list[-1]

for (sample in Zhou_sample_list) {
  
  sample_dir <- paste0("/external/rprshnas01/kcni/ychen/SCC-Bashwork/Temp/syn21670836/", sample, "/")
  assign(sample, CreateSeuratObject(counts = Read10X(data.dir = sample_dir), project = sample))
  
}

Zhou_assay_metadata <- read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21670836/snRNAseqAD_TREM2_assay_scRNAseq_metadata.csv")
Zhou_biospecimen_metadata <- read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21670836/snRNAseqAD_TREM2_biospecimen_metadata.csv")

Seu_zhou_obj <- merge(AD1, y = c(AD10, AD11, AD12, AD13, AD2, AD3, AD5, AD7,  AD8,  AD9, C1, C11, C12, C2, C3, C4, C5, C6, C7, C8, C9, P10, P11, P12, P13, P2, P3, P5, P6, P7, P9), 
                      add.cell.ids = Zhou_sample_list, project = "Zhou_et_al")

table(Seu_zhou_obj$orig.ident)

Zhou_meta_holder <- Seu_zhou_obj@meta.data
Zhou_meta_holder <- tibble::rownames_to_column(Zhou_meta_holder)
Zhou_meta_holder <- merge(Zhou_meta_holder, Zhou_assay_metadata, by.x = "orig.ident", by.y = "specimenID", all.x = T)
Zhou_meta_holder <- merge(Zhou_meta_holder, Zhou_biospecimen_metadata, by.x = "orig.ident", by.y = "specimenID", all.x = T)

Zhou_meta_holder$AD_Group <- stringr::str_extract(Zhou_meta_holder$orig.ident, "[A-Z]+")
table(Zhou_meta_holder$orig.ident, Zhou_meta_holder$AD_Group)
Zhou_meta_holder$AD_Group <- factor(Zhou_meta_holder$AD_Group, levels = c("C", "AD", "P"))

row.names(Zhou_meta_holder) <- Zhou_meta_holder$rowname
Zhou_meta_holder <- Zhou_meta_holder[,5:33]
Seu_zhou_obj <- AddMetaData(Seu_zhou_obj, metadata = Zhou_meta_holder)

Idents(Seu_zhou_obj) <- "AD_Group"

Seu_zhou_obj <- NormalizeData(Seu_zhou_obj, normalization.method = "LogNormalize", scale.factor = 1000000) #preprocess/normalize seurat object
Seu_zhou_obj[["RNA"]]@data[1:10,1:10] #see normalized data
Seu_zhou_obj[["RNA"]]@counts[1:10,1:10] #see counts for comparison

# add some additional/better metadata

Zhou_rosmap_metadata <- as.data.frame(read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/phenotype/ROSMAP_clinical.csv"))
meta_holder <- Seu_zhou_obj@meta.data
meta_holder <- meta_holder[,c("AD_Group", "individualID")]
meta_holder <- tibble::rownames_to_column(meta_holder)
meta_holder <- merge(meta_holder, Zhou_rosmap_metadata, by = "individualID", all.x = T, all.y = F)
row.names(meta_holder) <- meta_holder$rowname
meta_holder <- meta_holder[,4:20]
Seu_zhou_obj <- AddMetaData(Seu_zhou_obj, metadata = meta_holder)

zhou_subject_metadata <- Seu_zhou_obj@meta.data
zhou_subject_metadata <- zhou_subject_metadata[, c("AD_Group", "individualID", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")]
zhou_subject_metadata <- unique(zhou_subject_metadata)

table(zhou_subject_metadata$AD_Group, zhou_subject_metadata$braaksc)
table(zhou_subject_metadata$AD_Group, zhou_subject_metadata$ceradsc)
table(zhou_subject_metadata$AD_Group, zhou_subject_metadata$cogdx)
table(zhou_subject_metadata$AD_Group, zhou_subject_metadata$dcfdx_lv)

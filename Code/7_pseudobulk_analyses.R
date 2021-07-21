#### Packages #### 
# Pseudobulk analysis; following workflow from https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# Also checked with Shreejoy's code from https://github.com/stripathy/mathys_analysis/blob/main/pseudobulk_analysis.R
# Packages:

library(Seurat)
library(dplyr)
library(magrittr)
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(purrr)
library(DESeq2)

#### Pick/load test data ####

Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_test_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_test_object <- subset(Seu_test_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

#### Initial data pre-processing from Seurat ####

### load in markers/genes of interest, select and filter for ones in Seurat object

new_MTGnCgG_lfct2_results <- readr::read_csv("~/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2_results.csv")
pseudo_genes_of_interest <- new_MTGnCgG_lfct2_results$gene
pseudo_genes_of_interest <- intersect(pseudo_genes_of_interest, rownames(Seu_test_object))

#Seu_test_object <- subset(Seu_test_object, features = pseudo_genes_of_interest) #do not filter until later

### Get and check metadata dfs

pseudo_metadf <- as.data.frame(Seu_test_object@meta.data)
pseudo_metadf[, c("projid", "age_death", "msex", "braaksc", "ceradsc", "cogdx")] # see if there are these columns

#if the above failed
Projid_info <- readRDS("~/collabgit/AD_snRNAseq/data/ROSmaster.rds") #load ROS metadata
Projid_info <- Projid_info[,c("projid", "pathoAD", "gpath", "age_death", "msex", "braaksc", "ceradsc", "cogdx")]
pseudo_metadf <- merge(pseudo_metadf, Projid_info, by = "projid")

#set case vs controls
pseudo_metadf %<>% mutate(LOAD = case_when((braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                           (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                           TRUE ~ 'OTHER')) 

#for mathys only - MAYBE
#ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset
#ros_meta_small = ros_meta_small %>% dplyr::mutate(pathoAD = factor(pathoAD), msex = factor(msex))

#reset rownames (lost from merge)
rownames(pseudo_metadf) <- pseudo_metadf$TAG #for mathys

#set factors of interest
pseudo_metadf$subclass <- factor(pseudo_metadf$predicted.id)
pseudo_metadf$projid <- factor(pseudo_metadf$projid)
pseudo_metadf$LOAD <- factor(pseudo_metadf$LOAD)

#get counts df
pseudo_counts <- Seu_test_object@assays$RNA@counts

#no longer need Seurat object
remove(Seu_test_object)

#### Further data processing with SCE ####

#reoder metadata based on order of cells in count
pseudo_metadf <- pseudo_metadf[colnames(pseudo_counts),]
identical(colnames(pseudo_counts), rownames(pseudo_metadf)) #need true

#create SCE object
sce <- SingleCellExperiment(assays = list(counts = pseudo_counts), 
                            colData = pseudo_metadf)

### Explore SCE object

assays(sce) # Check the assays present
dim(counts(sce)) # Explore the raw counts for the dataset
counts(sce)[1:6, 1:6] # see some raw counts

dim(colData(sce)) # Explore the cellular metadata for the dataset
head(colData(sce)) # See some columns

### Aggregate data

# Get named vector of cluster names
kids <- purrr::set_names(levels(sce$subclass))
kids

# Get total number of clusters
nk <- length(kids)
nk

# Get named vector of sample names
sids <- purrr::set_names(levels(sce$projid))

# Get total number of samples 
ns <- length(sids)
ns

## Generate sample level metadata

# Determine the number of cells per sample
table(sce$projid)

# Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$projid))

# Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$projid)

# Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"subclass")
ei

## QC: Remove outlier cells

dim(sce)

# Calculate quality control (QC) metrics
sce <- calculateQCMetrics(sce)

# Get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(
  metric = sce$total_features_by_counts,
  nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
sce <- sce[, !sce$is_outlier] # nothing removed for mathys or cain
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ] # some genes removed for Mathys (~4k) and Cain (~7.3k)
dim(sce)

#### Aggregate data to sample level ####

# Subset metadata to only include the cluster and sample IDs to aggregate across
pseudo_groups <- colData(sce)[, c("subclass", "projid")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = pseudo_groups, fun = "sum") 

class(pb)
dim(pb)
pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                                    `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
                       lapply(function(u) 
                       set_colnames(t(u), 
                       stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb)

# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$subclass, sce$projid)

#### DESeq2 prep and QC ####

### format and subset data
# Get sample names for each of the cell type clusters

# prep. data.frame for plotting

get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(subclass = de_cluster_ids,
                    projid = de_samples)

gg_df <- left_join(gg_df, ei[, c("projid", "LOAD")]) 


metadata <- gg_df %>%
  dplyr::select(subclass, projid, LOAD) 

metadata     

# Generate vector of cluster IDs
clusters <- levels(metadata$subclass)
clusters

# Subset the metadata to only a select cell type
cluster_ind <- match("SST",clusters) #as example
cluster_metadata <- metadata[which(metadata$subclass == clusters[cluster_ind]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$projid
head(cluster_metadata)

# Subset the counts to only the given
counts <- pb[[clusters[cluster_ind]]]

cluster_counts <- counts[, which(colnames(counts) %in% rownames(cluster_metadata))]
cluster_counts <- as.data.frame(cluster_counts)

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))    

cluster_metadata = merge(cluster_metadata,  ei %>% select(projid, age_death, msex))
cluster_metadata$msex = factor(cluster_metadata$msex)
cluster_metadata$LOAD = factor(cluster_metadata$LOAD)
cluster_metadata[cluster_metadata$age_death == "90+","age_death"] <- 90
cluster_metadata$age_death <- as.numeric(cluster_metadata$age_death)
cluster_metadata$age_death = scale(cluster_metadata$age_death)

#cluster_metadata$pathoAD = factor(cluster_metadata$pathoAD)
rownames(cluster_metadata) <- cluster_metadata$projid

### Create object and QC

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ LOAD + msex + age_death)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "LOAD")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap::pheatmap(rld_cor, annotation = cluster_metadata[, c("LOAD"), drop=F])

#### DESeq2 run ####

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Output results of Wald test for contrast for AD vs control
levels(cluster_metadata$LOAD)[1]
levels(cluster_metadata$LOAD)[2]

contrast <- c("LOAD", levels(cluster_metadata$LOAD)[1], levels(cluster_metadata$LOAD)[2])
#contrast <- c("msex", levels(cluster_metadata$msex)[2], levels(cluster_metadata$msex)[1])


# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 res=res)

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()

# Check results output
res_tbl

# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

### ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- res_tbl %>%
  dplyr::arrange(pvalue) %>%
  dplyr::pull(gene) %>%
  head(n=20)

#top20_sig_genes = factor(c('GOLT1B', 'ATF6B', 'DDRGK1', 'TUBB2A', 'BEX2', 'ATPIF1', 'RASGEF1B', 'NGFRAP1', 'LINGO1', 'NTNG1'))

sst_marker_genes = new_MTGnCgG_lfct2_results %>% dplyr::filter(subclass == 'SST') %>% pull(gene)

top20_sig_norm <- data.frame(normalized_counts) %>%
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sst_marker_genes)

gathered_top20_sig <- top20_sig_norm %>%
  tidyr::gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

ei$samplename = ei$projid %>% make.names()
gathered_top20_sig <- inner_join(ei[, c("samplename", "LOAD" )], gathered_top20_sig, by = "samplename")

## plot using ggplot2
ggplot(gathered_top20_sig, aes(x = gene, 
                               y = normalized_counts, 
                               fill = LOAD)) +
  geom_boxplot() + 
  # geom_point(aes(x = gene, 
  #                y = normalized_counts, 
  #                color = pathoAD), 
  #            position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("SST Markers") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_tbl) +
  geom_point(aes(x = log2FoldChange, y = -log10(pvalue))) +
  ggtitle("Volcano plot of SST cells relative to control") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


#hahahah


#### Old pseudobulk analysis attempt by Sonny ####

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
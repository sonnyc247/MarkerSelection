## assign best subclasses from Allen Institute onto mathys

## ideally this would include some expr normalization step, like Seurat or cpm or something

# first load expression counts from mathys
library(Matrix)
library(magrittr)
library(tidyverse)

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

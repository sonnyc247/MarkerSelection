#### Subclass ####

library(tidyverse)

## read in and update sonny's markers with as many ensembl ids as possible

all_markers_ranked_orig = read_csv('/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_ITexpand_WL35IT_lfct11_minpct15_dup.csv')

# i noticed that some of sonny's markers are missing ensembl ids, so i'm getting those from the raw data below

# read in gene info off the SCC
aibs_gene_info = read_csv('/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_MTG_Smartseq_2018/human_MTG_2018-06-14_genes-rows.csv')
aibs_gene_info = aibs_gene_info %>% mutate(gene_name_in_r = make.names(gene))

# this adds back in entrez ids for all genes (assumes names haven't changed)
all_markers_ranked_w_entrez = left_join(all_markers_ranked_orig %>% dplyr::select(-entrez_id, -ensembl_id), 
                               aibs_gene_info %>% dplyr::select(-gene_name, -mouse_homologenes, -chromosome), by = c("gene" = "gene_name_in_r")) %>% 
  dplyr::select(-gene.y)

# next, i updated the ensembl ids
# weirdly, neither biomart nor annotation dbi seemed to have all of the gene symbols i was looking for, so i just got the stuff from hgnc's archives

# get hgnc gene mapping from website
hgnc_mapping = read_tsv(url('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'))
#hgnc_mapping = read_tsv('/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Inputs/hgnc_complete_set.txt')


# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
all_markers_ranked_w_ids = left_join(all_markers_ranked_w_entrez, 
                                     hgnc_mapping %>% dplyr::select(entrez_id, ensembl_gene_id) %>% dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(X1, gene, entrez_id, ensembl_id, everything())


## now, read in the seurat object for the all regions hodge data and then add the extra columns
Seu_AIBS_obj = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds')


# i decided to only calculate the extra fields for things that were calcualted after findallmakrers
use_markers = all_markers_ranked_w_ids %>% pull(gene)

gene_List <- rownames(Seu_AIBS_obj[["RNA"]]@data) %>% unlist

matching_inds = which(gene_List %in% use_markers)
use_genes = gene_List[matching_inds]



AIBS_lnCPM_InEx_ACC <- t(as.matrix(Seu_AIBS_obj[["RNA"]]@data[matching_inds, ])) %>% as.data.frame() #get transposed lnCPM matrix
# format count/expression matrix for group_dot_plot smooth running
library(tibble)
AIBS_lnCPM_InEx_ACC <- rownames_to_column(AIBS_lnCPM_InEx_ACC, var = 'sample_name') #get sample names as a column
colnames(AIBS_lnCPM_InEx_ACC)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(AIBS_lnCPM_InEx_ACC) <- AIBS_lnCPM_InEx_ACC$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot
gdp_anno$sample_name = row.names(gdp_anno)
gdp_anno$subclass_label_expanded_L35IT_label = gdp_anno$subclass_label_expanded_L35IT
gdp_anno$subclass_label_expanded_L35IT_id = gdp_anno$subclass_label_expanded_L35IT
gdp_anno$subclass_label_expanded_L35IT_color = 'grey'

new = full_join(AIBS_lnCPM_InEx_ACC, gdp_anno %>% dplyr::select(sample_name, subclass_label_expanded_L35IT))

# calculates averages of gene expression from every cell type using a trimmed mean - takes a while to run
cell_type_avgs = new %>% dplyr::select(-sample_name) %>% 
  group_by(subclass_label_expanded_L35IT) %>% 
  summarise_all(mean)

# this is a long data frame that merges subclass info with the cell type gene averages computed above
cell_type_avgs_long = cell_type_avgs %>% 
  pivot_longer(-subclass_label_expanded_L35IT, 
               names_to = 'gene', values_to = 'avg_expr')

# now we go through and calculate the top_log2FC feature, 
# we do this in two steps, we first calculate it for the condition where the subclass is most strongly expressing the gene

top_genes_per_type = cell_type_avgs_long %>% 
  group_by(gene) %>% 
  top_n(n = 2, wt = avg_expr) %>% # this gets the top two highest expressing subclasses expressing the gene
  arrange(-avg_expr) %>% # sorts them
  mutate(top_log2FC = log2(avg_expr[1]) - log2(avg_expr[2])) %>% # and then calculates top_log2FC
  top_n(n = 1, wt = avg_expr) %>% # gets the top highest expressing subclass
  dplyr::rename(subclass = subclass_label_expanded_L35IT)



# now this goes through all of the other subclasses to find their top_log2FC if they're not the highest expressing gene per subclass
df1 = top_genes_per_type %>% 
                  dplyr::select(gene, avg_expr) %>% 
                  dplyr::rename(top_avg_expr = avg_expr)
top_markers_second_best = inner_join(df1, cell_type_avgs_long %>% 
                                       dplyr::rename(subclass = subclass_label_expanded_L35IT )) %>% 
   mutate(top_log2FC = log2(avg_expr) - log2(top_avg_expr)) %>% dplyr::select(-top_avg_expr) %>% 
  filter(!top_log2FC == 0) # this last part removes subclasses for genes that are the highest expressing subclass per gene, as they're calculated in the first part
  
# this creates a new data frame which has the extra columns we wanted to calculate
extra_cols_df = bind_rows(top_genes_per_type, top_markers_second_best) %>% dplyr::select(gene, subclass, avg_expr, top_log2FC) %>% arrange(gene)
  
# this now is the final unranked data frame
final_markers_df_unranked = left_join(all_markers_ranked_w_ids, extra_cols_df)  %>% 
  dplyr::distinct(gene, ensembl_id, subclass, .keep_all = T)


# this generates the final list of markers, creates a few columns which are slightly different ways of ranking the genes
# shreejoy's favorite is the column bretigea_ranking_best, which is like the bretigea ranking plus adding a term to rank by descending pct.2
final_markers_df_ranked = final_markers_df_unranked %>% 
  dplyr::select(X1, gene, entrez_id, ensembl_id, everything()) %>% 
  group_by(subclass) %>% 
  mutate(bretigea_ranking = rank(desc((rank(avg_log2FC) + rank(top_log2FC) + rank(avg_expr))))) %>% # this is the implementation of ranking in the bretigea paper
  mutate(dan_ranking = rank(desc((rank(avg_log2FC) + rank(top_log2FC) )))  ) %>% # this is dan's suggestion - this strongly empahsizing specificity
  mutate(bretigea_ranking_best = rank(desc((rank(avg_log2FC) + rank(top_log2FC) + rank(avg_expr) + rank(desc(pct.2)))))  ) %>% # this seems to be a good balance of specificity and sensitivity
  ungroup(subclass) %>%
  arrange(subclass, bretigea_ranking_best) 

# this writes the final ranked markers to a file
write_csv(final_markers_df_ranked, file = '/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/Ranked_markers_noTrim_ALLReg_ITexpand_WL35IT_lfct11_minpct15_dup.csv')


### visualize a few cell type's markers using group dot plots

gdp_markers = final_markers_df_ranked %>% filter(subclass == 'L4 IT', !is.na(ensembl_id)) %>% distinct(gene, .keep_all = T) %>%
  slice_max(n = 20, order_by = -bretigea_ranking_best) %>% 
  pull(gene)

# devtools::install_github("AllenInstitute/scrattch.vis") # as needed
library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis



# do the plot
group_dot_plot(AIBS_lnCPM_InEx_ACC, 
               gdp_anno, 
               genes = gdp_markers,
               grouping = "subclass_label_expanded_L35IT", 
               log_scale = TRUE,
               font_size = 10,
               max_size = 20,
               rotate_counts = TRUE, 
               fill_stat = "tmean")

#

#### Class ####

library(tidyverse)

## read in and update sonny's markers with as many ensembl ids as possible

all_markers_ranked_orig = read_csv('/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/new_ALLReg_results_class_ITexpand_WL35IT_lfct11_minpct15_dup.csv')


# i noticed that some of sonny's markers are missing ensembl ids, so i'm getting those from the raw data below

# read in gene info off the SCC
aibs_gene_info = read_csv('/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_MTG_Smartseq_2018/human_MTG_2018-06-14_genes-rows.csv')
aibs_gene_info = aibs_gene_info %>% mutate(gene_name_in_r = make.names(gene))

# this adds back in entrez ids for all genes (assumes names haven't changed)
all_markers_ranked_w_entrez = left_join(all_markers_ranked_orig %>% dplyr::select(-entrez_id, -ensembl_id), 
                                        aibs_gene_info %>% dplyr::select(-gene_name, -mouse_homologenes, -chromosome), by = c("gene" = "gene_name_in_r")) %>% 
  dplyr::select(-gene.y)

# next, i updated the ensembl ids
# weirdly, neither biomart nor annotation dbi seemed to have all of the gene symbols i was looking for, so i just got the stuff from hgnc's archives

# get hgnc gene mapping from website
hgnc_mapping = read_tsv(url('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'))
#hgnc_mapping = read_tsv('/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Inputs/hgnc_complete_set.txt')


# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
all_markers_ranked_w_ids = left_join(all_markers_ranked_w_entrez, 
                                     hgnc_mapping %>% dplyr::select(entrez_id, ensembl_gene_id) %>% dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(X1, gene, entrez_id, ensembl_id, everything())


## now, read in the seurat object for the all regions hodge data and then add the extra columns
Seu_AIBS_obj = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds')


# i decided to only calculate the extra fields for things that were calcualted after findallmakrers
use_markers = all_markers_ranked_w_ids %>% pull(gene)

gene_List <- rownames(Seu_AIBS_obj[["RNA"]]@data) %>% unlist

matching_inds = which(gene_List %in% use_markers)
use_genes = gene_List[matching_inds]

AIBS_lnCPM_InEx_ACC <- t(as.matrix(Seu_AIBS_obj[["RNA"]]@data[matching_inds, ])) %>% as.data.frame() #get transposed lnCPM matrix
# format count/expression matrix for group_dot_plot smooth running
library(tibble)
AIBS_lnCPM_InEx_ACC <- rownames_to_column(AIBS_lnCPM_InEx_ACC, var = 'sample_name') #get sample names as a column
colnames(AIBS_lnCPM_InEx_ACC)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(AIBS_lnCPM_InEx_ACC) <- AIBS_lnCPM_InEx_ACC$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot
gdp_anno$sample_name = row.names(gdp_anno)
gdp_anno$subclass_label_expanded_L35IT_label = gdp_anno$subclass_label_expanded_L35IT
gdp_anno$subclass_label_expanded_L35IT_id = gdp_anno$subclass_label_expanded_L35IT
gdp_anno$subclass_label_expanded_L35IT_color = 'grey'

new = full_join(AIBS_lnCPM_InEx_ACC, gdp_anno %>% dplyr::select(sample_name, class_label))

# calculates averages of gene expression from every cell type using a trimmed mean - takes a while to run
cell_type_avgs = new %>% dplyr::select(-sample_name) %>% 
  group_by(class_label) %>% 
  summarise_all(mean)

# this is a long data frame that merges subclass info with the cell type gene averages computed above
cell_type_avgs_long = cell_type_avgs %>% 
  pivot_longer(-class_label, 
               names_to = 'gene', values_to = 'avg_expr')

# now we go through and calculate the top_log2FC feature, 
# we do this in two steps, we first calculate it for the condition where the subclass is most strongly expressing the gene

top_genes_per_type = cell_type_avgs_long %>% 
  group_by(gene) %>% 
  top_n(n = 2, wt = avg_expr) %>% # this gets the top two highest expressing subclasses expressing the gene
  arrange(-avg_expr) %>% # sorts them
  mutate(top_log2FC = log2(avg_expr[1]) - log2(avg_expr[2])) %>% # and then calculates top_log2FC
  top_n(n = 1, wt = avg_expr) %>% # gets the top highest expressing subclass
  dplyr::rename(class = class_label)



# now this goes through all of the other subclasses to find their top_log2FC if they're not the highest expressing gene per subclass
df1 = top_genes_per_type %>% 
  dplyr::select(gene, avg_expr) %>% 
  dplyr::rename(top_avg_expr = avg_expr)
top_markers_second_best = inner_join(df1, cell_type_avgs_long %>% 
                                       dplyr::rename(class = class_label)) %>% 
  mutate(top_log2FC = log2(avg_expr) - log2(top_avg_expr)) %>% dplyr::select(-top_avg_expr) %>% 
  filter(!top_log2FC == 0) # this last part removes subclasses for genes that are the highest expressing subclass per gene, as they're calculated in the first part

# this creates a new data frame which has the extra columns we wanted to calculate
extra_cols_df = bind_rows(top_genes_per_type, top_markers_second_best) %>% dplyr::select(gene, class, avg_expr, top_log2FC) %>% arrange(gene)

# this now is the final unranked data frame
final_markers_df_unranked = left_join(all_markers_ranked_w_ids, extra_cols_df)  %>% 
  dplyr::distinct(gene, ensembl_id, class, .keep_all = T)


# this generates the final list of markers, creates a few columns which are slightly different ways of ranking the genes
# shreejoy's favorite is the column bretigea_ranking_best, which is like the bretigea ranking plus adding a term to rank by descending pct.2
final_markers_df_ranked = final_markers_df_unranked %>% 
  dplyr::select(X1, gene, entrez_id, ensembl_id, everything()) %>% 
  group_by(class) %>% 
  mutate(bretigea_ranking = rank(desc((rank(avg_log2FC) + rank(top_log2FC) + rank(avg_expr))))) %>% # this is the implementation of ranking in the bretigea paper
  mutate(dan_ranking = rank(desc((rank(avg_log2FC) + rank(top_log2FC) )))  ) %>% # this is dan's suggestion - this strongly empahsizing specificity
  mutate(bretigea_ranking_best = rank(desc((rank(avg_log2FC) + rank(top_log2FC) + rank(avg_expr) + rank(desc(pct.2)))))  ) %>% # this seems to be a good balance of specificity and sensitivity
  ungroup(class) %>%
  arrange(class, bretigea_ranking_best) 

# this writes the final ranked markers to a file
write_csv(final_markers_df_ranked, file = '/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/All_hodge_regions/Ranked_class_markers_noTrim_ALLReg_ITexpand_WL35IT_lfct11_minpct15_dup.csv')


### visualize a few cell type's markers using group dot plots

gdp_markers = final_markers_df_ranked %>% filter(class == "Non-neuronal", !is.na(ensembl_id)) %>% distinct(gene, .keep_all = T) %>%
  slice_max(n = 20, order_by = -bretigea_ranking_best) %>% 
  pull(gene)

# devtools::install_github("AllenInstitute/scrattch.vis") # as needed
library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis



# do the plot
gdp_anno$class_id <- gdp_anno$class_label 

group_dot_plot(AIBS_lnCPM_InEx_ACC, 
               gdp_anno, 
               genes = gdp_markers,
               grouping = "class", 
               log_scale = TRUE,
               font_size = 10,
               max_size = 20,
               rotate_counts = TRUE, 
               fill_stat = "tmean")




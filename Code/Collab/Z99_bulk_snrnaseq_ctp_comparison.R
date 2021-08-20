library(tidyverse)
library(markerGeneProfile)
library(cowplot)
library(broom)
library(broom.mixed)
library(lme4)

## this is where i've estimated bulk rCTPs  - update this with your internal bulk MGPs from rosmap

estimations <-  mgpEstimate(
  exprData=ros_gene_mat,
  genes=new_marker_list_ensembl,
  geneColName='ensembl_id',
  outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
  geneTransform = NULL, # this is the default option for geneTransform
  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus=FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority=FALSE) 

# this turns the wide dataframe into a long dataframe with columns for projid, subclass, rCTP
mgp_df_big = estimations$estimates %>% as.data.frame() %>% tibble::rownames_to_column(var = "projid") %>% mutate(projid = as.numeric(projid)) %>% 
  pivot_longer(-projid, names_to = "subclass", values_to = 'rCTP')

# create a new column called dataset and label it bulk
mgp_df_big$dataset = 'bulk'

# create a new column called cell_type_proportion which is same as rCTP
mgp_df_big = mgp_df_big %>% mutate(cell_type_proportion = rCTP)


# ## get the snCTP data from sonny's repo
zhou_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/Proportions/MTGnCgG/MTGnCgG_20Pct/Filtered/zhou_cell_type_prop_df.csv'
mathys_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/Proportions/MTGnCgG/MTGnCgG_20Pct/Filtered/mathys_cell_type_prop_df.csv'
cain_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/Proportions/MTGnCgG/MTGnCgG_20Pct/Filtered/cain_cell_type_prop_df.csv'
# 
# 
# # # read in data frames from urls
zhou_ad = read_csv(url(zhou_ad_url))
mathys_ad = read_csv(url(mathys_ad_url))
cain_ad = read_csv(url(cain_ad_url))

# annotate datasets
zhou_ad$dataset = 'Zhou'
mathys_ad$dataset = 'Mathys'
cain_ad$dataset = 'Cain'
# need to remove age_death because they're weirdly characters 
zhou_ad %<>% dplyr::select(-age_death)
cain_ad %<>% dplyr::select(-age_death)


# bind all data frames together
ad_snrnaseq_df = bind_rows(mathys_ad, zhou_ad %>% mutate(projid = as.numeric(projid)), cain_ad, mgp_df_big)

ad_snrnaseq_df = ad_snrnaseq_df %>% select(-age_death, -msex, -pmi)

# join snCTP data frame with rosmaster (you need to have your own version of rosmaster for this - this just gets the metadata from rosmaster)
ad_snrnaseq_df_merged = right_join(rosmaster, ad_snrnaseq_df, by = 'projid')
ad_snrnaseq_df_merged$subclass = factor(ad_snrnaseq_df_merged$subclass %>% make.names())

## redefine AD cases and Controls using the same definitions from the resequencing paper

# based on definitions here: https://github.com/th1vairam/ampad-DiffExp/blob/bd9766224c2d1515586c9377db7e08a6cb62bcc9/gene_level_analysis/ROSMAP_geneLevel.Rmd#L157-L159
ad_snrnaseq_df_merged %<>% mutate(LOAD = case_when( (braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                                    (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                                    TRUE ~ 'OTHER')
) 

ad_snrnaseq_df_merged$dataset = factor(ad_snrnaseq_df_merged$dataset, levels = c('Mathys', 'Zhou', 'Cain', 'bulk'))
ad_snrnaseq_df_merged$LOAD = factor(ad_snrnaseq_df_merged$LOAD, levels = c('C', 'AD'))

# this finds duplicated subjects between bulk and snRNAseq
dup_subjects = ad_snrnaseq_df_merged %>% distinct(projid, dataset, .keep_all = T) %>% group_by(projid) %>% tally() %>% filter(n > 1) %>% pull(projid)

# get snCTPs for subjects that also have bulk data
df_wide = ad_snrnaseq_df_merged %>% filter(projid %in% dup_subjects, !dataset == 'bulk') %>%
  dplyr::select(projid, subclass, cell_type_proportion, dataset) %>% 
  pivot_wider(id_cols = c(projid, subclass), names_from = dataset, values_from = cell_type_proportion)

# create data frame with ctps for all datasets for projids with duplicated data
rCTP_snCTP_merged_df = full_join(df_wide, ad_snrnaseq_df_merged %>% filter(dataset == 'bulk', projid %in% dup_subjects) %>% 
                                   dplyr::select(projid, subclass, cell_type_proportion) %>% dplyr::rename(bulk = cell_type_proportion))
                                 
# estimate dataset by cell type correlation using paired bulk and snRNAseq samples
bulk_snrnaseq_corrs = rCTP_snCTP_merged_df %>% 
  pivot_longer(cols = Mathys:Cain, names_to = 'dataset', values_to = 'cell_type_proportion') %>% 
  group_by(dataset, subclass) %>% mutate(cell_type_cor = cor(bulk, cell_type_proportion, method = "spearman", use = "pairwise.complete.obs")) %>% 
  distinct(dataset, subclass, cell_type_cor)  %>% pivot_wider(names_from = dataset, values_from = cell_type_cor) %>% as.data.frame()

bulk_snrnaseq_corrs_long = bulk_snrnaseq_corrs %>% pivot_longer(cols = Mathys:Cain, names_to = 'dataset', values_to = 'bulk_sn_cor') 
bulk_snrnaseq_corrs_long$dataset = factor(bulk_snrnaseq_corrs_long$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
bulk_snrnaseq_corrs_long$subclass = factor(bulk_snrnaseq_corrs_long$subclass, levels = sort(make.names(new_cell_types)))


bulk_sn_cor_summary_plot = bulk_snrnaseq_corrs_long %>% full_join(.,beta_coefs_non_meta_df %>% dplyr::select(subclass, class)) %>% distinct() %>%
  ggplot(aes(x = subclass, y = bulk_sn_cor, fill = dataset, label = bulk_sn_cor)) + 
  geom_bar(stat = 'identity', position = "dodge") + 
  geom_text(aes(label=sprintf("%0.2f", bulk_sn_cor)), position = position_dodge(0.9), vjust = -.5) + 
  geom_hline(yintercept = 0) + 
  facet_grid(~class, scales = 'free_x', space = "free") + 
  ylab('Corr. bulk to sn proportions (Spearman)') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(c(-0.30, 0.75))


# create individual plots showing SST and IT associations
rCTP_snCTP_merged_df_long = left_join(rCTP_snCTP_merged_df %>% 
  pivot_longer(cols = Mathys:Cain, names_to = 'dataset', values_to = 'cell_type_proportion'), 
  ad_snrnaseq_df_merged %>% dplyr::select(projid, LOAD))
rCTP_snCTP_merged_df_long$dataset = factor(rCTP_snCTP_merged_df_long$dataset, levels = c('Mathys', 'Zhou', 'Cain'))

# sst and it plot
sst_it_sn_bulk_indiv_plot = rCTP_snCTP_merged_df_long %>% 
  filter(subclass %in% c("SST", "IT")) %>% ggplot(aes(x = bulk, y = cell_type_proportion, color = LOAD, group = 1)) + 
  geom_smooth(method = 'lm', se = F) + 
  geom_point() + facet_grid(rows = vars(subclass), cols = vars(dataset), scales = "free") + xlab('bulk proportion') + ylab('single nuc proportion') + 
  theme_cowplot()

# combined plot grid object for Supp Fig 4
plot_grid(bulk_sn_cor_summary_plot, sst_it_sn_bulk_indiv_plot, nrow = 2, rel_heights = c(2, 1.5))

library(tidyverse)
library(magrittr)
library(lme4)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(broom)
library(ggbeeswarm)

## load in the most up to date version of rosmaster - you'll have to change this
rosmaster = readRDS('data/ROSmaster.rds')
rosmaster$projid = rosmaster$projid %>% as.numeric()

## get sonny's subclass colors and naming
#subclass_meta_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Misc/Name_and_colour_scheme.csv'
#subclass_meta = read_csv(url(subclass_meta_url))
subclass_meta = read_csv('data/subclass_meta_info.csv')

## get the snCTP data from sonny's repo
zhou_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/zhou_cell_type_prop_df.csv'
mathys_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/mathys_cell_type_prop_df.csv'
cain_ad_url = 'https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/cain_cell_type_prop_df_noNA.csv'

# read in data frames from urls
zhou_ad = read_csv(url(zhou_ad_url))
mathys_ad = read_csv(url(mathys_ad_url))
cain_ad = read_csv(url(cain_ad_url))

# annotate datasets
zhou_ad$dataset = 'Zhou'
mathys_ad$dataset = 'Mathys'
cain_ad$dataset = 'Cain'
# need to remove age_death because they're weirdly characters 
zhou_ad %<>% select(-age_death)
cain_ad %<>% select(-age_death)

# bind all data frames together
ad_snrnaseq_df = bind_rows(mathys_ad, zhou_ad, cain_ad)

# join snCTP data frame with rosmaster
ad_snrnaseq_df$projid <- as.character(ad_snrnaseq_df$projid)
ad_snrnaseq_df_merged = right_join(rosmaster, ad_snrnaseq_df, by = 'projid')


## redefine AD cases and Controls using the same definitions from the resequencing paper

# based on definitions here: https://github.com/th1vairam/ampad-DiffExp/blob/bd9766224c2d1515586c9377db7e08a6cb62bcc9/gene_level_analysis/ROSMAP_geneLevel.Rmd#L157-L159
ad_snrnaseq_df_merged %<>% mutate(LOAD = case_when( (braaksc.x >= 4 & ceradsc.x <= 2 & cogdx.x == 4) ~ 'AD',
                                                   (braaksc.x <= 3 & ceradsc.x >= 3 & cogdx.x == 1) ~ 'C',
                                                   TRUE ~ 'OTHER')
                                                  ) 

cell_subclass_order = subclass_meta %>% pull(subclass)

# count up the number of folks meeting the case/Control criteria
ad_snrnaseq_df_merged %>% select(projid, dataset, gpath.y, pathoAD.x, braaksc.x, ceradsc.x, cogdx.x, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

ad_snrnaseq_df_merged$LOAD = factor(ad_snrnaseq_df_merged$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df_merged$dataset = factor(ad_snrnaseq_df_merged$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
ad_snrnaseq_df_merged$subclass = factor(ad_snrnaseq_df_merged$subclass, levels = cell_subclass_order)

# generate some plots showing that we didn't fuck up
sst_boxplots = ad_snrnaseq_df_merged %>% 
  filter(subclass == 'SST', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

# IT don't match bulk. oh well.
it_boxplots = ad_snrnaseq_df_merged %>% filter(subclass == 'IT', LOAD %in% c('C', 'AD')) %>% 
  distinct(projid, .keep_all = T) %>% ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('IT snCTP (%)') + 
  xlab('')

## calculate LOAD beta coefficients for snCTPs per dataset and each cell subclass

# this defines the list of unique cell types
cell_type_list = cell_subclass_order

# dropping L5 ET and VLMC because they're not present in each dataset
cell_type_list = cell_type_list[! cell_type_list %in% c('L5 ET', 'VLMC')] 

# this goes through each dataset and cell type 
# and fits a linear model and a beta coeff for LOAD with age, sex, and pmi as covariates
beta_coefs_non_meta_df = lapply(levels(ad_snrnaseq_df_merged$dataset), function(curr_dataset){
  print(curr_dataset)
  
    dataset_df = lapply(cell_type_list, function(curr_cell_type){
      print(curr_cell_type)
      df = ad_snrnaseq_df_merged %>% filter(LOAD %in% c('C', 'AD'))
      df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
      my_model = lm('scale(cell_type_proportion) ~ scale(age_death.x) + factor(msex.x) + scale(pmi) + LOAD ',
           data = df)
      model_df = tidy(my_model)
      model_df$dataset = curr_dataset
      model_df$subclass = curr_cell_type
      
      return(model_df) 
  }) %>% bind_rows()
    return(dataset_df)  
  #beta_out = 
}) %>% bind_rows()

# sets factor levels to match order we want
beta_coefs_non_meta_df$dataset = factor(beta_coefs_non_meta_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))

beta_coefs_non_meta_df = merge(beta_coefs_non_meta_df, subclass_meta, by.x = 'subclass', by.y = 'subclass')


# plots beta coefficients for LOAD across each dataset faceted by cell type
all_dataset_subclass_betas_plot = beta_coefs_non_meta_df %>% filter(term == 'LOADAD') %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_wrap(~class*subclass) + 
  ylab('LOAD (std. Beta)') + 
  xlab('')

# generate individual bar charts for sst and it - note that bars are 95% CIs, not standard errors
sst_dataset_betas_plot =  beta_coefs_non_meta_df %>% filter(term == 'LOADAD', subclass == 'SST') %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error) , width = .33) + 
  ylab('LOAD (std. Beta)') + 
  xlab('')

it_dataset_betas_plot =  beta_coefs_non_meta_df %>% filter(term == 'LOADAD', subclass == 'IT') %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error) , width = .33) + 
  ylab('LOAD (std. Beta)') + 
  xlab('')


## go through and perform the mega-analysis across datasets using a mixed effects model
cell_type_list = cell_subclass_order
cell_type_list = cell_type_list[! cell_type_list %in% c('L5 ET', 'VLMC')]

beta_coefs_meta_df = lapply(cell_type_list, function(curr_cell_type){
  print(curr_cell_type)
  
  # this inner loop just gets the dataset together
  dataset_df = lapply(levels(ad_snrnaseq_df_merged$dataset), function(curr_dataset){
    print(curr_dataset)
    df = ad_snrnaseq_df_merged %>% filter(LOAD %in% c('C', 'AD'))
    df$LOAD = factor(df$LOAD, levels = c('C', 'AD'))
    df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
    
    df$dataset = curr_dataset

    return(df) 
  }) %>% bind_rows()
  
  # this estimates the mixed effects model, one with LOAD and one without
  # using random effects for dataset and projid as there's a couple projid duplicates across the datasets (49 unique projids)
  
  model_w_load = lmer('scale(cell_type_proportion) ~ scale(age_death.x) + factor(msex.x) + scale(pmi.x) + LOAD + (1|dataset) + (1|projid)', 
            data = dataset_df )
  
  model_no_load = lmer('scale(cell_type_proportion) ~ scale(age_death.x) + factor(msex.x) + scale(pmi.x) + (1|dataset)  + (1|projid) ', 
            data = dataset_df)
  
  anova_ob = anova(model_w_load, model_no_load)
  pval_load = anova_ob$`Pr(>Chisq)`[2]
  
  lmer_df = tidy(model_w_load)
  lmer_df$subclass = curr_cell_type
  lmer_df$pval_load = pval_load
  return(lmer_df)  
  
}) %>% bind_rows()

# estimate FDRs and bonferroni pvalues
beta_coefs_load_meta_df = beta_coefs_meta_df %>% filter(term == 'LOADAD') 

beta_coefs_load_meta_df = beta_coefs_load_meta_df %>% mutate(fdr_load = p.adjust(pval_load, method = 'BH'),
                                                             bonf_p_load = p.adjust(pval_load, method = 'bonferroni'))

# merge data frame with subclass metadata data frame
beta_coefs_load_meta_df = merge(beta_coefs_load_meta_df, 
                           subclass_meta, by.x = 'subclass', by.y = 'subclass')
beta_coefs_load_meta_df$subclass = factor(beta_coefs_load_meta_df$subclass, levels = cell_type_list)



## generate a plot with the mega-analysis LOAD associations faceted by major cell class

mega_analysis_plot = beta_coefs_load_meta_df %>% 
  ggplot(aes(x = subclass, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error * 1.96, ymax = estimate + std.error * 1.96), width = .33) + 
  ylab('LOAD (std. Beta)') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  facet_grid(~class, scale = 'free_x', space = 'free_x') +
  xlab('')


## generate the final multi-panel figure using cowplot's plot_grid
top_plot = plot_grid(sst_boxplots, sst_dataset_betas_plot, 
                     it_boxplots, it_dataset_betas_plot, nrow = 2, ncol = 2, rel_widths = c(.6, .4), 
                     axis = 'l', align = 'v', labels = c('A', 'B', 'C', 'D'))

# new_plot = plot_grid(top_plot, bulk_vs_single_nuc_ctp_plot, nrow = 1, ncol = 2, rel_widths = c(1, .6), axis = 'b', align = 'v')

full_plot = plot_grid(top_plot, mega_analysis_plot, nrow = 2 , rel_heights = c(2, 1.5), 
                       axis = 'l', align = 'h', labels = c('', 'E'))


full_plot

ggsave('figures/snRNAseq_main_fig.png', device = 'png', plot = full_plot, width = 5, height = 7)

## write out csvs that will be used for supplementary tables showing the results of the beta coeffs, pvals, and FDRs for the 
# dataset level analyses and mega-analyses

# put together data frames with summary data
snCTP_dataset_level_stats = beta_coefs_non_meta_df %>% filter(term == 'LOADAD') %>% 
  dplyr::select(subclass, class, term, dataset, estimate:p.value) %>%
  arrange(class, subclass)

snCTP_mega_level_stats = beta_coefs_load_meta_df %>% filter(term == 'LOADAD') %>%   
  dplyr::select(subclass, class, term, estimate:bonf_p_load, -group) %>%
  arrange(class, subclass)

write_csv(x = snCTP_dataset_level_stats, path = "outputs/snCTP_dataset_level_stats.csv")
write_csv(x = snCTP_mega_level_stats, path = "outputs/snCTP_mega_level_stats.csv")



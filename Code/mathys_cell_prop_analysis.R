### from Shreejoy

library(Seurat)
library(tidyverse)
library(broom)

# readme - this script reads in metadata anootations from Mathys (after re-annotation by Sonny)
# cell proportions and then performs testing for which cell types are different in cases than controls

# script should work from anywhere on SCC

# define path to data -> this is to Sonny's reannotation of Mathys snRNAseq data using Hodge snRNAseq subclass labels
mathys_seurat_ob_path = '/external/rprshnas01/netdata_kcni/stlab/marker_genes/Seu_mathys_obj.rds'

# read datasets into environment
Seu_mathys_obj = readRDS(mathys_seurat_ob_path)

# pull mathys metadata with subclass annotations into data frame
mathys_meta_df = Seu_mathys_obj@meta.data %>% as.data.frame()

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS('/external/rprshnas01/netdata_kcni/dflab/data/rosmap/phenotype/ROSmaster.rds')

# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex)

# two samples are missing diagnoses - these are AD cases according to the mathys paper 
ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset
ros_meta_small = ros_meta_small %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))

# merge mathys cell level meta with subject case information
mathys_meta_df = merge(mathys_meta_df , ros_meta_small, by = 'projid')


mathys_meta_df = mathys_meta_df %>% rename(subclass = predicted.id)

# count up cell counts per subclass and total per subject
cell_type_counts = mathys_meta_df %>% group_by(projid, subclass) %>% summarize(cell_type_count = n())
tot_cell_counts = mathys_meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df = merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop = cell_type_count / tot_cell_counts, 
         cell_type_prop_se = sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

# merge cell proportions with subject level meta
mathys_cell_prop_meta_long = merge(cell_prop_df, ros_meta_small, by = 'projid')

# plot cell proportions per subclass by gpath (global pathology)
mathys_cell_prop_meta_long %>% ggplot(aes(x = gpath, y = cell_type_prop, color = pathoAD, group = 1)) + 
  geom_smooth(method = "lm") + 
  geom_linerange(aes(ymin = cell_type_prop - cell_type_prop_se, 
                     ymax = cell_type_prop + cell_type_prop_se), 
                 alpha = .5) + 
  geom_point(alpha =  1) +
  facet_wrap(~subclass, scales = 'free_y') + 
  ylab('Mathys cell type proportions')

# plot cell proportions per subclass by pathoAD (global pathology)
mathys_cell_prop_meta_long %>% ggplot(aes(x = pathoAD, y = cell_type_prop)) + 
  geom_boxplot() + 
  geom_point(alpha =  1) +
  facet_wrap(~subclass, scales = 'free_y') + 
  ylab('Mathys cell type proportions')

# perform linear modeling to assess cell proportion differences

test_cell_types = mathys_cell_prop_meta_long$subclass %>% unique 

lm_formula = 'scale(cell_type_prop) ~ pathoAD + age_death + msex' # modeling cell proportions as a fxn of pathoAD and age_death and msex

# im using an lapply loop here, but this can be changed out
cell_prop_stats = lapply(test_cell_types, function(subclass_name){
  print(subclass_name)
  
  m = lm(lm_formula, data = mathys_cell_prop_meta_long %>% filter(subclass == subclass_name)  )
  model_df = tidy(m)
  keep_df = model_df[2, ]
  keep_df[1, 6] = subclass_name
  return(keep_df)
  
  #lm('factor(pathoAD) ~ SST + age_death + msex ', data = mathys_cell_prop_wide_meta) %>% summary()
}) %>% bind_rows()
colnames(cell_prop_stats)[6] = 'subclass'

# cell_prop_stats is the thing that keeps track of effect sizes
cell_prop_stats

write_csv(cell_prop_stats, "cell_prop_stats.csv")

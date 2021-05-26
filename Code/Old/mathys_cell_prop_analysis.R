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
ros_meta = readRDS("/external/rprshnas01/public_datasets2/rosmap/phenotype/ROSmaster.rds")

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


#### Plot requests Feb 1, 2021 ####

### figure 2

selected_mathys_cell_prop_meta_long <- mathys_cell_prop_meta_long %>% filter(subclass %in% c("SST", "IT"))
selected_mathys_cell_prop_meta_long[selected_mathys_cell_prop_meta_long$subclass == "SST", "subclass"] <- "Inh_SST"
selected_mathys_cell_prop_meta_long[selected_mathys_cell_prop_meta_long$subclass == "IT", "subclass"] <- "Exc_IT"
selected_mathys_cell_prop_meta_long$pathoAD <- as.character(selected_mathys_cell_prop_meta_long$pathoAD)
selected_mathys_cell_prop_meta_long[selected_mathys_cell_prop_meta_long$pathoAD == "0", "pathoAD"] <- "Control"
selected_mathys_cell_prop_meta_long[selected_mathys_cell_prop_meta_long$pathoAD == "1", "pathoAD"] <- "AD"
selected_mathys_cell_prop_meta_long$pathoAD <- factor(selected_mathys_cell_prop_meta_long$pathoAD, levels = c("Control", "AD"))
selected_mathys_cell_prop_meta_long$cell_type_prop <- selected_mathys_cell_prop_meta_long$cell_type_prop * 100

p <- selected_mathys_cell_prop_meta_long %>% ggplot(aes(x = pathoAD, y = cell_type_prop, fill = pathoAD)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("white", "light blue")) +
  geom_point(alpha = 1) +
  facet_wrap(~subclass, scales = 'free_y') + 
  ylab('Single nucleus cell type proportion (%)') +
  xlab("Alzheimer's status") +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 14),
        text = element_text(size = 11), 
        axis.text = element_text(size = 9),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))

ggsave(plot = p, 
       width = 180, 
       dpi = 300, 
       units = "mm", 
       filename = "IT_SST_Box.jpg",
       device = "jpeg")

### figure 3

cell_prop_stats <- merge(cell_prop_stats, Name_and_colour_scheme, by.x = "subclass", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
cell_prop_stats$Our_label <- factor(cell_prop_stats$Our_label, levels = c("Inh_SST", "Inh_PVALB", "Inh_LAMP5", "Exc_L6b", "OPC", "Exc_L5/6 NP", "Inh_VIP", "Exc_L5/6 IT Car3", "Oligodendrocyte", "Exc_IT", "Astrocyte", "Exc_L6 CT", "VLMC", "Endothelial", "Pericyte", "Microglia", "Inh_PAX6"))

ggplot(cell_prop_stats, aes(x = Our_label)) + 
  geom_bar(stat = 'identity', width = 0.8, aes(y = estimate, fill = Our_label)) +
  geom_errorbar((aes(ymin = estimate - std.error, ymax = estimate + std.error))) +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_fill_manual(values=c("Astrocyte" = "#73ABFF",     
                             "Endothelial" = "#7F00FF",
                             "Exc_IT" = "#52FF26",          
                             "Exc_L5/6 NP" = "#006B99",     
                             "Inh_LAMP5" = "#FF7373",       
                             "Microglia" = "#FF26A8",       
                             "Oligodendrocyte" = "#311799", 
                             "OPC" = "#3D4BCC",            
                             "Pericyte" = "#992E8E",
                             "Inh_SST" = "#FFE500",
                             "Inh_VIP" = "#996517",
                             "VLMC" = "#B65CCC",
                             "Inh_PVALB" =	"#B6CC5C",
                             "Exc_L5/6 IT Car3" =	"#00CC14",
                             "Exc_L6 CT" =	"#459967",
                             "Inh_PAX6" = "#CC683D",
                             "Exc_L6b" = "#4DFFC9")) +
  facet_wrap(~AIBS_class_label, scales = "free_x") +
  xlab('Cell types') + 
  ylab('Cell type association with AD') +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 20, hjust = 0.95, vjust = 0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size = 12),
        text = element_text(size = 10), 
        axis.text = element_text(size = 6.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) 

ggsave(width = 180, 
       dpi = 300, 
       units = "mm", 
       filename = "AD_Beta_plot.jpg",
       device = "jpeg")

prediction_score_threshold_low = 0.5
prediction_score_threshold_high = 0.8

library(Seurat)

Seu_AIBS_obj = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds')

Seu_AIBS_obj <- FindVariableFeatures(Seu_AIBS_obj, selection.method = "vst", 
                                       nfeatures = 10000, verbose = FALSE) #need variable features for transferring

cain_seurat = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_cain_obj_update_11JUN21.rds')
cain_meta = cain_seurat@meta.data
cain_meta$subclass = factor(cain_meta$predicted.id.AllHodge_ExpSubclas)
cain_meta$prediction.score.max = cain_meta$prediction.score.max.AllHodge_ExpSubclas

cain_seurat = FindVariableFeatures(cain_seurat, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE) #need variable features for transferring


#transfer
tanchors <- FindTransferAnchors(reference = Seu_AIBS_obj, query = cain_seurat, dims = 1:30) # find anchors for mapping/transfering/integrating data
predictions <- TransferData(anchorset = tanchors, 
                            refdata = Seu_AIBS_obj$subclass_label_expanded_L35IT, #this is the key part of specifying what class/subclass/cluster labels you want to map onto
                            dims = 1:30)

cain_meta = full_join(cain_seurat@meta.data %>% dplyr::select(orig.ident:Patho_AD_assumed, Cog_assumed:individualID, 
                                                                pathoAD_Update_22MAY21:LOAD_Update_22MAY21), 
                        predictions %>% tibble::rownames_to_column(var = 'cell_name') %>% 
                          rename(subclass = predicted.id)) 




# count up cell counts per subclass and total per subject
cell_type_counts = cain_meta %>% group_by(projid, subclass, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(cell_type_count = n)

total_cell_count_per_individual = cain_meta %>% group_by(projid, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(total_cell_count_per_individual = n)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df = merge(total_cell_count_per_individual, cell_type_counts) %>% 
  mutate(cell_type_proportion = cell_type_count / total_cell_count_per_individual, 
         cell_type_prop_se = sqrt((cell_type_proportion * (1 - cell_type_proportion))/total_cell_count_per_individual))

cain_ad = cell_prop_df


## try same for zhou

zhou_seurat = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_zhou_obj_update_11JUN21.rds')

zhou_meta = zhou_seurat@meta.data
zhou_meta$subclass = factor(zhou_meta$predicted.id.AllHodge_ExpSubclas)
zhou_meta$prediction.score.max = zhou_meta$prediction.score.max.AllHodge_ExpSubclas

zhou_seurat = FindVariableFeatures(zhou_seurat, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE) #need variable features for transferring


#transfer
tanchors <- FindTransferAnchors(reference = Seu_AIBS_obj, query = zhou_seurat, dims = 1:30) # find anchors for mapping/transfering/integrating data
zhou_predictions <- TransferData(anchorset = tanchors, 
                            refdata = Seu_AIBS_obj$subclass_label_expanded_L35IT, #this is the key part of specifying what class/subclass/cluster labels you want to map onto
                            dims = 1:30)


zhou_meta = full_join(zhou_seurat@meta.data %>% 
                          dplyr::select(orig.ident:AD_Group, 
                                        projid:dcfdx_lv,
                                        pathoAD_Update_22MAY21:LOAD_Update_22MAY21) %>%
                          tibble::rownames_to_column(var = 'TAG'), 
                        zhou_predictions %>% tibble::rownames_to_column(var = 'TAG') %>% 
                          rename(subclass = predicted.id))

zhou_meta$subclass = factor(zhou_meta$subclass)
zhou_meta$projid = as.numeric(zhou_meta$projid)

# count up cell counts per subclass and total per subject
cell_type_counts = zhou_meta %>% group_by(projid, subclass, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(cell_type_count = n)

total_cell_count_per_individual = zhou_meta %>% group_by(projid, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(total_cell_count_per_individual = n)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df = merge(total_cell_count_per_individual, cell_type_counts) %>% 
  mutate(cell_type_proportion = cell_type_count / total_cell_count_per_individual, 
         cell_type_prop_se = sqrt((cell_type_proportion * (1 - cell_type_proportion))/total_cell_count_per_individual))

zhou_ad = cell_prop_df



mathys_seurat = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_mathys_obj_update_08Jun21.rds')
mathys_meta = mathys_seurat@meta.data
mathys_meta$subclass = factor(mathys_meta$predicted.id.AllHodge_ExpSubclas)
mathys_meta$prediction.score.max = mathys_meta$prediction.score.max.AllHodge_ExpSubclas


mathys_seurat = FindVariableFeatures(mathys_seurat, selection.method = "vst", 
                                   nfeatures = 3000, verbose = FALSE) #need variable features for transferring


#transfer
tanchors <- FindTransferAnchors(reference = Seu_AIBS_obj, query = mathys_seurat, dims = 1:30) # find anchors for mapping/transfering/integrating data
mathys_predictions <- TransferData(anchorset = tanchors, 
                            refdata = Seu_AIBS_obj$subclass_label_expanded_L35IT, #this is the key part of specifying what class/subclass/cluster labels you want to map onto
                            dims = 1:30)

mathys_meta = full_join(mathys_seurat@meta.data %>% dplyr::select(orig.ident:Subcluster, pathoAD_Update_22MAY21:LOAD_Update_22MAY21), 
                        mathys_predictions %>% tibble::rownames_to_column(var = 'TAG') %>% 
                          rename(subclass = predicted.id))

mathys_meta$subclass = factor(mathys_meta$subclass)

# count up cell counts per subclass and total per subject
cell_type_counts = mathys_meta %>% group_by(projid, subclass, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(cell_type_count = n)

total_cell_count_per_individual = mathys_meta %>% group_by(projid, .drop = F) %>% 
  mutate(qc_passing = case_when((subclass != 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                (subclass == 'L2/3 IT') & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                TRUE ~ F)) %>% 
  filter(qc_passing) %>% count(.drop = F) %>% 
  rename(total_cell_count_per_individual = n)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df = merge(total_cell_count_per_individual, cell_type_counts) %>% 
  mutate(cell_type_proportion = cell_type_count / total_cell_count_per_individual, 
         cell_type_prop_se = sqrt((cell_type_proportion * (1 - cell_type_proportion))/total_cell_count_per_individual))

mathys_ad = cell_prop_df



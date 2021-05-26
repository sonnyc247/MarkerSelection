#### packages ####

library(dplyr)
library(Seurat)
library(magrittr)

#### May 22 update, adding "LOAD" to all objects as well as CgG-mapped and M1-mapped identities ####

# Load seurat objects

Seu_update_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_update_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_update_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)

# read in rosmap metadata
ros_meta = readRDS('/external/rprshnas01/kcni/ychen/collabgit/AD_snRNAseq/data/ROSmaster.rds')

# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex, pmi, braaksc, ceradsc, cogdx)
ros_meta_small %<>% mutate(LOAD = case_when((braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                            (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                            TRUE ~ 'OTHER')) 
ros_meta_small = ros_meta_small %>% mutate(msex = factor(msex), LOAD = factor(LOAD), pathoAD = factor(pathoAD))
colnames(ros_meta_small)[2:10]
colnames(ros_meta_small)[2:10] <- paste0(colnames(ros_meta_small)[2:10], "_", "Update_22MAY21")
colnames(ros_meta_small)
remove(ros_meta)

# make mapping-update df

update_df <- merge(predictions_Mathys_CgG[,c("predicted.id", "prediction.score.max")], predictions_Mathys_M1[,c("predicted.id", "prediction.score.max")], by = "row.names") # Mathys
update_df <- merge(predictions_Cain_CgG[,c("predicted.id", "prediction.score.max")], predictions_Cain_M1[,c("predicted.id", "prediction.score.max")], by = "row.names") # Cain
update_df <- merge(predictions_Zhou_CgG[,c("predicted.id", "prediction.score.max")], predictions_Zhou_M1[,c("predicted.id", "prediction.score.max")], by = "row.names") # Zhou

table(update_df$predicted.id.x)
names(update_df)[2:5] <- c("predicted.id.CgG", "prediction.score.max.CgG", "predicted.id.M1", "prediction.score.max.M1")
table(update_df$predicted.id.M1)

# get projID from seurat object, finalize update df

id_holder <- Seu_update_object@meta.data
id_holder <- tibble::rownames_to_column(id_holder)
id_holder <- id_holder[,c("rowname", "projid")]
length(intersect(update_df$Row.names, id_holder$rowname)) == nrow(update_df)
update_df <- merge(id_holder, update_df, by.x = "rowname", by.y = "Row.names", all.x = T, all.y = T)
update_df <- merge(update_df, ros_meta_small, by = "projid", all.x = T, all.y = F)
rownames(update_df) <- update_df$rowname

# add metadata

test <- AddMetaData(object = Seu_update_object, metadata = update_df)
identical(test$projid, Seu_update_object$projid)
identical(unname(test$rowname), as.character(test$TAG)) #check for mathys
identical(unname(test$rowname), as.character(test$cell_name)) #check for cain
identical(unname(test$rowname), as.character(colnames(test))) #check for zhou
table(test$predicted.id.CgG, exclude = "ifany")
remove(test)
update_df <- update_df[,3:15]
Seu_update_object <- AddMetaData(object = Seu_update_object, metadata = update_df)

# save object

saveRDS(Seu_update_object, "~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds") # for mathys
saveRDS(Seu_update_object, "~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_22MAY21.rds") # for cain
saveRDS(Seu_update_object, "~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_22MAY21.rds") # for zhou

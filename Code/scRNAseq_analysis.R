#### All involved packages ####

library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(cowplot)
library(ggbeeswarm)
library(broom)

#
#### Data loading - Cain et al ####

### read in data
Cain_cell_metadata <- readr::read_csv("/external/rprshnas01/external_data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv")
Cain_gene_metadata <- readr::read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv")
Cain_COO_matrix <- readr::read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv.gz")

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

#### Transfer anchor - Cain et al ####

#prep work
Seu_cain_obj <- FindVariableFeatures(Seu_cain_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
Idents(new_Seu_AIBS_obj) <- "subclass_label" #so that we're mapping at the subclass level

length(intersect(VariableFeatures(Seu_cain_obj), VariableFeatures(new_Seu_AIBS_obj)))
length(intersect(VariableFeatures(Seu_cain_obj), rownames(new_Seu_AIBS_obj)))
length(intersect(rownames(Seu_cain_obj), VariableFeatures(new_Seu_AIBS_obj)))

Seu_cain_obj <- subset(Seu_cain_obj, subset = subtype == "None.NA", invert = TRUE)
Seu_cain_obj <- FindVariableFeatures(Seu_cain_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #check variable features for transferring

#transfer
tanchors <- FindTransferAnchors(reference = new_Seu_AIBS_obj, query = Seu_cain_obj, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = new_Seu_AIBS_obj$subclass_label, dims = 1:30)
Seu_cain_obj <- AddMetaData(Seu_cain_obj, metadata = predictions)

Confusion_matrix <- data.frame(unclass(table(Seu_cain_obj$predicted.id, Seu_cain_obj$subtype)))

length(intersect(VariableFeatures(new_Seu_AIBS_obj), tanchors@anchor.features))

       

#### Proportion analysis - Cain et al ####

# edit some metadata for groups

Seu_cain_obj$Cog_assumed <- Seu_cain_obj$U4
identical(Seu_cain_obj$Cog_assumed, Seu_cain_obj$U4)

Seu_cain_obj$Cog_assumed[Seu_cain_obj$Cog_assumed == "Cdx1"] <- "Cog1" #unify cog group
Seu_cain_obj$Cog_assumed[Seu_cain_obj$Cog_assumed == "Cdx4"] <- "Cog4"
unique(Seu_cain_obj$Cog_assumed)

Seu_cain_obj$group_assumed <- paste0(Seu_cain_obj$Cog_assumed, "_", Seu_cain_obj$Patho_AD_assumed)

# pull metadata with subclass annotations into data frame
cain_meta_df <- Seu_cain_obj@meta.data %>% as.data.frame()
cain_meta_df <- cain_meta_df %>% mutate(Patho_AD_assumed = factor(Patho_AD_assumed))
cain_meta_df <- cain_meta_df %>% mutate(Cog_assumed = factor(Cog_assumed))
cain_meta_df <- cain_meta_df %>% mutate(group_assumed = factor(group_assumed))
cain_meta_df <- cain_meta_df %>% rename(subclass = predicted.id)

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- cain_meta_df %>% group_by(specimenID, subtype) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$specimenID, cell_type_counts_initial$subtype))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$specimenID, "_", cell_type_counts_initial$subtype)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("specimenID", "subtype")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- cain_meta_df %>% group_by(specimenID) %>% summarize(tot_cell_counts = n())

table(cain_meta_df$subtype, cain_meta_df$group_assumed)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
cain_cell_prop_meta_long <- merge(cell_prop_df, unique(cain_meta_df[,c("specimenID", "Patho_AD_assumed", "Patho_AD_Maybe", "U4", "Cog_assumed", "group_assumed")]), by = "specimenID")

cain_cell_prop_meta_long %>% 
  filter(subclass == "SST") %>% 
  group_by(group_assumed, subclass) %>% 
  summarize(m = mean(cell_type_prop))

cain_cellsubtype_prop_meta_long %>% 
  filter(subtype == "Inh.Inh.SST") %>% 
  group_by(group_assumed, subtype) %>% 
  summarize(m = mean(cell_type_prop))

# add more metadata

label_holder <- unique(Seu_cain_obj@meta.data[, c("specimenID", "convertednames")])
meta_holder <- merge(cain_rosmap_metadata, label_holder, by.x = "Sample_ID",  by.y = "convertednames")
cain_cell_prop_meta_long <- merge(cain_cell_prop_meta_long, meta_holder, by = "specimenID", all.x = T)

# plot

ggplot(cain_cell_prop_meta_long, aes(x = group_assumed, y = cell_type_prop, fill = group_assumed)) + 
  geom_boxplot(aes(middle = mean(cell_type_prop))) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("red", "green", "light blue", "purple")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cain_cell_prop_meta_long[cain_cell_prop_meta_long$subclass == "SST",], aes(x = group_assumed, y = cell_type_prop, fill = group_assumed)) + 
  geom_boxplot(aes(middle = mean(cell_type_prop))) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("red", "green", "light blue", "purple")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cain_cellsubtype_prop_meta_long[cain_cellsubtype_prop_meta_long$subtype == "Inh.Inh.SST",], aes(x = group_assumed, y = cell_type_prop, fill = group_assumed)) + 
  geom_boxplot(aes(middle = mean(cell_type_prop))) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("red", "green", "light blue", "purple")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

# export df for plotting

export_holder <- cain_cell_prop_meta_long[,c(13,14,12,30,2:6,11,15:29)]
colnames(export_holder)[c(5,7:10)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error", "Assumed_grouping_from_ID")
write.csv(export_holder, "Cain_cell_type_prop_df.csv")


#### Data loading - Zhou et al ####

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

#### Transfer anchor - Zhou et al ####

#prep work
Seu_zhou_obj <- FindVariableFeatures(Seu_zhou_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
Idents(new_Seu_AIBS_obj) <- "subclass_label" #so that we're mapping at the subclass level

length(intersect(VariableFeatures(Seu_zhou_obj), VariableFeatures(new_Seu_AIBS_obj)))
length(intersect(VariableFeatures(Seu_zhou_obj), rownames(new_Seu_AIBS_obj)))
length(intersect(rownames(Seu_zhou_obj), VariableFeatures(new_Seu_AIBS_obj)))

#transfer
tanchors <- FindTransferAnchors(reference = new_Seu_AIBS_obj, query = Seu_zhou_obj, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = new_Seu_AIBS_obj$subclass_label, dims = 1:30)
Seu_zhou_obj <- AddMetaData(Seu_zhou_obj, metadata = predictions)

length(intersect(VariableFeatures(Seu_zhou_obj), tanchors@anchor.features))

#### Proportion analysis - Zhou et al ####

# subset for QC reason (don't want too few or too many genes detected)
Seu_zhou_obj_QCed <- subset(Seu_zhou_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# pull metadata with subclass annotations into data frame
zhou_meta_df <- Seu_zhou_obj_QCed@meta.data %>% as.data.frame()
zhou_meta_df <- zhou_meta_df %>% rename(subclass = predicted.id)

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- zhou_meta_df %>% group_by(orig.ident, subclass) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$orig.ident, cell_type_counts_initial$subclass))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$orig.ident, "_", cell_type_counts_initial$subclass)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("orig.ident", "subclass")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- zhou_meta_df %>% group_by(orig.ident) %>% summarize(tot_cell_counts = n())

table(zhou_meta_df$subclass, zhou_meta_df$orig.ident)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
zhou_cell_prop_meta_long <- merge(cell_prop_df, unique(zhou_meta_df[,c("orig.ident", "AD_Group")]), by = "orig.ident")

zhou_cell_prop_meta_long %>% 
  filter(subclass == "SST") %>% 
  group_by(AD_Group, subclass) %>% 
  summarize(m = mean(cell_type_prop))

# add some metadata

meta_holder <- Zhou_rosmap_metadata[Zhou_rosmap_metadata$individualID %in% unique(Seu_zhou_obj_QCed$individualID),]
label_holder <- unique(Seu_zhou_obj_QCed@meta.data[, c("orig.ident", "individualID")])
meta_holder <- merge(meta_holder, label_holder, by = "individualID")
zhou_cell_prop_meta_long <- merge(zhou_cell_prop_meta_long, meta_holder, by = "orig.ident", all.x = T)

# plot

ggplot(zhou_cell_prop_meta_long, aes(x = AD_Group, y = cell_type_prop, fill = AD_Group)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(zhou_cell_prop_meta_long[zhou_cell_prop_meta_long$subclass == "SST",], aes(x = AD_Group, y = cell_type_prop, fill = AD_Group)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(zhou_cell_prop_meta_long$subclass, zhou_cell_prop_meta_long$AD_Group)

# export df for plotting

export_holder <- zhou_cell_prop_meta_long[,c(9,10,1,8,2:6,7,11:25)]
colnames(export_holder)[c(3,5,7:10)] <- c("Sample_ID_from_files", "total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error", "Assumed_grouping_from_ID")
write.csv(export_holder, "Zhou_cell_type_prop_df.csv")

#### Proportion analysis - Mathys et al ####

# pull metadata with subclass annotations into data frame
mathys_meta_df <- Seu_mathys_obj@meta.data %>% as.data.frame()
mathys_meta_df <- mathys_meta_df %>% rename(subclass = predicted.id)

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS("/external/rprshnas01/external_data/rosmap/phenotype/ROSmaster.rds")

# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex)

# two samples are missing diagnoses - these are AD cases according to the mathys paper 
ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset
ros_meta_small = ros_meta_small %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))
remove(ros_meta)

# merge mathys cell level meta with subject case information
mathys_meta_df = merge(mathys_meta_df , ros_meta_small, by = 'projid')

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- mathys_meta_df %>% group_by(projid, subclass) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$projid, cell_type_counts_initial$subclass))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$projid, "_", cell_type_counts_initial$subclass)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("projid", "subclass")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- mathys_meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

table(mathys_meta_df$subclass, mathys_meta_df$projid)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
mathys_cell_prop_meta_long <- merge(cell_prop_df, unique(mathys_meta_df[,c("projid", "pathoAD", "gpath", "age_death", "msex")]), by = "projid")

mathys_cell_prop_meta_long %>% 
  filter(subclass == "SST") %>% 
  group_by(pathoAD, subclass) %>% 
  summarize(m = mean(cell_type_prop))

# plot

ggplot(mathys_cell_prop_meta_long, aes(x = pathoAD, y = cell_type_prop, fill = pathoAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(mathys_cell_prop_meta_long[mathys_cell_prop_meta_long$subclass == "SST",], aes(x = pathoAD, y = cell_type_prop, fill = pathoAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(mathys_cell_prop_meta_long$subclass, mathys_cell_prop_meta_long$pathoAD)

# export df for plotting

export_holder <- mathys_cell_prop_meta_long
colnames(export_holder)[c(2,4:6)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")
write.csv(export_holder, "mathys_cell_type_prop_df.csv")

#### Revisited Proportion analysis (filter IT) ####

# load object
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_prop_object <- subset(Seu_prop_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# filter bad-fit IT cells
ncol(Seu_prop_object)
table(Seu_prop_object$predicted.id) 
Seu_prop_object <- subset(Seu_prop_object, subset = predicted.id == "IT" & prediction.score.max < 0.8, invert = T) 
table(Seu_prop_object$predicted.id)

# pull metadata with subclass annotations into data frame
prop_meta_df <- Seu_prop_object@meta.data %>% as.data.frame()
prop_meta_df <- prop_meta_df %>% rename(subclass = predicted.id)

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS("~/collabgit/AD_snRNAseq/data/ROSmaster.rds")

# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex, braaksc, ceradsc, cogdx)

# FOR MATHYS: two samples are missing diagnoses - these are AD cases according to the mathys paper 
ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset

# format metadata and make new LOAD variable
ros_meta_small = ros_meta_small %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))
ros_meta_small %<>% mutate(LOAD = case_when((braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                            (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                            TRUE ~ 'OTHER')) 
remove(ros_meta) # no longer needed

# merge cell level meta with subject case information
prop_meta_df <- prop_meta_df[,c("projid", "subclass", "orig.ident", "nCount_RNA", "nFeature_RNA")] # for cain and zhou; to avoid duplicated columns
prop_meta_df = merge(prop_meta_df , ros_meta_small, by = 'projid')

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- prop_meta_df %>% group_by(projid, subclass) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$projid, cell_type_counts_initial$subclass))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$projid, "_", cell_type_counts_initial$subclass)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("projid", "subclass")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- prop_meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

table(prop_meta_df$subclass, prop_meta_df$projid)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
cell_prop_meta_long <- merge(cell_prop_df, unique(prop_meta_df[,c("projid", "pathoAD", "gpath", "age_death", "msex", "braaksc", "ceradsc", "cogdx", "LOAD")]), by = "projid")

cell_prop_meta_long %>% 
  filter(subclass == "IT") %>% 
  group_by(LOAD, subclass) %>% 
  summarize(m = mean(cell_type_prop))

# plot

ggplot(cell_prop_meta_long, aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cell_prop_meta_long[cell_prop_meta_long$subclass == "IT",], aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(cell_prop_meta_long$subclass, cell_prop_meta_long$LOAD)

# export df for plotting

export_holder <- cell_prop_meta_long
colnames(export_holder)[c(2,4:6)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")
write.csv(export_holder, "zhou_cell_type_prop_df.csv")


#### Further Revisited Proportion analysis (filter None.NA for Cain only) ####

# load object
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_prop_object <- subset(Seu_prop_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_prop_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_prop_object <- subset(Seu_prop_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# filter bad-fit IT cells
ncol(Seu_prop_object)
table(Seu_prop_object$predicted.id) 
Seu_prop_object <- subset(Seu_prop_object, subset = predicted.id == "IT" & prediction.score.max < 0.8, invert = T) 
table(Seu_prop_object$predicted.id)

# pull metadata with subclass annotations into data frame
prop_meta_df <- Seu_prop_object@meta.data %>% as.data.frame()
prop_meta_df <- prop_meta_df %>% rename(subclass = predicted.id)

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS("~/collabgit/AD_snRNAseq/data/ROSmaster.rds")

# select just the columns from ros master that we need
ros_meta_small = ros_meta %>% select(projid, pathoAD, gpath, age_death, msex, braaksc, ceradsc, cogdx)

# FOR MATHYS: two samples are missing diagnoses - these are AD cases according to the mathys paper 
ros_meta_small[is.na(ros_meta_small$pathoAD), 'pathoAD'] = 1 # 1 is the label for AD in this dataset

# format metadata and make new LOAD variable
ros_meta_small = ros_meta_small %>% mutate(pathoAD = factor(pathoAD), msex = factor(msex))
ros_meta_small %<>% mutate(LOAD = case_when((braaksc >= 4 & ceradsc <= 2 & cogdx == 4) ~ 'AD',
                                            (braaksc <= 3 & ceradsc >= 3 & cogdx == 1) ~ 'C',
                                            TRUE ~ 'OTHER')) 
remove(ros_meta) # no longer needed

# merge cell level meta with subject case information
prop_meta_df <- prop_meta_df[,c("projid", "subclass", "orig.ident", "nCount_RNA", "nFeature_RNA")] # for cain and zhou; to avoid duplicated columns
prop_meta_df = merge(prop_meta_df , ros_meta_small, by = 'projid')

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- prop_meta_df %>% group_by(projid, subclass) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$projid, cell_type_counts_initial$subclass))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$projid, "_", cell_type_counts_initial$subclass)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("projid", "subclass")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- prop_meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

table(prop_meta_df$subclass, prop_meta_df$projid)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
cell_prop_meta_long <- merge(cell_prop_df, unique(prop_meta_df[,c("projid", "pathoAD", "gpath", "age_death", "msex", "braaksc", "ceradsc", "cogdx", "LOAD")]), by = "projid")

cell_prop_meta_long[cell_prop_meta_long$subclass == "IT",] %>% 
  group_by(LOAD, subclass) %>% 
  summarize(m = mean(cell_type_prop))

# plot

ggplot(cell_prop_meta_long, aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cell_prop_meta_long[cell_prop_meta_long$subclass == "IT",], aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(cell_prop_meta_long$subclass, cell_prop_meta_long$LOAD)

# export df for plotting

export_holder <- cell_prop_meta_long
colnames(export_holder)[c(2,4:6)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")
write.csv(export_holder, "cain_cell_type_prop_df_noNA.csv")



#### Revisited mapping (CgG only) ####

# load objects
Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

#prep work
unique(Seu_ref_object$outlier_call)
table(Seu_ref_object$NeuN_Region)
Seu_ref_object
Seu_ref_object <- subset(Seu_ref_object, subset = NeuN_Region == "MTG_Neuronal", invert = TRUE)
table(Seu_ref_object$NeuN_Region)
Seu_ref_object

Idents(Seu_ref_object) <- "subclass_label" #so that we're mapping at the subclass level
Seu_ref_object <- FindVariableFeatures(Seu_ref_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
table(Idents(Seu_ref_object))

Seu_map_object <- FindVariableFeatures(Seu_map_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferrin

length(intersect(VariableFeatures(Seu_map_object), VariableFeatures(Seu_ref_object)))
length(intersect(VariableFeatures(Seu_map_object), rownames(Seu_ref_object)))
length(intersect(rownames(Seu_map_object), VariableFeatures(Seu_ref_object)))

#transfer
tanchors <- FindTransferAnchors(reference = Seu_ref_object, query = Seu_map_object, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = Seu_ref_object$subclass_label, dims = 1:30)

#### Revisited mapping (M1 AIBS data) ####

# make and prep M1 object

metadata <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/Bakken_et_al_2020/metadata.csv")
matrix <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/Bakken_et_al_2020/matrix.csv")

row.names(metadata) <- metadata$sample_name
row.names(matrix) <- matrix$sample_name
count_matrix <- matrix[,-1]
remove(matrix)

Seu_M1AIBS_obj <- CreateSeuratObject(counts = t(count_matrix), meta.data = metadata)

test <- Seu_M1AIBS_obj@meta.data #check that seurat object was created properly - with metadata
identical(metadata, test[,4:42])
remove(test)
identical(table(Seu_M1AIBS_obj$subclass_label), table(metadata$subclass_label))
remove(metadata)

test <- Seu_M1AIBS_obj[["RNA"]]@counts #check counts in Seurat object too
test <- as.data.frame(test[1:100,1:100])
test2 <- count_matrix[1:100,1:100]
test2 <- t(test2)
test2 <- as.data.frame(test2)
identical(rownames(test), rownames(test2))
test2 <- sapply(test2, as.numeric)
test2 <- as.data.frame(test2)
rownames(test2) <- rownames(test)
identical(test2, test)
remove(test)
remove(test2)
remove(count_matrix)

# normalize data

Seu_M1AIBS_obj <- NormalizeData(Seu_M1AIBS_obj, normalization.method = "LogNormalize", scale.factor = 1000000) #preprocess/normalize seurat object
Seu_M1AIBS_obj[["RNA"]]@data[1:20,1:40] #see normalized data
Seu_M1AIBS_obj[["RNA"]]@counts[1:20,1:40] #see counts for comparison

# check for outliers and other metrics

table(Seu_M1AIBS_obj$outlier_call, Seu_M1AIBS_obj$outlier_type, exclude = NULL)
table(Seu_M1AIBS_obj$region_label)
Idents(Seu_M1AIBS_obj) <- "subclass_label" #so that we're mapping/working at the subclass level
table(Idents(Seu_M1AIBS_obj))
table(Seu_M1AIBS_obj$subclass_label)

# find variable features

head(Idents(Seu_M1AIBS_obj)) #check that we're mapping at/using the subclass level - in case it affects variable features (it probably doesn't)
Seu_M1AIBS_obj <- FindVariableFeatures(Seu_M1AIBS_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
length(Seu_M1AIBS_obj@assays$RNA@var.features)

# save Seu object as rds

saveRDS(Seu_M1AIBS_obj, "~/git/Ex_Env_Storage/MarkerSelection/Seu_M1AIBS_obj.rds")

# load Seu objects

Seu_ref_object <- Seu_M1AIBS_obj
remove(Seu_M1AIBS_obj)
Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_M1AIBS_obj.rds")

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# prep as needed

Seu_ref_object <- FindVariableFeatures(Seu_ref_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
length(Seu_ref_object@assays$RNA@var.features) #check
Idents(Seu_ref_object) <- "subclass_label" #so that we're mapping at the subclass level
table(Idents(Seu_ref_object))

test <- Seu_map_object@assays$RNA@var.features
Seu_map_object <- FindVariableFeatures(Seu_map_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferrin
test2 <- Seu_map_object@assays$RNA@var.features
identical(test, test2)
remove(test)
remove(test2)

length(intersect(VariableFeatures(Seu_map_object), VariableFeatures(Seu_ref_object)))
length(intersect(VariableFeatures(Seu_map_object), rownames(Seu_ref_object)))
length(intersect(rownames(Seu_map_object), VariableFeatures(Seu_ref_object)))

#transfer
tanchors <- FindTransferAnchors(reference = Seu_ref_object, query = Seu_map_object, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = Seu_ref_object$subclass_label, dims = 1:30)



Seu_zhou_obj <- AddMetaData(Seu_zhou_obj, metadata = predictions)

length(intersect(VariableFeatures(Seu_zhou_obj), tanchors@anchor.features))

# export df for plotting

export_holder <- cell_prop_meta_long
colnames(export_holder)[c(2,4:6)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")
write.csv(export_holder, "cain_cell_type_prop_df_noNA.csv")

#### Confusion matrix cgg vs m1 ####

# qc thresh cut for IT cells

ncol(Seu_plot_object)
table(Seu_plot_object$predicted.id)
Seu_plot_object <- subset(Seu_plot_object, subset = predicted.id == "IT" & prediction.score.max < 0.8, invert = T) 
table(Seu_plot_object$predicted.id)

# get numbers for confusion matrix

df_for_CgGvM1_plot <- merge(predictions_Zhou_CgG, predictions_Zhou_M1, by = "row.names")
df_for_CgGvM1_plot <- merge(predictions_Mathys_CgG[!(predictions_Mathys_CgG$prediction.score.max < 0.8
                                                 & predictions_Mathys_CgG$predicted.id == "IT"),], 
                                                 predictions_Mathys_M1, 
                                                 by = "row.names",
                                                 all.x = T,
                                                 all.y = F)

confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$predicted.id.y, df_for_CgGvM1_plot$predicted.id.x)) 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

#heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

mathys_M1pct <- ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("CgG-mapped subclass") + 
  ylab("M1-mapped subclass") +
  labs(fill = "% of M1- \n mapped cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

plot_grid(mathys_M1pct,
          cain_M1pct,
          zhou_M1pct,
          labels = c("Mathys", "Cain", "Zhou"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

### export plot (most recently plotted)

ggsave(width = 180,
       height = 400,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "CgGvsM1_M1pct_IT8filter.pdf",
       device = "pdf")

#### Confusion matrix original vs m1 ####

# load data

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only

# get numbers for confusion matrix

metadata_for_plot <- Seu_map_object@meta.data

metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "Subcluster")] # for Mathys
identical(rownames(metadata_for_plot), rownames(predictions_Mathys_M1))
df_for_CgGvM1_plot <- merge(metadata_for_plot, predictions_Mathys_M1, by = "row.names")

table(metadata_for_plot$subtype) # for Cain
metadata_for_plot <- metadata_for_plot[,c("nCount_RNA", "nFeature_RNA", "subtype")] 
identical(rownames(metadata_for_plot), rownames(predictions_Cain_M1))
df_for_CgGvM1_plot <- merge(metadata_for_plot, predictions_Cain_M1, by = "row.names")
df_for_CgGvM1_plot <- merge(predictions_Cain_M1[!(predictions_Cain_M1$prediction.score.max < 0.8
                                                & predictions_Cain_M1$predicted.id == "L2/3 IT"),], 
                            metadata_for_plot, 
                            by = "row.names",
                            all.x = T,
                            all.y = F)


confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$Subcluster, df_for_CgGvM1_plot$predicted.id)) # for mathys 
confusion_martix_hold <- as.matrix(table(df_for_CgGvM1_plot$subtype, df_for_CgGvM1_plot$predicted.id)) # for cain 
confusion_martix_hold <- as.matrix(confusion_martix_hold)
confusion_martix_hold <- (confusion_martix_hold/rowSums(confusion_martix_hold))*100

# plot heatmap

display.brewer.all()
dev.off() #as needed to reset graphics

#heatmap(confusion_martix_hold, col=brewer.pal(9 ,"Blues"), Rowv=TRUE, Colv=TRUE)

pheatmap::pheatmap(confusion_martix_hold, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)
#pheatmap::pheatmap(confusion_martix_hold_filt, treeheight_row = 0, treeheight_col = 0, color = brewer.pal(9 ,"Blues"), cluster_rows = F, cluster_cols = F)

### plot in ggplot

melted_conf_mtx <- melt(t(confusion_martix_hold[nrow(confusion_martix_hold):1,]))

cain_M1pct_08 <- ggplot(melted_conf_mtx, aes(Var1,Var2, fill=value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = brewer.pal(9 ,"Blues")) +
  theme_classic() +
  xlab("Mapped subclass") + 
  ylab("Pre-mapping cell group") +
  labs(fill = "% of pre- \n mapping \n cell group") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        axis.title.y =  element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.spacing = element_blank(),
        panel.background = element_blank()) +
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,-10),
        axis.line.y = element_blank(),
        axis.line.x = element_blank()) 

plot_grid(mathys_M1pct,
          cain_M1pct,
          cain_M1pct_08,
          labels = c("Mathys", "Cain", "Cain L2/3 IT Filter .8"),
          label_size = 10,
          hjust = -0.15,
          align = "hv",
          ncol = 1)

### export plot (most recently plotted)

ggsave(width = 180,
       height = 550,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/Heatmap/",
       filename = "OrigvsM1.pdf",
       device = "pdf")

#### Proportion analysis revisited - All datasets - M1, CgG, unfiltered, filtered ####

### takes place after updating the saved seurat object in "Z9_updating_seurat_objects.R" 
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds")
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_22MAY21.rds") #load cain seurat object (instead)
Seu_prop_obj <- subset(Seu_prop_obj, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_22MAY21.rds") #load zhou seurat object (instead)
Seu_prop_obj <- subset(Seu_prop_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# pull metadata with subclass annotations into data frame
meta_df <- Seu_prop_obj@meta.data %>% as.data.frame()
remove(Seu_prop_obj)
meta_df_backup <- meta_df
meta_df <- meta_df_backup
meta_df <- meta_df[!(meta_df$predicted.id.CgG == "IT" & meta_df$prediction.score.max.CgG < 0.8),] # as appropriate
meta_df <- meta_df[!(meta_df$predicted.id.M1 == "L2/3 IT" & meta_df$prediction.score.max.M1 < 0.8),] # as appropriate
meta_df <- meta_df %>% rename(subclass = predicted.id.CgG) # as appropriate
meta_df <- meta_df %>% rename(subclass = predicted.id.M1) # as appropriate

# count up cell counts per subclass and total per subject
cell_type_counts_initial <- meta_df %>% group_by(projid, subclass) %>% summarize(cell_type_count = n())

cell_type_counts <- as.data.frame(table(cell_type_counts_initial$projid, cell_type_counts_initial$subclass))
cell_type_counts$bridger <- paste0(cell_type_counts$Var1, "_", cell_type_counts$Var2)
cell_type_counts_initial$bridger <- paste0(cell_type_counts_initial$projid, "_", cell_type_counts_initial$subclass)
cell_type_counts <- merge(cell_type_counts[,-3], cell_type_counts_initial[,3:4], by = "bridger", all.x = T, all.y = T)
cell_type_counts <- cell_type_counts[,2:4]
colnames(cell_type_counts)[1:2] <- c("projid", "subclass")
cell_type_counts$cell_type_count[is.na(cell_type_counts$cell_type_count)] <- 0

tot_cell_counts <- meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

table(meta_df$subclass, meta_df$projid)

# calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta
cell_prop_meta_long <- merge(cell_prop_df, unique(meta_df[,c("projid", "age_death_Update_22MAY21", "msex_Update_22MAY21", "pmi_Update_22MAY21", "LOAD_Update_22MAY21")]), by = "projid")
colnames(cell_prop_meta_long)[7:10] <- c("age_death", "msex", "pmi", "LOAD")

cell_prop_meta_long %>% 
  filter(subclass == "SST") %>% 
  group_by(LOAD, subclass) %>% 
  summarize(m = mean(cell_type_prop))

# plot

ggplot(cell_prop_meta_long, aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cell_prop_meta_long[cell_prop_meta_long$subclass == "SST",], aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot(aes(middle = mean(cell_type_prop))) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(cell_prop_meta_long$subclass, cell_prop_meta_long$LOAD)

# export df for plotting

export_holder <- cell_prop_meta_long
colnames(export_holder)[c(2,4:6)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")

write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/zhou_cell_type_prop_df.csv")
write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Unfiltered/zhou_cell_type_prop_df.csv")
write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Filtered/zhou_cell_type_prop_df.csv")
write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/zhou_cell_type_prop_df.csv")


#
#### lm of MGPs and for rosmap factors ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/CgG_Unfiltered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
mathys_cell_type_prop_df$dataset <- "Mathys"
zhou_cell_type_prop_df$dataset <- "Zhou"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) 

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, zhou_cell_type_prop_df, cain_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'Sst', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check subclasses and datasets
cell_type_list <- as.character(unique(ad_snrnaseq_df$subclass))
table(ad_snrnaseq_df$subclass, ad_snrnaseq_df$dataset)

cell_type_list = cell_type_list[! cell_type_list %in% c('L5/6 IT Car3')] #remove subclasses that are not in all datasets; CgG map
cell_type_list = cell_type_list[! cell_type_list %in% c('L5 IT', 'L5/6 NP', 'L6 IT', 'L6 IT Car3', 'L6b', 'Sst Chodl', 'VLMC', 'L5 ET')] #remove subclasses that are not in all datasets; M1 map

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD'), subclass %in% cell_type_list) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, "subclass"]
cell_type_list = cell_type_list[! cell_type_list %in% c('L6b')] #remove subclasses that are pure 0 in 1 dataset; CgG map
remove(test)

# calculate beta coefficients

beta_coefs_non_meta_df = lapply(levels(ad_snrnaseq_df$dataset), function(curr_dataset){
  print(curr_dataset)
  
  dataset_df = lapply(cell_type_list, function(curr_cell_type){
    print(curr_cell_type)
    df = ad_snrnaseq_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    
    return(model_df) 
  }) %>% bind_rows()
  return(dataset_df)  
  #beta_out = 
}) %>% bind_rows()

# for troubleshooting the above function (basically convert the above into for loops):

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = ad_snrnaseq_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$dataset == curr_dataset & df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
  }
  
}

### Plot results

# sets factor levels to match order we want

beta_coefs_non_meta_df$dataset = factor(beta_coefs_non_meta_df$dataset, levels = c('Mathys', 'Zhou', 'Cain'))
beta_coefs_non_meta_df = merge(beta_coefs_non_meta_df, subclass_meta, by.x = 'subclass', by.y = 'subclass') #cgg mapping
Seu_M1AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_M1AIBS_obj.rds")
subclass_meta_M1 <- unique(Seu_M1AIBS_obj@meta.data[c("subclass_label", "class_label", "class_color")])
remove(Seu_M1AIBS_obj)
colnames(subclass_meta_M1)[1:2] <- c("subclass", "class")
beta_coefs_non_meta_df = merge(beta_coefs_non_meta_df, subclass_meta_M1, by.x = 'subclass', by.y = 'subclass') #M1 mapping

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_wrap(~subclass) + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('')

### save results

beta_coefs_non_meta_df$mapping <- "M1"
beta_coefs_non_meta_df$filtering <- "Filtered"

beta_coefs_non_meta_df_CgG_Unfiltered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_CgG_Filtered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_M1_Unfiltered <- beta_coefs_non_meta_df
beta_coefs_non_meta_df_M1_Filtered <- beta_coefs_non_meta_df

### combine results

beta_coefs_non_meta_df_M1_Filtered$subclass <-toupper(beta_coefs_non_meta_df_M1_Filtered$subclass)
beta_coefs_non_meta_df_M1_Unfiltered$subclass <-toupper(beta_coefs_non_meta_df_M1_Unfiltered$subclass)

identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_M1_Unfiltered))
identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_M1_Filtered))
identical(colnames(beta_coefs_non_meta_df_CgG_Unfiltered), colnames(beta_coefs_non_meta_df_CgG_Filtered))

beta_coefs_non_meta_df_combined = bind_rows(beta_coefs_non_meta_df_CgG_Unfiltered, 
                                            beta_coefs_non_meta_df_CgG_Filtered,
                                            beta_coefs_non_meta_df_M1_Unfiltered,
                                            beta_coefs_non_meta_df_M1_Filtered)
beta_coefs_non_meta_df_combined[beta_coefs_non_meta_df_combined$class == "Non-neuronal", "class"] <- "Non-Neuronal"
beta_coefs_non_meta_df_combined$filtering = factor(beta_coefs_non_meta_df_combined$filtering, levels = c('Unfiltered', 'Filtered'))

### combined graph

beta_coefs_non_meta_df_combined %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate, fill = class_color)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_identity() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(subclass ~ mapping*filtering, scales = "free") + 
  #facet_wrap(~mapping*filtering, scales = "free") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('')


#### Original identity proportion analysis ####

### Load the datasets

Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds")
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_22MAY21.rds") #load cain seurat object (instead)
Seu_prop_obj <- subset(Seu_prop_obj, subset = subtype == "None.NA", invert = TRUE) # for cain object only

### pull metadata with subclass annotations into data frame; format metadata
meta_df <- Seu_prop_obj@meta.data %>% as.data.frame()
remove(Seu_prop_obj)
names(meta_df)

#for mathys
meta_df <- meta_df %>% rename(subclass = Subcluster, class = broad.cell.type)
meta_df <- meta_df %>% select(TAG, projid, class, subclass, age_death_Update_22MAY21, msex_Update_22MAY21, pmi_Update_22MAY21, LOAD_Update_22MAY21)
names(meta_df)[5:8] <- c("age_death", "msex", "pmi", "LOAD")

#for cain
meta_df <- meta_df %>% rename(subclass = subtype, class = broad_class)
meta_df <- meta_df %>% select(cell_name, projid, class, subclass, age_death_Update_22MAY21, msex_Update_22MAY21, pmi_Update_22MAY21, LOAD_Update_22MAY21)
names(meta_df)[5:8] <- c("age_death", "msex", "pmi", "LOAD")
    
### count up cell counts per subclass and total per subject
cell_type_counts <- as.data.frame(table(meta_df$subclass, meta_df$projid))
names(cell_type_counts) <- c("subclass", "projid", "cell_type_count")
tot_cell_counts <- meta_df %>% group_by(projid) %>% summarize(tot_cell_counts = n())

### calculate cell proportions and standard errors per subclass per subject
cell_prop_df <- merge(tot_cell_counts, cell_type_counts) %>% 
  mutate(cell_type_prop <- cell_type_count / tot_cell_counts, 
         cell_type_prop_se <- sqrt((cell_type_prop * (1 - cell_type_prop))/tot_cell_counts))

names(cell_prop_df)[5] <- "cell_type_prop"
names(cell_prop_df)[6] <- "cell_type_prop_se"

# merge cell proportions with subject level meta

cell_prop_meta_long <- merge(cell_prop_df, unique(meta_df[,c("projid", "age_death", "msex", "pmi", "LOAD")]), by = "projid")

cell_prop_meta_long %>% 
  filter(subclass == "Inh.Inh.ADARB2.LAMP5.1") %>% 
  group_by(LOAD, subclass) %>% 
  summarize(m = mean(cell_type_prop))

### plot

ggplot(cell_prop_meta_long, aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  facet_wrap(~ subclass, scales = "free") +
  xlab('Cell types') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

ggplot(cell_prop_meta_long[cell_prop_meta_long$subclass == "Inh.Inh.ADARB2.LAMP5.1",], aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
  geom_boxplot(aes(middle = mean(cell_type_prop))) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("white", "light blue", "dark blue")) +
  xlab('AD Group') + 
  ylab('Cell type proportion') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "white"),
        strip.text = element_text(colour = 'black'))

table(cell_prop_meta_long$subclass, cell_prop_meta_long$LOAD)

### export df 

export_holder <- merge(cell_prop_meta_long, unique(meta_df[,c("subclass", "class")]), by = "subclass")
export_holder <- export_holder[,c(2,3,1,11,4:10)]
colnames(export_holder)[c(2,5:7)] <- c("total_cell_count_per_individual", "cell_count_per_subclass_per_indiv", "cell_type_proportion", "cell_type_proportion_std_error")

write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/Original_Identities/mathys_cell_type_prop_df.csv")
write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/Original_Identities/cain_cell_type_prop_df.csv")

mathys_cell_type_prop_df <- export_holder
cain_cell_type_prop_df <- export_holder

### lm 

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
mathys_cell_type_prop_df$dataset <- "Mathys"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, cain_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Cain'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'In4', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD')) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, "subclass"]
remove(test)

# calculate beta coefficients

beta_coefs_non_meta_df_orig <- beta_coefs_non_meta_df[0,c(2:7,1)]

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  temp_df <- ad_snrnaseq_df[ad_snrnaseq_df$dataset == curr_dataset, ]
  cell_type_list <- unique(temp_df$subclass)
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = temp_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    beta_coefs_non_meta_df_orig <- bind_rows(beta_coefs_non_meta_df_orig, model_df)
  }
  
}

### Plot results

# sets factor levels to match order we want

beta_coefs_non_meta_df_orig$dataset = factor(beta_coefs_non_meta_df_orig$dataset, levels = c('Mathys', 'Cain'))
beta_coefs_non_meta_df_orig = merge(beta_coefs_non_meta_df_orig, unique(ad_snrnaseq_df[,c("subclass", "class")]), by.x = 'subclass', by.y = 'subclass') #cgg mapping

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df_orig %>% filter(term == 'LOADAD', class %in% c("Ex", "Exc")) %>% 
  ggplot(aes(x = subclass, y = estimate, fill = dataset)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("red", "blue")) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_wrap(~dataset, scales = "free_x", ncol = 1) + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

#
#### All mapped identity proportion analysis ####

### Load the datasets

cain_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/cain_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
mathys_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/mathys_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)
zhou_cell_type_prop_df <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/M1_Filtered/zhou_cell_type_prop_df.csv", row.names=1, stringsAsFactors=TRUE)

### lm 

### adjustments to and combining of imported datasets

cain_cell_type_prop_df$dataset <- "Cain"
zhou_cell_type_prop_df$dataset <- "Zhou"
mathys_cell_type_prop_df$dataset <- "Mathys"

identical(colnames(cain_cell_type_prop_df), colnames(mathys_cell_type_prop_df)) #checks before rbind all together
identical(colnames(cain_cell_type_prop_df), colnames(zhou_cell_type_prop_df)) #checks before rbind all together

ad_snrnaseq_df = bind_rows(mathys_cell_type_prop_df, cain_cell_type_prop_df, zhou_cell_type_prop_df)

ad_snrnaseq_df$LOAD = factor(ad_snrnaseq_df$LOAD, levels = c('C', 'AD', 'OTHER'))
ad_snrnaseq_df$dataset = factor(ad_snrnaseq_df$dataset, levels = c('Mathys', 'Cain', 'Zhou'))
ad_snrnaseq_df$msex = factor(ad_snrnaseq_df$msex)
str(ad_snrnaseq_df)

### some checks

# tally of individuals by case and dataset
ad_snrnaseq_df %>% select(projid, dataset, LOAD) %>%
  distinct(.keep_all = T) %>% group_by(dataset, LOAD) %>% tally()

# plot of data
ad_snrnaseq_df %>% 
  filter(subclass == 'Sst', LOAD %in% c('C', 'AD')) %>% 
  ggplot(aes(x = LOAD, y = cell_type_proportion * 100)) + 
  geom_boxplot(outlier.shape = NA, ) + 
  geom_quasirandom() + 
  facet_wrap(~dataset, scales = 'free_x') + 
  ylab('SST snCTP (%)') + 
  xlab('')

### calculate LOAD beta coefficients

# check for cases where all counts of a cell type = zero (mean would = exactly 0)
test <- ad_snrnaseq_df %>% select(cell_type_proportion, projid, subclass, dataset, LOAD) %>% filter(LOAD %in% c('C', 'AD')) %>%
  group_by(subclass, dataset, LOAD) %>% summarize(mean_prop = mean(cell_type_proportion, na.rm = TRUE))
test[test$mean_prop == 0, c("subclass", "dataset")]
remove(test)

# for troubleshooting the above function (basically convert the above into for loops):

beta_coefs_non_meta_df_ALL <- beta_coefs_non_meta_df[0,c(2:7,1)] # first time

for (curr_dataset in unique(ad_snrnaseq_df$dataset)) {
  print(curr_dataset)
  temp_df <- ad_snrnaseq_df[ad_snrnaseq_df$dataset == curr_dataset, ]
  cell_type_list <- unique(temp_df$subclass)
  
  #if (curr_dataset == "Zhou"){
  #cell_type_list <- cell_type_list[cell_type_list != "L6b"]
  #}  #CgG only
  
  for (curr_cell_type in cell_type_list){
    print(curr_cell_type)
    df = temp_df %>% filter(LOAD %in% c('C', 'AD'))
    df = df[df$subclass == curr_cell_type, ]
    my_model = lm('scale(cell_type_proportion) ~ scale(age_death) + factor(msex) + scale(pmi) + LOAD ',
                  data = df)
    model_df = tidy(my_model)
    model_df$dataset = curr_dataset
    model_df$subclass = curr_cell_type
    model_df$mapping = "M1" # as appropriate
    model_df$filtering = "Filtered" # as appropriate
    beta_coefs_non_meta_df_ALL <- bind_rows(beta_coefs_non_meta_df_ALL, model_df)
  }
  
}

beta_coefs_non_meta_df_ALL$mapping <-"CgG" # after first run
beta_coefs_non_meta_df_ALL$filtering <-"Filtered" # after first run

beta_coefs_non_meta_df_backup <- beta_coefs_non_meta_df_ALL # as desired
table(beta_coefs_non_meta_df_ALL$mapping, beta_coefs_non_meta_df_ALL$filtering, exclude = "ifany")

### Plot results

# final adjustments to metadata/factors

beta_coefs_non_meta_df_ALL$dataset = factor(beta_coefs_non_meta_df_ALL$dataset, levels = c('Mathys', 'Cain', "Zhou"))
beta_coefs_non_meta_df_ALL$filtering = factor(beta_coefs_non_meta_df_ALL$filtering, levels = c('Unfiltered', 'Filtered'))
subclass_meta_total <- bind_rows(subclass_meta, subclass_meta_M1)
subclass_meta_total[subclass_meta_total$class == "Non-neuronal", "class"] <- "Non-Neuronal"
beta_coefs_non_meta_df_ALL = merge(beta_coefs_non_meta_df_ALL, unique(subclass_meta_total[,c("subclass", "class")]), by.x = 'subclass', by.y = 'subclass') #cgg mapping

# plots beta coefficients for LOAD across each dataset faceted by cell type

beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = subclass, y = estimate, fill = filtering)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("blue", "cyan")) + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(dataset ~ mapping*filtering, scales = "free", space = "free") + 
  #facet_wrap(~dataset*mapping*filtering, scales = "free_x") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

beta_coefs_non_meta_df_ALL %>% filter(term == 'LOADAD', class == "Glutamatergic") %>% 
  ggplot(aes(x = dataset, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error) , width = .33) + 
  facet_grid(subclass ~ mapping*filtering, scales = "free", space = "free") + 
  #facet_wrap(~dataset*mapping*filtering, scales = "free_x") + 
  theme_classic() +
  ylab('LOAD (std. Beta)') + 
  xlab('') +
  theme(legend.position = "none") 

#save

write.csv(beta_coefs_non_meta_df_combined, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/LOADAD_lm_results/beta_coefs_non_meta_df_4in1_common_subclasses.csv")

#
#### Revisited mapping (Hodge all regions, expanded subclasses with L3/5 IT) ####

# load Seu objects

Seu_ref_object <- Seu_AIBS_obj

Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds.rds")

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_22MAY21.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_22MAY21.rds") #load zhou seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# prep as needed

table(Seu_ref_object$outlier_call, exclude = "ifany")
table(Seu_ref_object$subclass_label_expanded_L35IT, exclude = "ifany")

Seu_ref_object <- FindVariableFeatures(Seu_ref_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring
length(Seu_ref_object@assays$RNA@var.features) #check
Idents(Seu_ref_object) <- "subclass_label_expanded_L35IT" #so that we're mapping at the subclass level
table(Idents(Seu_ref_object))

test <- Seu_map_object@assays$RNA@var.features
Seu_map_object <- FindVariableFeatures(Seu_map_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferrin
test2 <- Seu_map_object@assays$RNA@var.features
identical(test, test2)
remove(test)
remove(test2)

length(intersect(VariableFeatures(Seu_map_object), VariableFeatures(Seu_ref_object)))
length(intersect(VariableFeatures(Seu_map_object), rownames(Seu_ref_object)))
length(intersect(rownames(Seu_map_object), VariableFeatures(Seu_ref_object)))

#transfer
tanchors <- FindTransferAnchors(reference = Seu_ref_object, query = Seu_map_object, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = Seu_ref_object$subclass_label_expanded_L35IT, dims = 1:30)

metada_to_add <- predictions[, c("predicted.id", "prediction.score.max")]
colnames(metada_to_add) <- paste0(colnames(metada_to_add), ".", "AllHodge_ExpSubclas")
table(metada_to_add$predicted.id.AllHodge_ExpSubclas)

Seu_map_object
Seu_map_object <- AddMetaData(Seu_map_object, metadata = metada_to_add)
table(Seu_map_object$predicted.id.AllHodge_ExpSubclas, exclude = "ifany")
table(Seu_map_object$predicted.id.CgG, Seu_map_object$predicted.id.AllHodge_ExpSubclas, exclude = "ifany")

length(intersect(VariableFeatures(Seu_map_object), tanchors@anchor.features))

# save updated seurat object

saveRDS(Seu_map_object, "~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_11JUN21.rds") 

#### Unused code from Shreejoy ####
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

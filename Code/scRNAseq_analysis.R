#### All involved packages ####

library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)


#### Data loading - Cain et al ####

### read in data
Cain_cell_metadata <- readr::read_csv("/external/rprshnas01/netdata_kcni/dflab/data/rosmap/rnaseq/syn21589957/ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv")
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

#transfer
tanchors <- FindTransferAnchors(reference = new_Seu_AIBS_obj, query = Seu_cain_obj, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = new_Seu_AIBS_obj$subclass_label, dims = 1:30)
Seu_cain_obj <- AddMetaData(Seu_cain_obj, metadata = predictions)

Confusion_matrix <- data.frame(unclass(table(Seu_cain_obj$predicted.id, Seu_cain_obj$subtype)))

length(intersect(VariableFeatures(new_Seu_AIBS_obj), tanchors@anchor.features))

       

#### Proportion analysis (a la Shreejoy) - Cain et al ####

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

#### Proportion analysis (a la Shreejoy) - Zhou et al ####

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

#### Proportion analysis (a la Shreejoy) - Mathys et al ####

# pull metadata with subclass annotations into data frame
mathys_meta_df <- Seu_mathys_obj@meta.data %>% as.data.frame()
mathys_meta_df <- mathys_meta_df %>% rename(subclass = predicted.id)

# read in rosmap metadata from dan's lab folder
ros_meta = readRDS("/external/rprshnas01/public_datasets2/rosmap/phenotype/ROSmaster.rds")

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
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

#### Initial Mathys et al ####

new_Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds")
Seu_mathys_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds")

Idents(new_Seu_AIBS_obj) <- "subclass_label"
Idents(Seu_mathys_obj)

new_Seu_AIBS_obj <- FindVariableFeatures(new_Seu_AIBS_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Seu_mathys_obj <- FindVariableFeatures(Seu_mathys_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

tanchors <- FindTransferAnchors(reference = new_Seu_AIBS_obj, query = Seu_mathys_obj, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = new_Seu_AIBS_obj$subclass_label, dims = 1:30)
Seu_mathys_obj <- AddMetaData(Seu_mathys_obj, metadata = predictions)

test <- data.frame(unclass(table(Seu_mathys_obj$predicted.id, Seu_mathys_obj$matched_group)))

Idents(Seu_mathys_obj) <- "predicted.id"

mathys_markers_mega_mast_anchor <- FindAllMarkers(Seu_mathys_obj, slot = "data", logfc.threshold = 0, min.pct = 0, only.pos = TRUE, return.thresh = 1, test.use = "MAST")

predicted_id_match <- data.frame(unclass(table(Seu_mathys_obj$predicted.id, Seu_mathys_obj$matched_group)))

mathys_markers_selected_mast_anchor <- FindAllMarkers(Seu_mathys_obj, slot = "data", features=intersect(new_CgG_results$gene, row.names(mathys_average)), logfc.threshold = 0, min.pct = 0, return.thresh = 1, test.use = "MAST")

# checks and plots

mathys_markers_mega_mast_anchor$markercheck <- paste0(mathys_markers_mega_mast_anchor$cluster, "_", mathys_markers_mega_mast_anchor$gene)
test <- new_CgG_results
test$markercheck <- paste0(test$subclass, "_", test$gene)
length(intersect(mathys_markers_mega_mast_anchor$markercheck, test$markercheck))

mathys_markers_mega_mast_anchor$volcano_group <- "Not significant"
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val < 0.05, "volcano_group"] <- "Independently significant"
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val_adj < 0.05, "volcano_group"] <- "Significant after correction"

for (cellgroup in unique(new_CgG_results$subclass)){
  temp_markers <- new_CgG_results[new_CgG_results$subclass == cellgroup, "gene"]
  mathys_markers_mega_mast_anchor[((mathys_markers_mega_mast_anchor$cluster == cellgroup) &
                                     (mathys_markers_mega_mast_anchor$gene %in% temp_markers))  
                                  , "volcano_group"] <- "Hodge marker"
} 

table(mathys_markers_mega_mast_anchor$volcano_group)

table(mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$volcano_group == "Hodge marker", "cluster"])
table(new_CgG_results[new_CgG_results$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster),"subclass"])

mathys_average <- AverageExpression(Seu_mathys_obj)
mathys_average <- mathys_average$RNA

intersect(test[(test$gene %in% row.names(mathys_average)) &
                 (test$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster))
               ,"markercheck"], mathys_markers_mega_mast_anchor$markercheck)

setdiff(test[(test$gene %in% row.names(mathys_average)) &
               (test$subclass %in% unique(mathys_markers_mega_mast_anchor$cluster)),"markercheck"], 
        mathys_markers_mega_mast_anchor$markercheck)

mathys_markers_mega_mast_anchor$p_val_adj_forplot <- log10(mathys_markers_mega_mast_anchor$p_val_adj)*(-1)
mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$p_val_adj_forplot == Inf, "p_val_adj_forplot"] <- 400

library(ggplot2)

ggplot(mathys_markers_mega_mast_anchor, aes(x=avg_logFC, y= p_val_adj_forplot, color=volcano_group)) + 
  geom_point(size = .25) +
  facet_wrap(~cluster, scales = "fixed")

ggplot(mathys_markers_mega_mast_anchor[mathys_markers_mega_mast_anchor$volcano_group == "Hodge marker",], aes(x=avg_logFC, y= p_val_adj_forplot, color=volcano_group)) + 
  geom_point(size = 1) +
  facet_wrap(~cluster, scales = "fixed") +
  theme(legend.position = "none")

Seu_mathys_obj$SST_or_Not <- Seu_mathys_obj$predicted.id
Seu_mathys_obj@meta.data[Seu_mathys_obj$SST_or_Not != "SST","SST_or_Not"] <- "Not_SST" #collapsing the non-sst groups
Idents(Seu_mathys_obj) <- "SST_or_Not"
VlnPlot(Seu_mathys_obj, features = "SST", slot = "data")
mathys_average_SST <- AverageExpression(Seu_mathys_obj)
mathys_average_SST <- mathys_average_SST$RNA

test <- t(mathys_average_SST["SST",])
test <- tibble::rownames_to_column(as.data.frame(test))
ggplot(data=test, aes(x=rowname, y=SST, fill=rowname)) + scale_fill_manual(values=c("white", "black")) +
  geom_bar(stat="identity", color = "black") + theme_bw() +theme(axis.title.x=element_blank()) + theme(legend.position = "none")

# gdp

library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis

gdp_plot <- t(as.data.frame(Seu_mathys_obj[["RNA"]]@data)) #get transposed lnCPM matrix
gdp_plot <- t(as.data.frame(new_Seu_AIBS_obj[["RNA"]]@data)) #get transposed lnCPM matrix

# format count/expression matrix for group_dot_plot smooth running
library(tibble)
gdp_plot <- rownames_to_column(as.data.frame(gdp_plot)) #get sample names as a column
colnames(gdp_plot)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_plot) <- gdp_plot$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(Seu_mathys_obj@meta.data) #create metadata copy for group_dot_plot
gdp_anno <- as.data.frame(new_Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot

colnames(gdp_anno)[4] <- "sample_name"
#gdp_anno$subclass_label <- gdp_anno$predicted.id
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"
gdp_anno[gdp_anno$subclass_label == "SST", "subclass_color"] <- "red"

gdp_anno$it2_query_subclass_id <- gdp_anno$it2_query_subclass_label #prepare "_id" needed for plotting
gdp_anno$it2_query_subclass_color <- gdp_anno$class_color #prepare "_color" needed for plotting; copy "class_color" for now
gdp_anno$it2_query_subclass_label <- factor(gdp_anno$it2_query_subclass_label, levels = c("Inh_SST", "Inh_VIP", "Inh_PVALB", "Inh_LAMP5", "Inh_PAX6", "Pyramidal", "Astro_FGFR3", "Oligo",  "OPC_MYT1", "Micro_C1QC", "Peri_MUSTN1", "Endo_CLDN5", "VLMC_CYP1B1")) #setting factor order

# set which genes to plot

gdp_markers <- new_CgG_results[new_CgG_results$subclass == "SST", "gene"] #get SST markers 
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

gdp_markers <- new_CgG_results[new_CgG_results$subclass == "SST", "gene"] #get SST markers
gdp_markers <- intersect(gdp_markers, row.names(mathys_average))
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

gdp_markers <- ACC_results[ACC_results$adapted_cluster_name == "PVALB", "gene"] #get PVALB markers 
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

# do the plot

gdp_anno_test <- merge(gdp_anno, Name_and_colour_scheme, by.x = "subclass_label", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
gdp_anno_test$old_labels <- gdp_anno_test$subclass_label
gdp_anno_test$subclass_label <- gdp_anno_test$Our_label
gdp_anno_test$subclass_label <- factor(gdp_anno_test$subclass_label, levels = c("Exc_IT",
                                                                                "Exc_L5 ET",
                                                                                "Exc_L5/6 IT Car3",
                                                                                "Exc_L5/6 NP",
                                                                                "Exc_L6 CT",  
                                                                                "Exc_L6b",
                                                                                "Inh_SST",
                                                                                "Inh_LAMP5",
                                                                                "Inh_PAX6",
                                                                                "Inh_PVALB",
                                                                                "Inh_VIP",
                                                                                "Astrocyte",
                                                                                "Endothelial",      
                                                                                "Microglia",
                                                                                "Oligodendrocyte",
                                                                                "OPC",
                                                                                "Pericyte",
                                                                                "VLMC"))

group_dot_plot(gdp_plot, 
               gdp_anno_test, 
               genes = gdp_markers, 
               grouping = "subclass", 
               log_scale = TRUE,
               font_size = 14,
               max_size = 29,
               rotate_counts = TRUE)

ggsave(dpi = 300, 
       limitsize = F,
       filename = "Hodge_SST_gdp.jpg",
       device = "jpeg")

ggsave(filename = "L5_6_NP_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)

### systemic gdp ###

#volcano

mathys_markers_selected_mast_anchor$volcano_group <- "Not significant"
mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$p_val < 0.05, "volcano_group"] <- "Significant"

for (cellgroup in unique(new_CgG_results$subclass)){
  temp_markers <- new_CgG_results[new_CgG_results$subclass == cellgroup, "gene"]
  mathys_markers_selected_mast_anchor[((mathys_markers_selected_mast_anchor$cluster == cellgroup) &
                                         (mathys_markers_selected_mast_anchor$gene %in% temp_markers))  
                                      , "volcano_group"] <- "Hodge marker"
} 

table(mathys_markers_selected_mast_anchor$volcano_group)

mathys_markers_selected_mast_anchor$p_val_forplot <- log10(mathys_markers_selected_mast_anchor$p_val)*(-1)
mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$p_val_forplot == Inf, "p_val_forplot"] <- 400

groups_of_interest <- c("SST", "VIP", "Oligodendrocyte", "Microglia")
groups_of_interest <- c("VLMC", "Pericyte", "Endothelial", "Microglia")

ggplot(mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster %in% groups_of_interest,], aes(x=avg_logFC, y= p_val_forplot, color=volcano_group)) + 
  geom_point(size = .75) +
  facet_wrap(~cluster, scales = "fixed") +
  scale_color_manual(values=c("red", "black", "orange")) +
  theme_bw() 

ggplot(mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster %in% groups_of_interest
                                           & mathys_markers_selected_mast_anchor$volcano_group == "Hodge marker",], aes(x=avg_logFC, y= p_val_forplot, color=volcano_group)) + 
  geom_point(size = .75) +
  facet_wrap(~cluster, scales = "fixed") +
  scale_color_manual(values=c("red", "black", "orange")) +
  theme_bw() +
  theme(legend.position = "none")

### log fc plot

mathys_markers_selected_mast_anchor$clust_gene <- paste0(mathys_markers_selected_mast_anchor$cluster,"_",mathys_markers_selected_mast_anchor$gene)
test2 <- new_CgG_results
test2$clust_gene <- paste0(test2$subclass,"_",test2$gene)
test <- merge(test2, mathys_markers_selected_mast_anchor, by = "clust_gene")

names(test)[c(8,14)] <- c("Hodge_logFC", "Mathys_logFC")

ggplot(test[test$cluster %in% groups_of_interest,], aes(x=Mathys_logFC, y= Hodge_logFC)) + 
  geom_point(size = .75) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~cluster, scales = "fixed") +
  theme_bw() 

test3 <- mathys_markers_selected_mast_anchor[mathys_markers_selected_mast_anchor$cluster == "SST",]

library(metap)

for (cellgroup in unique(test$subclass)){
  print(cellgroup)
  print(sumlog(test[test$subclass == cellgroup,"p_val"]))
}

source('https://raw.githubusercontent.com/yaowuliu/ACAT/master/R/ACAT.R')

for (cellgroup in unique(test$subclass)){
  print(cellgroup)
  print(ACAT(test[test$subclass == cellgroup,"p_val"]))
}

### More formal proportion analysis - closer to what is used for initial Zhou and Cain

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


#### Initial Cain et al ####

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

### proportion analysis

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
#### Initial Zhou et al ####

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

### Proportion analysis

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

#### Revisited Proportion analysis (filter IT based on prediction scores) ####

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


#### Further Revisited Proportion analysis (filter None.NA for Cain only), otherwise same as above ####

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



#### Revisited mapping (Hodge CgG neurons only) ####

# load objects
Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds")

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

#### Revisited mapping (AIBS M1 dataset reference) ####

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

#### Proportion analysis revisited - All datasets - M1 & CgG, IT prediction score unfiltered &filtered ####

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
#### Original identity proportion (snCTP) analysis ####

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

#### Analyses with Hodge all regions, expanded subclasses with L3/5 IT ####

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

### Proportion analysis

### load objects
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_08Jun21.rds")
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_11JUN21.rds") #load cain seurat object (instead)
Seu_prop_obj <- subset(Seu_prop_obj, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_prop_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_11JUN21.rds") #load zhou seurat object (instead)
Seu_prop_obj <- subset(Seu_prop_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# pull metadata with subclass annotations into data frame
meta_df <- Seu_prop_obj@meta.data %>% as.data.frame()
remove(Seu_prop_obj)
#meta_df_backup <- meta_df
#meta_df <- meta_df_backup
#meta_df <- meta_df[!(meta_df$predicted.id.CgG == "IT" & meta_df$prediction.score.max.CgG < 0.8),] # as appropriate
#meta_df <- meta_df[!(meta_df$predicted.id.M1 == "L2/3 IT" & meta_df$prediction.score.max.M1 < 0.8),] # as appropriate
meta_df <- meta_df %>% rename(subclass = predicted.id.AllHodge_ExpSubclas) # as appropriate

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

write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/AllHodge_Unfiltered/cain_cell_type_prop_df.csv")

#

#### Revisited mapping - CgG and MTG, original AIBS subclasses,20% var features, two-threshold QC ####

# load Seu objects

Seu_ref_object <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds")

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_08Jun21.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_11JUN21.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_11JUN21.rds") #load zhou seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

# prep as needed

table(Seu_ref_object$NeuN_Region, exclude = "ifany")
Seu_ref_object
Seu_ref_object <- subset(Seu_ref_object, 
                         subset = NeuN_Region %in% c("A1C_Neuronal",
                                                     "M1lm_Neuronal",
                                                     "M1ul_Neuronal",
                                                     "S1lm_Neuronal",
                                                     "S1ul_Neuronal",
                                                     "V1C_Neuronal"), 
                         invert = TRUE)
table(Seu_ref_object$NeuN_Region, exclude = "ifany")
Seu_ref_object

table(Seu_ref_object$outlier_call, exclude = "ifany")

table(Seu_ref_object$subclass_label, exclude = "ifany")

# run mapping/predictions

Seu_ref_object <- FindVariableFeatures(Seu_ref_object, selection.method = "vst", nfeatures = round(nrow(Seu_ref_object)/5), verbose = FALSE) #need variable features for transferring
length(Seu_ref_object@assays$RNA@var.features) #check

Seu_map_object <- FindVariableFeatures(Seu_map_object, selection.method = "vst", nfeatures = round(nrow(Seu_map_object)/5), verbose = FALSE) #need variable features for transferrin
length(Seu_map_object@assays$RNA@var.features) #check

length(intersect(VariableFeatures(Seu_map_object), VariableFeatures(Seu_ref_object)))
length(intersect(VariableFeatures(Seu_map_object), rownames(Seu_ref_object)))
length(intersect(rownames(Seu_map_object), VariableFeatures(Seu_ref_object)))

# transfer + predict

tanchors <- FindTransferAnchors(reference = Seu_ref_object, query = Seu_map_object, dims = 1:30)
predictions <- TransferData(anchorset = tanchors, refdata = Seu_ref_object$subclass_label, dims = 1:30)

# pass-fail mapping QC

prediction_score_threshold_low <- 0.5
prediction_score_threshold_high <- 0.8
bad_cell_type <- "IT"

predictions <- predictions %>% mutate(qc_passing = case_when((predicted.id != bad_cell_type) & (prediction.score.max > prediction_score_threshold_low) ~ T,
                                                             (predicted.id == bad_cell_type) & (prediction.score.max > prediction_score_threshold_high) ~ T,
                                                             TRUE ~ F)) 

# save the predictions df

write.csv(predictions, "~/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Mapping_full/predictions_Zhou_MTGnCgG_20PctMostVar.csv")

# add metadata to Seu Object

metada_to_add <- predictions[, c("predicted.id", "prediction.score.max", "qc_passing")]
colnames(metada_to_add) <- paste0(colnames(metada_to_add), ".", "MTGnCgG_20Pct")
table(metada_to_add$predicted.id.MTGnCgG_20Pct)

Seu_map_object
Seu_map_object <- AddMetaData(Seu_map_object, metadata = metada_to_add)
table(Seu_map_object$predicted.id.MTGnCgG_20Pct, exclude = "ifany")
table(Seu_map_object$predicted.id, Seu_map_object$predicted.id.MTGnCgG_20Pct, exclude = "ifany")

length(intersect(VariableFeatures(Seu_map_object), tanchors@anchor.features))

# save updated seurat object

saveRDS(Seu_map_object, "~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_27JUL21.rds") 

### Quantify

# pull metadata with subclass annotations into data frame
meta_df <- Seu_map_object@meta.data %>% as.data.frame()
remove(Seu_map_object)
meta_df_backup <- meta_df
meta_df <- meta_df[meta_df$qc_passing.MTGnCgG_20Pct == T,] # as appropriate
meta_df <- meta_df %>% rename(subclass = predicted.id.MTGnCgG_20Pct) # as appropriate

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

ggplot(cell_prop_meta_long[cell_prop_meta_long$subclass == "IT",], aes(x = LOAD, y = cell_type_prop, fill = LOAD)) + 
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

write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/MTGnCgG_20Pct/Unfiltered/zhou_cell_type_prop_df.csv")
write.csv(export_holder, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/sn_cell_type_proportions/MTGnCgG_20Pct/Filtered/zhou_cell_type_prop_df.csv")


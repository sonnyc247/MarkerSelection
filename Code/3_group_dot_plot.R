#### packages ####

library(Seurat)
library(tibble)
library(scrattch.vis)
options(stringsAsFactors = T) # following https://github.com/AllenInstitute/scrattch.vis
library(cowplot)

#### template/pilot/main code ####

### load data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds") #load hodge (mtg and cgg neurons) seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj_update_11JUN21.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

### extract and format data for plotting

#cpm (can modify for counts)
gene_list <- new_CgG_results_class$gene # get list of genes
gene_list <- intersect(gene_list, row.names(Seu_plot_object)) # intersect genes with Seurat object for plotting
gdp_data <- GetAssayData(Seu_plot_object, slot = "data")[gene_list,] # get lnCPM matrix
gdp_data <- as.matrix(gdp_data) # convert to non-sparse matrix
gdp_data <- t(gdp_data) #get transposed lnCPM matrix
gdp_data <- rownames_to_column(as.data.frame(gdp_data)) #get sample names as a column
colnames(gdp_data)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_data) <- gdp_data$sample_name #reset df rownames as sample names as well, just in case

#metadata/cell annotations
gdp_anno <- as.data.frame(Seu_plot_object@meta.data) #create metadata copy for group_dot_plot
gdp_anno <- rownames_to_column(as.data.frame(gdp_anno)) #get sample names as a column
colnames(gdp_anno)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_anno) <- gdp_anno$sample_name #reset df rownames as sample names as well, just in case

#the resolution we want to plot

gdp_anno$subclass_label <- gdp_anno$predicted.id.AllHodge_ExpSubclas #as needed/appropriate
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"
#gdp_anno[gdp_anno$subclass_label == "SST", "subclass_color"] <- "red" #as needed

# set which genes to plot

group_to_gdp <- "Glutamatergic"
gdp_markers <- new_CgG_results_class[new_CgG_results_class$class == group_to_gdp, "gene"] #get markers 
gdp_markers <- intersect(gdp_markers$gene, colnames(gdp_data))
#gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

#overlap_markers_SST <- gdp_markers
#overlap_markers_IT <- gdp_markers

# final adjustments of plotting df

#gdp_plot <- merge(gdp_anno, Name_and_colour_scheme, by.x = "subclass_label", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
#gdp_plot$old_labels <- gdp_plot$subclass_label
#gdp_plot$subclass_label <- gdp_plot$Our_label
#gdp_plot$subclass_label <- factor(gdp_plot$subclass_label, levels = c("Exc_IT",
#                                                                      "Exc_L5 ET",
#                                                                      "Exc_L5/6 IT Car3",
#                                                                      "Exc_L5/6 NP",
#                                                                      "Exc_L6 CT",  
#                                                                      "Exc_L6b",
#                                                                      "Inh_SST",
#                                                                      "Inh_LAMP5",
#                                                                      "Inh_PAX6",
#                                                                      "Inh_PVALB",
#                                                                      "Inh_VIP",
#                                                                      "Astrocyte",
#                                                                      "Endothelial",      
#                                                                      "Microglia",
#                                                                      "Oligodendrocyte",
#                                                                      "OPC",
#                                                                      "Pericyte",
#                                                                      "VLMC"))

#gdp_plot$subclass_label <- gdp_plot$old_labels

# do the plot

gdp_gluta <- group_dot_plot(gdp_data, 
               gdp_anno, 
               genes = gdp_markers, 
               grouping = "subclass", 
               log_scale = TRUE,
               font_size = 7,
               max_size = 10,
               rotate_counts = TRUE)

# combine plots as appropriate/needed

plot_grid(gdp_gluta,
          gdp_gaba, 
          labels = c("Gluta", "Gaba"),
          align = "hv",
          ncol = 1)

# export plot (most recently plotted)

ggsave(width = 180, 
       height = 250,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/GDPs/Final/",
       filename = "group_dot_plot.pdf",
       device = "pdf")

#### systemic hodge gdp ####

library(scrattch.vis)
options(stringsAsFactors = F) # following https://github.com/AllenInstitute/scrattch.vis

gdp_plot <- t(as.data.frame(new_Seu_AIBS_obj[["RNA"]]@data)) #get transposed lnCPM matrix

# format count/expression matrix for group_dot_plot smooth running
library(tibble)
gdp_plot <- rownames_to_column(as.data.frame(gdp_plot)) #get sample names as a column
colnames(gdp_plot)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_plot) <- gdp_plot$sample_name #reset df rownames as sample names as well, just in case

# further adding and tweeking data in metadata dataframe to suite group_dot_plot
gdp_anno <- as.data.frame(new_Seu_AIBS_obj@meta.data) #create metadata copy for group_dot_plot

#colnames(gdp_anno)[4] <- "sample_name"
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"

# set which genes to plot

for (cellgroup in names(humanMarkersCommon)[-5]) {
        gdp_anno[gdp_anno$subclass_label == cellgroup, "subclass_color"] <- "red"
        temp_markers <- unlist(humanMarkersCommon[cellgroup])
        gdp_markers <- new_CgG_results_pairdata[new_CgG_results_pairdata$gene %in% temp_markers, c("gene", "roc_power")] #get group markers 
        gdp_markers <- gdp_markers[order(gdp_markers$roc_power, decreasing = T),]
        
        group_dot_plot(gdp_plot, 
                       gdp_anno, 
                       genes = gdp_markers$gene, 
                       grouping = "subclass", 
                       log_scale = TRUE,
                       font_size = 10,
                       max_size = 20,
                       rotate_counts = TRUE)
        
        if (cellgroup == "L5/6 NP"){
                ggsave(filename = "L5_6_NP_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)
        } else if (cellgroup == "L5/6 IT Car3"){
                ggsave(filename = "L5_6_IT_Car3_gdp.pdf", path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)  
        } else{
                ggsave(filename = paste0(cellgroup,"_gdp.pdf"), path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/GDPs", height = (7 + (nrow(gdp_markers))/3), limitsize =  F)  
        }
        
        gdp_anno[gdp_anno$subclass_label == cellgroup, "subclass_color"] <- "white"
        
}
#### for publication ####

### load data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds") #load hodge (mtg and cgg neurons) seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj_update_27JUL21.rds") #load cain seurat object (instead)
Seu_plot_object <- subset(Seu_plot_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only

marker_reference <- read.csv("/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_usen.csv")
marker_reference <- marker_reference[,-c(1,3)]

### extract and format data for plotting

#cpm (can modify for counts)
gene_list <- as.character(marker_reference[marker_reference$subclass %in% c("SST", "IT") &
                              marker_reference$used == T, "gene"]) # get list of genes
#gene_list <- intersect(gene_list, row.names(Seu_plot_object)) # intersect genes with Seurat object for plotting
#gene_list <- gene_list[-3] # remove CORT which isn't in Cain
#gene_list[8] <- "LY86-AS1"
gdp_data <- GetAssayData(Seu_plot_object, slot = "data")[gene_list,] # get lnCPM matrix
gdp_data <- as.matrix(gdp_data) # convert to non-sparse matrix
gdp_data <- t(gdp_data) #get transposed lnCPM matrix
gdp_data <- rownames_to_column(as.data.frame(gdp_data)) #get sample names as a column
colnames(gdp_data)[1] <- "sample_names" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_data) <- gdp_data$sample_name #reset df rownames as sample names as well, just in case
#colnames(gdp_data)[9] <- "LY86.AS1"

#metadata/cell annotations
gdp_anno <- as.data.frame(Seu_plot_object@meta.data) #create metadata copy for group_dot_plot
gdp_anno <- rownames_to_column(as.data.frame(gdp_anno)) #get sample names as a column
colnames(gdp_anno)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_anno) <- gdp_anno$sample_name #reset df rownames as sample names as well, just in case

#the resolution we want to plot

gdp_anno$subclass_label <- gdp_anno$subclass_label
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"
#gdp_anno[gdp_anno$subclass_label == "SST", "subclass_color"] <- "red" #as needed

gdp_anno$subclass_label <- gdp_anno$predicted.id.MTGnCgG_20Pct
gdp_anno$subclass_id <- gdp_anno$predicted.id.MTGnCgG_20Pct
gdp_anno$subclass_color <- "white"


# set which genes to plot

group_to_gdp <- "SST"
gdp_markers <- marker_reference[marker_reference$subclass == group_to_gdp &
                                marker_reference$used == T, c("gene", "avg_logFC")] #get markers 
gdp_markers <- as.character(gdp_markers[order(gdp_markers$avg_logFC, decreasing = T),"gene"])
gdp_markers <- intersect(gdp_markers, colnames(gdp_data))
#gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

#overlap_markers_SST <- gdp_markers
#overlap_markers_IT <- gdp_markers

# final adjustments of plotting df

#gdp_plot <- merge(gdp_anno, Name_and_colour_scheme, by.x = "subclass_label", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
#gdp_plot$old_labels <- gdp_plot$subclass_label
#gdp_plot$subclass_label <- gdp_plot$Our_label
#gdp_plot$subclass_label <- factor(gdp_plot$subclass_label, levels = c("Exc_IT",
#                                                                      "Exc_L5 ET",
#                                                                      "Exc_L5/6 IT Car3",
#                                                                      "Exc_L5/6 NP",
#                                                                      "Exc_L6 CT",  
#                                                                      "Exc_L6b",
#                                                                      "Inh_SST",
#                                                                      "Inh_LAMP5",
#                                                                      "Inh_PAX6",
#                                                                      "Inh_PVALB",
#                                                                      "Inh_VIP",
#                                                                      "Astrocyte",
#                                                                      "Endothelial",      
#                                                                      "Microglia",
#                                                                      "Oligodendrocyte",
#                                                                      "OPC",
#                                                                      "Pericyte",
#                                                                      "VLMC"))

#gdp_plot$subclass_label <- gdp_plot$old_labels

# do the plot

gdp_SST_AIBS <- group_dot_plot(gdp_data, 
                               gdp_anno,
                               genes = gdp_markers,
                               grouping = "subclass",
                               log_scale = TRUE,
                               font_size = 8,
                               max_size = 20,
                               rotate_counts = TRUE)

# combine plots as appropriate/needed

plot_grid(gdp_SST,
          gdp_IT, 
          labels = c("A", "B"),
          align = "hv",
          ncol = 1)

# export plot (most recently plotted)

ggsave(width = 180, 
       height = 250,
       dpi = 300, 
       units = "mm", 
       limitsize = F,
       path = "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/Figures/GDPs/Final/",
       filename = "Cain_group_dot_plot.pdf",
       device = "pdf")

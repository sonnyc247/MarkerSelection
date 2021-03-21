### packages

library(Seurat)
library(tibble)
library(scrattch.vis)
options(stringsAsFactors = T) # following https://github.com/AllenInstitute/scrattch.vis
library(cowplot)

### load data

Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/new_Seu_AIBS_obj.rds") #load hodge (mtg and cgg neurons) seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj.rds") #load mathys seurat object
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_plot_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_plot_object <- subset(Seu_test_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

### extract and format data for plotting

#cpm (can modify for counts)
gdp_data <- t(as.data.frame(Seu_plot_object[["RNA"]]@data)) #get transposed lnCPM matrix
gdp_data <- rownames_to_column(as.data.frame(gdp_data)) #get sample names as a column
colnames(gdp_data)[1] <- "sample_name" #change column name of sample names (to match AIBS vignette on group_dot_plot)
rownames(gdp_data) <- gdp_data$sample_name #reset df rownames as sample names as well, just in case

#metadata/cell annotations
gdp_anno <- as.data.frame(Seu_plot_object@meta.data) #create metadata copy for group_dot_plot
colnames(gdp_anno)[4] <- "sample_name"
gdp_anno$subclass_label <- gdp_anno$predicted.id #as needed/appropriate
gdp_anno$subclass_id <- gdp_anno$subclass_label
gdp_anno$subclass_color <- "white"
#gdp_anno[gdp_anno$subclass_label == "SST", "subclass_color"] <- "red" #as needed

# set which genes to plot

group_to_gdp <- "IT"
gdp_markers <- Result_df_MTGandCgG_lfct2.0[Result_df_MTGandCgG_lfct2.0$subclass == group_to_gdp, "gene"] #get markers 
gdp_markers <- intersect(gdp_markers, colnames(gdp_data))
gdp_markers <- sort(gdp_markers[!is.na(gdp_markers)]) #remove NA, sort alphabetically

overlap_markers_SST <- gdp_markers
overlap_markers_IT <- gdp_markers

# final adjustments of plotting df

gdp_plot <- merge(gdp_anno, Name_and_colour_scheme, by.x = "subclass_label", by.y = "AIBS_subclass_label", all.x = T, all.y = F)
gdp_plot$old_labels <- gdp_plot$subclass_label
gdp_plot$subclass_label <- gdp_plot$Our_label
gdp_plot$subclass_label <- factor(gdp_plot$subclass_label, levels = c("Exc_IT",
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

gdp_plot$subclass_label <- gdp_plot$old_labels

# do the plot

gdp_mathys_IT <- group_dot_plot(gdp_data, 
               gdp_plot, 
               genes = overlap_markers_IT, 
               grouping = "subclass", 
               log_scale = TRUE,
               font_size = 7,
               max_size = 10,
               rotate_counts = TRUE)

# combine plots as appropriate/needed

plot_grid(gdp_hodge_IT,
          gdp_mathys_IT, 
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
       filename = "IT_gdp_h250.pdf",
       device = "pdf")

### THIS IS AN INCOMPLETE/UNREFINED FIRST-ATTEMPT AT GROUP-TO-GROUP CELL IDENTITY MAPPING ###
## THINGS TO DO INCLUDE: UPDATE/CONSOLIDATE/MAXIMIZE INTERSECTING GENE NAMES; MORE THOUGHTFUL SELECTING OF GENES/MARKERS TO USE TO MAP ##

#### packages ####

library(Seurat)

#### load data ####

Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds")

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_08Jun21.rds") #load mathys seurat object
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_cain_obj.rds") #load cain seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = subtype == "None.NA", invert = TRUE) # for cain object only
Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_zhou_obj.rds") #load zhou seurat object (instead)
Seu_map_object <- subset(Seu_map_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500) #for zhou object only

#### get averages ####

# for reference dataset
Idents(Seu_ref_object) <- "subclass_label_expanded_L35IT"

ref_avg_expression <- AverageExpression(Seu_ref_object, slot = "data")
ref_avg_expression <- ref_avg_expression$RNA

remove(Seu_ref_object)

# for mapping dataset
Idents(Seu_map_object) <- "Subcluster"

map_avg_expression <- AverageExpression(Seu_map_object, slot = "data")
map_avg_expression <- map_avg_expression$RNA

remove(Seu_map_object)

#### select intersecting genes ####

common_genes <- intersect(row.names(map_avg_expression), 
                          row.names(ref_avg_expression))

ref_avg_exp_common <- ref_avg_expression[common_genes,]
map_avg_expression <- map_avg_expression[common_genes,]

#### get correlations ####

group_cor_results <- cor(x = ref_avg_exp_common,
                         y = map_avg_expression)

#### visualize ####

heatmap(group_cor_results)

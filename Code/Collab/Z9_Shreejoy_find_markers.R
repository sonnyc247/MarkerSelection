#### get markers from all of hodge with new IT groups ####

### initial setup

library(Seurat)
library(MAST)
library(dplyr)

Seu_AIBS_obj <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj.rds")
table(Seu_AIBS_obj$outlier_call, exclude = "ifany") #check for outliers
table(Seu_AIBS_obj$NeuN_Region) #check our data composition

### make new IT groups

Seu_AIBS_obj$subclass_label_expanded <- Seu_AIBS_obj$subclass_label #make new variable, starting from subclasses

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 062",
                                                                                     "Neuron 063",
                                                                                     "Neuron 064",
                                                                                     "Neuron 065",
                                                                                     "Neuron 066",
                                                                                     "Neuron 067",
                                                                                     "Neuron 068",
                                                                                     "Neuron 069",
                                                                                     "Neuron 070",
                                                                                     "Neuron 071",
                                                                                     "Neuron 072")] <- "L2/3 IT" #set L2/3 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 073",
                                                                                     "Neuron 074",
                                                                                     "Neuron 075",
                                                                                     "Neuron 076")] <- "L6 IT" #set L6 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 079",
                                                                                     "Neuron 080",
                                                                                     "Neuron 081",
                                                                                     "Neuron 082",
                                                                                     "Neuron 083",
                                                                                     "Neuron 084",
                                                                                     "Neuron 085",
                                                                                     "Neuron 086",
                                                                                     "Neuron 087")] <- "L5 IT" #set L5 IT

Seu_AIBS_obj$subclass_label_expanded[Seu_AIBS_obj$cell_type_designation_label %in% c("Neuron 077",
                                                                                     "Neuron 078")] <- "Other IT" #set Leftover IT

table(Seu_AIBS_obj$subclass_label_expanded)
Idents(Seu_AIBS_obj) <- "subclass_label_expanded"
table(Idents(Seu_AIBS_obj)) #double check what subclasses we have and that they're set as active identity

### subset as needed ###

Idents(Seu_AIBS_obj) <- "NeuN_Region" #we identify samples by "NeuN_Region", which we made befire in 1_Loading_and_selecting_data.Rmd
Seu_AIBS_obj <- subset(Seu_AIBS_obj, subset = NeuN_Region %in% c("MTG_Neuronal", 
                                                                 "V1C_Neuronal", 
                                                                 "M1lm_Neuronal", 
                                                                 "S1ul_Neuronal", 
                                                                 "S1lm_Neuronal", 
                                                                 "M1ul_Neuronal", 
                                                                 "A1C_Neuronal"), invert = TRUE) #remove all non-CgG neurons

table(Seu_AIBS_obj$subclass_label_expanded) #double check what subclasses we have
Seu_AIBS_obj <- subset(Seu_AIBS_obj, subset = subclass_label_expanded == "L4 IT", invert = TRUE) #remove L4 IT for n = 1

### find markers

Idents(Seu_AIBS_obj) <- "subclass_label_expanded" #assign proper labels

new_AIBS_markers_mast_expIT_ALL <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 2, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #find markers
new_AIBS_markers_roc_expIT_ALL <- FindAllMarkers(Seu_AIBS_obj, slot = "data", logfc.threshold = 2, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #find markers using roc

### find markers for one group (vs all others as one gorup together, NOT vs a specific one other group)

marker_results <- FindMarkers(Seu_AIBS_obj, ident.1 = "SST", logfc.threshold = 2, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") # for mast
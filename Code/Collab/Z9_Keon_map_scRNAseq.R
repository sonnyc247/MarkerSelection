Seu_ref_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_AIBS_obj_update_07JUN21.rds.rds") # get seurat reference object

Seu_map_object <- readRDS("~/git/Ex_Env_Storage/MarkerSelection/Seu_mathys_obj_update_22MAY21.rds") # get the sc RNAseq object you want somehow


# prep as needed

Seu_ref_object <- FindVariableFeatures(Seu_ref_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #need variable features for transferring

Seu_map_object <- FindVariableFeatures(Seu_map_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) #may not nood, but better to have


#transfer
tanchors <- FindTransferAnchors(reference = Seu_ref_object, query = Seu_map_object, dims = 1:30) # find anchors for mapping/transfering/integrating data
predictions <- TransferData(anchorset = tanchors, 
                            refdata = Seu_ref_object$subclass_label_expanded_L35IT, #this is the key part of specifying what class/subclass/cluster labels you want to map onto
                            dims = 1:30)

View(predictions) #notice the biga** dataframe, you may not want to keep/use all of it
metada_to_add <- predictions[, c("predicted.id", "prediction.score.max")] #these two columns are probably what you want the most
colnames(metada_to_add) <- paste0(colnames(metada_to_add), ".", "AllHodge_ExpSubclas") #rename columns as desired

Seu_map_object <- AddMetaData(Seu_map_object, metadata = metada_to_add) # add new metadata (in this case, the mapped ids)

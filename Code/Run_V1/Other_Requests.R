### try just finding markers for inhibitory vs excitatory (vs non-neuronal cells)

Seu_AIBS_obj_ACC_it2 <- subset(Seu_AIBS_obj_ACC, idents = c("Inh_ADARB2", "Inh_LHX6"), invert = TRUE) #we want all cells except for these two groups, to be comparable with previous analyses; but can re-include these cells if desired, since they are inhibitory cells, after all

unique(Idents(Seu_AIBS_obj_ACC_it2)) #see identity

Seu_AIBS_obj_ACC_it2$InExGroups <- Idents(Seu_AIBS_obj_ACC_it2) #make new variable/grouping factor of inhibitory vs excitatory

levels(Seu_AIBS_obj_ACC_it2$InExGroups)[levels(Seu_AIBS_obj_ACC_it2$InExGroups)=="Inh_LAMP5"] <- "Inhibitory" #collapsing the inh groups
levels(Seu_AIBS_obj_ACC_it2$InExGroups)[levels(Seu_AIBS_obj_ACC_it2$InExGroups)=="Inh_PVALB"] <- "Inhibitory" #collapsing the inh groups
levels(Seu_AIBS_obj_ACC_it2$InExGroups)[levels(Seu_AIBS_obj_ACC_it2$InExGroups)=="Inh_SST"] <- "Inhibitory" #collapsing the inh groups
levels(Seu_AIBS_obj_ACC_it2$InExGroups)[levels(Seu_AIBS_obj_ACC_it2$InExGroups)=="Inh_VIP"] <- "Inhibitory" #collapsing the inh groups
levels(Seu_AIBS_obj_ACC_it2$InExGroups)[levels(Seu_AIBS_obj_ACC_it2$InExGroups)=="Inh_PAX6"] <- "Inhibitory" #collapsing the inh groups

unique(Seu_AIBS_obj_ACC_it2$InExGroups) #see groups
Idents(Seu_AIBS_obj_ACC_it2) <- "InExGroups" #assign Inhibitory vs Excitatory as grouping variables/active identities

InEx_ACC_AIBS_markers <- FindAllMarkers(Seu_AIBS_obj_ACC_it2, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #same thresholds as before, using MAST
InEx_ACC_AIBS_markers_roc <- FindAllMarkers(Seu_AIBS_obj_ACC_it2, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #same thresholds as before, using roc

### try just finding markers for neurons (vs non-neuronal cells)

Seu_AIBS_obj_ACC_it2$NeuNonnGroups <- Idents(Seu_AIBS_obj_ACC_it2) #make new variable/grouping factor of neurons

levels(Seu_AIBS_obj_ACC_it2$NeuNonnGroups)[levels(Seu_AIBS_obj_ACC_it2$NeuNonnGroups)=="Inhibitory"] <- "Neuron" #collapsing the inh groups
levels(Seu_AIBS_obj_ACC_it2$NeuNonnGroups)[levels(Seu_AIBS_obj_ACC_it2$NeuNonnGroups)=="Pyramidal"] <- "Neuron" #collapsing the excitatory/pyramidal groups

unique(Seu_AIBS_obj_ACC_it2$NeuNonnGroups) #see groups
Idents(Seu_AIBS_obj_ACC_it2) <- "NeuNonnGroups" #assign these group designations as active identity

NeuNonn_ACC_AIBS_markers <- FindAllMarkers(Seu_AIBS_obj_ACC_it2, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #same thresholds as before, using MAST
NeuNonn_ACC_AIBS_markers_roc <- FindAllMarkers(Seu_AIBS_obj_ACC_it2, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #same thresholds as before, using roc

### some quick comparisons/checks

Seu_AIBS_markers_multicomp <- InEx_ACC_AIBS_markers_roc # to compare vs ACC roc
Seu_AIBS_markers_multicomp <- NeuNonn_ACC_AIBS_markers_roc # to compare vs MTG roc

for (cluster_name in unique(InEx_ACC_AIBS_markers$cluster)) { #for each group (subclass, in this case)
  
  ref_num <- length(unique(Seu_AIBS_markers_multicomp[Seu_AIBS_markers_multicomp$cluster == cluster_name, "gene"])) #number of marker genes in reference
  res_num <- length(InEx_ACC_AIBS_markers[InEx_ACC_AIBS_markers$cluster == cluster_name, "gene"]) #number of marker genes in results
  inter_num <- length(intersect(unique(Seu_AIBS_markers_multicomp[Seu_AIBS_markers_multicomp$cluster == cluster_name, "gene"]), InEx_ACC_AIBS_markers[InEx_ACC_AIBS_markers$cluster == cluster_name, "gene"])) #number of overlapping marker genes
  
  print(c(cluster_name, ref_num, res_num, inter_num)) #print the results
  
  print("") #break
  
  print(intersect(Seu_AIBS_markers_multicomp[Seu_AIBS_markers_multicomp$cluster == cluster_name, "gene"] ,InEx_ACC_AIBS_markers[InEx_ACC_AIBS_markers$cluster == cluster_name, "gene"])) #print list of overlapping marker genes, if desired
  
  print("") #break
  
}

duplicated_markers <- Seu_AIBS_markers_multicomp[duplicated(Seu_AIBS_markers_multicomp$gene),"gene"] #get list of duplicated markers from the roc test
Seu_AIBS_markers_multicomp[Seu_AIBS_markers_multicomp$gene %in% duplicated_markers,] #see the rows with duplicated markers
roc_df <- Seu_AIBS_markers_multicomp[!(Seu_AIBS_markers_multicomp$gene %in% duplicated_markers),] #remove duplicated markers
Result_df <- merge(roc_df[,c(7,6,4,5,2,1,3)], NeuNonn_ACC_AIBS_markers[,c(1,5,7)], by = "gene", all.x = TRUE, all.y = FALSE) #combine the marker df; may want to double check the indices/order
Result_df <- Result_df[,c(1:5,8,9,6,7)] #reorder columns; again, careful to run this once + may want to check this is the order you want again
colnames(Result_df)[c(5:9)] <- c("avg_logFC","MAST_p_val","MAST_p_val_adj", "roc_myAUC", "roc_power") #rename columns to clarify which test generated which statistic

### finally, we want to add entrez IDs to the output/result df
# our mgi-to-entrez ID conversion comes from a file from AIBS http://celltypes.brain-map.org/api/v2/well_known_file_download/694416044

human_MTG_2018_06_14_genes_rows <- read_csv("~/git/MarkerSelection/Data/Inputs/human_MTG_2018-06-14_genes-rows.csv") #loading the conversion data
Result_df <- merge(Result_df, human_MTG_2018_06_14_genes_rows[,c(1,3)], by = "gene") #add entrez ids
Result_df <- Result_df[,c(1,10,2:9)] #reorder columns

ACC_InEx_results <- Result_df #store the results; as appropriate
ACC_NeuNonn_results <- Result_df #store the results; store the results; as appropriate

write.csv(ACC_results, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/ACC_results.csv") #save results/output
write.csv(MTG_results, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/MTG_results.csv") #save results/output
write.csv(ACC_InEx_results, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/ACC_InEx_results.csv") #save results/output
write.csv(ACC_NeuNonn_results, "/external/rprshnas01/kcni/ychen/git/MarkerSelection/Data/Outputs/ACC_NeuNonn_results.csv") #save results/output

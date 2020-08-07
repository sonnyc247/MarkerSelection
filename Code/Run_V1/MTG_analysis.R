### This script replicates Loading_AIBS_Data and Seurat_implementation for MTG instead of ACC; there should be no steps/code here not part of the more generalizable/better documented steps in the notebooks

# Reset filtered metadata
library(tidyr) #needed for the separate function
tome_sample_meta_filtered <- separate(data = tome_sample_meta, col = cell_type_alias_label, into = c("split_class", "split_layer", "split_subclass", "split_type"), sep = " ") #this splits the one column of cell identifiers into cell class, layer, subclass, and type
tome_sample_meta_filtered <- as.data.frame(tome_sample_meta_filtered) #fully convert metada df to df structure (else will get tibble warning message with the next line of code)
row.names(tome_sample_meta_filtered) <- tome_sample_meta_filtered$sample_name #sometimes the workflow requires rownames to be the sample names

# Filter metadata
tome_sample_meta_filtered <- tome_sample_meta_filtered[!(tome_sample_meta_filtered$class_label == "Exclude"),] # this is a merged-grouping of all Outlier and Donor split_class values, which we remove from the dataset
tome_sample_meta_filtered <- tome_sample_meta_filtered[!(tome_sample_meta_filtered$class_label!="Non-neuronal" & tome_sample_meta_filtered$region_label!="MTG"),] # for our purposes/demonstration, we remove all neurons not from the MTG
query_sample_name_list <- tome_sample_meta_filtered$sample_name #extracting sample names from our filtered dataframe

# Load the data
library(scrattch.io)
h5closeAll() #need to run, to prevent red wall of text when reading tome data (run this every time before every use of read_tome_sample_data)
AIBS_Rawcount_InEx_ACC <- AIBS_Rawcount_InEx #store the old count matrix
AIBS_Rawcount_InEx <- read_tome_sample_data(tome, query_sample_name_list, regions = "both", units = "counts", transform = "none", format = "data.frame") #this generates the actual count dataframe in long form, we are getting Introns and Exons (regions = "both"), raw counts without any transformation, in a df format
row.names(AIBS_Rawcount_InEx) <- AIBS_Rawcount_InEx$gene_name #set row names
AIBS_Rawcount_InEx <- AIBS_Rawcount_InEx[,2:length(colnames(AIBS_Rawcount_InEx))] #remove gene name column from df

# Into Seurat
Seu_AIBS_obj_ACC <- Seu_AIBS_obj #save a copy of previous Seurat object (key: this object still has "Inh_LHX6" and "Inh_ADARB2" cells; will want to remove to replicate previous analysis/compare with MTG analysis below)
remove(Seu_AIBS_obj_ACC_it2) #simplify variables/objects' redundancy

tome_sample_meta_filtered$query_subclass <- paste(tome_sample_meta_filtered$split_class, tome_sample_meta_filtered$split_subclass, sep="_") #generate a final set of labels for subclass by combining class labels and subclass gene labels
Seu_AIBS_obj <- CreateSeuratObject(counts = AIBS_Rawcount_InEx, meta.data = tome_sample_meta_filtered) #after confirming our grouping variable, this creates the Seurat object
Seu_AIBS_obj <- NormalizeData(Seu_AIBS_obj, normalization.method = "LogNormalize", scale.factor = 1000000) #log normalize data

Idents(Seu_AIBS_obj) <- "query_subclass" #apply/assign our identity labels of interest to cells
Seu_AIBS_obj$it2_query_subclass <- Idents(Seu_AIBS_obj) #making a new identity variable; starting by copying the previous active idents (which was "query_subclass" at the time)

levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Oligo_MOBP"] <- "Oligo" #collapsing the oligo groups
levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Oligo_OPALIN"] <- "Oligo" #collapsing the oligo groups

levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Exc_RORB"] <- "Pyramidal" #collapsing the exc groups
levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Exc_LINC00507"] <- "Pyramidal" #collapsing the exc groups
levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Exc_FEZF2"] <- "Pyramidal" #collapsing the exc groups
levels(Seu_AIBS_obj$it2_query_subclass)[levels(Seu_AIBS_obj$it2_query_subclass)=="Exc_THEMIS"] <- "Pyramidal" #collapsing the exc groups

Idents(Seu_AIBS_obj) <- "it2_query_subclass" #setting the active identity to the new set of identity variables

Seu_AIBS_obj_it2 <- subset(Seu_AIBS_obj, idents = c("Inh_ADARB2", "Inh_LHX6"), invert = TRUE) #we want all cells except for these two groups
Seu_AIBS_obj <- Seu_AIBS_obj_it2 #simplify variables/objects' redundancy
remove(Seu_AIBS_obj_it2) #simplify variables/objects' redundancy

Seu_AIBS_markers_MTG <- FindAllMarkers(Seu_AIBS_obj, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "MAST") #use same thresholds as ACC for consistency/comparability; use MAST test
Seu_AIBS_markers_MTGroc <- FindAllMarkers(Seu_AIBS_obj, logfc.threshold = 2.5, min.pct = .35, only.pos = TRUE, return.thresh = .05, test.use = "roc") #use same thresholds as ASS for consistency/comparability; use roc test

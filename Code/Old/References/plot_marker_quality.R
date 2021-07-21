### this is from https://github.com/stripathy/human_marker_genes/blob/master/scripts/plot_marker_quality.R

# plot marker quality
library(scrattch.vis)
library(markerGeneProfile)
library(cowplot)
library(ggplot2)
library(homologene)
library(parallel)
library(Matrix)
library(Signac)

# read in pre-saved human snucseq gene expression file
human_exon_intron_cpm = readRDS(file = '~/allen_human/data/human_exon_intron_cpm.rds')
allen_human_data_dir = '~/allen_human/data-raw/allenHuman/'
allen_human_data_output_dir = '~/allen_human/data/'


humanGenes = read_csv(paste0(allen_human_data_dir, 'human_MTG_2018-06-14_genes-rows.csv'))
rownames(human_exon_intron_cpm) = humanGenes$gene

humanMetaJoined = readRDS(
  file= paste0(allen_human_data_output_dir, 'humanMetaJoined.rds'))


allen_human_subclass_names = humanMetaJoined$new_subclass %>% levels
allen_human_subclass_names = allen_human_subclass_names[allen_human_subclass_names != 'NA.']

subclass_cluster_markers = readRDS(file= paste0(allen_human_data_output_dir, 'derived_human_subclass_markers.rds'))
data("mouseMarkerGenes")

ogan_markers = readRDS(file = '/external/rprshnas01/netdata_kcni/stlab/marker_genes/human_markers_quick_sel.rds')


plotMarkers = function(gene_list, mouse_genes = F){
  new_gene_list = gene_list
  if (mouse_genes){
    new_gene_list =  mouse2human(new_gene_list) %>% pull(humanGene)
  }
  intersecting_genes = intersect(new_gene_list[1:min(length(new_gene_list), 50)] %>% make.names(), rownames(human_exon_intron_cpm))
  # use_genes = use_genes

  # use_genes = new_cell_markers
  data_df <- cbind(sample_name = colnames(human_exon_intron_cpm),
                   as.data.frame(as.matrix(t(human_exon_intron_cpm[intersecting_genes,]))))
  colnames(data_df) = colnames(data_df) %>% make.names()



  plot = group_dot_plot(data_df,
                        humanMetaJoined,
                        genes = intersecting_genes,
                        grouping = "subclass",
                        log_scale = FALSE,
                        font_size = 8,
                        rotate_counts = F)
  return(plot)

}

plot_base_dir = '~/human_marker_genes/plots/neuroexpresso_markers/'
# plot markers for all MGP cell types
lapply(names(mouseMarkerGenes$Cortex), function(marker_list){
  cell_type_name = marker_list
  p1 = plotMarkers(mouseMarkerGenes$Cortex[marker_list] %>% unlist %>% as.character() , mouse_genes = T)
  fn = paste0('~/human_marker_genes/plots/neuroexpresso_markers/', cell_type_name, '.png')
  save_plot(filename = fn, plot = p1, base_width = 8)
  dev.off()
})

plot_base_dir = '~/human_marker_genes/plots/human_derived_markers/'

# plot markers for all new AIBS-data based cell types
mclapply(allen_human_subclass_names, function(marker_list){
  cell_type_name = marker_list
  p1 = plotMarkers(subclass_cluster_markers[[cell_type_name]])
  fn = paste0('~/human_marker_genes/plots/human_derived_markers/', cell_type_name, '.png')
  save_plot(filename = fn, plot = p1, base_width = 8)
  dev.off()
}, mc.cores = 20)

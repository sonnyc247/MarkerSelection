# Tripathy lab gene-expression marker selection, validation, and analysis pipeline

This repo serves as an introduction to and step-by-step tutorial for marker gene identification from RNAseq data in the Tripathy lab at the Krembil Centre for Neuroinformatics in the Centre for Addiction and Mental Health (CAMH). 

Currently, it also holds documentation and code used in marker selection, validation, and related analyses as part of Consens et al. (WIP) and Hunter et al. (WIP).

This document and repo will consist of the following sections:

1) Loading data
2) Finding markers
3) Validating markers
4) Other related analyses

## Data loading 

This step is largely demonstrated by the documentation and code in https://github.com/sonnyc247/MarkerSelection/blob/master/Code/1_Loading_and_selecting_data.Rmd.

The ultimate goal of this step is to obtain count matrices and associated metadata, which may involve some work. 

### Count matrix

We desire a count matrix with the format of samples as columns and genes as rows - each cell of the matrix therefore represents the number of RNAseq counts associated with a given gene and expressed in a given sample. This matrix can be in sparse format, but columns still need to represent samples while rows represent genes. 

[Insert image of example count matrix]

Sometimes, this is as simple as downloading a count matrix from a dataset and simply importing it into R. This is largely the case with the example above, where we downloaded and loaded the count matrix ("Gene expression matrix", which links to "matrix.csv" as of March 31, 2021) from https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq. 

However, sometimes you do not get a ready-to-use count matrix. In some cases, we need to start from raw fastq files when there is no count matrix available. For this case, we would need to align the reads and generate a count matrix from scratch by following the steps and/or rationale outlined in https://github.com/sonnyc247/PSQ_Pipeline. Other times, 

# Tripathy lab gene-expression marker selection, validation, and analysis pipeline

This repo serves as an introduction to and step-by-step tutorial for marker gene identification from RNAseq data in the Tripathy lab at the Krembil Centre for Neuroinformatics in the Centre for Addiction and Mental Health (CAMH). 

Currently, it also holds documentation and code used in marker selection, validation, and related analyses as part of Consens et al. (WIP) and Hunter et al. (WIP).

This document and repo will consist of the following sections:

1) Loading data
2) Finding markers
3) Validating markers
4) Other related analyses

## <ins>1. Loading data</ins>

<ins>Code:</ins> 

This step is largely documented and demonstrated by the documentation and code in https://github.com/sonnyc247/MarkerSelection/blob/master/Code/1_Loading_and_selecting_data.Rmd.

The ultimate goal of this step is to obtain count matrices and associated metadata, which may involve some work. 

### Count matrix

We desire a count matrix with the format of samples (such as individual biological cells) as columns and genes as rows - each value in the matrix therefore represents the number of RNAseq counts detected/transcribed from a given gene and in a given sample. This matrix can be in sparse format, but columns still need to represent samples while rows still represent genes. 

[Insert image of example count matrix]

Sometimes, this is as simple as downloading a count matrix from a dataset and simply importing it into R. This is largely the case with the demo code above, where we downloaded and loaded the count matrix ("Gene expression matrix" on the website, which links to "matrix.csv" as of March 31, 2021) from https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq. 

However, sometimes you do not get a ready-to-use count matrix. In some cases, we need to start from raw fastq files when there is no count matrix available. For this case, we would need to align the RNAseq sequences/reads and generate a count matrix from scratch by following the steps and rationale outlined in https://github.com/sonnyc247/PSQ_Pipeline. Other times, we get gene expression data stored in other formats, such as 10x genomics count files. To see an example of how we handle this case, look at https://github.com/sonnyc247/MarkerSelection/blob/master/Code/scRNAseq_analysis.R, particularly the code under the section "Data loading - Zhou et al", which takes place after placing the files for each sample in their own folder and removing the sample name/prefixes from said files. Other situations may require yet other solutions.

### Metadata

Metadata files are usually already formatted in our desired format - tables with samples organized as rows and features (such as batch number, patient ID, sex, age, etc.) as columns. The usual challenge in this aspect of data loading is finding complete and relevant metadata. Sometimes, the metadata provided are not comprehensive enough for the analyses we want. This is usually resolved by looking further for more metadata of the same dataset from other sources, potentially by directly contacting study authors of publications with publicly available datasets. 

At the end of the data loading step, you have one count matrix and one metadata table in R, both containing all the samples you will ever potentially want in your study (i.e: the greatest set of data, which you can later subset, rather than coming back to this step and retrieving more samples) as well as all relevant features. The number of samples and the set of sample identifiers (such as cell IDs) need to be exactly the same between the count matrix and the metadata. Feel free to remove all samples you are sure you do not need in your analysis, to save memory and processing time downstream.

## <ins>2. Finding markers</ins>

<ins>Code:</ins> 

This step is largely documented and demonstrated in https://github.com/sonnyc247/MarkerSelection/blob/master/Code/2_Find_DE_genes_in_Seurat.Rmd.

<ins>References:</ins>

1. We use Seurat to find markers (https://satijalab.org/seurat/), which is a well-documented and developed toolkit and R package for working with RNAseq data. 
2. In particular, we currently use Seurat V3 (https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8), though marker selection and differential gene expression testing is not specific to V3. 
3. For the stats behind our differential gene expression testing, we often use MAST (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5).

We mostly find markers using Seurat's FindAllMarkers function. This is a function that performs differential gene expression testing between each cell grouping (or "Identity" in Seurat) and all other cells. For example, if a dataset has 10 cell types, FindAllMarkers will return genes that are statistically differentially expressed in cells of cell type 1 vs cells in types 2-10, in cell type 2 vs types 1 and 3-10 combined, and so on. A gene which is specifically and differentially expressed (often enriched) in a cell type is considered a marker.  There are several options and considerations when using FindAllMarkers, the full list of which can be seen in its documentation: https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/FindAllMarkers.

In our lab, we often adjust:

1. logfc.threshold - how big of a difference between cell groups do we want; we've often used 2 to 2.5
2. min.pct - for a gene to be tested, it must be detected in at least this fraction of cells of either group being compared (e.g: in cell type 1 or types 2-10 combined, when looking for cell type 1 markers); we have often used 0.35   
3. test.use - which statistical test should we use for differential gene expression testing; we've often used MAST for single-cell-esq data, though we have used "roc" as well  
4. only.pos - do we want only positive markers (do we want our cell type to be identified by the presence/enrichment of a gene's RNA transcripts, or do we also include genes that are consistently absent or reduced in our cell types of interest)
5. return.thresh - what false discovery rate or other confidence cut-off do we want for our results; we sometimes use a very non-selective cut-off to let us see all uncorrected p-values in screening-type analyses

We often adjust these parameters depending on the dataset and the purpose of use. For example, a bad-quality single-cell RNAseq dataset may have a lot of zeros - many genes may be mistakenly missed/not detected. In these datasets, we might loosen the "min.pct" cut-off to allow more genes to be tested, even if they're only detected in a small fraction of cells. We might set "only.pos" as TRUE depending on how we want to use the markers. For example, we would not suggest to wet-lab experimenters to look for a cell via immunohistochemistry for a negative marker. Cell type 1 might consistently not express Gene A, and that's indeed a characteristic of cell type 1, but it would be difficult to design a molecular probe to detect the <ins>absense</ins> of Gene A to then find a "cell type 1" cell.

### Settings for some common use cases

For finding markers to use in markerGeneProfile (https://github.com/PavlidisLab/markerGeneProfile) as in Consens et al. and Hunter et al., we have often used the following settings for finding markers from single-nucleus RNAseq data provided by the Allen Institute for Brain Science:

1. logfc.threshold = 2 to 2.5   
2. min.pct = 0.35 
3. test.use = "MAST"; = "roc" (we do two runs and accept markers that turn up in at least one set of results)  
4. only.pos = "TRUE"
5. return.thresh set to default / not specified

After getting these results, we filter out genes that serve as markers/are differentially enriched in multiple cell types.

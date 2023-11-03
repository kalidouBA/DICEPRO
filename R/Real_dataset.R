#'@title GSE107011 dataset
#'
#'@description RNA-Seq profiling of 29 immune cell types and peripheral blood mononuclear cells
#'
#'@docType data
#
#'@keywords data
#'
#'@name GSE107011
#'
#'@usage data(GSE107011)
#'
#'@details Total RNA of 29 immune cell types (from 4 individuals) and peripheral blood mononuclear cells (PBMCs, from 13 individuals)
#'was extracted for gene expression profiling. The 13 PBMCs samples were processed with both microarray and RNA-Seq platforms.
#'
#'@details We performed RNA-Seq transcriptome profiling on 29 immune cell types consituting peripheral blood mononuclear cells (PBMCs)
#'sorted from 4 Singaporean-Chinese individuals (S4 cohort). We also performed RNA-Seq and microarray transcriptome profiling of PBMCs
#'from an extended cohort of 13 individuals (S13 cohort). The data was used first to characterize the transcriptomic signatures and
#'relationships among the 29 immune cell types. Then we explored the difference in mRNA composition in terms of transcripts proportions
#'and abundance. Lastly, we performed deep deconvolution for both microarray and RNA-Seq technologies.
#'Flow cytometry measurements of 6 cell populations:
#'\enumerate{
#' \item{Monocytes}
#' \item{CD4Tcells}
#' \item{Bcells}
#' \item{CD8Tcells}
#' \item{NK}
#' \item{Unknown}
#'}
"GSE107011"

#'@title GSE107572 dataset
#'
#'@description Blood-derived immune-cell mixtures from 9 healthy donors
#'
#'@docType data
#
#'@keywords data
#'
#'@name GSE107572
#'
#'@usage data(GSE107572)
#'
#'@details Total RNA of 29 immune cell types (from 4 individuals) and peripheral blood mononuclear cells (PBMCs, from 13 individuals)
#'was extracted for gene expression profiling. The 13 PBMCs samples were processed with both microarray and RNA-Seq platforms.
#'
#'@details We introduce quanTIseq, a method to quantify the tumor immune contexture, determined by the type and density of tumor-infiltrating
#'immune cells. quanTIseq is based on a novel deconvolution algorithm for RNA sequencing data that was validated with independent data sets.
#'Complementing the deconvolution output with image data from tissue slides enables in silico multiplexed immunodetection and provides an
#'efficient method for the immunophenotyping of a large number of tumor samples.
#'Flow cytometry measurements of 8 cell populations:
#'\enumerate{
#' \item{B.cells}
#' \item{T.cells.CD4}
#' \item{T.cells.CD8}
#' \item{Tregs}
#' \item{Neutrophils}
#' \item{NK.cells}
#' \item{Dendritic.cells}
#' \item{Monocytes}
#'}
"GSE107572"

#'@title GSE127813 dataset
#'
#'@description High-throughput tissue dissection and cell purification with digital cytometry
#'
#'@docType data
#
#'@keywords data
#'
#'@name GSE127813
#'
#'@usage data(GSE127813)
#'
#'@details Whole blood samples were collected from 12 healthy adult subjects and immediately processed to enumerate leukocyte composition by FACS
#'using an FDA-approved in vitro diagnostic test (IVD Multitest 6-color TBNK, Becton Dickinson) and automated hematology analyzer for blood leukocyte
#'differential counts (Sysmex XE-2100) in a CLIA hematology lab setting (Stanford Clinical Laboratories). Aliquots from the same whole blood samples
#'were stored in PAXgene blood RNA tubes (Qiagen) for subsequent RNA extraction and RNA-Seq library preparation.
#'
#'@details Tissue composition is a major determinant of phenotypic variation and a key factor influencing disease outcomes. Although scRNA-Seq
#'has emerged as a powerful technique for characterizing cellular heterogeneity, it is currently impractical for large sample cohorts and cannot
#'be applied to fixed specimens collected as part of routine clinical care. To overcome these challenges, we extended Cell type Identification
#'By Estimating Relative Subsets Of RNA Transcripts (CIBERSORT) into a new platform for in silico cytometry. Our approach enables the simultaneous
#'inference of cell type abundance and cell type-specific gene expression profiles (GEPs) from bulk tissue transcriptomes. The utility of this
#'integrated framework, called CIBERSORTx, is demonstrated in multiple tumor types, including melanoma, where single cell reference profiles are
#'used to dissect primary clinical specimens, revealing cell type-specific signatures of driver mutations and immunotherapy response. We anticipate
#'that digital cytometry will augment single cell profiling efforts, enabling cost-effective, high throughput tissue characterization without the
#'need for antibodies, disaggregation, or viable cells.
#'Flow cytometry measurements of 8 cell populations:
#'\enumerate{
#' \item{B cells}
#' \item{T cells CD4}
#' \item{T cells CD8}
#' \item{T cells}
#' \item{Neutrophils}
#' \item{NK cells}
#' \item{Lymphocytes}
#' \item{Monocytes}
#'}
"GSE127813"

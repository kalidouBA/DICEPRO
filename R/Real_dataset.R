#'@title REAL_DATASET dataset
#'
#'@description
#'\code{GSE107011} : RNA-Seq profiling of 29 immune cell types and peripheral blood mononuclear cells.
#'\code{GSE107572} : Blood-derived immune-cell mixtures from 9 healthy donors
#'\code{GSE127813} : High-throughput tissue dissection and cell purification with digital cytometry
#'
#'@docType data
#
#'@keywords data
#'
#'@name REAL_DATASET
#'@aliases GSE107011 GSE107572 GSE127813
#'
#'@usage data(REAL_DATASET)
#'
#'
#'@details We tested the methods on a set consisting of four different real-world bulk RNA-seq samples obtained using RNA-Seq technology available on GEO with accessions \code{GSE107011}, \code{GSE107572} and \code{GSE127813}.
#'\code{GSE107011} is a bulk RNA-Seq dataset from 127 total samples of which 12 PBMC samples were selected for this analysis. It contains expression profiles of 29 blood cell subtypes, of which are mainly analyzed
#'B, CD4+ T, CD8+ T, NK and Monocytes.
#'\code{GSE107572} contains preprocessed RNA-seq data in Transcripts Per Million (TPM) units from nine colorectal cancer patients. These data also has flow cytometry estimates for the
#'immune subpopulation corresponding to cells B, CD4+ T, CD8+ T, NK, Monocyte, Dendritic, Neutrophils and Tregs.
#'\code{GSE127813} dataset were collected from 12 healthy adult subjects to enumerate leukocyte composition by FACS using an in vitro diagnostic test and automated hematology analyzer for
#' blood leukocyte differential counts. The fractions of ground truth of  Neutrophils, Lymphocytes and Monocytes  are compared with the corresponding estimated fractions obtained using methods.
#'
#'Flow cytometry measurements of 5, 8, and 8 cell populations for \code{GSE107011}, \code{GSE107572}, \code{GSE127813} datasets respectively:
#'\enumerate{
#' \item{B cells}
#' \item{T cells CD4}
#' \item{T cells CD8}
#' \item{Tregs}
#' \item{Neutrophils}
#' \item{NK cells}
#' \item{Dendritic cells}
#' \item{Monocytes}
#' \item{Lymphocytes}
#'}
#

#'@format The data are composed of 3 objects:
#'\describe{
#'  \item{\code{GSE107011}: }{a \code{list} Total RNA of 29 immune cell types (from 4 individuals) and peripheral blood mononuclear cells (PBMCs, from 13 individuals) was extracted for gene expression profiling.}
#'  \item{\code{GSE107572}: }{a \code{list} Blood-derived immune-cell mixtures from 9 healthy donors}
#'  \item{\code{GSE127813}: }{a \code{list} Whole blood samples were collected from 12 healthy adult subjects and immediately processed to enumerate leukocyte composition by FACS using an FDA-approved in vitro diagnostic test.}
#'}
#'
#'
#'@references {
#' \code{GSE107011}: Xu W, Monaco G, Wong EH, Tan WLW et al. Mapping of γ/δ T cells reveals Vδ2+ T cells resistance to senescence. EBioMedicine 2019 Jan;39:44-58. PMID: 30528453
#' \code{GSE107572}: Finotello, F., Mayer, C., Plattner, C. et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Med 11, 34 (2019). https://doi.org/10.1186/s13073-019-0638-6
#' \code{GSE127813}: Newman AM, Steen CB, Liu CL, Gentles AJ et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 2019 Jul;37(7):773-782. PMID: 31061481
#'}

"GSE107011"
"GSE107572"
"GSE127813"
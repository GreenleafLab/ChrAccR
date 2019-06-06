#' ChrAccR: Analyzing chromatin accessibility data in R.
#'
#' Tools for analyzing chromatin accessibility data in R. Currently supports ATAC-seq and NOMe-seq data analysis.
#'
#' @import Matrix
#' @import methods
#' @importFrom stats median quantile cor cov var sd ecdf fisher.test dist as.dist hclust as.hclust convolve na.omit as.formula
#' @importFrom utils packageVersion read.table write.table combn data
#' @importFrom data.table data.table as.data.table fread
#' @import GenomicRanges
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom S4Vectors queryHits subjectHits queryLength DataFrame elementNROWS
#' @import GenomicAlignments
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlengths seqlengths<- seqnames genome genome<- organism sortSeqlevels
#' @import SummarizedExperiment
#' @import ggplot2
#' @import muLogR
#' @import muRtools
#' @docType package
#' @name ChrAccR
NULL

# # avoid NOTEs in R CMD CHECK
# utils::globalVariables(c(
# 	":=", ".", "..count..", ".N", ".SD",
# 	"V1", "V2"
# ))

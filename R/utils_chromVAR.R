################################################################################
# Utilities and helper functions for chromVAR
################################################################################
#' getGenomeObject
#'
#' retrieve the appropriate \code{BSgenome} for an assembly string
#'
#' @param assembly     string specifying the assembly
#' @param adjChrNames  should the prefix "chr" be added to main chromosomes if not already present and chrMT be renamed to chrM?
#' @return \code{BSgenome} object
#' @export
getGenomeObject <- function(assembly, adjChrNames=TRUE){
	mainREnum <- "^([1-9][0-9]?|[XYM]|MT)$"
	if (is.element(assembly, c("hg19"))){
		require(BSgenome.Hsapiens.UCSC.hg19)
		res <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
	} else if (is.element(assembly, c("GRCh37", "GRCh37_chr"))){
		require(BSgenome.Hsapiens.1000genomes.hs37d5)
		res <- BSgenome.Hsapiens.1000genomes.hs37d5
	} else if (is.element(assembly, c("hg38", "hg38_chr"))){
		require(BSgenome.Hsapiens.UCSC.hg38)
		res <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
	} else if (is.element(assembly, c("GRCh38", "GRCh38_chr"))){
		require(BSgenome.Hsapiens.NCBI.GRCh38)
		res <- BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens
	} else if (is.element(assembly, c("mm9"))){
		require(BSgenome.Mmusculus.UCSC.mm9)
		res <- BSgenome.Mmusculus.UCSC.mm9::Mmusculus
	} else if (is.element(assembly, c("mm10"))){
		require(BSgenome.Mmusculus.UCSC.mm10)
		res <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
	} else {
		stop(paste0("Unknown assembly:", assembly))
	}
	if (adjChrNames){
		prep <- grepl(mainREnum, seqnames(res))
		seqnames(res)[prep] <- paste0("chr", seqnames(res)[prep])
		seqnames(res)[seqnames(res)=="chrMT"] <- "chrM"
	}
	return(res)
}

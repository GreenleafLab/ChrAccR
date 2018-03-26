################################################################################
# Utilities and helper functions for chromVAR
################################################################################
#' getGenomeObject
#'
#' retrieve the appropriate \code{BSgenome} for an assembly string
#'
#' @param assembly     string specifying the assembly
#' @return \code{BSgenome} object
#' @export
getGenomeObject <- function(assembly){
	if (is.element(assembly, c("hg19"))){
		require(BSgenome.Hsapiens.UCSC.hg19)
		res <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
	} else if (is.element(assembly, c("GRCh37", "GRCh37_chr"))){
		require(BSgenome.Hsapiens.1000genomes.hs37d5)
		res <- BSgenome.Hsapiens.1000genomes.hs37d5
	} else if (is.element(assembly, c("hg38", "GRCh38", "hg38_chr", "GRCh38_chr"))){
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
	return(res)
}

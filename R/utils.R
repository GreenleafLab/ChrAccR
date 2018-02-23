#' cleanMem
#' 
#' clean the system mory by invoking the garbage collector
#' @param iter.gc number of times to invoke the garbage collector
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @export
cleanMem <- function(iter.gc=1L){
	if (getConfigElement("cleanMem")){
		for (i in 1:iter.gc){gc()}
	}
	invisible(NULL)
}

#' isCanonicalChrom
#' 
#' for a character string of chromosome names, determine if it is a canonical chromosome
#' (i.e. not not ChrUn*, *_random, ...)
#' @param ss character vector of chromosome names
#' @return logical vector stating whether the given chromosome names correspond to canonical chromosomes
#' @author Fabian Mueller
#' @export
isCanonicalChrom <- function(ss){
	re <- "^(chr)?([1-9][0-9]?|[XYM]|MT)$"
	return(grepl(re, ss))
}

#' prepareMotifmatchr
#' 
#' prepare objects for a \code{motifmatchr} analysis
#' @param genome character string specifying genome assembly
#' @param motifs either a character string (currently only "jaspar" is supported) or an object containing PWMs
#'               that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} object)
#' @return a list containing objects to be used as arguments for \code{motifmatchr}
#' @author Fabian Mueller
#' @noRd
prepareMotifmatchr <- function(genome, motifs){
	require(motifmatchr)
	res <- list()

	# get the species name and the genome sequence object based on the object
	spec <- NULL
	genomeObj <- NULL
	if (is.element(genome, c("mm9", "mm10"))){
		spec <- "Mus musculus"
		genomePkg <- paste0("BSgenome.Mmusculus.UCSC.", genome)
		require(genomePkg, character.only=TRUE)
		genomeObj <- get(genomePkg)
	} else if (is.element(genome, c("hg18", "hg19", "hg38"))){
		spec <- "Homo sapiens"
		genomePkg <- paste0("BSgenome.Hsapiens.UCSC.", genome)
		require(genomePkg, character.only=TRUE)
		genomeObj <- get(genomePkg)
	} else {
		logger.error(c("Unsupported genome:", genome))
	}

	# get the motif PWMs
	if (is.character(motifs)){
		if (motifs=="jaspar"){
			motifs <- getJasparMotifs(species=spec)
		} else {
			logger.error(c("Unsupported motifs:", motifs))
		}
	}
	res[["genome"]] <- genomeObj
	res[["motifs"]] <- motifs
	return(res)
}

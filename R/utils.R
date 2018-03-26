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
#' @param motifs either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'               that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object)
#' @return a list containing objects to be used as arguments for \code{motifmatchr}
#' @author Fabian Mueller
#' @noRd
prepareMotifmatchr <- function(genome, motifs){
	require(motifmatchr)
	res <- list()

	# get the species name and the genome sequence object based on the object
	spec <- NULL
	genomeObj <- genome
	if (!is.element("BSgenome", class(genomeObj))){
		genomeObj <- getGenomeObject(genome)
	}

	# get the motif PWMs
	motifL <- TFBSTools::PWMatrixList()
	if (is.character(motifs)){
		if (is.element("jaspar", motifs)){
			require(chromVAR)
			motifL <- c(motifL, getJasparMotifs(species=spec))
		}
		if (is.element("homer", motifs)){
			require(chromVARmotifs)
			data("homer_pwms")
			motifL <- c(motifL, homer_pwms)
		}
		if (is.element("encode", motifs)){
			require(chromVARmotifs)
			data("encode_pwms")
			motifL <- c(motifL, encode_pwms)
		}
		if (is.element("cisbp", motifs)){
			require(chromVARmotifs)
			if (spec == "Mus musculus"){
				data("mouse_pwms_v2")
				motifL <- c(motifL, mouse_pwms_v2)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v2")
				motifL <- c(motifL, human_pwms_v2)
			} else {
				logger.warning(c("Could not find cisBP annotation for species", spec))
			}
		} 
		if (length(motifL) < 1) {
			logger.error(c("No motifs were loaded. Unsupported motifs (?) :", motifs))
		}	
	} else if (is.element("PWMatrixList", class(motifs)) || is.element("PFMatrixList", class(motifs))) {
		motifL <- motifs
	} else {
		logger.error(c("unsupported value for motifs:", motifs))
	}	
	res[["genome"]] <- genomeObj
	res[["motifs"]] <- motifL
	return(res)
}

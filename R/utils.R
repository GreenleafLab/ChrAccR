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

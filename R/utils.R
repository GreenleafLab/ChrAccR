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

#' fastDelayedArrayToMatrix
#' 
#' faster subsetting by index of DelayedArrays via linear indexing.
#' Code taken from: https://github.com/Bioconductor/DelayedArray/issues/13
#' @param X \code{DelayedArray}
#' @param i row indices
#' @param j column indices
#' @return a regular matrix object representing the indexed submatrix
#' @author Fabian Mueller
#' @noRd
fastDelayedArrayToMatrix <- function(X, i=NULL, j=NULL){
	M <- X
	if (!is.null(i) || !is.null(j)){
		linIdx <- DelayedArray:::to_linear_index(list(i, j), dim(X))
		M <- matrix(X[linIdx], ncol=ncol(X))
	} else {
		M <- as.matrix(X)
	}
	return(M)
}

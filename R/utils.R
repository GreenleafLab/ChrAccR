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
#' @param i row indices (or logical)
#' @param j column indices (or logical)
#' @return a regular matrix object representing the indexed submatrix
#' @author Fabian Mueller
#' @noRd
fastDelayedArrayToMatrix <- function(X, i=NULL, j=NULL){
	M <- X
	if (!is.null(i)){
		if (is.logical(i)){
			if (length(i)!=nrow(X)) stop("Invalid index i (logical)")
			i <- which(i)
		}
		if (is.character(i)) stop("Invalid index i (character indexing not supported in fastDelayedArrayToMatrix)")
		if (any(i > nrow(X) | i < 1)) stop("Invalid index i")
	}
	if (!is.null(j)) {
		if (is.logical(j)){
			if (length(j)!=ncol(X)) stop("Invalid index j (logical)")
			j <- which(j)
		}
		if (is.character(j)) stop("Invalid index j (character indexing not supported in fastDelayedArrayToMatrix)")
		if (any(j > ncol(X) | j < 1)) stop("Invalid index j")
	}

	if (!is.null(i) || !is.null(j)){
		linIdx <- DelayedArray:::to_linear_index(list(i, j), dim(X))
		nc <- ncol(X)
		if (!is.null(j)) nc <- length(j)
		# linear indexing with more than .Machine$integer.max indices causes trouble
		# so does any element in the linear index that exceeds .Machine$integer.max
		# --> avoid the ploblem using MEMORY INEFFICIENT coercion to regular matrix
		if (any(linIdx >= .Machine$integer.max) || length(linIdx >= .Machine$integer.max)){
			logger.warning("Linear index exceeds INT_MAX --> coercing to regular matrix [fastDelayedArrayToMatrix]")
			M <- as.matrix(M)
			if (!is.null(i)) M <- M[i,,drop=FALSE]
			if (!is.null(j)) M <- M[,j,drop=FALSE]
		} else {
			M <- matrix(X[linIdx], ncol=nc)
		}
		rnames <- rownames(X)
		if (!is.null(i) && !is.null(rnames)) rnames <- rnames[i]
		rownames(M) <- rnames
		cnames <- colnames(X)
		if (!is.null(j) && !is.null(cnames)) cnames <- cnames[j]
		colnames(M) <- cnames
	} else {
		M <- as.matrix(X)
	}
	return(M)
}

#' getGroupsFromTable
#'
#' Retrieve groupings given a table containing some categorical columns
#'
#' @param tt              table to retrieve groupings for
#' @param cols         (Optional) predefined column names (in the form of a \code{character} vector) or indices (an
#'                        \code{integer} vector) to consider. All other columns in the annotation table will be ignored.
#' @param minGrpSize      Minimum number of items required to form a group in comparison
#' @param maxGrpCount     Maximum number of groups to be considered
#' @return List of groupings. Each element corresponds to a categorical column in the table and contains the row indices for each
#'         category
#'
#' @author Fabian Mueller
#' @export
getGroupsFromTable <- function(tt, cols=NULL, minGrpSize=2, maxGrpCount=nrow(tt)-1) {

	if (nrow(tt) < 2) logger.error("Required a table with at least 2 rows")

	if (!is.null(cols)) tt <- tt[,cols,drop=FALSE]
	if (ncol(tt) < 1) logger.error("Required a table with at least 1 column")

	idxVec <- 1:nrow(tt)
	
	res <- list()
	for (j in 1:ncol(tt)){
		cname <- colnames(tt)[j]

		vv <- tt[,j]
		rr <- tapply(idxVec, vv, identity)

		rr <- rr[sapply(rr, length) > 0] # ignore levels that are missing in the dataset
		rr <- rr[sapply(rr, function(x){!any(is.na(x))})] # ignore levels that are missing in the dataset
		passesMinSize <- sapply(rr, length) >= minGrpSize
		if (length(rr) > 1 && length(rr) <= maxGrpCount && all(passesMinSize)){
			res[[cname]] <- rr
		}
	}
	return(res)
}

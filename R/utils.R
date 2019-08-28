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

#' findOrderedNames
#' 
#' find the first occurrence of a name in a vector of strings
#' @param ss character vector in which the name should be found
#' @param orderedNames vector of names that will be queried. This method will go through them one by one and find the first occurrence
#'           in the order of the orderedNames provided
#' @param exact should only be exact matches be reported
#' @param ignore.case should casing be ignored
#' @return the string that matches the first occurrence in the order of \code{orderedNames}. Returns \code{NA} if no match is found.
#' @author Fabian Mueller
#' @export
findOrderedNames <- function(x, orderedNames, exact=TRUE, ignore.case=FALSE){
	rr <- NA
	for (cn in orderedNames){
		pat <- gsub(".", "\\.", cn, fixed=TRUE)
		if (exact) pat <- paste0("^", pat, "$")
		foundCn <- grep(pat, x, value=TRUE, ignore.case=ignore.case)[1]
		if (!is.na(foundCn)) {
			rr <- foundCn
			break
		}
	}
	return(rr)
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
		if (any(linIdx >= .Machine$integer.max) || length(linIdx) >= (.Machine$integer.max/2 - 1) || prod(dim(X)) >= .Machine$integer.max){
			logger.warning("Linear index or matrix dimensions exceed INT_MAX --> coercing to regular matrix [fastDelayedArrayToMatrix]")
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

#' safeMatrixStats
#' 
#' Compute matrix statistics selecting the appropriate function depending on the
#' matrix class of the input (supports sparse matrices and DelayedArrays)
#'
#' @param X	input matrix
#' @param statFun statistic. E.g. \code{"rowSums"}, \code{"colSums"}, \code{"rowMeans"}, \code{"colMeans"}, ...
#' @param ... arguments passed on to the matrix stats function. E.g. \code{na.rm}.
#' @return result of the corresponding matrix statistic
#' @author Fabian Mueller
#' @export
#}
safeMatrixStats <- function(X, statFun="rowSums", ...){
  cl <- class(X)
  pkg <- attr(cl, "package")
  if (is.character(pkg)) {
  	if (pkg=="Matrix") {
  		statFun <- paste0("Matrix::", statFun)
  	} else if (is.element(pkg, c("DelayedArray", "HDF5Array"))){
  		statFun <- paste0("BiocGenerics::", statFun)
  	}
  }
  # print(statFun)
  statFun <- eval(parse(text=statFun))
  return(statFun(X, ...))
}

#' custom_cicero_cds
#' 
#' More performant version of \code{cicero::make_cicero_cds}.
#' Allows for more parameter customization and doesn't run into memory trouble as easily.
#' @author hpliner, Jeffrey Granja
#' @noRd
custom_cicero_cds <- function(
		cds,
		reduced_coordinates,
		k=50,
		max_knn_iterations = 5000,
		summary_stats = NULL,
		size_factor_normalize = TRUE,
		silent = FALSE
) {
	require(dplyr)
	assertthat::assert_that(is(cds, "cell_data_set"))
	assertthat::assert_that(is.data.frame(reduced_coordinates) | is.matrix(reduced_coordinates))
	assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates), nrow(pData(cds))))
	assertthat::assert_that(setequal(row.names(reduced_coordinates), colnames(cds)))
	assertthat::assert_that(assertthat::is.count(k) & k > 1)
	assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
	if(!is.null(summary_stats)) {
		assertthat::assert_that(all(summary_stats %in% names(pData(cds))),
		                        msg = paste("One of your summary_stats is missing",
		                                    "from your pData table. Either add a",
		                                    "column with the name in",
		                                    "summary_stats, or remove the name",
		                                    "from the summary_stats parameter.",
		                                    collapse = " "))
		assertthat::assert_that(sum(vapply(summary_stats, function(x) {
		  !(is(pData(cds)[,x], "numeric") | is(pData(cds)[,x], "integer"))}, 1)) == 0,
		                        msg = paste("All columns in summary_stats must be",
		                                    "of class numeric or integer.",
		                                    collapse = " "))
	}
	assertthat::assert_that(is.logical(size_factor_normalize))
	assertthat::assert_that(is.logical(silent))
	reduced_coordinates <- as.data.frame(reduced_coordinates)
	reduced_coordinates <- reduced_coordinates[colnames(cds),]

	start <- Sys.time()
	# Create a k-nearest neighbors map
	message("\nFNN k-nearest search...")
	nn_map <- FNN::knn.index(reduced_coordinates, k=(k-1)) # no data.frame wrapper
	row.names(nn_map) <- row.names(reduced_coordinates)
	nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))

	good_choices <- seq_len(nrow(nn_map))
	choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
	chosen <- good_choices[choice]
	good_choices <- good_choices[good_choices != good_choices[choice]]
	it <- 0
	k2 <- k * 2 # Compute once

	# function for sapply
	get_shared <- function(other, this_choice) {
		k2 - length(union(cell_sample[other,], this_choice))
	}

	while (length(good_choices) > 0 & it < max_knn_iterations) { # slow
		if(!silent && (it %% 500 == 0)) message(sprintf("%s of %s iterations", it, max_knn_iterations))
		it <- it + 1
		choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
		new_chosen <- c(chosen, good_choices[choice])
		good_choices <- good_choices[good_choices != good_choices[choice]]
		cell_sample <- nn_map[new_chosen,]
		others <- seq_len(nrow(cell_sample) - 1)
		this_choice <- cell_sample[nrow(cell_sample),]
		shared <- sapply(others, get_shared, this_choice = this_choice)

		if (max(shared) < .9 * k) {
		  chosen <- new_chosen
		}
	}
	if (!silent) message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
	cell_sample <- nn_map[chosen,]
	cell_sample_map <- lapply(seq_len(nrow(cell_sample)), function(x) rownames(reduced_coordinates)[cell_sample[x,]]) %>% Reduce("rbind",.) %>% data.frame
	rownames(cell_sample_map) <- rownames(cell_sample)
	if(!silent) {
		# Only need this slow step if !silent
		combs <- combn(nrow(cell_sample), 2)
		combs <- combs[,sample(seq_len(ncol(combs)),min(ncol(combs),10^6))] #sample 1 M because this really doesnt matter
		shared <- apply(combs, 2, function(x) {  #slow
		  k2 - length(unique(as.vector(cell_sample[x,])))
		})

		message(paste0("\nOverlap QC metrics:\nCells per bin: ", k,
		               "\nMaximum shared cells bin-bin: ", max(shared),
		               "\nMean shared cells bin-bin: ", mean(shared),
		               "\nMedian shared cells bin-bin: ", median(shared)))

		if (mean(shared)/k > .1) warning("On average, more than 10% of cells are shared between paired bins.")
	}
	if (!silent) message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
	message("\nMaking aggregated scATAC Matrix...")
	exprs_old <- exprs(cds)

	new_exprs <- matrix(0, nrow = nrow(cell_sample), ncol = nrow(exprs_old))
	for(x in seq_len(nrow(cell_sample))){
		if(!silent && (x %% 500 == 0)){
		    message(sprintf("%s of %s iterations : %s minutes since start", x, nrow(cell_sample), round(difftime(Sys.time(),start,units="mins"),1)))
		}
		new_exprs[x,] <- Matrix::rowSums(exprs_old[,cell_sample[x,]])
	}
	remove(exprs_old)
	if (!silent) message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
	message("\nMaking aggregated CDS...")
	pdata <- pData(cds)
	new_pcols <- "agg_cell"
	if(!is.null(summary_stats)) { 
		new_pcols <- c(new_pcols, paste0("mean_",summary_stats)) 
	}

	new_pdata <- plyr::adply(cell_sample,1, function(x) {
		sub <- pdata[x,]
		df_l <- list()
		df_l["temp"] <- 1
		for (att in summary_stats) {
		  df_l[paste0("mean_", att)] <- mean(sub[,att])
		}
		data.frame(df_l)
	})

	new_pdata$agg_cell <- paste("agg", chosen, sep="")
	new_pdata <- new_pdata[,new_pcols, drop = FALSE] # fixes order, drops X1 and temp

	row.names(new_pdata) <- new_pdata$agg_cell
	row.names(new_exprs) <- new_pdata$agg_cell
	new_exprs <- as.matrix(t(new_exprs))

	fdf <- monocle3::fData(cds)
	new_pdata$temp <- NULL

	cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(new_exprs,
		cell_metadata = new_pdata,
		gene_metadata = fdf
	))

	cicero_cds <- monocle3::detect_genes(cicero_cds, min_expr = .1)
	cicero_cds <- monocle3::estimate_size_factors(cicero_cds)

	#cicero_cds <- suppressWarnings(BiocGenerics::estimateDispersions(cicero_cds))
	if (any(!c("chr", "bp1", "bp2") %in% names(monocle3::fData(cicero_cds)))) {
		monocle3::fData(cicero_cds)$chr <- NULL
		monocle3::fData(cicero_cds)$bp1 <- NULL
		monocle3::fData(cicero_cds)$bp2 <- NULL
		monocle3::fData(cicero_cds) <- cbind(
			monocle3::fData(cicero_cds),
			cicero::df_for_coords(row.names(monocle3::fData(cicero_cds)))
		)
	}
	if (!silent) message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))
	if (size_factor_normalize) {
		message("\nSize factor normalization...")
		cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(Matrix::t(Matrix::t(monocle3::exprs(cicero_cds))/monocle3::pData(cicero_cds)$Size_Factor), 
			cell_metadata = monocle3::pData(cicero_cds), 
			gene_metadata = monocle3::fData(cicero_cds)
		))
	}
	if (!silent) message(sprintf("%s minutes since start", round(difftime(Sys.time(),start,units="mins"),1)))

	# return(list(ciceroCDS = cicero_cds, knnMap = cell_sample_map))
	return(cicero_cds)
}

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
#' @param x character vector in which the name should be found
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

#' findNearestGeneForGr
#'
#' get gene annotation for a \code{GRanges} object by linking to the nearest gene
#'
#' @param gr    \code{GRanges} object
#' @param geneGr  gene annotation from which to pull the annotation. Can be \code{NULL} for
#'              automatic retrieval of annotation. Must be named or have a gene name column
#'              in the metadata
#' @param maxDist maximum distance for matching to nearest gene
#' @return \code{data.frame} containing information on the nearest gene for each element in \code{gr}
#'
#' @export
findNearestGeneForGr <- function(gr, geneGr=NULL, maxDist=1e5){
	if (is.null(geneGr)){
		genAss <- genome(gr)[1]
		annoPkg <- getChrAccRAnnotationPackage(genAss)
		if (is.null(annoPkg)) logger.error("Annotation package needed")
		geneGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
	}
	if (is.null(names(geneGr))){
		emd <- elementMetadata(geneGr)
		nameCol <- findOrderedNames(colnames(emd), c("gene_name", "genename", "gene_id", "geneid"), exact=TRUE, ignore.case=TRUE)
		if (!is.na(nameCol)){
			names(geneGr) <- emd[,nameCol]
		} else {
			logger.error("geneGr must have names")
		}
	}

	tssGr <- promoters(geneGr, upstream=0, downstream=1) #get the TSS coordinate
	dd <- distanceToNearest(gr, tssGr, ignore.strand=TRUE, select="arbitrary")
	dd <- dd[mcols(dd)[,"distance"] <= maxDist,] # remove too far matches

	res <- data.frame(
		gene_name=rep(as.character(NA), length(gr)),
		dist_to_tss=rep(as.integer(NA), length(gr)),
		stringsAsFactors=FALSE
	)
	geneGr.sub <- geneGr[subjectHits(dd)]
	geneStrand <- as.character(strand(geneGr.sub))
	emd <- elementMetadata(geneGr.sub)
	# assign negative distances for elements that are downstream of the gene
	dist.signed <- mcols(dd)[,"distance"]
	coord.q <- start(resize(gr[queryHits(dd)], width=1, fix="center", ignore.strand=TRUE))
	# isNeg.q <- strand(gr[queryHits(dd)])=="-"
	coord.s <- start(tssGr[subjectHits(dd)])
	isNeq.s <- geneStrand=="-"
	isDownstream <- (dist.signed > 0) & ((!isNeq.s & (coord.q > coord.s)) | (isNeq.s & (coord.q < coord.s)))
	dist.signed[isDownstream] <- -dist.signed[isDownstream]

	res[queryHits(dd),"gene_name"]   <- names(geneGr.sub)
	res[queryHits(dd),"dist_to_tss"] <- dist.signed
	res[queryHits(dd),"gene_chrom"]  <- as.character(seqnames(geneGr.sub))
	res[queryHits(dd),"gene_chromStart"]  <- start(geneGr.sub)
	res[queryHits(dd),"gene_chromEnd"]  <- end(geneGr.sub)
	res[queryHits(dd),"gene_strand"]  <- geneStrand

	return(res)
}

#' rowZscores
#' 
#' Performs z-score normalization on the rows of a matrix. (Basically a wrapper around \code{matrixStats})
#' @param X input matrix
#' @param na.rm should NAs be omitted?
#' @return z-score normalized matrix
#' @author Fabian Mueller
#' @export
rowZscores <- function(X, na.rm=FALSE){
	(X - rowMeans(X, na.rm=na.rm)) / matrixStats::rowSds(X, na.rm=na.rm)
}
#' colZscores
#' 
#' Performs z-score normalization on the columns of a matrix. (Basically a wrapper around \code{matrixStats})
#' @param X input matrix
#' @param na.rm should NAs be omitted?
#' @return z-score normalized matrix
#' @author Fabian Mueller
#' @export
colZscores <- function(X, na.rm=FALSE){
	t((t(X) - colMeans(X, na.rm=na.rm)) / matrixStats::colSds(X, na.rm=na.rm))
}

#' fastDelayedArrayToMatrix
#' 
#' [not needed anymore according to the DelayedArray developers]
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
		allNA <- all(is.na(vv))
		if (allNA){
			logger.warning(c("Detected only NAs for sample annotation column:", cname))
		} else {
			rr <- tapply(idxVec, vv, identity)

			rr <- rr[sapply(rr, length) > 0] # ignore levels that are missing in the dataset
			rr <- rr[sapply(rr, function(x){!any(is.na(x))})] # ignore levels that are missing in the dataset
			passesMinSize <- sapply(rr, length) >= minGrpSize 
			if (length(rr) > 1 && length(rr) <= maxGrpCount && !all(passesMinSize)){
				notpassed <- names(which(passesMinSize == FALSE)) #remove comparasions that doesn't pass min size
                rr<- within(rr, rm(notpassed)) 
                rr <- rr[- which(names(rr) %in% notpassed)]
			}
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

#' getRBF
#' 
#' retrieve a gaussian radial basis function for weighting
#'
#' @param sigma   decay parameter (shape parameter (epsilon) of a gaussian RBF= 1/(sqrt(2)*sigma))
#' @param cutoff  distance cutoff. Everything beyond will be 0
#' @param ymin    y assymptote, i.e. the assymptotic minimum of the function
#' @return parametrized function
#' @author Fabian Mueller
#' @noRd
getRBF <- function(sigma=10000, cutoff=Inf, ymin=0){
	hassym <- 0
	if (ymin > 0) hassym <- log(1 - ymin)
	rbf <- function(x){
		ifelse(is.na(x) | abs(x)>cutoff, 0, exp(-x^2/(2*sigma^2) + hassym) + ymin)
	}
	return(rbf)
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
	require(monocle3)
	assertthat::assert_that(is(cds, "cell_data_set"))
	assertthat::assert_that(is.data.frame(reduced_coordinates) | is.matrix(reduced_coordinates))
	assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates), nrow(monocle3::pData(cds))))
	assertthat::assert_that(setequal(row.names(reduced_coordinates), colnames(cds)))
	assertthat::assert_that(assertthat::is.count(k) & k > 1)
	assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
	if(!is.null(summary_stats)) {
		assertthat::assert_that(all(summary_stats %in% names(monocle3::pData(cds))),
		                        msg = paste("One of your summary_stats is missing",
		                                    "from your pData table. Either add a",
		                                    "column with the name in",
		                                    "summary_stats, or remove the name",
		                                    "from the summary_stats parameter.",
		                                    collapse = " "))
		assertthat::assert_that(sum(vapply(summary_stats, function(x) {
		  !(is(monocle3::pData(cds)[,x], "numeric") | is(monocle3::pData(cds)[,x], "integer"))}, 1)) == 0,
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
	pdata <- monocle3::pData(cds)
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


#' projectMatrix_UMAP
#' 
#' given a (count) matrix and dimension reduction result, return the projected UMAP coordinates
#' in the embedding space
#'
#' @param X       matrix to be projected (features X samples)
#' @param umapObj dimension reduction result as returned by \code{\link{dimRed_UMAP}}
#' @param binarize binarize the counts before projecting
#' @param addPcCoord also add PC coordinates to the resulting matrix
#' @return Projected UMAP coordinates
#' @author Fabian Mueller
#' @export
projectMatrix_UMAP <- function(X, umapObj, binarize=TRUE, addPcCoord=FALSE){
	if (!is.element(class(umapObj), c("DimRed_UMAP_sc", "iterativeLSIResultSc"))){
		logger.error("Invalid dimension reduction object")
	}
	if (nrow(X) != length(umapObj$regionGr)){
		logger.error("Incompatible number of regions [nrow(X)]")
	}
	if (binarize) X <- X > 0

	tf <- Matrix::t(Matrix::t(X) / Matrix::colSums(X)) #term frequency
	tfidf <- tf * umapObj$idfBase # inverse document frequency
	tfidf <- as.matrix(tfidf)

	pcaCoord_proj <- (t(tfidf) %*% attr(umapObj$pcaCoord, "SVD_U"))[, umapObj$pcs]
	colnames(pcaCoord_proj) <- paste0('PC', 1:ncol(pcaCoord_proj))
	umapCoord_proj <- uwot::umap_transform(pcaCoord_proj, umapObj$umapRes)
	rownames(umapCoord_proj) <- rownames(pcaCoord_proj)
	colnames(umapCoord_proj) <- colnames(umapObj$umapCoord)
	if (addPcCoord){
		umapCoord_proj <- cbind(umapCoord_proj, pcaCoord_proj)
	}
	return(umapCoord_proj)
}

#' smoothMagic
#' 
#' [!EXPERIMENTAL!] smooth a matrix using an adaptation of the MAGIC algorithm (doi:10.1016/j.cell.2018.05.061)
#'
#' @param X       matrix to be smoothed (data points [cells] represent rows and features [genes] represent columns)
#' @param X_knn   matrix to be used for finding nearest neighbors. Must have the same rows (cells) as X
#' @param k       number of nearest neighbors to compute
#' @param ka      The adaptive kernel width (sigma) will be set to the distance to the \code{ka}th nearest neighbor for each cell. Thus the kernel is wider in sparse areas and smaller in dense areas.
#' @param td      diffusion time (number of exponentiations of the kernel matrix)
#' @return list containing the smoothed matrix and the kernel matrix
#' @author Jeff Granja, Fabian Mueller
#' @export
#' @noRd
smoothMagic <- function(X=NULL, X_knn=NULL, k=15, ka=ceiling(k/4), td=3){
	# X <- t(SummarizedExperiment::assay(geneAct))
	# X_knn <- dro$pcaCoord[,dro$pcs]
	require(Matrix)
	require(FNN)
	if (is.null(X_knn)) X_knn <- X
	if (is.null(X_knn)) logger.error("Expected matrix to compute nearest neighbors")
	if (nrow(X_knn) != nrow(X)) logger.error("Incompatible row numbers of X and X_knn")

	Ncells <- nrow(X)
	logger.status(c("Computing", k, "- nearest neighbors ..."))
	knnRes <- FNN::get.knnx(X_knn, X_knn, k=k)
	dd <- knnRes$nn.dist
	idx <- knnRes$nn.index

	if (ka > 0){
		logger.info(paste0("Normalizing distances to the ", ka, "-th nearest neighbor for each data point"))
		dd <- dd / dd[,ka]
	}
	logger.status("Computing affinity/kernel matrix ...")
	A <- Matrix::sparseMatrix(rep(seq_len(Ncells), k), c(idx), x=c(dd), dims=c(Ncells, Ncells))
	# make the matrix symmetric
	A <- A + Matrix::t(A)
	# compute kernel from (normalized) distances
	A@x <- exp(-A@x)
	# normalize the kernel to represent probabilities
	M <- A / Matrix::rowSums(A)
	logger.status("Diffusing kernel ...")
	Mt <- M
	if (td > 1){
		for(i in seq_len(td-1)) Mt <- Mt %*% M
	}

	if (is.null(X)){
		Xs <- NULL
	} else {
		logger.status("Smoothing ...")
		Xs <- Mt %*% X
		rownames(Xs) <- rownames(X)
		# logger.status("Rescaling ...")
	}
	
	res <- list(
		Xs=Xs,
		kernelM=Mt
	)
	class(res) <- "SmoothMagicMatrix"
	return(res)
}

#' findGroupMarkers
#' 
#' [EXPERIMENTAL] Detect group markers by testing distribution of values corresponding to a group (foreground) to the corresponding values from a background set.
#' The background set can be matched by finding nearest neighbors to the foreground.
#'
#' @param X       matrix with values to be differentially tested for each group (features X samples)
#' @param grouping Group assignments. Factor vector or a vector that can be coerced to a factor. The length must match the number of samples (columns of X)
#' @param bgRatio For each group the number of samples that will be drawn as background set. This should be a factor, i.e. a value of 2 means that the background set for each group will be twice the size of the foreground set
#' @param matchFrac fraction of samples from the background set that are drawn from the pool of nearest neighbors of the foreground. This way the degree of matching between background and foreground can be adjusted.
#' @param X_knn   matrix to be used for finding nearest neighbors. Must have the same columns (samples) as X. If \code{NULL} (default), \code{X} will be used
#' @param knn_k   k to be used for nearest-neighbor matching
#' @param testMethod test to be used. Currently valid options are \code{"wilcoxon"} and \code{"ttest"}
#' @param subsampleFrac Fraction of cells to consider for each comparison. If the size of the dataset exceeds
#'                this number, a random subset of cells will be considered for each comparison.
#' @param eps     pseudocount to add before calculating log fold-changes
#' @return list containing test statistics for each group and details on the foreground/background samples drawn for each group
#' @author Fabian Mueller
#' @export
#' @noRd
findGroupMarkers <- function(X, grouping, bgRatio=1.0, matchFrac=1.0, X_knn=NULL, knn_k=100, testMethod="wilcoxon", subsampleFrac=1.0, eps=1){
	require(FNN)
	if (!is.factor(grouping)){
		grouping <- factor(grouping)
	}
	if (length(grouping) != ncol(X)){
		logger.error("Dimensions of grouping and matrix do not match")
	}
	if (is.element(testMethod, c("ttest", "wilcoxon")) && !is.matrix(X)){
		logger.warning("t-test and Wilcoxon test require a matrix object as input. --> X will be converted to dense")
	}

	if (is.null(X_knn)) X_knn <- X
	if (is.null(X_knn)) logger.error("Expected matrix to compute nearest neighbors")
	if (ncol(X_knn) != ncol(X)) logger.error("Incompatible column numbers of X and X_knn")
	X_knn <- t(X_knn)

	diffFun <- NULL
	if (testMethod=="ttest"){
		require(matrixTests)
		diffFun <- function(X, idx1, idx2){
			X1 <- X[,idx1]
			if (!is.matrix(X1)) X1 <- as.matrix(X1)
			X2 <- X[,idx2]
			if (!is.matrix(X2)) X2 <- as.matrix(X2)
			testRes <- matrixTests::row_t_welch(X1, X2)
			rr <- data.frame(
				pval  = testRes$pvalue,
				pval_fdr = p.adjust(testRes$pvalue, method="fdr"),
				stat  = testRes$statistic,
				mean1 = testRes$mean.x,
				mean2 = testRes$mean.y,
				var1  = testRes$var.x,
				var2  = testRes$var.y,
				n1    = testRes$obs.x,
				n2    = testRes$obs.y,
				df    = testRes$df,
				conf.low = testRes$conf.low,
				conf.high = testRes$conf.high
			)
			if (!is.null(rownames(X))) rr[,"name"] <- rownames(X)
			return(rr)
		}
	} else if (testMethod=="wilcoxon"){
		require(matrixTests)
		require(matrixStats)
		diffFun <- function(X, idx1, idx2){
			X1 <- X[,idx1]
			if (!is.matrix(X1)) X1 <- as.matrix(X1)
			X2 <- X[,idx2]
			if (!is.matrix(X2)) X2 <- as.matrix(X2)
			testRes <- matrixTests::row_wilcoxon_twosample(X1, X2)
			rr <- data.frame(
				pval  = testRes$pvalue,
				pval_fdr = p.adjust(testRes$pvalue, method="fdr"),
				stat  = testRes$statistic,
				mean1 = rowMeans(X1, na.rm=TRUE),
				mean2 = rowMeans(X2, na.rm=TRUE),
				var1  = rowVars(X1, na.rm=TRUE),
				var2  = rowVars(X2, na.rm=TRUE),
				n1    = testRes$obs.x,
				n2    = testRes$obs.y
			)
			if (!is.null(rownames(X))) rr[,"name"] <- rownames(X)
			return(rr)
		}
	} else {
		logger.error(c("Unknown differential test method:", testMethod))
	}
	allIdx <- seq_len(ncol(X))
	gIdxL <- tapply(seq_along(grouping), grouping, c)
	nGrps <- nlevels(grouping)

	logger.start("Sampling background")
		# sample background using nearest neighbors and random drawing
		idxL <- lapply(1:nGrps, FUN=function(i){
			g <- levels(grouping)[i]
			gIdx <- gIdxL[[g]]
			logger.status(c("Group:", g, paste0("(", i, " of ", nGrps, ")")))
			logger.info(c("Group size:", length(gIdx)))
			oogIdx <- allIdx[-gIdxL[[g]]] # out-of-group indices
			gN <- length(gIdx)
			if (subsampleFrac > 0 && subsampleFrac < 1){
				nSample <- ceiling(subsampleFrac * gN)
				logger.info(paste0("Subsampling to ", nSample, " of ", gN, " cells (", round(100*nSample/gN, 2), "%)"))
				gIdx <- sort(sample(gIdx, nSample))
				gN <- length(gIdx)
			}
			bgIdx <- c()
			bgN <- min(length(oogIdx), ceiling(bgRatio * gN))
			# background samples from nearest neighbor matches
			match_k <- 0
			if (matchFrac > 0) match_k <- max(0, min(bgN, ceiling(matchFrac * bgRatio * gN)))
			k <- min(knn_k, match_k)
			if (match_k > 0){
				knnRes <- FNN::get.knnx(X_knn[oogIdx,], X_knn[gIdx,], k=k)
				matchIdx <- oogIdx[unique(sort(as.vector(knnRes$nn.index)))]
				if (length(matchIdx) > match_k) matchIdx <- sort(sample(matchIdx, match_k))
				bgIdx <- sort(union(bgIdx, matchIdx))
			}
			if (length(bgIdx) > 0) logger.info(c("Matched samples:", length(bgIdx)))
			# fill remaining samples with random samples
			rand_k <- max(0, bgN - length(bgIdx))
			if (rand_k > 0){
				ss <- setdiff(oogIdx, bgIdx)
				rand_k <- min(rand_k, length(ss))
				rIdx <- sort(sample(ss, rand_k))
				bgIdx <- sort(union(bgIdx, rIdx))
				logger.info(c("Random samples:", rand_k))
			}
			return(list(group=gIdx, background=bgIdx))
		})
		names(idxL) <- levels(grouping)
	logger.completed()

	logger.start("Differential test")
		testL <- lapply(1:nGrps, FUN=function(i){
			g <- levels(grouping)[i]
			logger.status(c("Group:", g, paste0("(", i, " of ", nGrps, ")")))
			rr <- diffFun(X, idxL[[g]]$group, idxL[[g]]$background)
			rr[,"group"] <- factor(g, levels=levels(grouping))
			rr[,"log2FC"] <- log2((rr[,"mean1"] + eps)/(rr[,"mean2"] + eps))
			return(rr)
		})
		names(testL) <- levels(grouping)

	logger.completed()

	res <- list(
		testRes = testL,
		grouping = list(
			grouping = grouping,
			groupIdx = idxL
		)
	)
	class(res) <- "groupMarkers"
	return(res)
}

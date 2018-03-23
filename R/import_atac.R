#-------------------------------------------------------------------------------
#' DsATAC.snakeATAC
#' 
#' Create a DsATAC dataset from multiple input files output by snakeATAC
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param filePrefixCol column name specifying the file prefix for each sample in the sample annotation table
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param type         input data type. Currently only "insBed" (insertion beds) and "bam" (aligned reads) are supported
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object. Only relevant when \code{type=="insBam"}.
#' @param pairedEnd    is the input data paired-end? Only relevant when \code{type=="insBam"}.
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, regionSets=NULL, sampleIdCol=filePrefixCol, type="insBed", keepInsertionInfo=TRUE, pairedEnd=TRUE){
	if (!is.element(type, c("bam", "insBam", "insBed"))){
		logger.error(c("Unsupported import type:", type))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (is.element(type, c("bam", "insBam"))){
		require(GenomicAlignments)
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], ".noMT.filtered.deduped.bam"))
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	} else if (type=="insBed"){
		require(rtracklayer)
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], ".insertions.bed.gz"))
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	}

	if (!all(file.exists(inputFns))){
		missingSamples <- sampleIds[!file.exists(inputFns)]
		logger.error(c("Missing input files for samples:", paste(missingSamples, collapse=", ")))
	}
	if (is.null(regionSets)){
		logger.info(c("Using default region sets:", paste("tiling200bp", collapse=", ")))
		regionSets <- GRangesList(
			tiling200bp=getTilingRegions(genome, width=200L, onlyMainChrs=TRUE)
		)
	}
	if (length(names(regionSets)) < 1){
		logger.error("Region sets must be named")
	}
	if (length(regionSets) < 1){
		logger.warning("No region sets specified")
	}
	logger.start("Creating DsATAC object")
		obj <- DsATAC(sampleAnnot, genome)
		for (rt in names(regionSets)){
			logger.info(c("Including region set:", rt))
			obj <- regionAggregation(obj, regionSets[[rt]], rt, signal=NULL, dropEmpty=FALSE)
		}
	logger.completed()

	logger.start("Reading ATAC data")
		nSamples <- length(sampleIds)
		if (type=="bam"){
			# print(inputFns)
			obj <- addCountDataFromBam(obj, inputFns)
		} else if (type=="insBam"){
			logger.start(c("Adding insertion data from bam"))
				obj <- addInsertionDataFromBam(obj, inputFns, pairedEnd=pairedEnd)
			logger.completed()
			logger.start(c("Summarizing region counts"))
				obj <- addCountDataFromGRL(obj, getInsertionSites(obj))
			logger.completed()
			# optionally remove insertion information to save space
			if (!keepInsertionInfo){
				obj@insertions <- list()
			}
		} else if (type=="insBed"){
			for (i in seq_along(sampleIds)){
				sid <- sampleIds[i]
				logger.status(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
				fn <- inputFns[sid]
				grl <- GRangesList(import(fn, format = "BED"))
				grl[[1]] <- setGenomeProps(grl[[1]], genome, onlyMainChrs=TRUE)
				names(grl) <- sid
				obj <- addCountDataFromGRL(obj, grl)
				rm(grl)
			}
		}
	logger.completed()
	
	return(obj)
}

#-------------------------------------------------------------------------------
#' getNonOverlappingByScore
#' 
#' Retrieve the set of non-verlapping regions by iteratively picking the region with maximum score for
#' each set of consecutively overlapping regions
#' @param gr           \code{GRanges} object
#' @param scoreCol     name of the column to be used as scor in the \code{elementMetadata} of the \code{gr} object
#' @return \code{GRanges} object containing non-overlapping regions
#' @author Fabian Mueller
#' @noRd
getNonOverlappingByScore <- function(gr, scoreCol="score"){
	gr.rem <- gr

	res <- GRanges()
	seqlevels(res) <- seqlevels(gr)
	seqlengths(res) <- seqlengths(gr)
	genome(res) <- genome(gr)
	i <- 0
	while (length(gr.rem) > 0){
		i <- i + 1
		# logger.status(c("iteration", i)) #DEBUG
		scs <- elementMetadata(gr.rem)[,scoreCol]
		gr.merged <- reduce(gr.rem, min.gapwidth=0L, with.revmap=TRUE, ignore.strand=TRUE)
		# maxScoreIdx <- sapply(gr.merged, FUN=function(x){
		# 	idx <- elementMetadata(x)[,"revmap"][[1]]
		# 	return(idx[which.max(scs[idx])])
		# }) #too slow
		maxScoreIdx <- sapply(elementMetadata(gr.merged)[,"revmap"], FUN=function(idx){
			return(idx[which.max(scs[idx])])
		})
		bait <- gr.rem[maxScoreIdx]
		res <- c(res, bait)
		gr.rem <- gr.rem[!overlapsAny(gr.rem, bait, ignore.strand=TRUE)]
		# logger.info(c(length(gr.rem), "regions left")) #DEBUG
	}
	return(res)
}
#-------------------------------------------------------------------------------
# #' filterUnreplicatedRegions
# #' 
# #' Remove regions that are not covered in a sufficient number of members
# #' @param grl          \code{GRangesList} object containing one element for each group member
# #' @param percReq      percentile required to keep it.
# #'                     E.g. a value of 1 (default) means that all members of a group are required to contain that region in order
# #'                     to keep it.
# #' @return \code{GRangesList} object with regions with insufficient coverage removed
# #' @author Fabian Mueller
# #' @noRd
# filterUnreplicatedRegions <- function(grl, percReq=1.0){
# 	nReq <- as.integer(ceiling(percReq * length(grl)))
# 	grlu <- unlist(grl, use.names=FALSE)
# 	elementMetadata(grlu)[,".sIdx"] <- as.integer(NA) # add a column with the sample index
# 	elementMetadata(grlu)[,".rIdx"] <- as.integer(NA) # add a column with the region index
# 	grl <- relist(grlu, grl)
# 	rm(grlu)
# 	for (i in seq_along(grl)){
# 		elementMetadata(grl[[i]])[,".sIdx"] <- i 
# 		elementMetadata(grl[[i]])[,".rIdx"] <- seq_along(grl[[i]]) 
# 	}
# 	gr <- unlist(grl)
# 	grm <- reduce(gr, min.gapwidth=0L, with.revmap=TRUE, ignore.strand=TRUE) #merge overlapping regions
# 	revmap <- mcols(grm)$revmap
# 	# gr.relist <- relist(gr[unlist(revmap)], revmap)

# 	sidxL <- relist(mcols(gr)[unlist(revmap), ".sIdx"], revmap) # list of sample indices per merged region
# 	keepm <- sapply(sidxL, FUN=function(x){
# 		res <- length(unique(unlist(x))) >= nReq
# 		return(res)
# 	})
# 	nClusters <- length(sidxL)
# 	nFiltered <- sum(!keepm)
# 	logger.info(c("Removed", nFiltered, "of", nClusters, paste0("(", round(nFiltered/nClusters*100, 2), "%)"), "region clusters covered less than", nReq, "times"))
# 	grm <- grm[keepm]
# 	idx.all <- sort(unique(unlist(mcols(grm)$revmap))) #vector of all indices
# 	gr <- gr[idx.all] # keep only regions which are in the filtered clusters

# 	# map filtered regions back to the original indices in the input GRangesList
# 	for (i in seq_along(grl)) {
# 		keep.region.idx <- sort(elementMetadata(gr)[elementMetadata(gr)[,".sIdx"]==i, ".rIdx"])
# 		grl[[i]] <- grl[[i]][keep.region.idx]
# 	}

# 	#remove the added columns
# 	grlu <- unlist(grl, use.names=FALSE)
# 	elementMetadata(grlu)[,".sIdx"] <- NULL
# 	elementMetadata(grlu)[,".rIdx"] <- NULL
# 	grl <- relist(grlu, grl)

# 	return(grl)
# }
#-------------------------------------------------------------------------------
#' getPeakSet.snakeATAC
#' 
#' Retrieve a consensus set of ATAC peaks from the snakeATAC pipline run
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param filePrefixCol column name specifying the file prefix for each sample in the sample annotation table
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param type         input data type. Currently only "summits_no_fw" (non-overlapping, fixed-width peaks deduced from summits)
#' @param unifWidth    width of the peaks if the results have uniform peak lengths
#' @param replicateCol column name specifying the replicate group for cross-checking coverage across replicates
#' @param replicatePercReq percentile of replicates in a group required to contain a peak in order to keep it.
#'                     E.g. a value of 1 (default) means that all replicates in a group are required to contain that peak in order
#'                     to keep it.
#' @param keepOvInfo   keep annotation columns in the elementMetadata of the results specifying whether a consensus peak overlaps with a
#'                     peak in each sample
#' @return \code{GRanges} object containing consensus peak set
#' @author Fabian Mueller
#' @export
getPeakSet.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, sampleIdCol=filePrefixCol, type="summits_no_fw", unifWidth=500L, replicateCol=NA, replicatePercReq=1.0, keepOvInfo=FALSE){
	if (!is.element(type, c("summits_no_fw"))){
		logger.error(c("Unsupported import type:", type))
	}
	if (replicatePercReq > 1 || replicatePercReq < 0){
		logger.error(c("Invalid value for replicatePercReq. Must be in [0,1]"))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (type=="summits_no_fw"){
		require(rtracklayer)
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], "_summits.bed"))
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	}
	if (!all(file.exists(inputFns))){
		missingSamples <- sampleIds[!file.exists(inputFns)]
		logger.error(c("Missing input files for samples:", paste(missingSamples, collapse=", ")))
	}

	peakFun <- NULL
	res <- NULL
	if (type=="summits_no_fw"){
		peakFun <- function(fn){
			rr <- import(fn, format="BED")
			rr <- setGenomeProps(rr, genome, onlyMainChrs=TRUE)
			rr <- rr[isCanonicalChrom(as.character(seqnames(rr)))]
			# scale scores to their percentiles
			scs <- elementMetadata(rr)[,"score"]
			elementMetadata(rr)[,"score_norm"] <- ecdf(scs)(scs)
			elementMetadata(rr)[,"sampleId"] <- sid
			#extend peaks around the summit
			rr <- trim(promoters(rr, upstream=ceiling(unifWidth/2), downstream=ceiling(unifWidth/2)+1)) #extend each summit on each side by half the width
			rr <- rr[width(rr)==median(width(rr))] #remove too short regions which might have been trimmed
			return(rr)
		}
	}

	for (sid in sampleIds){
		logger.status(c("Reading peak summits from sample:", sid))
		peakSet.cur <- peakFun(inputFns[sid])
		
		#add coverage info for all samples
		elementMetadata(peakSet.cur)[,paste0(".ov.", sampleIds)] <- FALSE # as.logical(NA)

		if (is.null(res)){
			#remove overlapping peaks in initial sample based on normalized scores
			peakSet.cur <- getNonOverlappingByScore(peakSet.cur, scoreCol="score_norm")
			#initialize peak set with all peaks from the first sample
			elementMetadata(peakSet.cur)[,paste0(".ov.", sid)] <- TRUE
			res <- peakSet.cur
		} else {
			# add new peaks and remove the overlapping ones by taking the peaks with the best score
			res <- getNonOverlappingByScore(c(res, peakSet.cur), scoreCol="score_norm")
			elementMetadata(res)[,paste0(".ov.", sid)] <- overlapsAny(res, peakSet.cur, ignore.strand=TRUE)
		}
	}

	if (is.element(replicateCol, colnames(sampleAnnot)) && replicatePercReq>0){
		logger.start("Accounting for peak reproducibility across replicates")
			groupF <- factor(sampleAnnot[,replicateCol])
			gRepMat <- matrix(as.logical(NA), nrow=length(res), ncol=nlevels(groupF))
			colnames(gRepMat) <- levels(groupF)
			for (gg in levels(groupF)){
				sidsRepl <- sampleIds[groupF==gg]
				ovMat <- as.matrix(elementMetadata(res)[,paste0(".ov.", sidsRepl)])
				nReq <- as.integer(ceiling(replicatePercReq * length(sidsRepl)))
				gRepMat[,gg] <- rowSums(ovMat) >= nReq
			}
			keep <- rowAnys(gRepMat)
			nTotal <- length(res)
			nRem <- nTotal - sum(keep)
			logger.info(c("Removed", nRem, "of", nTotal, "peaks", paste0("(",round(nRem/nTotal*100, 2),"%)"), "because they were not reproduced by more than", round(replicatePercReq*100, 2), "% of samples in any replicate group"))
			res <- res[keep]
		logger.completed()
	}
	if (!keepOvInfo){
		#remove the helper columns
		for (sid in sampleIds){
			elementMetadata(res)[,paste0(".ov.", sid)] <- NULL
		}
	}
	#sort
	res <- sortSeqlevels(res)
	res <- sort(res)
	return(res)
}

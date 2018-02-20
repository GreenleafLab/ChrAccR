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
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, regionSets=NULL, sampleIdCol=filePrefixCol, type="insBed"){
	if (!is.element(type, c("bam", "insBed"))){
		logger.error(c("Unsupported import type:", type))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (type=="bam"){
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
		if (type=="bam"){
			# print(inputFns)
			obj <- addCountDataFromBam(obj, inputFns)
		} else if (type=="insBed"){
			for (sid in sampleIds){
				logger.status(c("Importing sample:", sid))
				fn <- inputFns[sid]
				grl <- GRangesList(import(fn, format = "BED", genome=genome))
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
		maxScoreIdx <- sapply(gr.merged, FUN=function(x){
			idx <- elementMetadata(x)[,"revmap"][[1]]
			return(idx[which.max(scs[idx])])
		})
		bait <- gr.rem[maxScoreIdx]
		res <- c(res, bait)
		gr.rem <- gr.rem[!overlapsAny(gr.rem, bait, ignore.strand=TRUE)]
		logger.info(c(length(gr.rem), "regions left")) #DEBUG
	}
	return(res)
}

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
#' @return \code{GRanges} object containing consensus peak set
#' @author Fabian Mueller
#' @export
getPeakSet.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, sampleIdCol=filePrefixCol, type="summits_no_fw", unifWidth=500L){
	if (!is.element(type, c("summits_no_fw"))){
		logger.error(c("Unsupported import type:", type))
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

	res <- NULL
	if (type=="summits_no_fw"){
		grl <- list()
		for (sid in sampleIds){
			logger.status(c("Reading peak summits from sample:", sid))
			fn <- inputFns[sid]
			gr.cur <- import(fn, format="BED", genome=genome)
			gr.cur <- gr.cur[isCanonicalChrom(as.character(seqnames(gr.cur)))]
			# scale scores to their percentiles
			scs <- elementMetadata(gr.cur)[,"score"]
			elementMetadata(gr.cur)[,"score_norm"] <- ecdf(scs)(scs)
			elementMetadata(gr.cur)[,"sampleId"] <- sid
			grl[[sid]] <- gr.cur
		}
		grl <- GRangesList(grl)
		gr <- unlist(grl)
		gr <- trim(promoters(gr, upstream=ceiling(unifWidth/2), downstream=ceiling(unifWidth/2)+1)) #extend each summit on each side by 250bp
		gr <- gr[width(gr)==median(width(gr))] #remove too short regions which might have been trimmed

		logger.start("Retrieving consensus peaks")
		res <- getNonOverlappingByScore(gr, scoreCol="score_norm")
		logger.completed()
		#sort
		res <- sortSeqlevels(res)
		res <- sort(res)
	}
	return(res)
}

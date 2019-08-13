#-------------------------------------------------------------------------------
#' DsATAC.snakeATAC
#' 
#' Create a DsATAC dataset from multiple input files output by snakeATAC
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param filePrefixCol column name specifying the file prefix for each sample in the sample annotation table. If \code{dataDir} is not empty (i.e. not \code{""})
#'                     filenames are assumed to be relative to that directory and  a corresponding filename suffix will be appended
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located. If it is the empty character (\code{""}; default) it is assumed that \code{filePrefixCol} specifies the full path to the input files
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param type         input data type. Currently only "insBed" (insertion beds), "insBed" (insertion info inferred from bam files (aligned reads); default) and "bam" (aligned reads) are supported
#' @param diskDump     should large data objects (count matrices, fragment data, ...) be disk-backed to save main memory
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object. Only relevant when \code{type=="insBam"}.
#' @param bySample     process sample-by-sample to save memory (currently only has an effect for \code{type=="insBam"})
#' @param pairedEnd    is the input data paired-end? Only relevant when \code{type=="insBam"}.
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir="", regionSets=NULL, sampleIdCol=filePrefixCol, type="insBam", diskDump=FALSE, keepInsertionInfo=TRUE, bySample=FALSE, pairedEnd=TRUE){
	if (diskDump){
		if (!requireNamespace("DelayedArray") || !requireNamespace("HDF5Array")) logger.error(c("Could not load dependency: DelayedArray, HDF5Array"))
	}
	if (!is.element(type, c("bam", "insBam", "insBed"))){
		logger.error(c("Unsupported import type:", type))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (is.element(type, c("bam", "insBam"))){
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], ".noMT.filtered.deduped.bam"))
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	} else if (type=="insBed"){
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], ".insertions.bed.gz"))
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	} else if (type=="covBigWig"){
		if (nchar(dataDir) > 0){
			inputFns <- file.path(dataDir, paste0(sampleAnnot[,filePrefixCol], ".100bp_coverage.bw"))
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
		obj <- DsATAC(sampleAnnot, genome, diskDump=diskDump, diskDump.fragments=keepInsertionInfo)
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
			if (bySample){
				if (obj@diskDump && obj@diskDump.fragments){
					# when disk-dumping, writeHDF5Array can be painfully slow for large DelayedArray objects with many indices/operations that need to be realized before saving.
					# this is a workaround, iterating over each sample and directly writing to disk-based HDF5 matrix
					logger.start(c("Adding insertion data from fragment files"))
						for (i in seq_along(inputFns)){
							sid <- names(inputFns)[i]
							logger.status(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
							obj <- addInsertionDataFromBam(obj, inputFns[i], pairedEnd=pairedEnd, .diskDump=obj@diskDump.fragments)
						}
					logger.completed()
					# logger.start("[DEBUG] tmp saving DsATAC object")
					# 	saveDsAcc(obj, file.path("/scratch/users/muellerf/temp", getHashString("DsATAC_tmp")))
					# logger.completed()
					logger.start("Agregating region count data")
						rebe <- DelayedArray::getRealizationBackend() # store previous realization backend setting to be able to reset it later
						DelayedArray::setRealizationBackend("HDF5Array")
						rTypes <- getRegionTypes(obj)
						rSinkL <- lapply(rTypes, FUN=function(rt){
							rSink <- DelayedArray::RealizationSink(as.integer(c(getNRegions(obj, rt), nSamples)), type="integer")
							return(list(
								sink=rSink,
								grid=DelayedArray::colGrid(rSink, ncol=1L)
							))
						})
						names(rSinkL) <- rTypes
						
						for (j in seq_along(sampleIds)){
							sid <- sampleIds[j]
							logger.start(c("Summarizing counts for sample", ":", sid, paste0("(", j, " of ", nSamples, ")")))
								tmpDs <- addCountDataFromGRL(obj, getInsertionSites(obj, samples=sid))
								for (rt in rTypes){
									cm <- getCounts(tmpDs, rt, j=j, asMatrix=TRUE)
									DelayedArray::write_block(rSinkL[[rt]]$sink, rSinkL[[rt]]$grid[[j]], cm)
								}
								rm(tmpDs)
							logger.completed()
						}
						for (rt in rTypes){
							DelayedArray::close(rSinkL[[rt]]$sink)
							obj@counts[[rt]] <- as(rSinkL[[rt]]$sink, "DelayedArray")
							colnames(obj@counts[[rt]]) <- sampleIds
						}
						DelayedArray::setRealizationBackend(rebe) # reset realization backend to previous setting
					logger.completed()
				} else {
					logger.start(c("Adding insertion and region count data from bam"))
						for (i in seq_along(inputFns)){
							sid <- names(inputFns)[i]
							logger.start(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
								obj <- addInsertionDataFromBam(obj, inputFns[i], pairedEnd=pairedEnd, .diskDump=obj@diskDump.fragments)
								obj <- addCountDataFromGRL(obj, getInsertionSites(obj, samples=sid))
								# optionally remove insertion information to save space
								if (!keepInsertionInfo){
									obj@fragments <- list()
								}
							logger.completed()
						}
					logger.completed()
				}
			} else {
				logger.start(c("Adding insertion data from bam"))
					obj <- addInsertionDataFromBam(obj, inputFns, pairedEnd=pairedEnd, .diskDump=obj@diskDump.fragments)
				logger.completed()
				logger.start(c("Summarizing region counts"))
					obj <- addCountDataFromGRL(obj, getInsertionSites(obj))
				logger.completed()
				# optionally remove insertion information to save space
				if (!keepInsertionInfo){
					obj@fragments <- list()
				}
			}
		} else if (type=="insBed"){
			for (i in seq_along(sampleIds)){
				sid <- sampleIds[i]
				logger.status(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
				fn <- inputFns[sid]
				grl <- GRangesList(rtracklayer::import(fn, format = "BED"))
				grl[[1]] <- setGenomeProps(grl[[1]], genome, onlyMainChrs=TRUE)
				names(grl) <- sid
				obj <- addCountDataFromGRL(obj, grl)
				rm(grl)
			}
		} else if (type=="covBigWig"){
			# TODO: check if this actually works
			for (i in seq_along(sampleIds)){
				sid <- sampleIds[i]
				logger.status(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
				fn <- inputFns[sid]
				grl <- GRangesList(rtracklayer::import(fn, format = "BigWig"))
				grl[[1]] <- setGenomeProps(grl[[1]], genome, onlyMainChrs=TRUE)
				names(grl) <- sid
				obj <- addSignalDataFromGRL(obj, grl, aggrFun=function(x){mean(x, na.rm=TRUE)})
				rm(grl)
			}
		}
	logger.completed()
	
	return(obj)
}

#-------------------------------------------------------------------------------
#' DsATAC.bam
#' 
#' Create a DsATAC dataset from multiple input bam files
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param bamFiles     either a character vector of the same length as sampleAnnot has rows, specifying the file paths of the bam files for each
#'                     sample or a single character string specifying the column name in \code{sampleAnnot} where the file paths can be found
#' @param genome       genome assembly
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name in the sample annotation table containing unique sample identifiers. If \code{NULL} (default), the function will look for a column that contains the word "sample"
#' @param diskDump     should large data objects (count matrices, fragment data, ...) be disk-backed to save main memory
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object. Only relevant when \code{type=="insBam"}.
#' @param pairedEnd    is the input data paired-end? Only relevant when \code{type=="insBam"}.
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' # download and unzip the dataset
#' datasetUrl <- "https://s3.amazonaws.com/muellerf/data/ChrAccR/data/tutorial/tcells.zip"
#' downFn <- "tcells.zip"
#' download.file(datasetUrl, downFn)
#' unzip(downFn, exdir=".")
#' # prepare the sample annotation table
#' sampleAnnotFn <- file.path("tcells", "samples.tsv")
#' bamDir <- file.path("tcells", "bam")
#' sampleAnnot <- read.table(sampleAnnotFn, sep="\t", header=TRUE, stringsAsFactors=FALSE)
#' # add a column that ChrAccR can use to find the correct bam file for each sample
#' sampleAnnot[,"bamFilenameFull"] <- file.path(bamDir, sampleAnnot[,"bamFilename"])
#' # prepare the dataset
#' dsa_fromBam <- DsATAC.bam(sampleAnnot, "bamFilenameFull", "hg38", regionSets=NULL, sampleIdCol="sampleId")
#' }
DsATAC.bam <- function(sampleAnnot, bamFiles, genome, regionSets=NULL, sampleIdCol=NULL, diskDump=FALSE, keepInsertionInfo=TRUE, pairedEnd=TRUE){
	if (!is.character(bamFiles)) logger.error("Invalid value for bamFiles. Expected character")
	if (length(bamFiles)==1 && is.element(bamFiles, colnames(sampleAnnot))){
		bamFnCol <- bamFiles
	} else {
		if (length(bamFiles)!=nrow(sampleAnnot)) logger.error("Invalid value for bamFiles. must match the number of rows in sampleAnnot")
		bamFnCol <- ".bamPath"
		sampleAnnot[,bamFnCol] <- bamFiles
	}
	if (!is.null(sampleIdCol) && !(is.character(sampleIdCol) && length(sampleIdCol)==1)) logger.error("Invalid value for sampleIdCol Expected character or NULL")
	if (is.null(sampleIdCol)){
		sampleIdCol <- grep("sample", colnames(sampleAnnot), value=TRUE)
		if (length(sampleIdCol) != 1) logger.error("Could not uniquely determine sample id column from sample annotation column names")
		logger.info(c("Automatically determined column '", sampleIdCol, "' as sample identifier column"))
	}
	DsATAC.snakeATAC(sampleAnnot, bamFnCol, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, type="insBam", diskDump=diskDump, keepInsertionInfo=keepInsertionInfo, bySample=TRUE, pairedEnd=pairedEnd)
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
#' @param replicateCol column name specifying the replicate group for cross-checking coverage across replicates
#' @param replicatePercReq percentile of replicates in a group required to contain a peak in order to keep it.
#'                     E.g. a value of 1 (default) means that all replicates in a group are required to contain that peak in order
#'                     to keep it.
#' @param replicateConsSelect if set, the peak set will also be checked for consistency, i.e. in order to retain a peak
#'                     it has to be consistently be present or absent in each replicate group (as specified in \code{replicatePercReq} percent of samples)
#' @param keepOvInfo   keep annotation columns in the elementMetadata of the results specifying whether a consensus peak overlaps with a
#'                     peak in each sample
#' @return \code{GRanges} object containing consensus peak set
#' @author Fabian Mueller
#' @export
getPeakSet.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, sampleIdCol=filePrefixCol, type="summits_no_fw", unifWidth=500L, replicateCol=NA, replicatePercReq=1.0, replicateConsSelect=FALSE, keepOvInfo=FALSE){
	if (!is.element(type, c("summits_no_fw", "summits_filt_no_fw"))){
		logger.error(c("Unsupported import type:", type))
	}
	if (replicatePercReq > 1 || replicatePercReq < 0){
		logger.error(c("Invalid value for replicatePercReq. Must be in [0,1]"))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds

	inputFns <- c()
	if (is.element(type, c("summits_no_fw", "summits_filt_no_fw"))){
		if (nchar(dataDir) > 0){
			fns <- NULL
			if (type=="summits_no_fw"){
				fns <- paste0(sampleAnnot[,filePrefixCol], "_summits.bed")
			} else if (type=="summits_filt_no_fw"){
				fns <- paste0(sampleAnnot[,filePrefixCol], "_peaks.narrowPeak")
			}
			inputFns <- file.path(dataDir, fns)
		} else {
			inputFns <- sampleAnnot[,filePrefixCol]
		}
		names(inputFns) <- sampleIds
	}
	if (!is.character(inputFns)) inputFns <- as.character(inputFns)
	if (!all(file.exists(inputFns))){
		missingSamples <- sampleIds[!file.exists(inputFns)]
		logger.error(c("Missing input files for samples:", paste(missingSamples, collapse=", ")))
	}

	peakFun <- NULL
	res <- NULL
	if (type=="summits_no_fw"){
		peakFun <- function(fn, sid){
			rr <- rtracklayer::import(fn, format="BED")
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
	} else if (type=="summits_filt_no_fw"){
		peakFun <- function(fn, sid){
			rr <- readMACS2peakFile(fn)
			rr <- setGenomeProps(rr, genome, onlyMainChrs=TRUE)
			rr <- rr[isCanonicalChrom(as.character(seqnames(rr)))]
			elementMetadata(rr)[,"calledPeakStart"] <- start(rr)
			elementMetadata(rr)[,"calledPeakEnd"] <- end(rr)
			start(rr) <- end(rr) <- elementMetadata(rr)[,"summit"]

			scs <- elementMetadata(rr)[,"negLog10qval"]
			elementMetadata(rr)[,"score_norm"] <- ecdf(scs)(scs)
			elementMetadata(rr)[,"sampleId"] <- sid

			# only retain peaks with q-value < 0.01
			rr <- rr[elementMetadata(rr)[,"negLog10qval"] > -log10(0.01)]

			rr <- trim(promoters(rr, upstream=ceiling(unifWidth/2), downstream=ceiling(unifWidth/2)+1)) #extend each summit on each side by half the width
			rr <- rr[width(rr)==median(width(rr))] #remove too short regions which might have been trimmed
			return(rr)
		}
	}

	logger.start("Reading peak sets")
		peakGrl <- lapply(sampleIds, FUN=function(sid){
			logger.status(c("sample:", sid))
			return(peakFun(inputFns[sid], sid))
		})
		names(peakGrl) <- sampleIds
	logger.completed()

	grps <- NULL
	if (is.element(replicateCol, colnames(sampleAnnot)) && replicatePercReq>0){
		grps <- sampleAnnot[,replicateCol]
	}

	res <- getConsensusPeakSet(peakGrl, mode="no_by_score", grouping=grps, groupAgreePerc=replicatePercReq, groupConsSelect=replicateConsSelect, scoreCol="score_norm", keepOvInfo=keepOvInfo)
	return(res)
}
#-------------------------------------------------------------------------------
#' getSampleMetrics.snakeATAC
#' 
#' Retrieve sample summary statistics from the output a snakeATAC pipline run
#' @param sampleAnnot  data.frame specifying the sample annotation table. Must have valid rownames corresponding to the sample ids used in
#'                     the snakeAtac filenames
#' @param snakeDir     snakeATAC base directory (where the files are located)
#' @param withPeaks    flag indicating whether to output peak statistics
#' @return \code{data.frame} containing sample summary statistics. the original sample annotation table will be appended to the summary output
#' @author Fabian Mueller
#' @export
getSampleMetrics.snakeATAC <- function(sampleAnnot, snakeDir, withPeaks=TRUE){
	if (length(rownames(sampleAnnot)) != nrow(sampleAnnot)){
		logger.error(c("sampleAnnot does not contain valid rownames"))
	}
	if (all(rownames(sampleAnnot)==as.character(1:nrow(sampleAnnot)))){
		logger.warning(c("sampleAnnot has numerical rownames. Should it?"))
	}

	#compiled countTables - disjoint
	cct <- readTab(file.path(snakeDir, "output", "bams", "qc", "compiled_counts.disjoint.txt"))
	colnames(cct) <- c("sampleId", paste0("reads_", colnames(cct)[2:5]))
	rownames(cct) <- cct$sampleId
	#initialize metrics summary table
	sampleMetrics <- cct
	unknownSamples <- setdiff(cct$sampleId, rownames(sampleAnnot))
	if (length(unknownSamples) > 0){
		logger.error(c("The following sample names (ids) could not be found in the sample annotation:", paste(unknownSamples, collapse=", ")))
	}

	#compiled countTables - stepwise
	cct <- readTab(file.path(snakeDir, "output", "bams", "qc", "compiled_counts.txt"))
	colnames(cct) <- c("sampleId", paste0("reads_", colnames(cct)[2:5]))
	rownames(cct) <- cct$sampleId
	sampleMetrics <- data.frame(
		sampleMetrics,
		cct[sampleMetrics$sampleId,2:ncol(cct)],
		stringsAsFactors=FALSE
	)

	#compiled countTables - percent
	cct <- readTab(file.path(snakeDir, "output", "bams", "qc", "compiled_counts.fraction.txt"))
	colnames(cct) <- c("sampleId", paste0("readFrac_", colnames(cct)[2:4]))
	rownames(cct) <- cct$sampleId
	sampleMetrics <- data.frame(
		sampleMetrics,
		cct[sampleMetrics$sampleId,2:ncol(cct)],
		stringsAsFactors=FALSE
	)

	#get flagstats numbers
	cct <- readTab(file.path(snakeDir, "output", "bams", "qc", "compiled_flagstats.txt"))
	colnames(cct) <- c("sampleId", paste0("flagstat_", colnames(cct)[2:ncol(cct)]))
	rownames(cct) <- cct$sampleId
	sampleMetrics <- data.frame(
		sampleMetrics,
		cct[sampleMetrics$sampleId,2:ncol(cct)],
		stringsAsFactors=FALSE
	)

	# TSS enrichment scores
	fn <- file.path(snakeDir, "output", "coverage_data", "compiled_TSS.txt")
	if (file.exists(fn)){
		cct <- readTab(fn)
		colnames(cct) <- c("sampleId", "tssEnrichment")
		rownames(cct) <- cct$sampleId
		sampleMetrics <- data.frame(
			sampleMetrics,
			tssEnrichment=cct[sampleMetrics$sampleId,"tssEnrichment"],
			stringsAsFactors=FALSE
		)
	}

	if (withPeaks){
		#get the number of called peaks
		peakDir <- file.path(snakeDir, "output", "peaks")
		numPeaks <- sapply(sampleMetrics$sampleId, FUN=function(sid){
			R.utils::countLines(file.path(peakDir, paste0(sid, "_summits.bed")))
		})
		sampleMetrics[,"numPeaks"] <- numPeaks
	}

	#add sample annotation
	sampleMetrics <- data.frame(
		sampleMetrics,
		sampleAnnot[sampleMetrics$sampleId,],
		stringsAsFactors=FALSE
	)
	return(sampleMetrics)
}

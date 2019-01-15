#-------------------------------------------------------------------------------
#' DsATAC.cellranger
#' 
#' Create a DsATAC dataset from multiple input files output by 10x cellranger
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param sampleDirPrefixCol column name specifying the directory prefix for each sample in the sample annotation table
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param diskDump     should large data objects (count matrices, fragment data, ...) be disk-backed to save main memory
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object.
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.cellranger <- function(sampleAnnot, sampleDirPrefixCol, genome, dataDir="", regionSets=NULL, sampleIdCol=sampleDirPrefixCol, diskDump=TRUE, keepInsertionInfo=TRUE){
	if (!is.element(type, c("bam", "insBam", "insBed"))){
		logger.error(c("Unsupported import type:", type))
	}

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds


	if (nchar(dataDir) > 0){
		sampleDirs <- file.path(dataDir, paste0(sampleAnnot[,sampleDirPrefixCol]))
	} else {
		sampleDirs <- sampleAnnot[,sampleDirPrefixCol]
	}
	names(sampleDirs) <- sampleIds

	if (!all(dir.exists(sampleDirs))){
		missingSamples <- sampleIds[!file.exists(sampleDirs)]
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

	# TODO: Read cells for each sample and construct joined sample annotation table

	# TODO: add sample metadata from summary.json or summary.csv

	# TODO: add cell metadata from singlecell.csv

	# TODO: optionally merge peak sets and add to region sets (peaks.bed)

	# TODO: check if this does not become too large for (tens to hundreds of) thousands of cells per sample
	#       do TENxMatrix-class objects from the 'HDF5Array' package help?
	#       do sparse matrices from the 'Matrix' package help?

	logger.start("Creating DsATAC object")
		obj <- DsATAC(sampleAnnot, genome, diskDump=diskDump)
		for (rt in names(regionSets)){
			logger.info(c("Including region set:", rt))
			obj <- regionAggregation(obj, regionSets[[rt]], rt, signal=NULL, dropEmpty=FALSE)
		}
	logger.completed()

	# TODO: adapt to read from fragment files for each sample
	logger.start("Reading ATAC data")
		nSamples <- length(sampleIds)
		logger.start(c("Adding insertion and region count data from bam"))
			for (i in seq_along(sampleDirs)){
				sid <- names(sampleDirs)[i]
				logger.start(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
					obj <- addInsertionDataFromBam(obj, sampleDirs[i], pairedEnd=pairedEnd, .diskDump=obj@diskDump)
					obj <- addCountDataFromGRL(obj, getInsertionSites(obj, samples=sid))
					# optionally remove insertion information to save space
					if (!keepInsertionInfo){
						obj@fragments <- list()
					}
				logger.completed()
			}
		logger.completed()
	logger.completed()
	
	return(obj)
}


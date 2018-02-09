#-------------------------------------------------------------------------------
#' DsATAC.bw
#' 
#' Create a DsATAC dataset from multiple input files in bigwig format
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param filePrefixCol column name specifying the file prefix for each sample in the sample annotation table
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param type         input data type. Currently only "bam" is supported
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.snakeATAC <- function(sampleAnnot, filePrefixCol, genome, dataDir, regionSets=NULL, sampleIdCol=filePrefixCol, type="bam"){
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
	#TODO:convenience: check if files exists and have the correct format before importing
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

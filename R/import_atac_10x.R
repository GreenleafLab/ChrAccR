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
DsATAC.cellranger <- function(sampleAnnot, sampleDirPrefixCol, genome, dataDir="", regionSets=NULL, addPeakRegions=TRUE, sampleIdCol=sampleDirPrefixCol, diskDump=FALSE, keepInsertionInfo=TRUE){

	# dataDir <- "~/myscratch/scATAC"
	# sampleAnnot <- data.frame(
	# 	sampleId=c("sample1", "sample2"),
	# 	sampleSubDir=c("scATAC_test", "scATAC_test"),
	# 	stringsAsFactors = FALSE
	# )
	# sampleDirPrefixCol <- "sampleSubDir"
	# sampleIdCol <- "sampleId"
	# genome <- "hg38"

	sampleIds <- as.character(sampleAnnot[,sampleIdCol])
	rownames(sampleAnnot) <- sampleIds


	if (nchar(dataDir) > 0){
		sampleDirs <- file.path(dataDir, paste0(sampleAnnot[,sampleDirPrefixCol]))
	} else {
		sampleDirs <- sampleAnnot[,sampleDirPrefixCol]
	}
	sampleDirs <- file.path(sampleDirs, "outs")
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
	cellAnnot <- do.call("rbind", lapply(1:nrow(sampleAnnot), FUN=function(i){
		sid <- sampleIds[i]
		sa <- readTab(file.path(sampleDirs[sid], "singlecell.csv"), sep=",")
		# trueCellIdx <- !(sa[, "cell_id"] %in% c("None")) & sa[, "is__cell_barcode"]==1
		trueCellIdx <- sa[, "is__cell_barcode"]==1
		sa <- sa[trueCellIdx, ]
		sa.qc <- readTab(file.path(sampleDirs[sid], "summary.csv"), sep=",")
		# alternatively, think about reading from the JSON, which contains a bit of additional info
		# sa.qc <- jsonlite::fromJSON(file.path(sampleDirs[sid], "summary.json"))
		colnames(sa.qc) <- paste0(".sampleQC.", colnames(sa.qc))
		sa.sample <- sampleAnnot[sid,]
		rownames(sa.sample) <- NULL
		sa <- data.frame(
			cellId=paste(sid, sa[,"cell_id"], sep=""),
			sa.sample,
			sa.qc,
			sa,
			stringsAsFactors = FALSE
		)
		
		return(sa)
	}))
	rownames(cellAnnot) <- cellAnnot[,"cellId"]

	# TODO: optionally merge peak sets and add to region sets (peaks.bed)
	if (addPeakRegions){

	}

	# TODO: check if this does not become too large for (tens to hundreds of) thousands of cells per sample
	#       do TENxMatrix-class objects from the 'HDF5Array' package help?
	#       do sparseMatrix objects from the 'Matrix' package help?

	logger.start("Creating DsATAC object")
		obj <- ChrAccR:::DsATAC(cellAnnot, genome, diskDump=diskDump, sparseCounts=TRUE)
		for (rt in names(regionSets)){
			logger.info(c("Including region set:", rt))
			obj <- regionAggregation(obj, regionSets[[rt]], rt, signal=NULL, dropEmpty=FALSE)
		}
	logger.completed()

	# TODO: adapt to read from fragment files for each sample
	logger.start("Reading ATAC data")
		nSamples <- length(sampleIds)
		logger.start(c("Adding insertion and region count data from fragment info"))
			for (i in seq_along(sampleDirs)){
				sid <- names(sampleDirs)[i]
				logger.start(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
					sampleIdx <- cellAnnot[,sampleIdCol]==sid
					barcode2cellId <- cellAnnot[sampleIdx,"cellId"]
					names(barcode2cellId) <- cellAnnot[sampleIdx,"barcode"]

					logger.start("Preparing fragment data")
						fragGr <- readTab(file.path(sampleDirs[i], "fragments.tsv.gz"), header=FALSE)
						colnames(fragGr) <- c("chrom", "chromStart", "chromEnd", "barcode", "duplicateCount")
						fragGr <- df2granges(fragGr, chrom.col=1L, start.col=2L, end.col=3L, strand.col=NULL, coord.format="B1RI", assembly=obj@genome, doSort=TRUE, adjNumChromNames=TRUE)
						fragGr <- fragGr[elementMetadata(fragGr)[,"barcode"] %in% names(barcode2cellId)] # only take into account fragments that can be mapped to cells
						fragGrl <- split(fragGr, elementMetadata(fragGr)[,"barcode"])
						names(fragGrl) <- barcode2cellId[names(fragGrl)]
					logger.completed()

					


					for (cid in names(fragGrl)){
						fgr <- fragGrl[[cid]]
						if (diskDump){
							fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
							saveRDS(fgr, fn)
							fgr <- fn
						}
						obj@fragments[[cid]] <- fgr
						obj <- ChrAccR:::addCountDataFromGRL(obj, getInsertionSites(obj, samples=cid))
						# TODO: check if this is fast enough. If not, try removing coercing to GRangesList in getInsertionSites() which could be a bottleneck

						# optionally remove insertion information to save space
						if (!keepInsertionInfo){
							obj@fragments <- list()
						}
					}
				logger.completed()
			}
		logger.completed()
	logger.completed()
	
	return(obj)
}


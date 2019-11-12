#-------------------------------------------------------------------------------
#' DsATACsc.fragments
#' 
#' Create a DsATACsc dataset from multiple input fragment files
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param fragmentFiles vector of fragment files or the column name in the sample annotation table containing thse file names.
#'                     fragment files must be tab-separated with columns "chrom", "chromStart", "chromEnd", "barcode" and "duplicateCount"
#'                     and must not contain a header line
#' @param genome       genome assembly
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param cellAnnot    (optional) annotation table of all cells in the dataset. Must contain a \code{'cellId'} and \code{'cellBarcode'} columns.
#' @param diskDump     should large data objects (count matrices, fragment data, ...) be disk-backed to save main memory
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object.
#' @return \code{\linkS4class{DsATACsc}} object
#' @author Fabian Mueller
#' @export
DsATACsc.fragments <- function(sampleAnnot, fragmentFiles, genome, regionSets=NULL, sampleIdCol=NULL, cellAnnot=NULL, diskDump=FALSE, keepInsertionInfo=FALSE){
	if (!is.character(fragmentFiles)) logger.error("Invalid value for fragmentFiles. Expected character")
	if (length(fragmentFiles)==1 && is.element(fragmentFiles, colnames(sampleAnnot))){
		fragmentFiles <- sampleAnnot[,fragmentFiles]
	} else {
		if (length(fragmentFiles)!=nrow(sampleAnnot)) logger.error("Invalid value for fragmentFiles. must match the number of rows in sampleAnnot")
	}
	if (!is.null(sampleIdCol) && !(is.character(sampleIdCol) && length(sampleIdCol)==1)) logger.error("Invalid value for sampleIdCol Expected character or NULL")
	if (is.null(sampleIdCol)){
		sampleIdCol <- grep("sample", colnames(sampleAnnot), value=TRUE)
		if (length(sampleIdCol) != 1) logger.error("Could not uniquely determine sample id column from sample annotation column names")
		logger.info(c("Automatically determined column '", sampleIdCol, "' as sample identifier column"))
	}
	sampleIds <- sampleAnnot[,sampleIdCol]
	nSamples <- length(sampleIds)

	if (is.null(regionSets)){
		logger.info(c("Using default region sets:", paste("tiling200bp", collapse=", ")))
		regionSets <- GRangesList(
			tiling200bp=getTilingRegions(genome, width=200L, onlyMainChrs=TRUE)
		)
	}
	
	if (is.null(cellAnnot)){
		logger.start("Preparing cell annotation")
			cellAnnot <- do.call("rbind", lapply(1:nSamples, FUN=function(i){
				sid <- sampleIds[i]
				logger.status(c("Sample:", sid, paste0("(", i, " of ", nSamples, ")")))
				fragTab <- readTab(fragmentFiles[i], header=FALSE)
				colnames(fragTab) <- c("chrom", "chromStart", "chromEnd", "barcode", "duplicateCount")
				fragCounts <- table(fragTab[,"barcode"])
				sa.sample <- sampleAnnot[sid,]
				rownames(sa.sample) <- NULL
				sa <- data.frame(
					.sampleId=sid,
					cellId=paste(sid, names(fragCounts), sep="_"),
					cellBarcode=names(fragCounts),
					sa.sample,
					nFrags=fragCounts,
					stringsAsFactors=FALSE
				)
				return(sa)
			}))
		logger.completed()
	} else {
		if (!is.element("cellId", colnames(cellAnnot))){
			logger.error("Invalid parameter: cellAnnot. Expected table with column 'cellId'")
		}
	}
	logger.start("Creating DsATAC object")
		obj <- DsATACsc(cellAnnot, genome, diskDump=diskDump, diskDump.fragments=keepInsertionInfo, sparseCounts=TRUE)
		for (rt in names(regionSets)){
			logger.info(c("Including region set:", rt))
			obj <- regionAggregation(obj, regionSets[[rt]], rt, signal=NULL, dropEmpty=FALSE)
		}
		if (obj@diskDump.fragments && .hasSlot(obj, "diskDump.fragments.nSamplesPerFile")){
			obj@diskDump.fragments.nSamplesPerFile <- 500L
		}
	logger.completed()

	logger.start("Reading scATAC data (fragment data)")
		for (i in seq_along(fragmentFiles)){
			sid <- sampleIds[i]
			logger.start(c("Importing sample", ":", sid, paste0("(", i, " of ", nSamples, ")")))
				sampleIdx <- cellAnnot[,sampleIdCol]==sid
				barcode2cellId <- cellAnnot[sampleIdx,"cellId"]
				names(barcode2cellId) <- cellAnnot[sampleIdx,"cellBarcode"]

				logger.start("Preparing fragment data")
					fragGr <- readTab(fragmentFiles[i], header=FALSE)
					colnames(fragGr) <- c("chrom", "chromStart", "chromEnd", "barcode", "duplicateCount")
					fragGr <- df2granges(fragGr, chrom.col=1L, start.col=2L, end.col=3L, strand.col=NULL, coord.format="B1RI", assembly=obj@genome, doSort=TRUE, adjNumChromNames=TRUE)
					fragGr <- fragGr[elementMetadata(fragGr)[,"barcode"] %in% names(barcode2cellId)] # only take into account fragments that can be mapped to cells
					if (length(fragGr) < 2) logger.error("Too few fragments corresponding to actual cells")
					fragGrl <- split(fragGr, elementMetadata(fragGr)[,"barcode"])
					names(fragGrl) <- barcode2cellId[names(fragGrl)]
				logger.completed()

				logger.start("Preparing insertion data")
					insGrl <- lapply(fragGrl, getInsertionSitesFromFragmentGr)
				logger.completed()
				logger.start("Summarizing count data")
					obj <- addCountDataFromGRL(obj, insGrl)
				logger.completed()
				
				if (keepInsertionInfo) {
					chunkedFragmentFiles <- obj@diskDump.fragments && .hasSlot(obj, "diskDump.fragments.nSamplesPerFile") && obj@diskDump.fragments.nSamplesPerFile > 1
					logger.start("Adding fragment data to data structure")
						if (chunkedFragmentFiles){
							nCells <- length(fragGrl)
							chunkL <- split(1:nCells, rep(1:ceiling(nCells/obj@diskDump.fragments.nSamplesPerFile), each=obj@diskDump.fragments.nSamplesPerFile)[1:nCells])
							names(chunkL) <- NULL
							for (k in 1:length(chunkL)){
								iis <- chunkL[[k]]
								cids <- names(fragGrl)[iis]
								fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
								saveRDS(fragGrl[iis], fn)
								obj@fragments[cids] <- rep(list(fn), length(iis))
							}
						} else {
							for (cid in names(fragGrl)){
								fgr <- fragGrl[[cid]]
								if (obj@diskDump.fragments){
									fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
									saveRDS(fgr, fn)
									fgr <- fn
								}
								obj@fragments[[cid]] <- fgr
							}
						}
					logger.completed()
				}
				
			logger.completed()
		}
		# just to make sure: reorder fragment list
		if (keepInsertionInfo) {
			obj@fragments <- obj@fragments[getSamples(obj)]
		}
	logger.completed()
	return(obj)
}

#-------------------------------------------------------------------------------
#' DsATAC.cellranger
#' 
#' Create a DsATAC dataset from multiple input files output by 10x cellranger
#' @param sampleAnnot  data.frame specifying the sample annotation table
#' @param sampleDirPrefixCol column name specifying the directory prefix for each sample in the sample annotation table
#' @param genome       genome assembly
#' @param dataDir      directory where the files are located
#' @param regionSets   a list of GRanges objects which contain region sets over which count data will be aggregated
#' @param addPeakRegions should a merged set of peaks be created as one of the region sets (merged, non-overlapping peaks of width=500bp from the peaks of individual samples)
#' @param sampleIdCol  column name or index in the sample annotation table containing unique sample identifiers
#' @param diskDump     should large data objects (count matrices, fragment data, ...) be disk-backed to save main memory
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object.
#' @return \code{\linkS4class{DsATAC}} object
#' @author Fabian Mueller
#' @export
DsATAC.cellranger <- function(sampleAnnot, sampleDirPrefixCol, genome, dataDir="", regionSets=NULL, addPeakRegions=TRUE, sampleIdCol=sampleDirPrefixCol, diskDump=FALSE, keepInsertionInfo=FALSE){
	# sampleAnnot <- data.frame(
	# 	sampleId=c("sample1", "sample2"),
	# 	sampleSubDir=c("scATAC_test", "scATAC_test"),
	# 	stringsAsFactors = FALSE
	# )
	# ds <- DsATAC.cellranger(sampleAnnot, "sampleSubDir", "hg38", dataDir="~/myscratch/scATAC", regionSets=NULL, addPeakRegions=TRUE, sampleIdCol="sampleId", diskDump=FALSE, keepInsertionInfo=TRUE)

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

	# for reading sample metadata from JSON
	json2df <- function(fn){
		ll <- jsonlite::fromJSON(fn)
		nr <- length(ll[[1]])
		ll <- lapply(ll, FUN=function(x){
			if (is.null(x)) return(NA)
			if (length(x)!=nr) logger.error("Could not read data.frame from JSON: not all elements have the same length")
			return(x)
		})
		df <- data.frame(ll, stringsAsFactors = FALSE)
		return(df)
	}

	# TODO: Read cells for each sample and construct joined sample annotation table
	cellAnnotL <- lapply(1:nrow(sampleAnnot), FUN=function(i){
		sid <- sampleIds[i]
		sa <- readTab(file.path(sampleDirs[sid], "singlecell.csv"), sep=",")
		# trueCellIdx <- !(sa[, "cell_id"] %in% c("None")) & sa[, "is__cell_barcode"]==1
		trueCellIdx <- sa[, "is__cell_barcode"]==1
		sa <- sa[trueCellIdx, ]
		colnames(sa) <- paste0(".CR.cellQC.", colnames(sa))
		# sa.qc <- readTab(file.path(sampleDirs[sid], "summary.csv"), sep=",") #reading from csv, which contains only a subset of information
		sa.qc <- json2df(file.path(sampleDirs[sid], "summary.json"))
		colnames(sa.qc) <- paste0(".CR.sampleQC.", colnames(sa.qc))
		sa.sample <- sampleAnnot[sid,]
		rownames(sa.sample) <- NULL
		sa <- data.frame(
			.sampleId=sid,
			cellId=paste(sid, sa[,".CR.cellQC.cell_id"], sep=""),
			sa.sample,
			sa.qc,
			sa,
			stringsAsFactors = FALSE
		)
		
		return(sa)
	})
	colnames.union <- c()
	colnames.intersect <- colnames(cellAnnotL[[1]])
	for (i in seq_along(cellAnnotL)){
		colnames.union <- union(colnames.union, colnames(cellAnnotL[[i]]))
		colnames.intersect <- intersect(colnames.intersect, colnames(cellAnnotL[[i]]))
	}
	if (!all(colnames.union %in% colnames.intersect)){
		logger.warning(c("The following columns could not be found in all CellRanger summary files and will be discarded from the annotation:", paste(setdiff(colnames.union, colnames.intersect),collapse=", ")))
	}

	cellAnnot <- do.call("rbind", lapply(cellAnnotL, FUN=function(x){x[,colnames.intersect]}))
	rownames(cellAnnot) <- cellAnnot[,"cellId"]
	cellAnnot[,"cellBarcode"] <- cellAnnot[,".CR.cellQC.barcode"]

	# optionally merge peak sets and add to region sets (peaks.bed)
	if (addPeakRegions){
		unifWidth=501L
		logger.start("Retrieving (consensus) non-ovelapping, fixed-width peak set")
			peakSet <- NULL
			for (i in seq_along(sampleDirs)){
				sid <- names(sampleDirs)[i]
				pGr <- rtracklayer::import(file.path(sampleDirs[i], "peaks.bed"), format="BED")
				pGr <- setGenomeProps(pGr, genome, onlyMainChrs=TRUE)
				pGr <- trim(resize(pGr, width=unifWidth, fix="center", ignore.strand=TRUE))
				pGr <- pGr[width(pGr)==median(width(pGr))] #remove too short regions which might have been trimmed
				elementMetadata(pGr)[,"dummyScore"] <- 1L # dummy score coloumn. All identical to 1. Results in the first peak of an overlap being selected

				if (is.null(peakSet)){
					#initialize peak set with all peaks from the first sample
					#remove overlapping peaks in initial sample based on the dummy score
					peakSet.cur <- getNonOverlappingByScore(pGr, scoreCol="dummyScore")
					# elementMetadata(peakSet.cur)[,paste0(".ov.", sid)] <- TRUE
					peakSet <- peakSet.cur
				} else {
					# add new peaks and remove the overlapping ones by taking the peaks with the best score
					peakSet <- getNonOverlappingByScore(c(peakSet, pGr), scoreCol="dummyScore")
					# elementMetadata(peakSet)[,paste0(".ov.", sid)] <- overlapsAny(peakSet, peakSet.cur, ignore.strand=TRUE)
				}
			}
			peakSet <- sortSeqlevels(peakSet)
			peakSet <- sort(peakSet)
			regionSets <- c(regionSets, list(.peaks=peakSet))
		logger.completed()
	}

	fragmentFiles <- file.path(sampleDirs, "fragments.tsv.gz")
	obj <- DsATACsc.fragments(sampleAnnot, fragmentFiles, genome=genome, regionSets=regionSets, sampleIdCol=sampleIdCol, cellAnnot=cellAnnot, diskDump=diskDump, keepInsertionInfo=keepInsertionInfo)
	
	return(obj)
}


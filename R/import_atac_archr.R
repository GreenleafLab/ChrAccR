#-------------------------------------------------------------------------------
#' DsATACsc.archr
#' 
#' Create a DsATACsc dataset from an \code{ArchR} project
#' @param ap	\code{ArchR} project object
#' @param keepInsertionInfo flag indicating whether to maintain the insertion information in the resulting object.
#' @param diskDump.fragments Keep fragment coordinates stored on disk rather than in main memory. This saves memory, but increases runtime and I/O.
#' @return \code{\linkS4class{DsATACsc}} object
#' @author Fabian Mueller
#' @export
DsATACsc.archr <- function(ap, keepInsertionInfo=FALSE, diskDump.fragments=keepInsertionInfo){
	require(ArchR)
	if (class(ap)!="ArchRProject") logger.error("Invalid input: expected ArchRProject")
	afs <- ArchR::getArrowFiles(ap)
	if (!all(file.exists(afs))) logger.error("Could not find all arrow files")

	gg <- ArchR::getGenome(ap)
	ga <- ArchR::getGenomeAnnotation(ap)
	genomeAss <- ""
	if (grepl("Hsapiens.*hg38", gg)){
		genomeAss <- "hg38"
	} else {
		logger.error(c("unsupported genome:", gg))
	}

	# check genome compatibility
	logger.status("Checking genome compatibility")
	sls <- muRtools::getSeqlengths4assembly(genomeAss, onlyMainChrs=TRUE, adjChrNames=TRUE)
	chromSizes_archr <- ga$chromSizes
	chromNames_archr <- as.character(seqnames(chromSizes_archr))

	chromSizes_archr <- muRtools::setGenomeProps(chromSizes_archr, genomeAss, dropUnknownChrs=TRUE, onlyMainChrs=TRUE, adjChrNames=TRUE, silent=FALSE)
	chromNames <- as.character(seqnames(chromSizes_archr))
	if (!all(chromNames_archr %in% chromNames)){
		logger.error(c("Could not find the following chromosomes in the annotation:", paste(chromNames_archr[!(chromNames_archr %in% chromNames)], collapse=",")))
	}

	sls_archr <- end(chromSizes_archr)
	names(sls_archr) <- chromNames

	if (any(sls_archr != sls[chromNames])){
		logger.error("Incompatible chromosome lengths detected")
	}

	logger.start("Preparing cell annotation")
		cellAnnot <- as.data.frame(ArchR::getCellColData(ap))
		# rename the annotation columns to "archr."
		colnames(cellAnnot) <- paste0("archr.", colnames(cellAnnot)) 

		cellIds_archr <- rownames(cellAnnot)
		cellIds <- gsub("#", "_", cellIds_archr)
		rownames(cellAnnot) <- cellIds
		cellBcs <- gsub("(.+)#(.+)$", "\\2", cellIds_archr)

		# map special column names to annotation
		cellAnnot[,"cellId"] <- cellIds
		cellAnnot[,"cellBarcode"] <- cellBcs
		cellAnnot[,"archr.cellId"] <- cellIds_archr
		cn_map <- c(
			"archr.TSSEnrichment"=".tssEnrichment",
			"archr.Sample"=".sampleId",
			"archr.nFrags"="nFrags"
		)
		idx <- colnames(cellAnnot) %in% names(cn_map)
		colnames(cellAnnot)[idx] <- cn_map[colnames(cellAnnot)[idx]]

		sampleIds <- ArchR::getSampleNames(ap)
	logger.completed()

	blacklistGr <- ga$blacklist
	getTileGr <- function(tileSize){
		# tGr <- muRtools::getTilingRegions(genomeAss, width=tileSize, onlyMainChrs=TRUE, adjChrNames=TRUE)
		# the above line is of by one basepair. 
		# make sure that the intervals start at *0 and end at *9 (except the first range of each chromosome (starts at 1) and the last range (ends at length))
		tGr <- do.call("c", lapply(chromNames, FUN=function(chrom){
			ss <- seq(0, sls[chrom], by=tileSize) # [-1] : exclude the 0
			ss[1] <- 1
			ee <- c(ss[2:length(ss)] - 1, sls[chrom])
			names(ee) <- NULL
			gr <- GRanges(seqnames=chrom, ranges=IRanges(start=ss, end=ee))
			gr <- muRtools::setGenomeProps(gr, genomeAss, dropUnknownChrs=TRUE, onlyMainChrs=TRUE, adjChrNames=TRUE)
			return(gr)
		}))
		# # sanity check
		# nTiles <- trunc(sls[chromNames] / tileSize) + 1
		# nTiles2 <- elementNROWS(split(tGr, seqnames(tGr)))[chromNames]
		# if (any(nTiles != nTiles2)) logger.error("Incompatible tile numbers")
		return(tGr)
	}
	logger.start("Preparing region sets")
		regionSets <- list()
		tileGr <- NULL
		tileDf <- ArchR::.getFeatureDF(afs, "TileMatrix")
		if (!is.null(tileDf) && class(tileDf)=="DataFrame"){
			logger.status("Preparing tiling regions")
			x <- tileDf[,"start"]
			tileW <- median(x[2:length(x)]-x[1:(length(x)-1)])
			logger.info(c("Determined tiling width as", tileW, "bp"))
			tileGr <- getTileGr(tileW)
			regionSets[[paste0("tiling", tileW, "bp")]] <- tileGr
		}
		geneAnnot <- ArchR::getGeneAnnotation(ap)
		if (!is.null(geneAnnot) && all(c("TSS", "genes") %in% names(geneAnnot))){
			logger.status("Preparing gene regions")
			regionSets[["tssWindow"]] <- promoters(geneAnnot$TSS, upstream=100, downstream=100)
			regionSets[["promoter"]] <- promoters(geneAnnot$genes, upstream=1500, downstream=500)
		}
		peakGr <- ArchR::getPeakSet(ap)
		if (!is.null(peakGr) && class(peakGr) == "GRanges"){
			logger.status("Preparing peak regions")
			regionSets[["archr.peaks"]] <- peakGr
		}
		redDim <- ArchR::getReducedDims(ap, reducedDims="IterativeLSI", returnMatrix=FALSE)
		if (is.element("LSIFeatures", names(redDim))){
			logger.status("Preparing iterativeLSI feature regions")
			featDf <- redDim$LSIFeatures
			if (class(featDf)=="DataFrame" && !is.null(tileGr)){
				logger.info("Assuming iterativeLSI has been applied to tiling windows")
				tileIdxL <- tapply(featDf[,"idx"], featDf[,"seqnames"], c)
				tileGrl <- split(tileGr, seqnames(tileGr))
				cns <- intersect(names(tileGrl), names(tileIdxL))
				gr <- do.call("c", lapply(cns, FUN=function(chrom){
					tileGrl[[chrom]][tileIdxL[[chrom]]]
				}))
				regionSets[["archr.itlsi.feats"]] <- gr
			}
		}
		for (rt in names(regionSets)){
			regionSets[[rt]] <- muRtools::setGenomeProps(regionSets[[rt]], genomeAss, dropUnknownChrs=TRUE, onlyMainChrs=TRUE, adjChrNames=TRUE, silent=TRUE)
		}
	logger.completed()


	logger.start("Creating DsATAC object")
		obj <- DsATACsc(cellAnnot, genomeAss, diskDump=FALSE, diskDump.fragments=diskDump.fragments, sparseCounts=TRUE)
		for (rt in names(regionSets)){
			logger.info(c("Including region set:", rt))
			obj <- regionAggregation(obj, regionSets[[rt]], rt, signal=NULL, dropEmpty=FALSE)
		}
		if (obj@diskDump.fragments && .hasSlot(obj, "diskDump.fragments.nSamplesPerFile")){
			obj@diskDump.fragments.nSamplesPerFile <- 500L
		}
	logger.completed()

	logger.start("Preparing fragments")
		fragGrl <- lapply(sampleIds, FUN=function(sid){
			logger.status(c("Reading fragments for sample:", sid))
			gr <- ArchR::getFragmentsFromArrow(afs[sid], cellNames=ArchR::getCellNames(ap), verbose=FALSE)
			gr <- muRtools::setGenomeProps(gr, genomeAss, dropUnknownChrs=TRUE, onlyMainChrs=TRUE, adjChrNames=TRUE)
			cids <- elementMetadata(gr)[,"RG"]
			cids <- gsub("#", "_", cids)
			elementMetadata(gr) <- NULL
			elementMetadata(gr)[,"sampleId"] <- sid
			elementMetadata(gr)[,"cellId"] <- cids
			return(split(gr, elementMetadata(gr)[,"cellId"]))
		})
		logger.status("Combining samples ...")
		fragGrl <- do.call("c", fragGrl)
		if (!all(cellIds %in% names(fragGrl))) logger.error("Could not find all cells in fragment files")
		fragGrl <- fragGrl[cellIds]

		logger.start("Preparing insertion data")
			insGrl <- lapply(fragGrl, getInsertionSitesFromFragmentGr)
		logger.completed()
		logger.start("Summarizing count data")
			obj <- addCountDataFromGRL(obj, insGrl)
		logger.completed()
		
		if (keepInsertionInfo) {
			chunkedFragmentFiles <- obj@diskDump.fragments && .hasSlot(obj, "diskDump.fragments.nSamplesPerFile") && obj@diskDump.fragments.nSamplesPerFile > 1
			logger.start("Adding fragment data to data structure")
				nCells <- length(fragGrl)
				if (chunkedFragmentFiles){
					chunkL <- split(1:nCells, rep(1:ceiling(nCells/obj@diskDump.fragments.nSamplesPerFile), each=obj@diskDump.fragments.nSamplesPerFile)[1:nCells])
					names(chunkL) <- NULL
					for (k in 1:length(chunkL)){
						iis <- chunkL[[k]]
						cids <- names(fragGrl)[iis]
						fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
						saveRDS(fragGrl[iis], fn, compress=TRUE)
						obj@fragments[cids] <- rep(list(fn), length(iis))
					}
				} else {
					cids <- names(fragGrl)
					if (obj@diskDump.fragments){
						obj@fragments[cids] <- lapply(fragGrl, FUN=function(x){
							fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
							saveRDS(x, fn, compress=TRUE)
							return(fn)
						})
					} else {
						obj@fragments[cids] <- fragGrl
					}
				}
			logger.completed()
		}
	logger.completed()

	return(obj)
}


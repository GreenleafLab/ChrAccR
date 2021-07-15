#' DsATAC
#'
#' A class for storing ATAC-seq accessibility data
#' inherits from \code{\linkS4class{DsAcc}}
#' 
#' @section Slots:
#' \describe{
#'   \item{\code{fragments}}{
#'		\code{GRanges} object storing sequencing fragments. Alternativily pointers to files in which this data isf stored
#'      as R data object
#'   }
#'   \item{\code{counts}}{
#'		List of count matrices for each summarized region type (dimension: regions X samples).
#'      Depending on the settings for the slots \code{diskDump} and \code{sparseCounts}, the matrices are
#'      either (a) regular matrices, (b) \code{HDF5Array}/\code{DelayedArray} or (c) sparse matrices.
#'   }
#'   \item{\code{countTransform}}{
#'		list of character vectors specifying which transformations have been applied to the count matrices
#'   }
#'   \item{\code{sparseCounts}}{
#'		Flag indicating whether count data will be stored as sparse matrices rather than regular matrices
#'   }
#'   \item{\code{diskDump.fragments}}{
#'		Flag indicating whether fragment data will be kept on disk rather than in main memory.
#'   }
#' }
#' @name DsATAC-class
#' @rdname DsATAC-class
#' @author Fabian Mueller
#' @exportClass DsATAC
setClass("DsATAC",
	slots = list(
		fragments = "list",
		counts = "list",
		countTransform = "list",
		sparseCounts = "logical",
		diskDump.fragments = "logical",
		diskDump.fragments.nSamplesPerFile = "integer"
	),
	contains = "DsAcc",
	package = "ChrAccR"
)
setMethod("initialize","DsATAC",
	function(
		.Object,
		fragments,
		coord,
		counts,
		sampleAnnot,
		genome,
		diskDump,
		diskDump.fragments,
		diskDump.fragments.nSamplesPerFile,
		sparseCounts
	) {
		.Object@fragments  <- fragments
		.Object@coord       <- coord
		.Object@counts      <- counts
		.Object@countTransform <- rep(list(character(0)), length(.Object@counts))
		names(.Object@countTransform) <- names(.Object@counts)
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object@diskDump    <- diskDump
		.Object@diskDump.fragments <- diskDump.fragments
		.Object@diskDump.fragments.nSamplesPerFile <- diskDump.fragments.nSamplesPerFile
		.Object@sparseCounts <- sparseCounts
		.Object@pkgVersion  <- packageVersion("ChrAccR")
		.Object
	}
)

#' @param sampleAnnot \code{data.frame} object containing sample annotation
#' @param genome    character string containing genome assembly
#' @noRd
DsATAC <- function(sampleAnnot, genome, diskDump=FALSE, diskDump.fragments=TRUE, sparseCounts=FALSE){
	obj <- new("DsATAC",
		list(),
		list(),
		list(),
		sampleAnnot,
		genome,
		diskDump,
		diskDump.fragments,
		diskDump.fragments.nSamplesPerFile=1L,
		sparseCounts
	)
	return(obj)
}

################################################################################
# Getters
################################################################################
if (!isGeneric("getCounts")) {
	setGeneric(
		"getCounts",
		function(.object, ...) standardGeneric("getCounts"),
		signature=c(".object")
	)
}
#' getCounts-methods
#'
#' Return table of count values
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param type    character string specifying the region type
#' @param i       (optional) row (region) indices
#' @param j       (optional) column (sample) indices
#' @param asMatrix return a matrix object instead of the internal representation
#' @param naIsZero should \code{NA}s in the count matrix be considered 0 value (instead of unknown/missing)
#' @param allowSparseMatrix if \code{asMatrix}: allow for sparse matrices as returned data format
#' @return Matrix containing counts for
#'         each region and sample
#'
#' @rdname getCounts-DsATAC-method
#' @docType methods
#' @aliases getCounts
#' @aliases getCounts,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getCounts",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		i=NULL,
		j=NULL,
		asMatrix=TRUE,
		naIsZero=TRUE,
		allowSparseMatrix=FALSE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		res <- .object@counts[[type]]

		# # OLD
		# if (.object@diskDump){
		# 	if (!is.null(i) || !is.null(j) || asMatrix){
		# 		# DelayedArray can have serious performance issues, when indexing is not done efficiently
		# 		# --> workaround-function: fastDelayedArrayToMatrix()
		# 		res <- fastDelayedArrayToMatrix(res, i=i, j=j) # not needed any more (?): issue should be fixed starting with HDF5Array version 1.11.11 (https://github.com/Bioconductor/DelayedArray/issues/13)
		# 		cns <- colnames(res)
		# 		if (!asMatrix){
		# 			res <- as(res, "HDF5Array")
		# 			colnames(res) <- cns
		# 		}
		# 	}
		# } else {
		# 	if (!is.null(i)) res <- res[i,,drop=FALSE]
		# 	if (!is.null(j)) res <- res[,j,drop=FALSE]
		# 	if (asMatrix && !is.matrix(res) && !(.object@sparseCounts && allowSparseMatrix)){
		# 		res <- as.matrix(res)
		# 		if (!naIsZero && .object@sparseCounts){
		# 			res[res==0] <- NA
		# 		}
		# 	}
		# }
		# # ###

		if (!is.null(i)) res <- res[i,,drop=FALSE]
		if (!is.null(j)) res <- res[,j,drop=FALSE]
		if (asMatrix && !is.matrix(res)){
			if (.object@diskDump) {
				res <- as.matrix(res)
			} else if (.object@sparseCounts && !allowSparseMatrix){
				res <- as.matrix(res)
				if (!naIsZero){
					res[res==0] <- NA
				}
			}
		}

		if (naIsZero){
			res[is.na(res)] <- 0
		}
		return(res)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("getCountsSE")) {
	setGeneric(
		"getCountsSE",
		function(.object, ...) standardGeneric("getCountsSE"),
		signature=c(".object")
	)
}
#' getCountsSE-methods
#'
#' Return a \code{SummarizedExperiment} object of count values
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param type    character string specifying the region type
#' @param naIsZero should \code{NA}s in the count matrix be considered 0 value (instead of unknown/missing)
#' @return \code{SummarizedExperiment} containing counts for each region and sample
#'
#' @rdname getCountsSE-DsATAC-method
#' @docType methods
#' @aliases getCountsSE
#' @aliases getCountsSE,DsATAC-method
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' se <- getCountsSE(dsa, "IA_prog_peaks")
#' se
#' }
setMethod("getCountsSE",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		naIsZero=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		#count matrix
		cm <- ChrAccR::getCounts(.object, type, asMatrix=TRUE, naIsZero=naIsZero, allowSparseMatrix=TRUE)
		coords <- getCoord(.object, type)
		se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=cm), rowRanges=coords, colData=S4Vectors::DataFrame(getSampleAnnot(.object)))
		return(se)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getFragmentGr")) {
	setGeneric(
		"getFragmentGr",
		function(.object, ...) standardGeneric("getFragmentGr"),
		signature=c(".object")
	)
}
#' getFragmentGr-methods
#'
#' Return a \code{GRanges} object of fragment data for a given sample
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param sampleId sample identifier
#' @return \code{GRanges} object containing fragment data
#'
#' @rdname getFragmentGr-DsATAC-method
#' @docType methods
#' @aliases getFragmentGr
#' @aliases getFragmentGr,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getFragmentGr",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		sampleId
	) {
		if (!is.element(sampleId, getSamples(.object))) logger.error(c("Invalid samples:", sampleId))
		if (!is.element(sampleId, names(.object@fragments))) logger.error(c("Object does not contain insertion information for sample:", sampleId))
		fragGr <- .object@fragments[[sampleId]]
		if (is.character(fragGr)){
			if (file.exists(fragGr)) {
				fragGr <- readRDS(fragGr)
				if (is.list(fragGr) || is.element(class(fragGr), c("GRangesList", "CompressedGRangesList"))){
					fragGr <- fragGr[[sampleId]]
				}
			} else {
				logger.error(c("Could not load fragment data from file:", fragGr))
			}
		}
		return(fragGr)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getFragmentGrl")) {
	setGeneric(
		"getFragmentGrl",
		function(.object, ...) standardGeneric("getFragmentGrl"),
		signature=c(".object")
	)
}
#' getFragmentGrl-methods
#'
#' Return a list of \code{GRanges} objects of fragment data for a given set of samples
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param sampleIds sample identifiers
#' @param asGRangesList should a \code{GRangesList} object be returned instead of a regular \code{list}
#' @return A named list of \code{GRanges} objects containing fragment data
#'
#' @rdname getFragmentGrl-DsATAC-method
#' @docType methods
#' @aliases getFragmentGrl
#' @aliases getFragmentGrl,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getFragmentGrl",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		sampleIds,
		asGRangesList=FALSE
	) {
		if (!all(sampleIds %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(sampleIds, getSamples(.object)), collapse=";")))
		if (!all(sampleIds %in% names(.object@fragments))) logger.error(c("Object does not contain insertion information for samples:", paste(setdiff(sampleIds, names(.object@fragments)), collapse=";")))

		fragL <- .object@fragments[sampleIds]
		# load disk-dumped fragment GRanges into memory
		isDd <- sapply(fragL, is.character)
		if (any(isDd)){
			ddFns_all <- unlist(fragL[isDd])
			sidsDd <- names(fragL)[isDd] # sample names of disk dumped fragments
			fnToSampleL <- tapply(sidsDd, ddFns_all, c)
			ddFns <- unique(ddFns_all)
			# large GRanges list of all samples that have been disk-dumped
			ddGrl <- lapply(ddFns, FUN=function(fn){
				# logger.status(fn)
				if (!file.exists(fn)) logger.error(c("Could not load fragment data from file:", fn))
				rr <- readRDS(fn)
				if (is.element(class(rr), c("GRangesList", "CompressedGRangesList"))){
					rr <- as.list(rr)
				}
				if (!is.list(rr)) {
					rr <- list(rr)
					names(rr) <- fnToSampleL[[fn]]
				}
				return(rr)
			})
			names(ddGrl) <- NULL
			ddGrl <- do.call("c", ddGrl)
			
			fragL[sidsDd] <- ddGrl[sidsDd]
		}
		if (asGRangesList){
			fragL <- GenomicRanges::GRangesList(fragL) # can potentially take long and throw error: 'long vectors not supported yet: memory.c:3715'
		}
		return(fragL)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("getFragmentNum")) {
	setGeneric(
		"getFragmentNum",
		function(.object, ...) standardGeneric("getFragmentNum"),
		signature=c(".object")
	)
}
#' getFragmentNum-methods
#'
#' Return the number of fragments in the \code{\linkS4class{DsATAC}} object
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param sampleIds sample identifiers
#' @return a vector of fragment counts per sample
#'
#' @rdname getFragmentNum-DsATAC-method
#' @docType methods
#' @aliases getFragmentNum
#' @aliases getFragmentNum,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getFragmentNum",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		sampleIds=getSamples(.object)
	) {
		if (!all(sampleIds %in% getSamples(.object))) logger.error(c("Invalid sampleIds:", paste(setdiff(sampleIds, getSamples(.object)), collapse=", ")))
		fgrl <- getFragmentGrl(.object, getSamples(.object), asGRangesList=FALSE)
		res <- sapply(fgrl, length)
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getInsertionSites")) {
	setGeneric(
		"getInsertionSites",
		function(.object, ...) standardGeneric("getInsertionSites"),
		signature=c(".object")
	)
}
#' getInsertionSites-methods
#'
#' Return a list of insertion sites (Tn5 cut sites) for each sample
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param samples sample identifiers
#' @param asGRangesList should a \code{GRangesList} object be returned instead of a regular \code{list}
#' @return \code{list} or \code{GRangesList} containing Tn5 cut sites for each sample
#'
#' @rdname getInsertionSites-DsATAC-method
#' @docType methods
#' @aliases getInsertionSites
#' @aliases getInsertionSites,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getInsertionSites",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		samples=getSamples(.object),
		asGRangesList=FALSE
	) {
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		if (!all(samples %in% names(.object@fragments))) logger.error(c("Object does not contain insertion information for samples:", paste(setdiff(samples, names(.object@fragments)), collapse=", ")))
		fGrl <- getFragmentGrl(.object, samples, asGRangesList=FALSE)
		res <- lapply(fGrl, getInsertionSitesFromFragmentGr)
		names(res) <- samples
		if (asGRangesList){
			res <- GRangesList(res)
		}
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getCoverage")) {
	setGeneric(
		"getCoverage",
		function(.object, ...) standardGeneric("getCoverage"),
		signature=c(".object")
	)
}
#' getCoverage-methods
#'
#' Return a list of genome-wide coverage from insertion sites
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param samples sample identifiers
#' @return \code{list} of \code{Rle} objects of sample coverage tracks
#'
#' @rdname getCoverage-DsATAC-method
#' @docType methods
#' @aliases getCoverage
#' @aliases getCoverage,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getCoverage",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		samples=getSamples(.object)
	) {
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		insGrl <- getInsertionSites(.object, samples)
		sampleCovgRle <- lapply(samples, FUN=function(sid){
			# logger.status(c("Computing genome-wide coverage for sample", sid))
			return(GenomicRanges::coverage(insGrl[[sid]]))
		})
		names(sampleCovgRle) <- samples
		return(sampleCovgRle)
	}
)

################################################################################
# Display
################################################################################
setMethod("show","DsATAC",
	function(object) {
		ss <- getSamples(object)
		str.ss <- paste(getSamples(object), collapse=", ")
		if (length(ss) > 5) str.ss <- paste(c(getSamples(object)[1:5], "..."), collapse=", ")
		rts <- getRegionTypes(object)
		str.rts <- "no region types"
		if (length(rts) > 0) str.rts <- paste0(length(rts), " region types: ", paste(rts, collapse=", "))
		str.frags <- "no fragment data"
		if (length(object@fragments) > 0) str.frags <- paste0("fragment data for ", length(object@fragments), " samples")
		if (object@diskDump.fragments) str.frags <- paste0(str.frags, " [disk-backed]")
		str.disk <- "[in memory object]"
		if (object@diskDump) str.disk <- "[contains disk-backed data]"

		if (class(object)=="DsATACsc"){
			cat("Single-cell ")
		}
		cat("DsATAC chromatin accessibility dataset \n")
		cat(str.disk, "\n")
		cat("contains:\n")
		cat(" * ", length(ss), " samples: ", str.ss, " \n")
		cat(" * ", str.frags, " \n")
		cat(" * ", str.rts, " \n")
		if (length(rts) > 0) {
			for (rt in rts){
				cat(" *  * ", getNRegions(object, rt), "regions of type", rt, " \n")
			}
		}
	}
)

################################################################################
# Summary functions
################################################################################

# #-------------------------------------------------------------------------------
# if (!isGeneric("getRegionMapping")) {
# 	setGeneric(
# 		"getRegionMapping",
# 		function(.object, ...) standardGeneric("getRegionMapping"),
# 		signature=c(".object")
# 	)
# }
# #' getRegionMapping-methods
# #'
# #' Retrieve a mapping from regions to the indices of counts in the dataset
# #'
# #' @param .object \code{\linkS4class{DsATAC}} object
# #' @param type    character string specifying a name for the region type
# #' @return list containing vectors of indices of GCs for each region of the
# #'         specified type
# #' 
# #' @rdname getRegionMapping-DsATAC-method
# #' @docType methods
# #' @aliases getRegionMapping
# #' @aliases getRegionMapping,DsATAC-method
# #' @author Fabian Mueller
# #' @export
# setMethod("getRegionMapping",
# 	signature(
# 		.object="DsATAC"
# 	),
# 	function(
# 		.object,
# 		type
# 	) {
# 		baseGr <- getCoord(.object, type="counts")
# 		regGr  <- getCoord(.object, type=type)
# 		oo <- findOverlaps(baseGr, regGr)
# 		sl <- tapply(GenomicRanges::queryHits(oo), GenomicRanges::subjectHits(oo), c)

# 		res <- rep(list(integer(0)), getNRegions(.object, type=type))
# 		res[as.integer(names(sl))] <- sl
# 		return(res)
# 	}
# )
#-------------------------------------------------------------------------------
if (!isGeneric("regionAggregation")) {
	setGeneric(
		"regionAggregation",
		function(.object, ...) standardGeneric("regionAggregation"),
		signature=c(".object")
	)
}
#' regionAggregation-methods
#'
#' Aggregate signal counts across a set of regions
#'
#' @param .object   \code{\linkS4class{DsATAC}} object
#' @param regGr     \code{GRanges} object containing regions to summarize
#' @param type      character string specifying a name for the region type
#' @param signal    character string specifying a name for the region type for the signal to be aggregated
#'                  If it is \code{NULL} (default), the new region type will be initialized with NA values.
#'                  If it is \code{"insertions"} count data will be initialized from insertion sites (if 
#'                  fragment data is present in the object).
#' @param aggrFun   aggregation function for signal counts. Will only be used if \code{signal!="insertions"}
#'                  Currently \code{sum}, \code{mean} and \code{median} (default) are supported.
#' @param dropEmpty discard all regions with no observed signal counts
#' @param bySample  [only relevant if \code{signal=="insertions"}]. Process data sample-by-sample to save memory.
#' @param chunkSize [only relevant if \code{signal=="insertions" & !bySample}] number of samples to process per chunk (saves memory).
#' 				    If \code{NULL} or larger than the number of samples, only one chunk will be processed.
#' @return a new \code{\linkS4class{DsATAC}} object with aggregated signal counts per regions
#'
#' 
#' @rdname regionAggregation-DsATAC-method
#' @docType methods
## @aliases regionAggregation
#' @aliases regionAggregation,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("regionAggregation",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regGr,
		type,
		signal=NULL,
		aggrFun="median",
		dropEmpty=TRUE,
		bySample=TRUE,
		chunkSize=5000L
	) {
		if (!is.element(aggrFun, c("sum", "mean", "median"))){
			logger.error(c("Unknown signal count aggregation function:", aggrFun))
		}
		if (!is.null(signal)){
			if (!is.element(signal, c(getRegionTypes(.object), "insertions"))){
				logger.error(c("invalid signal region type", signal))
			}
			if (type==signal){
				logger.error(paste0("Region type ('", type, "') is not allowed to be the same as signal ('", signal, "')"))
			}
		}
		if (is.element(type, getRegionTypes(.object))){
			logger.warning(c("Overwriting aggregated region type:", type))
		}

		rsFun <- rowSums
		# sparse matrices
		if (.object@sparseCounts){
			rsFun <- Matrix::rowSums
		}
		# DelayedArray
		if (.object@diskDump){
			rsFun <- BiocGenerics::rowSums
		}

		#sort the regions
		coordOrder <- order(as.integer(seqnames(regGr)), start(regGr), end(regGr), as.integer(strand(regGr)))
		regGr <- regGr[coordOrder]

		.object@coord[[type]] <- regGr	

		#initialize empty count matrix for new region type
		sampleIds <- getSamples(.object)
		nSamples <- length(sampleIds)
		nRegs    <- length(regGr)
		# emptyVec <- rep(as.integer(NA), nRegs)
		# .object@counts[[type]] <- data.table(emptyVec)
		# for (i in 1:nSamples){
		# 	.object@counts[[type]][[i]] <- emptyVec
		# }
		if (.object@sparseCounts){
			.object@counts[[type]] <- Matrix::sparseMatrix(i=c(), j=c(), x=1, dims=c(nRegs,nSamples))
		} else {
			.object@counts[[type]] <- matrix(as.integer(NA), nrow=nRegs, ncol=nSamples)
			if (.object@diskDump) .object@counts[[type]] <- as(.object@counts[[type]], "HDF5Array")
		}
		
		colnames(.object@counts[[type]]) <- sampleIds
		.object@countTransform[[type]] <- character(0)

		#Aggregate signal
		doAggr <- FALSE
		if (!is.null(signal) && !is.element(signal, c("insertions")) && is.element(aggrFun, c("sum", "mean", "median"))){
			doAggr <- TRUE
			signalGr <- getCoord(.object, signal)
			if (genome(regGr)[1]!=genome(signalGr)[1]){
				logger.warning(c("Potentially incompatible genome assemblies (object:", genome(signalGr)[1], ", regGr:", genome(regGr)[1], ")"))
			}
			oo <- findOverlaps(signalGr, regGr)
			if (any(duplicated(queryHits(oo)))) logger.info("Some signals map to multiple regions")
			dtC <- data.table(getCounts(.object, signal, i=queryHits(oo)), mergedIndex=subjectHits(oo))

			if (aggrFun=="sum") {
				rr <- dtC[, lapply(.SD, sum, na.rm=TRUE), by=.(mergedIndex)]
			} else if (aggrFun=="mean") {
				rr <- dtC[, lapply(.SD, mean, na.rm=TRUE), by=.(mergedIndex)]
			} else if (aggrFun=="median") {
				rr <- dtC[, lapply(.SD, median, na.rm=TRUE), by=.(mergedIndex)]
			} else {
				logger.error(c("Unknown signal count aggregation function:", aggrFun))
			}
			# .object@counts[[type]][rr[["mergedIndex"]],] <- rr[,!"mergedIndex"]
			.object@counts[[type]][rr[["mergedIndex"]],] <- as.matrix(rr[,!"mergedIndex"])
			rm(dtC, rr); cleanMem() #clean-up
		}
		if (!is.null(signal) && signal=="insertions"){
			doAggr <- TRUE
			if (bySample){
				i <- 0
				for (sid in sampleIds){
					i <- i + 1
					doMsg <- nSamples < 500 || (i %% 500 == 0)
					if (doMsg) logger.status(c("Aggregating counts for sample", sid, paste0("(", i, " of ", nSamples, ")"), "..."))
					.object@counts[[type]][,sid] <- as.matrix(countOverlaps(regGr, getInsertionSites(.object, samples=sid)[[1]], ignore.strand=TRUE))
				}
			} else {
				chunkIdxL <- list(1:nSamples)
				if (!is.null(chunkSize) & chunkSize < nSamples){
					chunkIdxL <- split(1:nSamples, as.integer((1:nSamples - 1) / chunkSize))
				}
				scm <- NULL
				for (k in seq_along(chunkIdxL)){
					logger.start(c("Processing chunk", k, "of", length(chunkIdxL)))
						nSamples_cur <- length(chunkIdxL[[k]])
						sampleIds_cur <- sampleIds[chunkIdxL[[k]]]
						logger.status("Retrieving joined insertion sites ...")
						fGrl <- getFragmentGrl(.object, sampleIds_cur, asGRangesList=FALSE)
						nFrags <- sapply(fGrl, length)
						fGr <- do.call("c", unname(fGrl))
						# fGr <- unlist(fGrl, use.names=FALSE)
						elementMetadata(fGr)[,".sampleIdx"] <- rep(seq_along(nFrags), nFrags)
						igr <- getInsertionSitesFromFragmentGr(fGr)
						logger.status("computing overlaps with regions ...") #TODO: this step can take a lot of memory (reduce chunksize?)
						oo <- findOverlaps(regGr, igr, ignore.strand=TRUE)
						# convenient: when constructing sparse matrices, if (i,j) indices are repeated, the x-values are summed up
						scm_cur <- Matrix::sparseMatrix(
							i=queryHits(oo),
							j=elementMetadata(igr)[subjectHits(oo), ".sampleIdx"],
							x=rep(1, length(queryHits(oo))),
							dims=c(nRegs, nSamples_cur)
						)
						colnames(scm_cur) <- sampleIds_cur

						scm <- cbind(scm, scm_cur)
						rm(igr, oo, scm_cur); cleanMem()
					logger.completed()
				}
				logger.status("Adding data to object ...")
				if (.object@sparseCounts){
					.object@counts[[type]] <- scm
				} else {
					.object@counts[[type]][,] <- as.matrix(scm)
				}
			}
		}

		if (doAggr){
			logger.info(c("Aggregated signal counts across", nrow(.object@counts[[type]]), "regions"))
			if (.object@sparseCounts) {
				hasNas <- anyNA(.object@counts[[type]]) # this should not be TRUE
				if (prod(dim(.object@counts[[type]])) < .Machine$integer.max && hasNas) {
					hasValM <- !is.na(.object@counts[[type]]) # TODO: this can take a lot of memory. Use a more efficient method to detect NAs in sparse matrices
					hasValM <- hasValM & .object@counts[[type]] != 0
				} else {
					if (hasNas) logger.warning("Sparse matrix with NAs detected --> assuming NAs are 0")
					# avoid this error: "Cholmod error 'problem too large'"
					hasValM <- .object@counts[[type]] != 0
				}
			} else {
				hasValM <- !is.na(.object@counts[[type]]) & .object@counts[[type]] > 0
			}
			rows2keep <- rsFun(hasValM) > 0
			logger.info(c("  of which", sum(rows2keep), "regions contained signal counts"))
			#discard regions where all signal counts are unobserved
			if (dropEmpty){
				.object@coord[[type]]   <- .object@coord[[type]][rows2keep]
				.object@counts[[type]]  <- .object@counts[[type]][rows2keep,]
			}
		}
		if (.object@sparseCounts) .object@counts[[type]] <- drop0(.object@counts[[type]])

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("mergeSamples")) {
	setGeneric(
		"mergeSamples",
		function(.object, ...) standardGeneric("mergeSamples"),
		signature=c(".object")
	)
}
#' mergeSamples-methods
#'
#' Merge signal and insertion data across samples
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param mergeGroups  factor or character vector or column name in sample annotation table.
#'                Can alternatively be a (named) list containing sample indices or names
#'                for each group to merge.
#' @param countAggrFun aggregation function for signal counts.
#'                Currently \code{sum} (default), \code{mean} and \code{median} are supported.
#' @return a new \code{\linkS4class{DsATAC}} object with samples merged
#'
#' 
#' @rdname mergeSamples-DsATAC-method
#' @docType methods
#' @aliases mergeSamples
#' @aliases mergeSamples,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("mergeSamples",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		mergeGroups,
		countAggrFun="sum"
	) {
		inputObj <- .object
		if (!is.element(countAggrFun, c("sum", "mean", "median"))){
			logger.error(c("Unknown signal count aggregation function:", countAggrFun))
		}
		if (.object@sparseCounts && is.element(countAggrFun, c("mean", "median"))){
			logger.warning("Mean and median merging of samples can be slow due to conversion of sparse matrices")
		}
		sampleNames <- getSamples(.object)
		nSamples <- length(sampleNames)
		ph <- getSampleAnnot(.object)
		rownames(ph) <- sampleNames # not sure if really necessary
		if (is.character(mergeGroups) && length(mergeGroups) == 1 && is.element(mergeGroups, colnames(ph))){
			mergeGroups <- ph[,mergeGroups]
		}
		# index list supplied
		if (is.list(mergeGroups)){
			if (is.null(names(mergeGroups))) names(mergeGroups) <- paste0("g", seq_along(mergeGroups))
			# convert to integer indices
			mergeGroups <- lapply(mergeGroups, FUN=function(x){
				if (is.character(x)){
					return(match(x, sampleNames))
				} else if (is.logical(x)){
					return(which(x))
				} else {
					return(x)
				}
			})
			ph_cat <- do.call(rbind, lapply(names(mergeGroups), FUN=function(gn){
				data.frame(ph[mergeGroups[[gn]],,drop=FALSE], .mergeGroup=gn, stringsAsFactors=FALSE)
			}))
			phm <- muRtools::aggregateDf(ph_cat, ph_cat[,".mergeGroup"])
			phm <- phm[names(mergeGroups), ]
			mgL <- mergeGroups
		} else {
			if ((!is.factor(mergeGroups) && !is.character(mergeGroups)) || length(mergeGroups) != nrow(ph)){
				logger.error("Invalid merge groups")
			}
			if (is.factor(mergeGroups)) mergeGroups <- as.character(mergeGroups)

			phm <- muRtools::aggregateDf(ph, mergeGroups)
			mgL <- lapply(rownames(phm), FUN=function(mg){which(mergeGroups==mg)})
			names(mgL) <- rownames(phm)
		}

		.object@sampleAnnot <- phm

		#count data
		mergeFun <- NULL
		if(countAggrFun=="sum"){
			mergeFun <- function(X){rowSums(X, na.rm=TRUE)}
			if (.object@diskDump) {
				mergeFun <- function(X){BiocGenerics::rowSums(X, na.rm=TRUE)}
			} else if (.object@sparseCounts) {
				mergeFun <- function(X){Matrix::rowSums(X, na.rm=TRUE)}
			}
		} else if(countAggrFun=="mean"){
			mergeFun <- function(X){rowMeans(X, na.rm=TRUE)}
		} else if(countAggrFun=="median"){
			mergeFun <- function(X){matrixStats::rowMedians(X, na.rm=TRUE)}
		}
		regTypes <- getRegionTypes(.object)
		for (rt in regTypes){
			logger.status(paste0("Merging samples (region set: '", rt, "')..."))
			cm <- ChrAccR::getCounts(.object, rt, asMatrix=FALSE) #design question: should we set the parameter naIsZero==FALSE?
			cmm <- do.call("cbind", lapply(mgL, FUN=function(iis){
				if (.object@sparseCounts && is.element(countAggrFun, c("mean", "median"))){
					aggMat <- ChrAccR::getCounts(.object, rt, j=iis, asMatrix=TRUE) #design question: should we set the parameter naIsZero==FALSE?
				} else {
					aggMat <- cm[,iis,drop=FALSE]
				}
				return(mergeFun(aggMat))
			}))
			# .object@counts[[rt]] <- data.table(cmm)
			if (.object@sparseCounts) {
				.object@counts[[rt]] <- as(cmm, "sparseMatrix")
				.object@counts[[rt]] <- drop0(.object@counts[[rt]])
			} else {
				.object@counts[[rt]] <- cmm
			}
			if (.object@diskDump) .object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
		}

		#insertion data: concatenate GRanges objects
		if (length(.object@fragments) == nSamples){
			logger.start(paste0("Merging sample fragment data"))			
				.object@fragments <- lapply(1:length(mgL), FUN=function(i){
					iis <- mgL[[i]]
					logger.status(c("Merged sample", i, "of", length(mgL)))
					curSns <- sampleNames[iis]
					curFragL <- getFragmentGrl(inputObj, curSns, asGRangesList=FALSE)
					nFrags <- sapply(curFragL, length)
					catRes <- do.call("c", unname(curFragL))
					# catRes <- unlist(curFragL, use.names=FALSE)
					elementMetadata(catRes)[,".sample"] <- rep(curSns, nFrags)

					if (.object@diskDump.fragments) {
						# logger.status(c("... saving to disk"))
						fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext=".rds")
						saveRDS(catRes, fn, compress=TRUE)
						catRes <- fn
					}
					return(catRes)
				})
				names(.object@fragments) <- names(mgL)
				if (.object@diskDump.fragments){
					# currently merging samples does not support multiple samples per fragment file
					.object@diskDump.fragments.nSamplesPerFile <- 1L
					# TODO: support multiple samples per merged fragment file
				}
			logger.completed()
		}

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("removeFragmentData")) {
	setGeneric(
		"removeFragmentData",
		function(object) standardGeneric("removeFragmentData"),
		signature=c("object")
	)
}
#' removeFragmentData-methods
#'
#' Removes fragment data from \code{\linkS4class{DsATAC}} object (e.g. to save space)
#'
#' @param object	\code{\linkS4class{DsATAC}} object
#' @return the modified object (without fragment data)
#'
#' @rdname removeFragmentData-DsATAC-method
#' @docType methods
#' @aliases removeFragmentData
#' @aliases removeFragmentData,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("removeFragmentData",
	signature(
		object="DsATAC"
	),
	function(
		object
	) {
		object@fragments <- list()
		return(object)
	}
)
#-------------------------------------------------------------------------------
# EXPERIMENTAL
if (!isGeneric("undiskFragmentData")) {
	setGeneric(
		"undiskFragmentData",
		function(object) standardGeneric("undiskFragmentData"),
		signature=c("object")
	)
}
#' undiskFragmentData-methods
#'
#' converts disk-backed fragment data to in-memory fragment data
#'
#' @param object	\code{\linkS4class{DsATAC}} object
#' @return the modified object
#'
#' @rdname undiskFragmentData-DsATAC-method
#' @docType methods
#' @aliases undiskFragmentData
#' @aliases undiskFragmentData,DsATAC-method
#' @author Fabian Mueller
#' @noRd
setMethod("undiskFragmentData",
	signature(
		object="DsATAC"
	),
	function(
		object
	) {
		object@fragments <- getFragmentGrl(object, names(object@fragments))
		object@diskDump.fragments <- FALSE
		return(object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("addCountDataFromBam")) {
	setGeneric(
		"addCountDataFromBam",
		function(.object, ...) standardGeneric("addCountDataFromBam"),
		signature=c(".object")
	)
}
#' addCountDataFromBam-methods
#'
#' Add count data to DsATAC object based on bam files
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param fns     a named vector of bam file locations. The names must correspond to sample identifiers in the object
#' @return a new \code{\linkS4class{DsATAC}} object with read counts for each sample and region set
#'
#' @rdname addCountDataFromBam-DsATAC-method
#' @docType methods
#' @aliases addCountDataFromBam
#' @aliases addCountDataFromBam,DsATAC-method
#' @author Fabian Mueller
#' @noRd
setMethod("addCountDataFromBam",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		fns
	) {
		if (!requireNamespace("GenomicAlignments")) logger.error(c("Could not load dependency: GenomicAlignments"))
		# TODO: adjust for Tn5 insertion:
		# + strand: i + 4
		# - strand: i - 5
		# Buenrostro, et al. (2013). Nature Methods, 10(12), 1213-1218.
		# "For peak calling and footprinting, we adjusted the read start sites to represent the center of the transposon binding event. Previous descriptions of the Tn5 transposase show that the transposon binds as a dimer and inserts two adaptors separated by 9 bp (ref. 11). Therefore, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the - strand were offset -5bp"
		ResizeReads <- function(reads, width=1, fix="start", ...) {
			reads <- as(reads, "GRanges")
			stopifnot(all(strand(reads) != "*"))
			resize(reads, width=width, fix=fix, ...)
		}
		sids <- names(fns)
		if (length(sids)!=length(fns)){
			logger.error("The vector specifying filenames (fns) must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}
		rts <- getRegionTypes(.object)
		for (rt in rts){
			logger.status(c("Counting reads in region set:", rt))
			gr <- getCoord(.object, rt)
			ov.rse <- summarizeOverlaps(gr, fns, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, preprocess.reads=ResizeReads)
			.object@counts[[rt]][,sids] <- as.matrix(assays(ov.rse)$count)
			if (.object@sparseCounts) .object@counts[[rt]] <- drop0(.object@counts[[rt]])
			.object@countTransform[[rt]] <- character(0)
		}

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("addCountDataFromGRL")) {
	setGeneric(
		"addCountDataFromGRL",
		function(.object, ...) standardGeneric("addCountDataFromGRL"),
		signature=c(".object")
	)
}
#' addCountDataFromGRL-methods
#'
#' Add count data to DsATAC object based on a list of GRanges objects 
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param grl     NAMED \code{GRangesList} or NAMED list of \code{GRanges} objects. Names must correspond to sample ids in the object
#' @param silent  limit log messages
#' @return a new \code{\linkS4class{DsATAC}} object with read counts for each sample and region set
#'
#' @rdname addCountDataFromGRL-DsATAC-method
#' @docType methods
#' @aliases addCountDataFromGRL
#' @aliases addCountDataFromGRL,DsATAC-method
#' @author Fabian Mueller
#' @noRd
setMethod("addCountDataFromGRL",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		grl,
		silent=FALSE
	) {
		sids <- names(grl)
		if (length(sids)!=length(grl)){
			logger.error("The list of GRanges must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}
		if (class(grl)=="GRangesList") {
			nElem <- elementNROWS(grl)
			gr.c <- unlist(grl, use.names=FALSE)
		} else {
			nElem <- sapply(grl, length)
			gr.c <- do.call("c", unname(grl))
		}
		
		if (length(gr.c) < 1) logger.error("[addCountDataFromGRL] invalid GRL: Must be of length 1 or more")
		sampleIds <- rep(names(grl), times=nElem)
		sampleIds.cm <- getSamples(.object)
		rts <- getRegionTypes(.object)

		for (rt in rts){
			if (!silent) logger.status(c("Counting reads in region set:", rt))
			gr.ds <- getCoord(.object, rt)
			oo <- findOverlaps(gr.ds, gr.c, ignore.strand=TRUE)
			idxDt <- data.table::as.data.table(cbind(
				queryHits(oo), #row indices (regions) in count matrix
				match(sampleIds[subjectHits(oo)], sampleIds.cm) #column indices (samples) in count matrix
			))
			# count the number of occurrences between each index pair
			idxDt <- idxDt[,.N, by=names(idxDt)]
			# if you in the future encounter some weird error here (e.g. 'Error in `[.data.frame`(x, i, j) : undefined columns selected')
			# this is likely due to some weird namespace issue. --> keep data.table in the 'depends' section of DESCRIPTION
			idxM <- as.matrix(idxDt[,c(1,2)])
			if (.object@diskDump){
				# DelayedArray does not support indexing via matrix
				# --> workaround
				perColumnIdxIdx <- tapply(1:nrow(idxM),idxM[,2],c) #index in the index matrix for each column
				for (cis in names(perColumnIdxIdx)){
					ci <- as.integer(cis)
					ii <- perColumnIdxIdx[[cis]]
					.object@counts[[rt]][idxM[ii,1],ci] <- as.matrix(idxDt$N[ii])
				}
			} else {
				.object@counts[[rt]][idxM] <- idxDt$N
			}
			
			if (.object@sparseCounts) .object@counts[[rt]] <- drop0(.object@counts[[rt]])
		}
		return(.object)
	}
)
# # Old, slower method f(for reference)
# setMethod("addCountDataFromGRL",
# 	signature(
# 		.object="DsATAC"
# 	),
# 	function(
# 		.object,
# 		grl
# 	) {
# 		sids <- names(grl)
# 		if (length(sids)!=length(grl)){
# 			logger.error("The list of GRanges must be named")
# 		}
# 		if (!all(sids %in% getSamples(.object))){
# 			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
# 		}
# 		rts <- getRegionTypes(.object)
# 		for (rt in rts){
# 			logger.status(c("Counting reads in region set:", rt))
# 			gr.ds <- getCoord(.object, rt)
# 			for (sid in sids){
# 				# logger.status(c("sample:", sid))
# 				gr.c <- grl[[sid]]
# 				.object@counts[[rt]][,sid] <- as.matrix(countOverlaps(gr.ds, gr.c, ignore.strand=TRUE))
# 			}
# 			if (.object@sparseCounts) .object@counts[[rt]] <- drop0(.object@counts[[rt]])
# 		}

# 		return(.object)
# 	}
# )
#-------------------------------------------------------------------------------
if (!isGeneric("addSignalDataFromGRL")) {
	setGeneric(
		"addSignalDataFromGRL",
		function(.object, ...) standardGeneric("addSignalDataFromGRL"),
		signature=c(".object")
	)
}
#' addSignalDataFromGRL-methods
#'
#' Add signal data to DsATAC object based on a list of GRanges objects 
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param grl     NAMED \code{GRangesList} or NAMED list of \code{GRanges} objects. Names must correspond to sample ids in the object
#' @param aggrFun aggregation method
#' @return a new \code{\linkS4class{DsATAC}} object with counts/signal for each sample and region set
#'
#' @rdname addSignalDataFromGRL-DsATAC-method
#' @docType methods
#' @aliases addSignalDataFromGRL
#' @aliases addSignalDataFromGRL,DsATAC-method
#' @author Fabian Mueller
#' @noRd
setMethod("addSignalDataFromGRL",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		grl,
		aggrFun=function(x){mean(x, na.rm=TRUE)}
	) {
		sids <- names(grl)
		if (length(sids)!=length(grl)){
			logger.error("The list of GRanges must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}
		rts <- getRegionTypes(.object)
		for (rt in rts){
			logger.status(c("Aggregating scores in region set:", rt))
			gr.ds <- getCoord(.object, rt)
			for (sid in sids){
				gr.c <- grl[[sid]]
				scs <- elementMetadata(gr.c)[,"score"]
				if (is.null(scs)) {
					logger.error(c("GRanges must contain a score column (sample", sid, ")"))
				}
				# TODO: check for correctness
				oo <- findOverlaps(gr.ds, gr.c, ignore.strand=TRUE)
				.object@counts[[rt]][sort(unique(queryHits(oo))), sid] <- as.matrix(tapply(scs[subjectHits(oo)], queryHits(oo), aggrFun))
			}
			if (.object@sparseCounts) .object@counts[[rt]] <- drop0(.object@counts[[rt]])
		}

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("addInsertionDataFromBam")) {
	setGeneric(
		"addInsertionDataFromBam",
		function(.object, ...) standardGeneric("addInsertionDataFromBam"),
		signature=c(".object")
	)
}
#' addInsertionDataFromBam-methods
#'
#' Add insertion data to DsATAC object based on bam or fragment bed files
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param fns     a named vector of bam or bed file locations. The names must correspond to sample identifiers in the object
#' @param pairedEnd flag indicating whether the bam files are from paired-end sequencing
#' @param .diskDump for internal use only (should fragment data be stored as RDS instead of GRanges)
#' @param mode    either \code{"bam"} (default) for reading bam files or \code{"fragBed"} for fragment bed files
#' @return a new \code{\linkS4class{DsATAC}} object with fragments/insertions for each sample
#'
#' @rdname addInsertionDataFromBam-DsATAC-method
#' @docType methods
#' @aliases addInsertionDataFromBam
#' @aliases addInsertionDataFromBam,DsATAC-method
#' @author Fabian Mueller
#' @noRd
setMethod("addInsertionDataFromBam",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		fns,
		pairedEnd=TRUE,
		.diskDump=.object@diskDump.fragments,
		mode = "bam"
	) {
		offsetTn <- TRUE
		sids <- names(fns)
		if (length(sids)!=length(fns)){
			logger.error("The vector specifying filenames (fns) must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}

		# only relevant if fragments from multiple samples are written to the same file (chunkedFragmentFiles)
		chunkedFragmentFiles <- .object@diskDump.fragments && .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
		curChunkL <- list()
		curChunkFn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")

		for (sid in sids){
			logger.start(c("Reading insertion data for sample:", sid))
				if (!is.null(.object@fragments[[sid]])){
					logger.warning(c("Overwriting insertion data for sample:", sid))
				}
				fragGr <- NULL
				if (mode=="bam"){
					ga <- NULL
					if (pairedEnd){
						ga <- GenomicAlignments::readGAlignmentPairs(fns[sid], use.names=FALSE, param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isProperPair=TRUE)))
					} else {
						ga <- GenomicAlignments::readGAlignments(fns[sid], use.names=FALSE)
					}
					ga <- setGenomeProps(ga, .object@genome, onlyMainChrs=TRUE, silent=TRUE)
					fragGr <- getATACfragments(ga, offsetTn=TRUE)
				} else if (mode=="fragBed"){
					fragGr <- rtracklayer::import(fns[sid], format = "BED")
					fragGr <- setGenomeProps(fragGr, .object@genome, onlyMainChrs=TRUE, silent=TRUE)

					if (offsetTn){
						# shift inserts inward due to the Tn5 dimer offset:
						# --> +4 bp
						# -------------------------
						# -------------------------
						#                 -4 bp <--
						fragGr <- GenomicRanges::resize(fragGr, width(fragGr)-4, fix="end", ignore.strand=TRUE)
						fragGr <- GenomicRanges::resize(fragGr, width(fragGr)-4, fix="start", ignore.strand=TRUE)
					}
				}
				logger.info(c(length(fragGr), "fragments found"))
				if (.diskDump){
					if (chunkedFragmentFiles){
						curChunkL[[sid]] <- fragGr
						fragGr <- curChunkFn
						if (length(curChunkL) >= .object@diskDump.fragments.nSamplesPerFile){
							saveRDS(curChunkL, curChunkFn, compress=TRUE)
							# reset after writing chunk
							curChunkL <- list()
							curChunkFn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
						}
					} else {
						fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
						saveRDS(fragGr, fn, compress=TRUE)
						fragGr <- fn
					}
				}
				.object@fragments[[sid]] <- fragGr
			logger.completed()
		}
		if (.diskDump && length(curChunkL) > 0){
			saveRDS(curChunkL, curChunkFn, compress=TRUE)
		}

		return(.object)
	}
)
################################################################################
# Manipulating DsATAC objects
################################################################################
if (!isGeneric("removeRegions")) {
	setGeneric(
		"removeRegions",
		function(.object, ...) standardGeneric("removeRegions"),
		signature=c(".object")
	)
}
#' removeRegions-methods
#'
#' Remove the specified sites or regions from an object
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param indices a vector of indices of sites/regions to be removed. Can be numeric, integer or logical.
#' @param type    character string specifying a name for the region type (sefault: sites)
## @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsATAC}} object with sites/regions removed
#' 
#' @rdname removeRegions-DsATAC-method
#' @docType methods
## @aliases removeRegions
#' @aliases removeRegions,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("removeRegions",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		indices,
		type
		# reaggregate=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))

		if (!is.vector(indices) || !(is.numeric(indices) || is.logical(indices))){
			logger.error(c("Unsupported type for index vector"))
		}
		nRegs <- getNRegions(.object, type)
		inds2keep <- rep(TRUE, nRegs)
		if (is.numeric(indices)){
			if (any(indices > nRegs | indices < 1)) {
				logger.error(c("Invalid values in indices"))
			}
			inds2keep[indices] <- FALSE
		} else if (is.logical(indices)){
			inds2keep <- !indices
		}
		if (sum(inds2keep)>=nRegs){
			logger.info("Nothing to be done: keeping object as is")
			return(.object)
		}

		.object@coord[[type]] <- .object@coord[[type]][inds2keep]
		# .object@counts[[type]]  <- .object@counts[[type]][inds2keep,]
		.object@counts[[type]] <- ChrAccR::getCounts(.object, type, i=inds2keep, j=NULL, asMatrix=FALSE, naIsZero=FALSE)

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("removeRegionType")) {
	setGeneric(
		"removeRegionType",
		function(.object, ...) standardGeneric("removeRegionType"),
		signature=c(".object")
	)
}
#' removeRegionType-methods
#'
#' Remove the specified region type from an object
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param type    character string specifying a name for the region type (sefault: sites)
#' @return a new \code{\linkS4class{DsATAC}} object with the region type removed
#' 
#' @rdname removeRegionType-DsATAC-method
#' @docType methods
#' @aliases removeRegionType
#' @aliases removeRegionType,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("removeRegionType",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type
	) {
		rts <- getRegionTypes(.object)
		if (!is.element(type, rts)) logger.error(c("Unsupported region type:", type))
		rts <- setdiff(rts, type)

		.object@coord <- .object@coord[rts]
		.object@counts <- .object@counts[rts]
		.object@countTransform <- .object@countTransform[rts]

		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("removeSamples")) {
	setGeneric(
		"removeSamples",
		function(.object, ...) standardGeneric("removeSamples"),
		signature=c(".object")
	)
}
#' removeSamples-methods
#'
#' Remove samples from a \code{\linkS4class{DsATAC}} object
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param indices a vector of indices of samples to be removed. Can be numeric, integer or logical.
#' @return a new \code{\linkS4class{DsATAC}} object with sites/regions removed
#' 
#' @rdname removeSamples-DsATAC-method
#' @docType methods
#' @aliases removeSamples
#' @aliases removeSamples,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("removeSamples",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		indices
		# reaggregate=TRUE
	) {
		if (!is.vector(indices) || !(is.numeric(indices) || is.logical(indices))){
			logger.error(c("Unsupported type for index vector"))
		}
		nSamples <- length(getSamples(.object))
		inds2keep <- rep(TRUE, nSamples)
		if (is.numeric(indices)){
			if (any(indices > nSamples | indices < 1)) {
				logger.error(c("Invalid values in indices"))
			}
			inds2keep[indices] <- FALSE
		} else if (is.logical(indices)){
			inds2keep <- !indices
		}
		if (sum(inds2keep)>=nSamples){
			logger.info("Nothing to be done: keeping object as is")
			return(.object)
		}
		.object@sampleAnnot <- .object@sampleAnnot[inds2keep,]

		for (rt in getRegionTypes(.object)){
			# .object@counts[[rt]]  <- .object@counts[[rt]][,..inds2keep]
			.object@counts[[rt]]  <- .object@counts[[rt]][,inds2keep]
		}

		if (length(.object@fragments) == nSamples){
			.object@fragments <- .object@fragments[inds2keep]
			chunkedFragmentFiles <- .object@diskDump.fragments && .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
			if (chunkedFragmentFiles){
				logger.warning("Repackaging of chunked fragment files not supported yet. Use 'ChrAccR:::undiskFragmentData' if you want to save the dataset later.")
				# TODO: Make this possible. See code for 'filterChroms' for template
			}
		}
		
		return(.object)
	}
)
#-------------------------------------------------------------------------------
#' Subsetting DsATAC datasets by sample
#' 
#' NOTE: '[' operator for DsATAC does not reorder samples or deal with index multiplicity
#'
#' @param x DsATAC object
#' @param i sample names or indices
setMethod("[", "DsATAC",
	function(x, i){
		if (missing(i)) return(x)
		sampleIds <- getSamples(x)
		nSamples <- length(sampleIds)
		inds2keep <- rep(FALSE, nSamples)
		if (is.numeric(i)){
			if (any(i > nSamples | i < 1)) {
				logger.error(c("Invalid values in i"))
			}
			inds2keep[i] <- TRUE
		} else if (is.logical(i)){
			inds2keep <- i
		} else if (is.character(i)){
			if (!all(i %in% sampleIds)){
				logger.error(c("Could not find the following sample ids:", paste(setdiff(i, sampleIds), collapse=",")))
			}
			inds2keep <- sampleIds %in% i
		}
		if (sum(inds2keep)>=nSamples){
			logger.info("Nothing to be done: keeping object as is")
			return(x)
		}
		if (is.character(i) || is.numeric(i)){
			logger.info("NOTE: '[' operator for DsATAC does not reorder samples")
			if (any(duplicated(i))){
				logger.info("NOTE: '[' operator for DsATAC does not deal with index multiplicity")
			}
		}
		
		return(removeSamples(x, !inds2keep))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("transformCounts")) {
	setGeneric(
		"transformCounts",
		function(.object, ...) standardGeneric("transformCounts"),
		signature=c(".object")
	)
}
#' transformCounts-methods
#'
#' transform count data for an ATAC seq dataset
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param method  transformation method to be applied. Currently only 'log2', 'log10', 'quantile' (quantile normalization), 'percentile' (percentile normalization),'rankPerc' (rank percentile), 'vst' (DESeq2 Variance Stabilizing Transformation), 'batchCorrect' (limma batch effect removal), tf-idf', 'CPM' (counts per million), and 'RPKM' (RPKM normalization) are supported
#' @param regionTypes character vector specifying a name for the region type in which count data should be normalized(default: all region types)
#' @param ...    other arguments depending on the \code{method} used. For \code{'batchCorrect'} it should be arguments passed on to \code{limma::removeBatchEffect} (most importantly, the \code{batch} argument).
#' @return a new \code{\linkS4class{DsATAC}} object with normalized count data
#' 
#' @rdname transformCounts-DsATAC-method
#' @docType methods
#' @aliases transformCounts
#' @aliases transformCounts,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("transformCounts",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		method="quantile",
		regionTypes=getRegionTypes(.object),
		...
	) {
		if (!all(regionTypes %in% getRegionTypes(.object))){
			logger.error(c("Unsupported region type:", paste(setdiff(regionTypes, getRegionTypes(.object)), collapse=", ")))
		}
		if (!is.element(method, c("quantile", "percentile", "rankPerc", "log2", "log10", "RPKM", "CPM", "RPKMtss", "CPMtss", "vst", "batchCorrect", "tf-idf"))) logger.error(c("Unsupported normalization method type:", method))

		# choose appropriate functions for row and column functions
		# depending on the matrix type
		rsFun <- rowSums
		csFun <- colSums
		# sparse matrices
		if (.object@sparseCounts){
			rsFun <- Matrix::rowSums
			csFun <- Matrix::colSums
		}
		# DelayedArray
		if (.object@diskDump){
			rsFun <- BiocGenerics::rowSums
			csFun <- BiocGenerics::colSums
		}

		if (method == "quantile"){
			logger.start(c("Performing quantile normalization"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					# .object@counts[[rt]] <- data.table(normalize.quantiles(as.matrix(.object@counts[[rt]])))
					.object@counts[[rt]] <- preprocessCore::normalize.quantiles(ChrAccR::getCounts(.object, rt, asMatrix=TRUE)) #design question: should we set the parameter naIsZero==FALSE?
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("quantileNorm", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "rankPerc"){
			logger.start(c("Applying rank percentile transformation"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					.object@counts[[rt]] <- muRtools::normalizeRank(ChrAccR::getCounts(.object, rt, asMatrix=TRUE), out="percentile")
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("rankPercNorm", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "percentile"){
			logger.start(c("Applying percentile transformation"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					.object@counts[[rt]] <- muRtools::normalizePercentile(ChrAccR::getCounts(.object, rt, asMatrix=TRUE))
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("percentile", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (is.element(method, c("RPKM", "CPM", "RPKMtss", "CPMtss"))){
			doRpkm <- grepl("^RPKM", method)
			doTssBg <- grepl("tss$", method)
			if (doTssBg){
				annoPkg <- getChrAccRAnnotationPackage(.object@genome)
				if (is.null(annoPkg)) logger.error("Annotation package needed")
				tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
				tssGr <- suppressWarnings(trim(resize(tssGr, width=101L, fix="center", ignore.strand=TRUE)))
			}
			logger.start(c("Performing", method, "normalization"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					cm <- as.matrix(.object@counts[[rt]])
					# cm <- .object@counts[[rt]]
					regLen <- rep(1000L, getNRegions(.object, rt)) # for CPM, treat all regions as length 1000
					rGr <- NULL
					if (doRpkm || doTssBg) rGr <- getCoord(.object, rt)
					if (doRpkm) regLen <- width(rGr)
					sizeFac <- NULL
					if (doTssBg){
						ovTss <- overlapsAny(rGr, tssGr, ignore.strand=TRUE)
						sizeFac <- matrix(csFun(cm[ovTss,,drop=FALSE], na.rm=TRUE), ncol=ncol(cm), nrow=nrow(cm), byrow=TRUE)
					} else {
						sizeFac <- matrix(csFun(cm, na.rm=TRUE), ncol=ncol(cm), nrow=nrow(cm), byrow=TRUE)
					}
					# .object@counts[[rt]] <- data.table(cm/(regLen * sizeFac) * 1e3 * 1e6)
					.object@counts[[rt]] <- cm/(regLen * sizeFac) * 1e3 * 1e6
					if (.object@diskDump) .object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("RPKM", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (is.element(method, c("log2", "log10"))){
			c0 <- 1
			logFun <- log10
			if (method=="log2") logFun <- log2
			logger.start(c(method, "transforming counts"))
				if (.object@sparseCounts) logger.error("Generating sparse matrix for matrix with possible true zero entries (", method, ")")
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cm <- as.matrix(.object@counts[[rt]])
					cnames <- colnames(cm)
					idx <- !is.na(cm) & cm!=0
					cm[idx] <- logFun(cm[idx] + c0)
					.object@counts[[rt]] <- cm
					if (.object@diskDump) .object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c(method, .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "vst"){
			logger.start(c("Applying DESeq2 VST"))
				if (.object@sparseCounts) logger.warning("Generating sparse matrix for matrix with possible true zero entries (VST)")
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					dds <- DESeq2::DESeqDataSet(getCountsSE(.object, rt), design=~1)
					# .object@counts[[rt]] <- data.table(assay(DESeq2::vst(dds, blind=TRUE)))
					.object@counts[[rt]] <- assay(DESeq2::vst(dds, blind=TRUE))
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					.object@countTransform[[rt]] <- c("deseq.vst", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "batchCorrect"){
			logger.start(c("Applying batch effect correction"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					.object@counts[[rt]] <- limma::removeBatchEffect(.object@counts[[rt]], ...)
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					.object@countTransform[[rt]] <- c("batchCorrect", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "tf-idf"){
			# TF-IDF transformation as applied in LSI (Shendure lab) (Cusanovich, et al. (2018). A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 1-35)
			# Recycled some code from: https://github.com/shendurelab/mouse-atac
			logger.start(c("Applying TF-IDF transformation"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cm <- .object@counts[[rt]] > 0 #indicator matrix: are there any counts in that region
					cnames <- colnames(cm)
					if (class(cm)=="lgCMatrix"){
						tf <- Matrix::t(Matrix::t(cm) / Matrix::colSums(cm)) #term frequency
						idf <- tf * log(1 + ncol(cm) / Matrix::rowSums(cm)) # inverse document frequency
					} else {
						cm <- cm & !is.na(.object@counts[[rt]])
						tf <- t(t(cm) / csFun(cm)) #term frequency
						idf <- tf * log(1 + ncol(cm) / rsFun(cm)) # inverse document frequency
					}
					.object@counts[[rt]] <- idf
					if (.object@diskDump) .object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("TF-IDF", .object@countTransform[[rt]])
				}
			logger.completed()
		}
		return(.object)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("filterLowCovg")) {
	setGeneric(
		"filterLowCovg",
		function(.object, ...) standardGeneric("filterLowCovg"),
		signature=c(".object")
	)
}
#' filterLowCovg-methods
#'
#' Filter regions with low read counts
#'
#' @param .object     \code{\linkS4class{DsATAC}} object
#' @param thresh      regions with read counts below this threshold will be considered lowly covered regions (default: regions with fewer than 1 read will be discarded)
#' @param reqSamples  the percentile of samples required to meet or exceed the threshold in order for a region to be retained.
#'                    must be in the interval [0, 1) (default: 0.75 = 75 percent)
#' @param regionTypes character vector specifying the names of the region types to which filtering should be applied (default: all region types)
#' @return a new \code{\linkS4class{DsATAC}} object with low coverage regions removed
#' 
#' @rdname filterLowCovg-DsATAC-method
#' @docType methods
#' @aliases filterLowCovg
#' @aliases filterLowCovg,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("filterLowCovg",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		thresh=1L,
		reqSamples=0.75,
		regionTypes=getRegionTypes(.object)
	) {
		if (!all(regionTypes %in% getRegionTypes(.object))){
			logger.error(c("Unsupported region type:", paste(setdiff(regionTypes, getRegionTypes(.object)), collapse=", ")))
		}
		N <- length(getSamples(.object))
		numAllowed <- reqSamples
		if (numAllowed < 1 && numAllowed >=0){
			numAllowed <- as.integer(ceiling(numAllowed * N))
		}
		if (numAllowed > N || numAllowed < 1){
			logger.error(c("Invalid number of samples. Must be in the interval [0,1) or an integer in [1, n_samples]"))
		}
		percAllowed <- round(numAllowed/N, 2)
		logger.status(c("Removing regions with read counts lower than", thresh, "in more than", N-numAllowed, "samples", paste0("(", (1-percAllowed)*100,"%)")))
		for (rt in regionTypes){
			rsFun <- rowSums
			if (.object@sparseCounts) rsFun <- Matrix::rowSums
			rem <- rsFun(getCounts(.object, rt, naIsZero=TRUE, allowSparseMatrix=TRUE) >= thresh) < numAllowed
			nRem <- sum(rem)
			nRegs <- getNRegions(.object, rt)
			if (nRem >= nRegs) logger.error(c("Cannot remove all regions of type:", rt))
			if (nRem > 0){
				.object <- removeRegions(.object, rem, rt)
			}
			logger.info(c("Removed", nRem, "regions", paste0("(", round(nRem/nRegs, 4)*100, "%)"), "of type", rt))
		}
		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("filterChroms")) {
	setGeneric(
		"filterChroms",
		function(.object, ...) standardGeneric("filterChroms"),
		signature=c(".object")
	)
}
#' filterChroms-methods
#'
#' Filter out regions based on chromosome list
#'
#' @param .object     \code{\linkS4class{DsATAC}} object
#' @param exclChrom   vector of chromosome names to filter out
#' @return a new \code{\linkS4class{DsATAC}} object filtered for chromosomes
#' 
#' @rdname filterChroms-DsATAC-method
#' @docType methods
#' @aliases filterChroms
#' @aliases filterChroms,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("filterChroms",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		exclChrom=c("chrX", "chrY", "chrM")
	) {
		for (rt in getRegionTypes(.object)){
			isExclChrom <- as.character(seqnames(getCoord(.object, rt))) %in% exclChrom
			.object <- removeRegions(.object, isExclChrom, rt)
			logger.info(c("Removed", sum(isExclChrom), "regions", paste0("(", round(sum(isExclChrom)/length(isExclChrom), 4)*100, "%)"), "of type", rt))
		}
		if (length(.object@fragments) > 0){
			logger.status("Filtering fragment data")
			chunkedFragmentFiles <- .object@diskDump.fragments && .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
			if (chunkedFragmentFiles){
				isFile <- sapply(.object@fragments, is.character)
				if (!all(isFile)) logger.error("Expected all disk-dumped fragment files")
				fns <- unlist(.object@fragments)
				names(fns) <- names(.object@fragments)
				for (fn in unique(fns)){
					fragGrl <- readRDS(fn)
					sids <- names(fragGrl)
					# filter
					fragGrl_filt <- endoapply(fragGrl, FUN=function(x){
						idx <- !(as.character(seqnames(x)) %in% exclChrom)
						return(x[idx])
					})
					names(fragGrl_filt) <- sids
					# save the list object to disk
					fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
					saveRDS(fragGrl_filt, fn, compress=TRUE)
					# replace old filename references
					.object@fragments[sids] <- rep(list(fn), length(sids))
				}
			} else {
				for (sid in names(.object@fragments)){
					fragGr <- getFragmentGr(.object, sid)
					idx <- !(as.character(seqnames(fragGr)) %in% exclChrom)
					fragGr <- fragGr[idx]
					if (.object@diskDump.fragments){
						fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
						saveRDS(fragGr, fn, compress=TRUE)
						fragGr <- fn
					}
					.object@fragments[[sid]] <- fragGr
				}
			}
		}
		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("filterByGRanges")) {
	setGeneric(
		"filterByGRanges",
		function(.object, ...) standardGeneric("filterByGRanges"),
		signature=c(".object")
	)
}
#' filterByGRanges-methods
#'
#' Filter out regions based on a GRanges object
#'
#' @param .object     \code{\linkS4class{DsATAC}} object
#' @param gr   		  \code{GRanges} object used for filtering
#' @param method      character string specifying the method. Can be \code{"white"} to treat \code{gr} as a whitelist
#'                    (i.e. only regions and fragments overlapping with it are retained) or \code{"black"} (default) to treat it as a blacklist.
#' @return a new, filtered \code{\linkS4class{DsATAC}} object
#' 
#' @rdname filterByGRanges-DsATAC-method
#' @docType methods
#' @aliases filterByGRanges
#' @aliases filterByGRanges,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("filterByGRanges",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		gr,
		method="black"
	) {
		if (class(gr)!="GRanges") logger.error("Invalid gr object. Expected GRanges.")
		if (!is.element(method, c("white", "black"))) logger.error(c("invalid method for filtering:", method))

		for (rt in getRegionTypes(.object)){
			isExcl <- overlapsAny(getCoord(.object, rt), gr) # blacklist
			if (method=="white") isExcl <- !isExcl
			.object <- removeRegions(.object, isExcl, rt)
			logger.info(c("Removed", sum(isExcl), "regions", paste0("(", round(sum(isExcl)/length(isExcl), 4)*100, "%)"), "of type", rt))
		}
		if (length(.object@fragments) > 0){
			logger.status("Filtering fragment data")
			chunkedFragmentFiles <- .object@diskDump.fragments && .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
			if (chunkedFragmentFiles){
				isFile <- sapply(.object@fragments, is.character)
				if (!all(isFile)) logger.error("Expected all disk-dumped fragment files")
				fns <- unlist(.object@fragments)
				names(fns) <- names(.object@fragments)
				for (fn in unique(fns)){
					fragGrl <- readRDS(fn)
					sids <- names(fragGrl)
					# filter
					fragGrl_filt <- endoapply(fragGrl, FUN=function(x){
						idx <- overlapsAny(x, gr) # whitelist
						if (method=="black") idx <- !idx
						return(x[idx])
					})
					names(fragGrl_filt) <- sids
					# save the list object to disk
					fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
					saveRDS(fragGrl_filt, fn, compress=TRUE)
					# replace old filename references
					.object@fragments[sids] <- rep(list(fn), length(sids))
				}
			} else {
				for (sid in names(.object@fragments)){
					fragGr <- getFragmentGr(.object, sid)
					idx <- overlapsAny(fragGr, gr) # whitelist
					if (method=="black") idx <- !idx
					fragGr <- fragGr[idx]
					if (.object@diskDump.fragments){
						fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
						saveRDS(fragGr, fn, compress=TRUE)
						fragGr <- fn
					}
					.object@fragments[[sid]] <- fragGr
				}
			}
		}
		return(.object)
	}
)
################################################################################
# Analysis Utils
################################################################################
if (!isGeneric("regionSetCounts")) {
	setGeneric(
		"regionSetCounts",
		function(.object, ...) standardGeneric("regionSetCounts"),
		signature=c(".object")
	)
}
#' regionSetCounts-methods
#'
#' Overlap the insertion data with a list of region sets
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param rsl     \code{GRangesList} or NAMED list of \code{GRanges} objects. Each element corresponds to a region set for which the summary statistics are reported
#' @param bySample for internal use: iterate over samples (instead of retrieving one giant insertion list for all samples) in order to save memory (at the tradeoff of compute time)
#' @return a matrix of overlap counts for each region set and sample
#'
#' @rdname regionSetCounts-DsATAC-method
#' @docType methods
#' @aliases regionSetCounts
#' @aliases regionSetCounts,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("regionSetCounts",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		rsl,
		bySample=FALSE
	) {
		nr <- length(rsl)
		res <- matrix(as.numeric(NA), nrow=nr, ncol=length(getSamples(.object)))
		if (bySample){
			if (class(rsl)=="GRangesList"){
				nc <- elementNROWS(rsl)
				rslGr <- unlist(rsl, use.names=FALSE)
			} else {
				nc <- sapply(rsl, length)
				rslGr <- do.call("c", unname(rsl))
			}
			
			idx.rsl <- rep(1:nr, times=nc)

			res <- do.call("cbind", lapply(getSamples(.object), FUN=function(sid){
				logger.status(c("Summarizing fragment counts for sample", sid))
				insGr <- getInsertionSites(.object, sid)[[1]]
				ov <- countOverlaps(rslGr, insGr, ignore.strand=TRUE)
				rr <- tapply(ov, idx.rsl, sum)
				rr <- rr[as.character(1:nr)] # make sure the counts are returned in the same order as in rsl (not sure, if really necessary)
				return(rr)
			}))
			rownames(res) <- names(rsl)
			colnames(res) <- getSamples(.object)
		} else {
			# insGrl <- getInsertionSites(.object)
			res <- countPairwiseOverlaps(rsl, getInsertionSites(.object), ignore.strand=TRUE)
		}
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getInsertionKmerFreq")) {
	setGeneric(
		"getInsertionKmerFreq",
		function(.object, ...) standardGeneric("getInsertionKmerFreq"),
		signature=c(".object")
	)
}
#' getInsertionKmerFreq-methods
#'
#' compute kmer frequencies at insertion sites for each sample
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param samples   sample identifiers
#' @param k          length of the kmer
#' @param normGenome should the result be normalized by genome-wide kmer frequencies
#' @return a \code{matrix} containing kmer frequencies (one row for each kmer and one column for each sample in the dataset)
#' 
#' @rdname getInsertionKmerFreq-DsATAC-method
#' @docType methods
#' @aliases getInsertionKmerFreq
#' @aliases getInsertionKmerFreq,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getInsertionKmerFreq",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		samples=getSamples(.object),
		k=6,
		normGenome=FALSE
	) {
		if (!requireNamespace("Biostrings")) logger.error(c("Could not load dependency: Biostrings"))
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		go <- getGenomeObject(.object@genome)
		res <- do.call("cbind", lapply(samples, FUN=function(sid){
			logger.status(c("Preparing insertion kmer-frequencies for sample", sid))
			insGr <-  trim(resize(shift(getInsertionSites(.object, sid)[[1]], -ceiling(k/2)), width=k, fix="start", ignore.strand=TRUE))
			kmerFreq <- Biostrings::oligonucleotideFrequency(Biostrings::Views(go, insGr), width=k, simplify.as="collapsed")
			return(kmerFreq)
		}))
		colnames(res) <- samples
		if (normGenome) {
			logger.status(c("Normalizing using genome-wide kmer-frequencies"))
			kmerFreq.g <- Biostrings::oligonucleotideFrequency(Biostrings::Views(go, getGenomeGr(.object@genome, onlyMainChrs=TRUE)), width=k, simplify.as="collapsed")
			res <- res/kmerFreq.g
		}
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("aggregateRegionCounts")) {
	setGeneric(
		"aggregateRegionCounts",
		function(.object, ...) standardGeneric("aggregateRegionCounts"),
		signature=c(".object")
	)
}
#' aggregateRegionCounts-methods
#'
#' Agregate counts across a set of regions, e.g. for footprinting analysis
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionGr   \code{GRanges} object specifying the regions to aggregate over
#' @param samples    sample identifiers
#' @param countAggrFun aggration function to be used for summarizing the insertion counts at each position. Possible values include \code{"sum"}, \code{"mean"}, and \code{"median"}
#' @param norm       method used for normalizing the resulting position-wise counts.
#'                   Currently only \code{'tailMean'} is supported, which computes normalization factors as the mean signal in the tails of the window
#' @param normTailW  fraction of the region window to be used on each side of the window to be used for normalization if \code{norm} is one of \code{'tailMean'}
#' @param kmerBiasAdj compute Tn5 bias and use it to adjust the counts as in Corces, et al., Science, (2018)
#' @param k          length of the kmer to be used for sequence bias correction. Only relevant if \code{kmerBiasAdj==TRUE}.
#' @param sampleCovg to save compute time, a sample coverage track list (as computed by \code{getCoverage(.object)}) can be supplied. If not, it will be computed on the fly.
#' @param sampleKmerFreqM to save compute time, a matrix of sample kmer frequency at insertion sites (as computed by \code{getInsertionKmerFreq(.object, ...)}) can be supplied.
#'                   If not, it will be computed on the fly. Only relevant if \code{kmerBiasAdj==TRUE}.
#' @param regionKmerFreqM to save compute time, a matrix of region kmer frequencies (kmers X window width).
#'                   Must have the same number of rows as the specified (or computed) \code{sampleKmerFreqM} (kmers)
#'                   and the same number of columns as the window width (median width of \code{regionGr}).
#'                   Only relevant if \code{kmerBiasAdj==TRUE}.
#' @param silent     limit log messages
#' @return a \code{data.frame} containing position-wise counts (raw, normalized and optionally Tn5-bias-corrected) for each sample
#' 
#' @rdname aggregateRegionCounts-DsATAC-method
#' @docType methods
#' @aliases aggregateRegionCounts
#' @aliases aggregateRegionCounts,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("aggregateRegionCounts",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionGr,
		samples=getSamples(.object),
		countAggrFun="sum",
		norm="tailMean",
		normTailW=0.1,
		kmerBiasAdj=TRUE,
		k=6,
		sampleCovg=NULL,
		sampleKmerFreqM=NULL,
		regionKmerFreqM=NULL,
		silent=FALSE
	) {
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		if (!is.element(countAggrFun, c("sum", "mean", "median"))) logger.error(c("Invalid value for countAggrFun:", countAggrFun))

		ww <- width(regionGr)
		wm <- as.integer(median(ww))
		idx <- ww==wm
		if (!all(idx)){
			if (!silent) logger.warning(c("not all elements in GRanges have the same width. --> discarding", sum(!idx), "of", length(idx), "regions that do not."))
			regionGr <- regionGr[idx]
		}

		if (kmerBiasAdj && !is.null(regionKmerFreqM) && ncol(regionKmerFreqM)!=wm){
			logger.error(c("Invalid value for regionKmerFreqM. expected a window size of", wm, "rows, but input has", ncol(regionKmerFreqM)))
		}

		if (!is.element(norm, c("tailMean"))) logger.error("Invalid value for 'norm' (normalization method)")
		if (normTailW < 1 && normTailW>0){
			normTailW <- ceiling(wm*normTailW)
		} else if (is.integer(normTailW) && normTailW > 0 && normTailW <= wm/2){
			# do nothing: treat integer values as basepairs
			normTailW <- normTailW
		} else {
			logger.error("Invalid value for tail-normalization window")
		}
		if (is.null(sampleCovg)){
			if (!silent) logger.start("Computing sample coverage")
				sampleCovg <- getCoverage(.object, samples=samples)
			if (!silent) logger.completed()
		} else {
			if (!all(samples %in% names(sampleCovg))) logger.error("'sampleCovg' does not cover all samples")
		}
		kmerFreqM <- NULL
		go <- NULL
		tn5bias <- rep(as.numeric(NA), wm)
		if (kmerBiasAdj){
			go <- getGenomeObject(.object@genome)
			if (is.null(sampleKmerFreqM)) {
				logger.start("Computing sample kmer frequencies")
					sampleKmerFreqM <- getInsertionKmerFreq(.object, samples=samples, k=k, normGenome=TRUE)
				logger.completed()
			} else {
				if (!all(samples %in% colnames(sampleKmerFreqM))) logger.error("'sampleKmerFreqM' does not cover all samples")
			}
			if (is.null(regionKmerFreqM)) {
				logger.start("Computing region kmer frequencies")
					kmerFreqM <- do.call("cbind", lapply(0:(wm-1), FUN=function(i){
						if (i%%50 == 0) logger.status(paste0("i=",i))
						wGr <- trim(resize(GenomicRanges::shift(regionGr, i-ceiling(k/2)), width=k, fix="start", ignore.strand=TRUE))
						rr <- Biostrings::oligonucleotideFrequency(Biostrings::Views(go, wGr), width=k, simplify.as="collapsed")
						return(rr)
					}))
				logger.completed()
			} else {
				kmerFreqM <- regionKmerFreqM
			}
			if (!all(rownames(kmerFreqM)==rownames(sampleKmerFreqM))) logger.error("kmers in frequency matrices do not match")
		}
		# given a RleList object (cov) with genomic coverage (as computed by GenomicRanges::coverage) and a GRanges
		# object (gr) specifying the genomic locations of the features of interests (of uniform length)
		# returns the summed/piled-up coverages for each position across all elements in gr
		# adapted and optimized code from Jeff Granja
		fastFootprint <- function(cov, gr, aggrFun="sum"){
			int <- intersect(names(cov), unique(seqnames(gr)))
			cov <- cov[int]
			gr <- gr[which(as.character(seqnames(gr)) %in% int)]
			suppressWarnings(seqlengths(gr)[int] <- sapply(cov,length))
			gr <- trim(gr)
			w <- as.integer(median(width(gr)))
			gr <- gr[width(gr)==w,] #only select elements that have the same length
			grL <- split(gr, seqnames(gr))
			covM <- do.call("rbind", lapply(names(cov), function(chrom){
				v <- as.matrix(Biostrings::Views(cov[[chrom]], ranges(grL[[chrom]])))
				v[is.na(v)] <- 0 #too handle errors should not be needed
				# revert negative strand regions
				revIdx <- as.character(strand(grL[[chrom]]))=="-"
				if (any(revIdx)) v[revIdx,] <- v[revIdx, w:1]
				return(v)
			}))
			res <- NULL
			if (aggrFun=="sum"){
				res <- colSums(covM)
			} else if (aggrFun=="mean"){
				res <- colMeans(covM)
			} else if (aggrFun=="median"){
				res <- matrixStats::colMedians(covM)
			} else {
				logger.error(c("Unknown aggrFun in fastFootprint:", aggrFun))
			}
			return(res)
		}
		if (!silent) logger.start("Aggregating counts")
			countL <- lapply(samples, FUN=function(sid){
				if (!silent) logger.status(c("Sample:", sid))
				return(fastFootprint(sampleCovg[[sid]], regionGr, aggrFun=countAggrFun))
			})
			names(countL) <- samples
		if (!silent) logger.completed()
		res <- do.call("rbind", lapply(samples, FUN=function(sid){
			cs <- countL[[sid]]
			tn5Bias <- c()
			if (kmerBiasAdj) {
				tn5Bias <- as.vector(t(kmerFreqM) %*% matrix(sampleKmerFreqM[,sid]))
			}
			normFac <- as.numeric(NA)
			normFac.tn5 <- as.numeric(NA)
			if (norm=="tailMean"){
				normFac <- 1/mean(cs[c(1:normTailW, (wm-normTailW+1):wm)], na.rm=TRUE)
				if (kmerBiasAdj) {
					normFac.tn5 <- 1/mean(tn5Bias[c(1:normTailW, (wm-normTailW+1):wm)], na.rm=TRUE)
				}
			}
			sampleDf <- data.frame(
				sampleId=sid,
				pos=1:wm,
				count=cs,
				countNorm=cs*normFac,
				stringsAsFactors=FALSE
			)
			if (kmerBiasAdj){
				sampleDf <- cbind(sampleDf, data.frame(
					countsBiasCor=cs/tn5Bias,
					countNormBiasCor=(cs*normFac)/(tn5Bias*normFac.tn5),
					Tn5bias=tn5Bias,
					Tn5biasNorm=tn5Bias*normFac.tn5,
					stringsAsFactors=FALSE
				))
			}
			return(sampleDf)
		}))
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getMotifEnrichment")) {
	setGeneric(
		"getMotifEnrichment",
		function(.object, ...) standardGeneric("getMotifEnrichment"),
		signature=c(".object")
	)
}
#' getMotifEnrichment-methods
#'
#' Perform enrichment analysis for (TF) motifs of a query set of regions.
#' Fisher's Exact Test is employed to test the association of motif present in the
#' query set against the background of all regions of that type covered in the object
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param type       character string specifying the region type
#' @param idx        logical vector or indices of the same length as \code{length(getCoord(.object))}
#'                   specifies the query set
#' @param motifs     either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'                   that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object)
#' @return a \code{data.frame} summarizing Fisher's Exact Test enrichment statistics for each motif
#' 
#' @rdname getMotifEnrichment-DsATAC-method
#' @docType methods
#' @aliases getMotifEnrichment
#' @aliases getMotifEnrichment,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getMotifEnrichment",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		idx,
		motifs="jaspar"
	) {
		res <- NULL
		#count matrix
		coords <- getCoord(.object, type)
		if (is.logical(idx)){
			if (length(idx)!=length(coords)){
				logger.error(c("Could not get motif enrichment: the supplied index vector is logical, but its length does not match"))
			}
		} else if (is.numeric(idx)){
			idx <- as.integer(idx)
			if (any(idx < 0 | idx > length(coords))){
				logger.error(c("Could not get motif enrichment: the supplied index vector numeric but indices are out of range"))
			}
			#create logical vector
			tmp <- idx
			idx <- rep(FALSE, length(coords))
			idx[tmp] <- TRUE
			rm(tmp)
		} else {
			logger.error(c("Could not get motif enrichment: invalid index vector"))
		}

		# for motifmatchr
		mmObjs <- prepareMotifmatchr(.object@genome, motifs)
		#MAIN
		se <- getCountsSE(.object, type)
		mm <- motifmatchr::matchMotifs(mmObjs[["motifs"]], se, genome=mmObjs[["genome"]])
		regionMotifMatch <- as.matrix(motifmatchr::motifMatches(mm))

		motifNames <- colData(mm)[["name"]]
		names(motifNames) <- NULL
		if (!any(duplicated(motifNames))){
			colnames(regionMotifMatch) <- motifNames
		} else {
			motifNames <- colnames(regionMotifMatch)
		}

		res <- do.call("rbind", lapply(motifNames, FUN=function(mo){
			# hasMotif <- regionMotifMatch[,mo]
			# a <- sum(hasMotif & idx)
			# b <- sum(!hasMotif & idx)
			# c <- sum(hasMotif & !idx)
			# d <- sum(!hasMotif & !idx)
			# fr <- fisher.test(matrix(c(a,c,b,d), nrow=2, ncol=2), alternative="greater")
			if (length(unique(idx)) < 2 || length(unique(regionMotifMatch[,mo])) < 2){
				fr <- list(p.value=as.numeric(NA), estimate=as.numeric(NA))
			} else {
				fr <- fisher.test(x=idx, y=regionMotifMatch[,mo], alternative="greater")
			}
			df <- data.frame(
				pVal=fr$p.value,
				oddsRatio=fr$estimate
			)
			return(df)
		}))
		res[,"qVal"] <- rep(as.numeric(NA), nrow(res))
		if (!any(is.na(res[,"pVal"]))) {
			res[,"qVal"]  <- tryCatch(
				 qvalue::qvalue(res[,"pVal"])$qvalue,
				 error = function(ee) {
				 	if (ee$message == "missing or infinite values in inputs are not allowed"){
				 		return(qvalue::qvalue(res[,"pVal"], lambda=0)$qvalue)
				 	} else {
				 		logger.warning(c("Could not compute qvalues:", ee$message))
				 		return(rep(as.numeric(NA), nrow(res)))
				 	}
				}
			)
		}
		res[,"motif"] <- motifNames
		rownames(res) <- motifNames
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getChromVarDev")) {
	setGeneric(
		"getChromVarDev",
		function(.object, ...) standardGeneric("getChromVarDev"),
		signature=c(".object")
	)
}
#' getChromVarDev-methods
#'
#' Compute chromVar deviations
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param type       character string specifying the region type
#' @param motifs     either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'                   that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object)
#' @return Deviations object as returned by \code{chromVAR::computeDeviations}
#' 
#' @rdname getChromVarDev-DsATAC-method
#' @docType methods
#' @aliases getChromVarDev
#' @aliases getChromVarDev,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getChromVarDev",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		motifs="jaspar"
	) {
		res <- NULL

		countSe <- getCountsSE(.object, type, naIsZero=TRUE)
		genomeObj <- getGenomeObject(.object@genome)
		genome(countSe) <- BSgenome::providerVersion(genomeObj) # hack to override inconsistent naming of genome versions (e.g. hg38 and GRCh38)

		countSe <- chromVAR::addGCBias(countSe, genome=genomeObj)

		# for motifmatchr
		mmInput <- prepareMotifmatchr(genomeObj, motifs)
		mmObj <- motifmatchr::matchMotifs(mmInput[["motifs"]], countSe, genome=genomeObj)
		
		ridx <- safeMatrixStats(assay(countSe), "rowSums") > 0 & safeMatrixStats(assay(mmObj), "rowSums") > 0 # only consider regions where there is an actual motif match and counts
		res <- chromVAR::computeDeviations(object=countSe[ridx,], annotations=mmObj[ridx,])

		return(res)
	}
)

################################################################################
# Footprinting
################################################################################
if (!isGeneric("getMotifFootprints")) {
	setGeneric(
		"getMotifFootprints",
		function(.object, ...) standardGeneric("getMotifFootprints"),
		signature=c(".object")
	)
}
#' getMotifFootprints-methods
#'
#' Perform enrichment analysis for (TF) motif footprinting
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param motifNames character vector of motif names
#' @param samples sample identifiers
#' @param motifFlank number of base pairs flanking the motif on each side
#' @param type       (PLACEHOLDER ARGUMENT: NOT IMPLEMENTED YET) character string specifying the region type or \code{".genome"} (default) for genome-wide profiling
#' @param motifDb    either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'                   that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object) OR a list of \code{GRanges} objects specifying motif occurrences
#' @return a \code{list} of footprinting results with one element for each motif. Each motif's results contain summary data frames with aggregated counts
#'         across all motif occurrences and a \code{ggplot} object for plotting footprints
#' 
#' @rdname getMotifFootprints-DsATAC-method
#' @docType methods
#' @aliases getMotifFootprints
#' @aliases getMotifFootprints,DsATAC-method
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' motifNames <- c("MA1419.1_IRF4", "MA0139.1_CTCF", "MA0037.3_GATA3")
#' # motifNames <- grep("(IRF4|CTCF|GATA3)$", names(prepareMotifmatchr("hg38", "jaspar")$motifs), value=TRUE, ignore.case=TRUE) # alternative by searching
#' samples <- c("TeffNaive_U_1001", "TeffNaive_U_1002", "TeffMem_U_1001", "TeffMem_U_1002")
#' fps <- getMotifFootprints(dsa, motifNames, samples)
#' fps[["MA1419.1_IRF4"]]$plot
#' }
setMethod("getMotifFootprints",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		motifNames,
		samples=getSamples(.object),
		motifFlank=250L,
		type=".genome",
		motifDb="jaspar"
	) {
		universeGr <- NULL
		if (type!=".genome"){
			logger.warning("Footprinting: Using types other than '.genome' is experimental")
			universeGr <- getCoord(.object, type)
			elementMetadata(universeGr) <- NULL
		}
		logger.start("Finding motif occurrences")
			motifObj <- NULL
			motifKmerFreqML <- NULL # region window k-mer frequencies for bias correction
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (!is.null(annoPkg) & is.character(motifDb)){
				logger.info(c("Using annotation from package:", annoPkg))
				motifObj        <- get("getMotifAnnotation", asNamespace(annoPkg))(anno=motifDb, type="motifs")
				motifGrl        <- get("getMotifAnnotation", asNamespace(annoPkg))(anno=motifDb, type="motifOccGrl")
				if (type==".genome"){
					motifKmerFreqML <- get("getMotifAnnotation", asNamespace(annoPkg))(anno=motifDb, type="motifWindowKmerFreq")
				}
			} else if (is.list(motifDb) || any(grepl("GRanges",class(motifDb)))){
				logger.info("Provided motif occurrences as (GRanges) list")
				motifGrl <- motifDb
			} else {
				logger.info("Using motifmatchr")
				if (is.character(motifDb)){
					motifObj <- prepareMotifmatchr(.object@genome, motifDb)$motifs # currently only used for motif logo plotting. Could be omitted if that is not desired
				}
				motifGrl <- getMotifOccurrences(motifNames, motifDb=motifDb, genome=.object@genome)
			}

			motifLens <- sapply(motifNames, FUN=function(mn){
				as.integer(median(width(motifGrl[[mn]])))
			})
			names(motifLens) <- motifNames
		logger.completed()

		obj_filt <- .object
		if (!is.null(universeGr)){
			# remove region aggregates to save time and space (we just need insertions)
			for (rt in getRegionTypes(obj_filt)){
				obj_filt <- removeRegionType(obj_filt, rt)
			}
			logger.start(c("Filtering insertions by region type:", type))
				obj_filt <- filterByGRanges(obj_filt, universeGr, method="white")
			logger.completed()
		}

		logger.start("Computing sample coverage")
			sampleCovg <- getCoverage(obj_filt, samples=samples)
		logger.completed()
		logger.start("Computing sample insertion k-mer frequencies")
			sampleKmerFreqM <- getInsertionKmerFreq(obj_filt, samples=samples, k=6, normGenome=TRUE)
		logger.completed()

		res <- lapply(motifNames, FUN=function(mn){
			logger.start(c("Motif:", mn))
				motifCenGr <- motifGrl[[mn]]
				if (!is.null(universeGr)){
					idx <- overlapsAny(motifCenGr, universeGr, type="within", ignore.strand=TRUE)
					logger.info(paste0("Retaining ", sum(idx), " of ", length(idx), " (", round(100*sum(idx)/length(idx), 2), "%) motif occurrences overlapping with query region type"))
					motifCenGr <- motifCenGr[idx]
				}
				motifCenGr <- unique(trim(resize(motifCenGr, width=2*motifFlank+1, fix="center", ignore.strand=TRUE)))
				regKmerFreqM <- NULL # region window k-mer frequencies for bias correction
				if (length(motifKmerFreqML) > 0){
					regKmerFreqM <- motifKmerFreqML[[mn]]
				}

				footprintDf <- aggregateRegionCounts(obj_filt, motifCenGr, samples=samples, sampleCovg=sampleCovg, sampleKmerFreqM=sampleKmerFreqM, regionKmerFreqM=regKmerFreqM)
				footprintDf$pos <- footprintDf$pos - (motifFlank+1) # offset the position: aggregateRegionCounts returns positions in [0,regionWidth]

				motifLen <- motifLens[mn]
				
				# Definitions according to [Corces, Granja, et al. (2018). The chromatin accessibility landscape of primary human cancers. Science, 362(6413)]
				absPos <- abs(footprintDf[,"pos"])
				i_base <- ceiling(motifLen/2) + 5
				x_base <- footprintDf[absPos <= i_base, "countNormBiasCor"]
				# 10% trimmed mean
				x_base[x_base > quantile(x_base, 0.9) | x_base < quantile(x_base, 0.1)] <- NA
				fp_baseMean  <- mean(x_base, na.rm=TRUE)
				fp_flankMean <- mean(footprintDf[absPos > i_base & absPos <= 50, "countNormBiasCor"], na.rm=TRUE)
				fp_bgMean    <- mean(footprintDf[absPos >= 200, "countNormBiasCor"], na.rm=TRUE)

				# plot
				pp <- plotFootprint(footprintDf, pwm=motifObj[[mn]], colorMap=NULL)
				# motifUp <- -ceiling(motifLen/2) + 1
				# motifDown <- floor(motifLen/2)
				# ppm <- ggplot(footprintDf, aes(x=pos, y=countNormBiasCor, color=sampleId, group=sampleId, fill=sampleId)) + 
				# 	  annotate("rect", xmin=motifUp, xmax=motifDown, ymin=-Inf, ymax=Inf, fill="#d9d9d9") +
				#       geom_line() #+ geom_smooth(aes(group=cellType),size=2,alpha=0.15,method="gam",formula=y ~ s(x, bs = "cs"))

				# ppb <- ggplot(footprintDf, aes(x=pos, y=Tn5biasNorm, color=sampleId, group=sampleId, fill=sampleId)) + 
				# 	  annotate("rect", xmin=motifUp, xmax=motifDown, ymin=-Inf, ymax=Inf, fill="#d9d9d9") +
				#       geom_line()

				# # add sequence logo
				# logo <- hmSeqLogo(motifObj[[mn]], unit(0.5, "npc"), unit(0.5, "npc"), 1, 1, ic.scale=TRUE)
				# logoHeight <- (max(footprintDf[,"countNormBiasCor"], na.rm=TRUE) - min(footprintDf[,"countNormBiasCor"], na.rm=TRUE)) * 0.2
				# logoY <- max(footprintDf[,"countNormBiasCor"], na.rm=TRUE) - logoHeight
				# ppm <- ppm + cowplot::draw_plot(logo, 100, logoY, 100, logoHeight)

				# pp <- cowplot::plot_grid(ppm, ppb + theme(legend.position="none"), ncol=1, rel_heights=c(2, 1), align="v", axis="lr")

				rr <- list(
					footprintDf=footprintDf,
					plot=pp,
					motifLen=motifLen,
					nMotifs=length(motifCenGr),
					fpFlankAcc=log2(fp_flankMean/fp_bgMean),
					fpDepth=log2(fp_baseMean/fp_flankMean)
				)
			logger.completed()
			return(rr)
		})
		names(res) <- motifNames
		
		return(res)
	}
)


################################################################################
# Differential analysis
################################################################################
if (!isGeneric("getDESeq2Dataset")) {
	setGeneric(
		"getDESeq2Dataset",
		function(.object, ...) standardGeneric("getDESeq2Dataset"),
		signature=c(".object")
	)
}
#' getDESeq2Dataset-methods
#'
#' Retrieve a differential expression dataset computed with DESeq2
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionType character string specifying the region type
#' @param designCols column names in the sample annotation potentially used to create the design matrix
#' @param compTab    if design columns are not specified, you can specify a comparison table directly. These comparison tables can be obtained by \code{getComparisonTable(...)}
#' @param ...        parameters passed on to \code{DESeq2::DESeq}
#' @return \code{DESeqDataSet} as returned by \code{DESeq2::DESeq}
#' 
#' @rdname getDESeq2Dataset-DsATAC-method
#' @docType methods
#' @aliases getDESeq2Dataset
#' @aliases getDESeq2Dataset,DsATAC-method
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' dds <- getDESeq2Dataset(dsa, "IA_prog_peaks", designCols=c("donor", "stimulus", "cellType"))
#' dds
#' }
setMethod("getDESeq2Dataset",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionType,
		designCols=NULL,
		compTab=NULL,
		...
	) {
		if (!requireNamespace("DESeq2")) logger.error(c("Could not load dependency: DESeq2"))
		if (!is.element(regionType, getRegionTypes(.object))) logger.error(c("Unsupported region type:", regionType))
		# annotation data
		ph <- getSampleAnnot(.object)

		if (!is.null(compTab) && nrow(compTab) > 0){
			vsAllIdx <- compTab[,"grp1Name"]==".ALL" | compTab[,"grp2Name"]==".ALL"
			vsAllCols <- unique(compTab[vsAllIdx, "compCol"])
			designCols <- union(designCols, unique(compTab[!vsAllIdx,"compCol"]))
			if (length(vsAllCols) > 0){
				logger.status(c("Adding 1-vs-all comparisons"))
				for (i in which(vsAllIdx)){
					gn <- setdiff(c(compTab[i, "grp1Name"], compTab[i, "grp2Name"]), ".ALL")
					if (length(gn)!=1) logger.error(c("Invalid comparison:", compTab[i,"compName"]))
					cn <- compTab[i,"compCol"]
					addCn <- paste0(".", cn, ".", gn, ".vsAll")
					ph[, addCn] <- ifelse(as.character(ph[,cn])==gn, gn, ".ALL")
					designCols <- union(designCols, addCn)
				}
			}
		}

		if (!all(designCols %in% colnames(ph))){
			logger.error(c("Invalid design columns. The following column names are not contained in the sample annotation:", paste(setdiff(designCols, colnames(ph)), collapse=",")))
		}
		#count matrix
		cm <- ChrAccR::getCounts(.object, regionType, asMatrix=TRUE)
		
		#remove columns from the design that do not have multiple levels
		idx <- sapply(designCols, FUN=function(cc){
			length(unique(ph[,cc]))>1
		})
		if (sum(idx) < length(designCols)){
			logger.warning(c("The following design columns will not be considered because they do not have multiple levels:", paste(designCols[!idx], collapse=",")))
			designCols <- designCols[idx]
		}
		#remove columns from the design that do not have replicates
		idx <- sapply(designCols, FUN=function(cc){
			is.numeric(ph[,cc]) || all(table(ph[,cc]) > 1)
		})
		if (sum(idx) < length(designCols)){
			logger.warning(c("The following design columns will not be considered because they do not have replicates:", paste(designCols[!idx], collapse=",")))
			designCols <- designCols[idx]
		}
		designF <- as.formula(paste0("~", paste(designCols,collapse="+")))
		
		dds <- DESeq2::DESeqDataSetFromMatrix(countData=cm, colData=ph, design=designF)
		rowRanges(dds) <- getCoord(.object, regionType)
		dds <- DESeq2::DESeq(dds, ...)

		return(dds)
	}
)

if (!isGeneric("getDiffAcc")) {
	setGeneric(
		"getDiffAcc",
		function(.object, ...) standardGeneric("getDiffAcc"),
		signature=c(".object")
	)
}
#' getDiffAcc-methods
#'
#' Compute differential accessibility
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionType character string specifying the region type
#' @param comparisonCol column name in the sample annotation table to base the comparison on
#' @param grp1Name   name of the first group in the comparison. if not specified, it will be taken as the first factor level specified in the 
#'                   sample annotation table in \code{'comparisonCol'}.
#' @param grp2Name   name of the second group (reference) in the comparison. if not specified, it will be taken as the first factor level specified in the 
#'                   sample annotation table in \code{'comparisonCol'}.
#' @param adjustCols column names in the sample annotation potentially used to create the design matrix
#' @param method     Method for determining differential accessibility. Currently only \code{'DESeq2'} is supported
#' @param diffObj    optional differential object to avoid computing it for each comparison and thus reduce runtime
#' @return a \code{data.frame} containing differential accessibility information
#' 
#' @rdname getDiffAcc-DsATAC-method
#' @docType methods
#' @aliases getDiffAcc
#' @aliases getDiffAcc,DsATAC-method
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' daTab <- getDiffAcc(dsa, "IA_prog_peaks", "stimulus", grp1Name="S", grp2Name="U", adjustCols=c("cellType", "donor"))
#' }
setMethod("getDiffAcc",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionType,
		comparisonCol,
		grp1Name=NULL,
		grp2Name=NULL,
		adjustCols=character(0),
		method='DESeq2',
		diffObj=NULL
	) {
		if (!is.element(method, c("DESeq2"))) logger.error(c("Invalid method for calling differential accessibility:", method))
		if (!is.element(comparisonCol, colnames(getSampleAnnot(.object)))) logger.error(c("Comparison column not found in sample annotation:", comparisonCol))
		comparisonCol0 <- comparisonCol
		contrastF <- factor(getSampleAnnot(.object)[,comparisonCol])
		if (length(levels(contrastF)) < 2)  logger.error(c("Invalid comparison column. There should be at least 2 groups."))

		if (is.null(grp1Name)) grp1Name <- levels(contrastF)[1]
		if (is.null(grp2Name)) grp2Name <- levels(contrastF)[2]
		if (!is.element(grp1Name, c(levels(contrastF), ".ALL"))) logger.error(c("Invalid group name (1). Sample annotation has no samples associated with that group:", grp1Name))
		if (!is.element(grp2Name, c(levels(contrastF), ".ALL"))) logger.error(c("Invalid group name (2). Sample annotation has no samples associated with that group:", grp2Name))
		sidx.grp1 <- which(contrastF==grp1Name)
		if (grp1Name==".ALL") sidx.grp1 <- which(contrastF!=grp2Name)
		sidx.grp2 <- which(contrastF==grp2Name)
		if (grp2Name==".ALL") sidx.grp2 <- which(contrastF!=grp1Name)


		if (method=="DESeq2"){
			is1vsAll <- grp1Name==".ALL" || grp2Name==".ALL"
			if (is1vsAll){
				if (grp1Name==grp2Name) logger.error(".ALL vs .ALL is not legal")
				gn <- setdiff(c(grp1Name, grp2Name), ".ALL")
				comparisonCol <- paste0(".", comparisonCol, ".", gn, ".vsAll")
			}

			if (!is.null(diffObj)){
				if (!is.element("DESeqDataSet", class(diffObj))) logger.error(c("Invalid 'diffObj' for method", method))
				#check if the elements match
				gr <- getCoord(.object, regionType)
				if (length(gr)!=nrow(diffObj)) logger.error("dimensions of the 'diffObj' do not match.")
				if (!all(colnames(diffObj)==getSamples(.object))) logger.error("Invalid 'diffObj'. Sample names must match the DsATAC object.")
				dds <- diffObj
			} else {
				compTab <- NULL
				designCs <- adjustCols
				if (is1vsAll){
					compTab <- getComparisonTable(.object, compNames=paste0(grp1Name, " vs ", grp2Name,  " [", comparisonCol0, "]"))
				} else {
					designCs <- unique(c(designCs, comparisonCol))
				}
				dds <- getDESeq2Dataset(.object, regionType, designCs, compTab=compTab)				
			}

			diffRes <- DESeq2::results(dds, contrast=c(comparisonCol, grp1Name, grp2Name))
			dm <- data.frame(diffRes)
			rankMat <- cbind(
				# rank(-dm[,"baseMean"]), na.last="keep", ties.method="min"),
				rank(-abs(dm[,"log2FoldChange"]), na.last="keep", ties.method="min"),
				rank(dm[,"pvalue"], na.last="keep", ties.method="min")
			)
			dm[,"cRank"] <- matrixStats::rowMaxs(rankMat, na.rm=FALSE)
			# dm[,"cRank"] <- rowMaxs(rankMat, na.rm=TRUE)
			dm[!is.finite(dm[,"cRank"]),"cRank"] <- NA
			dm[,"cRank_rerank"] <- rank(dm[,"cRank"], na.last="keep", ties.method="min")

			l10fpkm <- log10(DESeq2::fpkm(dds, robust=TRUE)+1)
			grp1.m.l10fpkm <- rowMeans(l10fpkm[, sidx.grp1, drop=FALSE], na.rm=TRUE)
			grp2.m.l10fpkm <- rowMeans(l10fpkm[, sidx.grp2, drop=FALSE], na.rm=TRUE)
			vstCounts <- assay(DESeq2::vst(dds, blind=FALSE))
			grp1.m.vst <- rowMeans(vstCounts[, sidx.grp1, drop=FALSE], na.rm=TRUE)
			grp2.m.vst <- rowMeans(vstCounts[, sidx.grp2, drop=FALSE], na.rm=TRUE)

			res <- data.frame(
				log2BaseMean=log2(dm[,"baseMean"]),
				meanLog10FpkmGrp1=grp1.m.l10fpkm,
				meanLog10FpkmGrp2=grp2.m.l10fpkm,
				meanVstCountGrp1=grp1.m.vst,
				meanVstCountGrp2=grp2.m.vst,
				dm
			)
			# add group names to column names
			for (cn in c("meanLog10FpkmGrp", "meanVstCountGrp")){
				colnames(res)[colnames(res)==paste0(cn,"1")] <- paste0(cn, "1_", grp1Name)
				colnames(res)[colnames(res)==paste0(cn,"2")] <- paste0(cn, "2_", grp2Name)
			}
		}
		return(res)
	}
)
################################################################################
# Export
################################################################################
if (!isGeneric("exportCountTracks")) {
	setGeneric(
		"exportCountTracks",
		function(.object, ...) standardGeneric("exportCountTracks"),
		signature=c(".object")
	)
}
#' exportCountTracks-methods
#'
#' export count data as genome tracks (e.g. for visualization in the browser)
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param type       character string specifying the region type
#' @param outDir     output directory. Must be existing.
#' @param formats    browser format. Currently only bed and "igv" are supported
#' @param groupBy    a column in the sample annotation table to group by (the mean will be computed)
#' @return nothing of particular interest
#' 
#' @rdname exportCountTracks-DsATAC-method
#' @docType methods
#' @aliases exportCountTracks
#' @aliases exportCountTracks,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("exportCountTracks",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		outDir,
		formats=c("bed", "igv"),
		groupBy=NULL
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		if (!dir.exists(outDir)) logger.error(c("Output directory:", outDir, "does not exist."))
		#count matrix
		cm <- ChrAccR::getCounts(.object, type, asMatrix=TRUE)
		sampleNames <- getSamples(.object)
		coords <- getCoord(.object, type)
		if (!is.null(groupBy) && is.character(groupBy) && is.element(groupBy, colnames(getSampleAnnot(.object)))){
			grps <- factor(getSampleAnnot(.object)[,groupBy])
			cm <- do.call("cbind", tapply(getSamples(.object), grps, FUN=function(sids){
				rowMeans(cm[,sids, drop=FALSE], na.rm=TRUE)
			}, simplify=FALSE))
			colnames(cm) <- levels(grps)
			sampleNames <- levels(grps)
		}

		elementMetadata(coords) <- NULL
		elementMetadata(coords) <- cm

		if (is.element("bed", formats)){
			for (sn in sampleNames){
				fn <- file.path(outDir, paste0(sn, "_counts.bed"))
				granges2bed(coords, fn, score=elementMetadata(coords)[,sn], addAnnotCols=FALSE, colNames=FALSE, doSort=TRUE, bigBed=FALSE)
			}
		}
		if (is.element("igv", formats)){
			fn <- file.path(outDir, "DsATAC_counts.igv")
			granges2igv(coords, fn, addStrand=FALSE, addAnnotCols=TRUE, doSort=TRUE, toTDF=TRUE)
		}
		#TODO: work in progress
		invisible(NULL)
	}
)

################################################################################
# Data processing / inference
################################################################################
if (!isGeneric("callPeaks")) {
	setGeneric(
		"callPeaks",
		function(.object, ...) standardGeneric("callPeaks"),
		signature=c(".object")
	)
}
#' callPeaks-methods
#'
#' Performs peak calling based on insertion sites
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param samples sample identifiers for which peak calling is performed
#' @param method  peak calling method. Currently only \code{'macs2_summit_fw_no'} is supported. See details section.
#' @param methodOpts list of other options depending on the \code{'method'} parameter (see details section).
#' @return \code{GRangesList} of peak coordinates for each sample
#' 
#' @details
#' The following methods are currently supported
#' \describe{
#'    \item{\code{'macs2_summit_fw_no'}}{
#'		Fixed-width, non-overlapping peaks based on MACS2 summit calls: 
#'      1. Call peaks using system call to MACS2. You can specify the MACS2 executable in \code{methodOpts$macs2.exec}.
#' 		2. Identify peak summits
#' 		3. extend peak summits on each side by a number of basepairs (specified in \code{methodOpts$fixedWidth}; default: 250bp) to obtain unified peak widths
#' 		4. Find non-overlapping peaks by taking the peak with the best MACS2 score from each set of partially overlapping peaks
#'    }
#' }
#'
#' @rdname callPeaks-DsATAC-method
#' @docType methods
#' @aliases callPeaks
#' @aliases callPeaks,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("callPeaks",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		samples=getSamples(.object),
		method='macs2_summit_fw_no',
		methodOpts=list(
			macs2.exec="macs2",
			macs2.params=c(
				"--shift", "-75",
				"--extsize", "150",
				"-p", "0.01"
			),
			fixedWidth=250,
			genomeSizesFromObject=FALSE
		)
	) {
		if (!is.element(method, c("macs2_summit_fw_no"))) logger.error(c("Invalid 'method':", method))
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		if (!all(samples %in% names(.object@fragments))) logger.error(c("Object does not contain insertion information for samples:", paste(setdiff(samples, names(.object@fragments)), collapse=", ")))
		peakGrl <- NULL
		if (method=="macs2_summit_fw_no"){
			if (!is.element("macs2.exec", names(methodOpts))) logger.error("Invalid 'methodOps' for method 'macs2_summit_fw_no' (missing 'macs2.exec')")
			if (!is.element("macs2.params", names(methodOpts))) logger.error("Invalid 'methodOps' for method 'macs2_summit_fw_no' (missing 'macs2.params')")
			if (!is.element("fixedWidth", names(methodOpts))) logger.error("Invalid 'methodOps' for method 'macs2_summit_fw_no' (missing 'fixedWidth')")
			if (!is.element("genomeSizesFromObject", names(methodOpts))) logger.error("Invalid 'methodOps' for method 'macs2_summit_fw_no' (missing 'genomeSizesFromObject')")
			argV <- c(
				"--nomodel",
				"--call-summits",
				"--nolambda",
				"--keep-dup", "all",
				# "-B",  "--SPMR",
				methodOpts$macs2.params
			)
			genomeSizeArg <- ""
			if (methodOpts$genomeSizesFromObject){
				# # find the set of the largest tiling regions you can find in the object
				# rts <- gtools::mixedsort(grep("^(tiling|t[0-9]+)", getRegionTypes(.object), value=TRUE), decreasing=TRUE)
				# if (length(rts) < 1 || is.na(rts[1])) logger.error("No appropriate tiling region set found for determining genome size")
				# ggr <- getCoord(.object, rts[1])
				# genomeSizeArg <- as.character(sum(width(ggr))) 
				# logger.info(paste0("Using a genome size of ", genomeSizeArg, " bp (determined from region type: ", rts[1], ")"))

				fgr <- do.call("c", unname(getFragmentGrl(.object, samples, asGRangesList=FALSE)))
				fgr <- reduce(fgr, min.gapwidth=9L, ignore.strand=TRUE) #adjust gapwidth to account for Tn5 offsetting
				genomeSizeArg <- as.character(sum(width(fgr)))
				logger.info(paste0("Using a genome size of ", genomeSizeArg, " bp (determined from aggregating fragments)"))
			} else if (is.element(.object@genome, c("hg19", "hg38"))){
				genomeSizeArg <- "hs"
			} else if (is.element(.object@genome, c("mm9", "mm10"))){
				genomeSizeArg <- "mm"
			} else {
				logger.error(c("Unsupported genome for peak calling:", .object@genome))
			}

			nSamples <- length(samples)
			samplePrefixes <- sapply(samples, FUN=function(x){getHashString(pattern=x)})
			genomeAss <- .object@genome

			callDir <- tempdir()

			cmdr <- getConfigElement("muPipeR_cmdr")
			doCmdR <- !is.null(cmdr) && inherits(cmdr, "CommandR")
			if (doCmdR){
				callDir <- getConfigElement("tmpDir")
			}

			logger.start("Preparing insertion bed files")
				insFns <- sapply(seq_along(samples), FUN=function(i){
					sid <- samples[i]
					logger.status(c("Sample:", sid, paste0("(", i, " of ", nSamples,")")))
					fp <- samplePrefixes[i]
					insFn <- file.path(callDir, paste0(fp, "_ins.bed"))
					# logger.status(c("[DEBUG:] Retrieving insertion sites..."))
					insGr <- getInsertionSites(.object, sid)[[1]]
					# logger.status(c("[DEBUG:] Writing to temp file..."))
					coordOnly <- all(as.character(strand(insGr)) %in% c("*", "."))
					granges2bed(insGr, insFn, score=NULL, addAnnotCols=FALSE, colNames=FALSE, doSort=TRUE, coordOnly=coordOnly)
					return(insFn)
				})
				names(insFns) <- samples
			logger.completed()

			# variable environment (for parallel computation using lapplyExec)
			envList <- list(
				samples=samples,
				insFns=insFns,
				callDir=callDir,
				samplePrefixes=samplePrefixes,
				genomeAss=genomeAss,
				methodOpts=methodOpts,
				argV=argV,
				genomeSizeArg=genomeSizeArg
			)
			# function to be used for peak calling
			callPeakFun <- function(i){
				sid <- samples[i]
				logger.status(c("Calling peaks for sample:", sid, paste0("(", i, " of ", length(samples),")")))
				
				# peakFn <- file.path(callDir, paste0(fp, "_summits.bed"))
				peakFn <- file.path(callDir, paste0(samplePrefixes[i], "_peaks.narrowPeak"))

				# logger.status(c("[DEBUG:] Calling MACS2..."))
				aa <- c(
					"callpeak",
					"-g", genomeSizeArg,
					"--name", samplePrefixes[i],
					"--treatment", insFns[sid],
					"--outdir", callDir,
					"--format", "BED",
					argV
				)
				system2(methodOpts$macs2.exec, aa, wait=TRUE, stdout="", stderr="") # MACS2 log messages to console
				# system2(methodOpts$macs2.exec, aa, wait=TRUE, stdout=FALSE, stderr=FALSE) # suppress MACS2 log messages

				# logger.status(c("[DEBUG:] Reading MACS2 output..."))
				# peakGr <- rtracklayer::import(peakFn, format="BED")
				peakGr <- readMACS2peakFile(peakFn)
				peakGr <- setGenomeProps(peakGr, genomeAss, onlyMainChrs=TRUE, silent=TRUE)
				peakGr <- peakGr[isCanonicalChrom(as.character(seqnames(peakGr)))]
				elementMetadata(peakGr)[,"calledPeakStart"] <- start(peakGr)
				elementMetadata(peakGr)[,"calledPeakEnd"] <- end(peakGr)
				start(peakGr) <- end(peakGr) <- elementMetadata(peakGr)[,"summit"]
				# scale scores to their percentiles
				# scs <- elementMetadata(peakGr)[,"score"]
				scs <- elementMetadata(peakGr)[,"negLog10qval"]
				elementMetadata(peakGr)[,"score_norm"] <- ecdf(scs)(scs)
				elementMetadata(peakGr)[,"name"] <- gsub(paste0("^", samplePrefixes[i]), sid, elementMetadata(peakGr)[,"name"])#replace the hashstring in the name by just the sample id

				# logger.status(c("[DEBUG:] Extending summits..."))
				peakGr <- trim(promoters(peakGr, upstream=methodOpts$fixedWidth, downstream=methodOpts$fixedWidth+1)) #extend each summit
				peakGr <- peakGr[width(peakGr)==median(width(peakGr))] #remove too short regions which might have been trimmed
				# logger.status(c("[DEBUG:] Finding non-overlapping peaks..."))
				peakGr <- getNonOverlappingByScore(peakGr, scoreCol="score_norm")
				# peakGr <- ChrAccR:::getNonOverlappingByScore(peakGr, scoreCol="score_norm")
				peakGr <- peakGr[order(as.integer(seqnames(peakGr)),start(peakGr), end(peakGr), as.integer(strand(peakGr)))] #sort
				logger.info(c("Identified", length(peakGr), "peaks"))
				return(peakGr)
			}
			logger.start("Calling peaks using MACS2")
				if (doCmdR){
					peakGrl <- lapplyExec(
						cmdr,
						lapply(seq_along(samples), identity),
						callPeakFun,
						env=envList,
						Rexec="Rscript",
						name=paste0("chraccr_callPeaks")
					)
				} else {
					peakGrl <- lapply(seq_along(samples), callPeakFun)
				}
				names(peakGrl) <- samples

				# cleanup
				fns <- file.path(callDir, paste0(samplePrefixes, "*"))
				unlink(fns)
			logger.completed()
			peakGrl <- GRangesList(peakGrl)
		}
		return(peakGrl)
	}
)


################################################################################
# Plotting
################################################################################
if (!isGeneric("plotInsertSizeDistribution")) {
	setGeneric(
		"plotInsertSizeDistribution",
		function(.object, ...) standardGeneric("plotInsertSizeDistribution"),
		signature=c(".object")
	)
}
#' plotInsertSizeDistribution-methods
#'
#' Plot insert size distribution
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param sampleId   sample to be plotted
#' @return \code{ggplot} object containing insert size distribution plot
#' 
#' @rdname plotInsertSizeDistribution-DsATAC-method
#' @docType methods
#' @aliases plotInsertSizeDistribution
#' @aliases plotInsertSizeDistribution,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("plotInsertSizeDistribution",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		sampleId
	) {
		insSizeDf <- data.frame(
			insertionSize=width(getFragmentGr(.object, sampleId))
		)

		pp <- ggplot(insSizeDf) + aes(x=insertionSize, y=..count..) + geom_density(alpha=.2, fill="#8c1515") +
		      xlim(c(0,700)) + ylab("#fragments")
		return(pp)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getTssEnrichment")) {
	setGeneric(
		"getTssEnrichment",
		function(.object, ...) standardGeneric("getTssEnrichment"),
		signature=c(".object")
	)
}
#' getTssEnrichment-methods
#'
#' Get TSS enrichment data and plot
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param sampleId   sample to be plotted
#' @param tssGr      \code{GRanges} object containing TSS coordinates or NULL to get default set from annotation package
#' @param flank      number of bases flanking each TSS that will be added on each side
#' @param normTailW  number of bases on each side whose counts will be used to normalize the data
#' @param smoothW    radius of the window (in bp) that will be used to smooth the data, i.e. the total width of the
#'                   smoothing window will be twice that number
#' @param silent     limit log messages
#' @return a list containing TSS enrichment data and a \code{ggplot} object containing TSS enrichment plot
#' 
#' @rdname getTssEnrichment-DsATAC-method
#' @docType methods
#' @aliases getTssEnrichment
#' @aliases getTssEnrichment,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getTssEnrichment",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		sampleId,
		tssGr=NULL,
		flank=2000L,
		normTailW=100L,
		smoothW=25L,
		silent=FALSE
	) {
		if (is.null(tssGr)){
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (is.null(annoPkg)) logger.error("Annotation package needed")
			tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
		}
		if (!all(width(tssGr)==1)) logger.error("tssGr must be a GRanges object in which each element has width=1")

		#extend the window by the flanking and smoothing lengths
		tssGr <- unique(trim(resize(tssGr, width=2*(flank+smoothW)+1, fix="center", ignore.strand=TRUE)))
		# get (normalized) count data
		tssCountDf <- aggregateRegionCounts(.object, tssGr, samples=sampleId, countAggrFun="mean", kmerBiasAdj=FALSE, normTailW=normTailW+smoothW, silent=silent)
		# offset the position: aggregateRegionCounts returns positions in [0,regionWidth]
		tssCountDf$pos <- tssCountDf$pos - (flank+smoothW+1)

		# TSS enrichment: max normalized count within window 
		idx.main <- abs(tssCountDf[,"pos"])<=(flank-normTailW)
		tsse <- max(tssCountDf[idx.main,"countNorm"], na.rm=TRUE)
		tsse.s <- as.numeric(NA)

		countCol <- "countNorm"
		if (smoothW > 1){
			# smoothing: take the mean normalized count in the corresponding window
			smoothedCounts <- convolve(tssCountDf[,"countNorm"], rep(1, 2*smoothW+1), type="filter")/(smoothW*2+1)
			tsse.s <- max(smoothedCounts[idx.main], na.rm=TRUE)
			tssCountDf <- tssCountDf[abs(tssCountDf[,"pos"])<=flank,] # revert the extension by the smoothing window
			tssCountDf[,"countNormSmoothed"] <- smoothedCounts
			countCol <- "countNormSmoothed"
		}

		pp <- ggplot(tssCountDf) + aes_string(x="pos", y="countNorm") + geom_vline(xintercept=c(0), color="#4d4f53") + geom_point(color="#969696") +
			  geom_line(aes_string(x="pos", y=countCol), color="#8c1515", size=2)

		res <- list(
			countDf=tssCountDf,
			tssEnrichment=tsse,
			tssEnrichment.smoothed=tsse.s,
			plot=pp
		)
		return(res)
	}
)

################################################################################
# Single-cell/ high-sample number analyses
################################################################################
if (!isGeneric("getTssEnrichmentBatch")) {
	setGeneric(
		"getTssEnrichmentBatch",
		function(.object, ...) standardGeneric("getTssEnrichmentBatch"),
		signature=c(".object")
	)
}
#' getTssEnrichmentBatch-methods
#'
#' Get TSS enrichment data and plot
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param tssGr      \code{GRanges} object containing TSS coordinates or NULL to get default set from annotation package
#' @param sampleIds  sampleIds for which TSS enrichment should be computed
#' @param tssW       size of the core TSS window
#' @param flank      number of bases flanking each TSS that will be added on each side
#' @param normTailW  number of bases on each side whose counts will be used to normalize the data
#' @param smoothW    diameter of the window (in bp) that will be used to smooth the data
#' @return a list containing TSS enrichment data
#' 
#' @rdname getTssEnrichmentBatch-DsATAC-method
#' @docType methods
#' @aliases getTssEnrichmentBatch
#' @aliases getTssEnrichmentBatch,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getTssEnrichmentBatch",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		tssGr=NULL,
		sampleIds=getSamples(.object),
		tssW=201L,
		flank=2000L,
		normTailW=200L,
		smoothW=51L
	) {
		if (is.null(tssGr)){
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (is.null(annoPkg)) logger.error("Annotation package needed")
			tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
		}
		if (!all(width(tssGr)==1)) logger.error("tssGr must be a GRanges object in which each element has width=1")

		#extend the window by the flanking lengths
		L <- 2*flank+1
		tsswGr <- suppressWarnings(trim(resize(tssGr, width=L, fix="center", ignore.strand=TRUE)))

		nRegs <- length(tssGr)
		if (length(tsswGr) != nRegs) logger.error("length error: tss window set")

		nSamples <- length(sampleIds)
		logger.status("Retrieving fragment data ...")
		fGr <- getFragmentGrl(.object, sampleIds, asGRangesList=FALSE)
		# logger.status("[DEBUG] ...(1)...")
		nFrags <- sapply(fGr, length)
		# subsample fragments if they exceed the threshold (to avoid long computation and long vector errors)
		# doSubsample <- .Machine$integer.max/2-1
		doSubsample <- 1e9
		if (sum(nFrags) > doSubsample){
			prop <- doSubsample/sum(nFrags)
			logger.warning(c("Subsampling to", round(prop*100, 2), "% of fragments"))
			nFrags <- floor(nFrags*prop)
			fGr <- lapply(1:length(fGr), FUN=function(i){
				fGr[[i]][sort(sample(1:length(fGr[[i]]), nFrags[i]))]
			})
		}
		# logger.status("[DEBUG] ...(2)...")
		fGr <- do.call("c", unname(fGr))
		# logger.status("[DEBUG] ...(3)...")
		# fGr <- unlist(fGr, use.names=FALSE)
		elementMetadata(fGr)[,".sampleIdx"] <- rep(seq_along(nFrags), nFrags)
		logger.status("Retrieving joined insertion sites ...")
		igr <- getInsertionSitesFromFragmentGr(fGr)
		rm(fGr) # clean-up
		cleanMem() # TODO: can this be removed?
		logger.status("computing overlaps with TSS regions ...")
		oo <- findOverlaps(tsswGr, igr, ignore.strand=TRUE)
		tssGro <- tssGr[queryHits(oo)]
		iGro <- igr[subjectHits(oo)]
		rm(igr) # clean-up
		cleanMem() # TODO: can this be removed?
		# dd <- -grSignedDistance(tssGro, iGro)
		dd <- start(iGro) - start(tssGro)
		negIdx <- as.logical(strand(tssGro)=="-")
		dd[negIdx] <- -dd[negIdx]
		sampleIdx <- elementMetadata(iGro)[,".sampleIdx"]

		# tssCut <- floor(tssW/2)
		# nTss <- tapply(seq_along(iGro), sampleIdx, FUN=function(idx){
		# 	idx <- idx[abs(dd[idx]) <= tssCut]
		# 	length(unique(queryHits(oo)[idx]))
		# })
		# bgCut <- flank - normTailW

		logger.status("computing TSS profiles ...")
		countV <- rep(0L, 2*flank+1)
		names(countV) <- as.character(-flank:flank)
		profileMat <- do.call("rbind", tapply(dd, sampleIdx, FUN=function(x){
			posCounts <- table(x)
			res <- countV
			res[names(posCounts)] <- posCounts
			res
		}))
		if (!all(rownames(profileMat)==as.character(1:nSamples))) {
			logger.warning("could not retrieve profiles for all samples. NAs will be returned")
			missingIdxc <- setdiff(as.character(1:nSamples), rownames(profileMat))
			missM <- matrix(as.integer(NA), nrow=length(missingIdxc), ncol=ncol(profileMat))
			rownames(missM) <- missingIdxc
			colnames(missM) <- colnames(profileMat)
			profileMat <- rbind(profileMat, missM)
			profileMat <- profileMat[as.character(1:nSamples),]
		}
		if (ncol(profileMat)!=L) logger.error("Computing profile matrix failed: incorrect dimensions")
		rownames(profileMat) <- sampleIds
		colnames(profileMat) <- NULL
		# profileMat <- profileMat + 1 # pseudo-counts on aggregate matrix
		bgIdx <- c(1:normTailW, (L-normTailW+1):L) # indices of the regions to use as background
		tssWidx <- (flank-floor(tssW/2)):(flank+floor(tssW/2))+1 # indices corresponding to window around TSS

		logger.status("normalizing TSS profiles ...")
		bgMeans <- rowMeans(profileMat[,bgIdx], na.rm=TRUE)
		normFactors <- pmax(bgMeans, 0.5) # low count samples/cells (no effect if pseudocounts are added)

		profileMatNorm <- profileMat/normFactors # divide rows by row normalization factor

		profileMatSmoothed <- t(apply(profileMatNorm, 1, FUN=function(x){
			zoo::rollmean(x, smoothW, fill=1)
			# convolve(c(fillV, x, fillV), rep(1, smoothW), type="filter")/smoothW
		}))

		profileMatNorm[!is.finite(profileMatNorm)] <- NA
		profileMatSmoothed[!is.finite(profileMatSmoothed)] <- NA
		tsse <- matrixStats::rowMaxs(profileMatNorm[,tssWidx], na.rm=TRUE)
		tsse[!is.finite(tsse)] <- NA
		tsse.s <- matrixStats::rowMaxs(profileMatSmoothed[,tssWidx], na.rm=TRUE)
		tsse.s[!is.finite(tsse.s)] <- NA

		res <- list(
			profileMatNorm=t(profileMatNorm),
			profileMatSmoothed=t(profileMatSmoothed),
			pos=(-flank):flank,
			tssEnrichment=tsse,
			tssEnrichment.smoothed=tsse.s
		)
		return(res)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("getQuickTssEnrichment")) {
	setGeneric(
		"getQuickTssEnrichment",
		function(.object, ...) standardGeneric("getQuickTssEnrichment"),
		signature=c(".object")
	)
}
#' getQuickTssEnrichment-methods
#'
#' [Experimental] Quick, heuristic version of TSS enrichment to just get scores for each
#' sample in the dataset. Useful, e.g. for single cells.
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param tssGr      \code{GRanges} object containing TSS coordinates
#' @param sampleIds  sampleIds for which TSS enrichment should be computed
#' @param tssW		 width to consider arount the TSS
#' @param distBg     number of bases flanking each TSS that will be added on each side
#' @return a vector of TSS enrichment values for each sample/cell in the dataset
#' 
#' @details
#' Computes TSS enrichment
#' as the ratio of total insertion sites at a window (of width \code{tssW} bp) directly at the TSS and 2 background regions
#' symmetrically located (\code{distBg} bp) upstream and downstream of the TSS
#' 
#' @rdname getQuickTssEnrichment-DsATAC-method
#' @docType methods
#' @aliases getQuickTssEnrichment
#' @aliases getQuickTssEnrichment,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getQuickTssEnrichment",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		tssGr=NULL,
		sampleIds=getSamples(.object),
		tssW=201L,
		distBg=1900L
	) {
		if (is.null(tssGr)){
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (is.null(annoPkg)) logger.error("Annotation package needed")
			tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
		}
		if (!all(width(tssGr)==1)) logger.error("tssGr must be a GRanges object in which each element has width=1")

		# background windows
		# bgLeftGr <- suppressWarnings(trim(resize(GenomicRanges::shift(tssGr, -distBg), width=floor(tssW/2), fix="center", ignore.strand=TRUE)))
		# bgRightGr <- suppressWarnings(trim(resize(GenomicRanges::shift(tssGr, distBg), width=floor(tssW/2), fix="center", ignore.strand=TRUE)))
		bgLeftGr <- suppressWarnings(trim(resize(GenomicRanges::shift(tssGr, -distBg), width=tssW, fix="center", ignore.strand=TRUE)))
		bgRightGr <- suppressWarnings(trim(resize(GenomicRanges::shift(tssGr, distBg), width=tssW, fix="center", ignore.strand=TRUE)))
		#extend the window by the flanking and smoothing lengths
		tsswGr <- suppressWarnings(trim(resize(tssGr, width=tssW, fix="center", ignore.strand=TRUE)))

		nRegs <- length(tssGr)
		if (length(tsswGr) != nRegs)    logger.error("length error: tss window set")
		if (length(bgLeftGr) != nRegs)  logger.error("length error: left background set")
		if (length(bgRightGr) != nRegs) logger.error("length error: right background set")

		nSamples <- length(sampleIds)
		logger.status("Retrieving joined insertion sites ...")
		fGr <- getFragmentGrl(.object, sampleIds, asGRangesList=FALSE)
		nFrags <- sapply(fGr, length)
		fGr <- do.call("c", unname(fGr))
		# fGr <- unlist(fGr, use.names=FALSE)
		elementMetadata(fGr)[,".sampleIdx"] <- rep(seq_along(nFrags), nFrags)
		igr <- getInsertionSitesFromFragmentGr(fGr)
		logger.status("computing overlaps with TSS regions ...")
		oo <- findOverlaps(tsswGr, igr, ignore.strand=TRUE)
		# convenient: when constructing sparse matrices, if (i,j) indices are repeated, the x-values are summed up
		cm_tss <- Matrix::sparseMatrix(
			i=queryHits(oo),
			j=elementMetadata(igr)[subjectHits(oo), ".sampleIdx"],
			x=rep(1, length(queryHits(oo))),
			dims=c(nRegs, nSamples)
		)
		colnames(cm_tss) <- sampleIds

		logger.status("computing overlaps with background flanking regions ...")
		ool <- findOverlaps(bgLeftGr, igr, ignore.strand=TRUE)
		oor <- findOverlaps(bgRightGr, igr, ignore.strand=TRUE)
		cm_bg <- Matrix::sparseMatrix(
			i=c(queryHits(ool), queryHits(oor)),
			j=elementMetadata(igr)[c(subjectHits(ool), subjectHits(oor)), ".sampleIdx"],
			x=rep(1, length(queryHits(ool))+length(queryHits(oor))),
			dims=c(nRegs, nSamples)
		)
		colnames(cm_bg) <- sampleIds

		# tsse <- Matrix::colSums(cm_tss, na.rm=TRUE) / Matrix::colSums(cm_bg, na.rm=TRUE)
		# normalize per basepair		
		normCounts_tss <- Matrix::colSums(cm_tss, na.rm=TRUE) / tssW
		normCounts_bg  <- Matrix::colSums(cm_bg, na.rm=TRUE) / (2*tssW) # 2 background regions per TSS
		tsse <- normCounts_tss / pmax(normCounts_bg, 1)
		return(tsse)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("getMonocleCellDataSet")) {
	setGeneric(
		"getMonocleCellDataSet",
		function(.object, ...) standardGeneric("getMonocleCellDataSet"),
		signature=c(".object")
	)
}
#' getMonocleCellDataSet-methods
#'
#' Obtain \code{cell_data_set} object for analysis using the \code{monocle3} package
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionType name of the region type to be exported
#' @param binarize   should the counts be binarized
#' @return a \code{cell_data_set} object containing the counts for the specified region type
#' 
#' @rdname getMonocleCellDataSet-DsATAC-method
#' @docType methods
#' @aliases getMonocleCellDataSet
#' @aliases getMonocleCellDataSet,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getMonocleCellDataSet",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionType,
		binarize=TRUE
	) {
		if (!is.element(regionType, getRegionTypes(.object))) logger.error(c("Unsupported region type:", regionType))

		regGr <- getCoord(.object, regionType)
		regDf <- data.frame(
			row.names = paste(seqnames(regGr),start(regGr),end(regGr),sep="_"),
			site_name = paste(seqnames(regGr),start(regGr),end(regGr),sep="_"),
			chr = gsub("chr", "", as.character(seqnames(regGr))),
			bp1 = start(regGr),
			bp2 = end(regGr)
		)
		names(regGr) <- regDf$site_name
		cm <- ChrAccR::getCounts(.object, regionType, allowSparseMatrix=TRUE)

		if (binarize){
			if (.object@sparseCounts){
				cm@x[cm@x > 0] <- 1
			} else {
				cm[cm > 0] <- 1
			}
		}
		cdsObj <- suppressWarnings(monocle3::new_cell_data_set(
			cm,
			cell_metadata = getSampleAnnot(.object),
			gene_metadata = regDf
		))
		monocle3::fData(cdsObj)$site_name <- as.character(monocle3::fData(cdsObj)$site_name)
		monocle3::fData(cdsObj)$chr <- as.character(monocle3::fData(cdsObj)$chr)
		monocle3::fData(cdsObj)$bp1 <- as.numeric(monocle3::fData(cdsObj)$bp1)
		monocle3::fData(cdsObj)$bp2 <- as.numeric(monocle3::fData(cdsObj)$bp2)

		cdsObj <- monocle3::detect_genes(cdsObj)
		cdsObj <- monocle3::estimate_size_factors(cdsObj)

		return(cdsObj)
	}
)

#-------------------------------------------------------------------------------

if (!isGeneric("getCiceroGeneActivities")) {
	setGeneric(
		"getCiceroGeneActivities",
		function(.object, ...) standardGeneric("getCiceroGeneActivities"),
		signature=c(".object")
	)
}
#' getCiceroGeneActivities-methods
#'
#' Obtain Cicero gene activities
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionType region type of regions that will be linked to the promoter (typical some sort of peak annotation)
#' @param promoterGr \code{GRanges} object of promoter coordinates
#' @param maxDist    maximum distance to consider for region-region interactions
#' @param corCutOff  cutoff of correlation coefficients (Pearson) to consider for region-region interactions
#' @param dimRedCoord matrix of reduced dimension coordinates. must have coordinates for all samples/cells in the dataset
#' @param knn.k      parameter k for Cicero's k-nearest-neighbor method
#' @return an \code{SummarizedExperiment} object containing gene activities for all cells/samples in the dataset
#' 
#' @rdname getCiceroGeneActivities-DsATAC-method
#' @docType methods
#' @aliases getCiceroGeneActivities
#' @aliases getCiceroGeneActivities,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getCiceroGeneActivities",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionType,
		promoterGr=NULL,
		maxDist=250000L,
		corCutOff=0.35,
		dimRedCoord=NULL,
		knn.k=50
	) {
		if (!is.element(regionType, getRegionTypes(.object))) logger.error(c("Unsupported region type:", regionType))
		if (is.null(promoterGr)){
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (is.null(annoPkg)) logger.error("Annotation package needed")
			promoterGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="promoterGr")
		}
		if (is.null(names(promoterGr))){
			emd <- elementMetadata(promoterGr)
			if (is.element("gene_name", colnames(emd))){
				names(promoterGr) <- emd[,"gene_name"]
			} else {
				logger.error("promoterGr must have names")
			}
		}

		logger.start("Creating Cicero object")
			cdsObj <- getMonocleCellDataSet(.object, regionType, binarize=TRUE)
			if (is.null(dimRedCoord)){
				logger.error("Automatic dimension reduction not implemented yet. Must supply dimension reduction coordinates")
			} else {
				if (!all(colnames(cdsObj) %in% rownames(dimRedCoord))) logger.error("all samples/cells must be contained in dimRedCoord object (which must have valid rownames)")
				# ciceroObj <- cicero::make_cicero_cds(cdsObj, k=knn.k, reduced_coordinates=dimRedCoord[colnames(cdsObj),])
				ciceroObj <- custom_cicero_cds(cdsObj, k=knn.k, reduced_coordinates=dimRedCoord[colnames(cdsObj),], silent=TRUE)
			}
		logger.completed()

		logger.start("Computing grouped correlations")
			gr <- getCoord(.object, regionType) # should be in the same order as fData(ciceroObj)
			names(gr) <- monocle3::fData(ciceroObj)$site_name
			
			oo <- suppressWarnings(as.matrix(findOverlaps(resize(resize(gr, 1, "center"), 2*maxDist + 1, "center"), resize(gr, 1, "center"), ignore.strand=TRUE)))
			oo <- data.frame(i=matrixStats::rowMins(oo), j=matrixStats::rowMaxs(oo))
			oo <- oo[!duplicated(oo) & oo[,"i"]!=oo[,"j"],]
			# fast(ish) correlations
			X <- monocle3::exprs(ciceroObj)
			rownames(X) <- NULL
			colnames(X) <- NULL
			mu <- safeMatrixStats(X, "rowMeans", na.rm=TRUE)

			logger.start("Computing Pearson correlation coefficients")
				# chunked processing
				chkSize <- 1e5
				chkStarts <- seq(1, nrow(oo), by=chkSize)
				chkEnds <- nrow(oo)
				nChunks <- length(chkStarts)
				if (nChunks > 1) chkEnds <- c(chkStarts[2:length(chkStarts)]-1, chkEnds)
				cc <- do.call("c", lapply(1:nChunks, FUN=function(i){
					logger.status(c("Chunk", i, "of", nChunks))
					oo_cur <- oo[chkStarts[i]:chkEnds[i], ]
					# centered submatrices
					Xc1 <- as.matrix(X[oo_cur[,1],,drop=FALSE] - mu[oo_cur[,1]])
					Xc2 <- as.matrix(X[oo_cur[,2],,drop=FALSE] - mu[oo_cur[,2]])
					return(rowSums(Xc1*Xc2, na.rm=TRUE) / (sqrt(rowSums(Xc1^2, na.rm=TRUE)) * sqrt(rowSums(Xc2^2, na.rm=TRUE))))
				}))
			logger.completed()
		logger.completed()
		logger.start("Computing gene activities")
			conns <- data.frame(
				Peak1 = names(gr[oo[,"i"]]),
				Peak2 = names(gr[oo[,"j"]]),
				coaccess = cc,
				stringsAsFactors = FALSE
			)
			
			promoterDf <- data.frame(chromosome=seqnames(promoterGr),start=start(promoterGr),end=end(promoterGr), gene=names(promoterGr))
			cdsObj <- cicero::annotate_cds_by_site(cdsObj, promoterDf)

			nSites <- safeMatrixStats(monocle3::exprs(cdsObj), "colSums")
			names(nSites) <- row.names(monocle3::pData(cdsObj))
			unnorm_ga <- cicero::build_gene_activity_matrix(cdsObj, conns, coaccess_cutoff=corCutOff)
			naCells <- safeMatrixStats(unnorm_ga, "colSums")==0
			if (sum(naCells) > 0) {
				logger.warning(c("detected", sum(naCells), "cells/samples with no gene activities"))
				ciceroGA <- Matrix::sparseMatrix(i=c(), j=c(), x=1, dims=c(nrow(unnorm_ga),ncol(unnorm_ga)), dimnames=list(rownames(unnorm_ga), colnames(unnorm_ga)))
				rr <- cicero::normalize_gene_activities(unnorm_ga[,!naCells], nSites[!naCells])
				ciceroGA[,colnames(rr)] <- rr
			} else {
				ciceroGA <- cicero::normalize_gene_activities(unnorm_ga, nSites)
			}
		logger.completed()
		
		seCicero <- SummarizedExperiment::SummarizedExperiment(
			assays = S4Vectors::SimpleList(gA = ciceroGA),
			rowRanges = promoterGr[rownames(ciceroGA)],
			colData = getSampleAnnot(.object)
		)
		S4Vectors::metadata(seCicero)$method <- "cicero"
		S4Vectors::metadata(seCicero)$params <- list(
			maxDist = maxDist,
			corCutOff = corCutOff,
			knn.k = knn.k
		)
		return(seCicero)
	}
)


if (!isGeneric("getRBFGeneActivities")) {
	setGeneric(
		"getRBFGeneActivities",
		function(.object, ...) standardGeneric("getRBFGeneActivities"),
		signature=c(".object")
	)
}
#' getRBFGeneActivities-methods
#'
#' [EXPERIMENTAL] Obtain gene activities by weighting counts using a Gaussian radial basis function (RBF)
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param regionType region type of regions that will be linked to the promoter (typical some sort of peak annotation)
#' @param tssGr      \code{GRanges} object of TSS coordinates or gene body coordinates
#' @param maxDist    maximum distance to consider for associating a region to a TSS
#' @param sigma      decay parameter of the radial basis function (shape parameter (epsilon) of a gaussian RBF= 1/(sqrt(2)*sigma))
#' @param minWeight  weight assymptote, i.e. the assymptotic minimum of the RBF
#' @param binarize   binarize counts before weighting
#' @param useGeneBody include gene bodies in the analysis.
#'                   I.e. counts within the gene body will receive full weight
#'                   instead of downweighting starting from the TSS
#' @param normCounts normalize gene activity "counts". Can be \code{"none"} (no normalization),
#'                   \code{"total"} (normalize by total counts per cell) or 
#' 					 \code{"weighted"} (default; normalize by distance-weighted counts per cell)
#' @return an \code{SummarizedExperiment} object containing gene activities for all cells/samples in the dataset
#' 
#' @rdname getRBFGeneActivities-DsATAC-method
#' @docType methods
#' @aliases getRBFGeneActivities
#' @aliases getRBFGeneActivities,DsATAC-method
#' @author Fabian Mueller
#' @export
#' @noRd
setMethod("getRBFGeneActivities",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		regionType,
		tssGr=NULL,
		maxDist=100000L,
		sigma=10000,
		minWeight=0.25,
		binarize=FALSE,
		useGeneBody=FALSE,
		normCounts="weighted"
	) {
		if (!is.element(regionType, getRegionTypes(.object))) logger.error(c("Unsupported region type:", regionType))
		if (!is.element(normCounts, c("none", "total", "weighted"))) logger.error("Invalid normCounts argument")
		if (is.null(tssGr)){
			annoPkg <- getChrAccRAnnotationPackage(.object@genome)
			if (is.null(annoPkg)) logger.error("Annotation package needed")
			if (useGeneBody){
				tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="geneGr")
			} else {
				tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno="gencode_coding", type="tssGr")
			}
		}
		if (is.null(names(tssGr))){
			emd <- elementMetadata(tssGr)
			if (is.element("gene_name", colnames(emd))){
				names(tssGr) <- emd[,"gene_name"]
			} else {
				logger.error("tssGr must have names")
			}
		}
		logger.start("Computing RBF weights")
			rGr <- getCoord(.object, regionType) # region coords		
			tz <- -1 # true zero replacement
			hassym <- 0
			if (minWeight > 0) hassym <- log(1 - minWeight)
			weightM <- lapply(1:length(tssGr), FUN=function(i){
				# print(i)
				rr <- GenomicRanges::distance(tssGr[i], rGr, ignore.strand=TRUE)
				rr[rr==0] <- tz # true zero distances (regions overlapping TSS)
				rr[is.na(rr) | abs(rr)>maxDist] <- 0
				rr <- as(rr, "sparseMatrix")
				idxTz <- rr@x==tz # which entries are true zeroes
				rr@x <- exp(-rr@x^2/(2*sigma^2) + hassym) + minWeight # apply RBF weights
				rr@x[idxTz] <- 1 # reset true zero weights
				return(rr)
			})
			weightM <- t(do.call("cbind", weightM))
			rownames(weightM) <- names(tssGr)
		logger.completed()

		logger.start("Computing activity scores")
			cm <- ChrAccR::getCounts(.object, regionType, allowSparseMatrix=TRUE)
			if (!is(cm, 'sparseMatrix')) cm <- as(cm, "sparseMatrix")
			if (binarize) cm@x[cm@x > 0] <- 1
			# gaM <- weightM %*% cm
			gaM <- as.matrix(weightM %*% cm) # multiply as sparse matrices: much faster

			scaleFac <- rep(1, ncol(gaM))
			if (normCounts == "total"){
				scaleFac <- 1/Matrix::colSums(cm, na.rm=TRUE) # normalize by total counts
			} else if (normCounts == "weighted"){
				scaleFac <- 1/colSums(gaM, na.rm=TRUE) # normalize by weighted counts
			}
			gaM <- t(t(gaM) * scaleFac) # R operates column-wise while reusing the scale factors
			
			se <- SummarizedExperiment::SummarizedExperiment(
				assays = S4Vectors::SimpleList(gA = gaM),
				rowRanges = tssGr,
				colData = getSampleAnnot(.object)
			)
			S4Vectors::metadata(se)$method <- "RBF"
			S4Vectors::metadata(se)$params <- list(
				maxDist = maxDist,
				sigma = sigma,
				minWeight = minWeight,
				binarize = binarize
			)
		logger.completed()
		return(se)
	}
)

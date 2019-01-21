#' DsATAC
#'
#' A class for storing NOMe accessibility data
#'
#' @section Slots:
#' \describe{
#'   \item{\code{covg}}{
#'		List of GC read coverage for sites and summarized regions.
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'    \item{\code{\link{samples,DsATAC-method}}}{
#'      Retrieve a vector of sample IDs for the dataset
#'    }
#' }
#'
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
		diskDump.fragments = "logical"
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
#' @return \code{data.table} or \code{matrix} containing counts for
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
		naIsZero=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		res <- .object@counts[[type]]
		if (!is.null(i)) res <- res[i,,drop=FALSE]
		if (!is.null(j)) res <- res[,j,drop=FALSE]
		if (asMatrix && !is.matrix(res)){
			res <- as.matrix(res)
			if (!naIsZero && .object@sparseCounts){
				res[res==0] <- NA
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
setMethod("getCountsSE",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type,
		naIsZero=TRUE
	) {
		require(SummarizedExperiment)
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		#count matrix
		cm <- ChrAccR::getCounts(.object, type, asMatrix=TRUE, naIsZero=naIsZero)
		coords <- getCoord(.object, type)
		se <- SummarizedExperiment(assays=list(counts=cm), rowRanges=coords, colData=DataFrame(getSampleAnnot(.object)))
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
#' Return a a \code{GRanges} object of fragment data} for a given sample
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
			} else {
				logger.error(c("Could not load fragment data from file:", fragGr))
			}
		}
		return(fragGr)
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
#' @return \code{GRangesList} containing Tn5 cut sites for each sample
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
		samples=getSamples(.object)
	) {
		if (!all(samples %in% getSamples(.object))) logger.error(c("Invalid samples:", paste(setdiff(samples, getSamples(.object)), collapse=", ")))
		if (!all(samples %in% names(.object@fragments))) logger.error(c("Object does not contain insertion information for samples:", paste(setdiff(samples, names(.object@fragments)), collapse=", ")))
		res <- list()
		for (sid in samples){
			res[[sid]] <- getInsertionSitesFromFragmentGr(getFragmentGr(.object, sid))
		}
		return(GRangesList(res))
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
# 		sl <- tapply(queryHits(oo), subjectHits(oo), c)

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
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param regGr   \code{GRanges} object containing regions to summarize
#' @param type    character string specifying a name for the region type
#' @param signal  character string specifying a name for the region type for the signal to be aggregated
#'                if it is \code{NULL} (default), the new region type will be initialized with NA values
#' @param aggrFun aggregation function for signal counts.
#'                Currently \code{sum}, \code{mean} and \code{median} (default) are supported.
#' @param dropEmpty discard all regions with no observed signal counts
#' @return a new \code{\linkS4class{DsATAC}} object with aggregated signal counts per regions
#'
#' 
#' @rdname regionAggregation-DsATAC-method
#' @docType methods
#' @aliases regionAggregation
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
		dropEmpty=TRUE
	) {
		if (!is.element(aggrFun, c("sum", "mean", "median"))){
			logger.error(c("Unknown signal count aggregation function:", aggrFun))
		}
		if (!is.null(signal)){
			if (!is.element(signal, c(getRegionTypes(.object)), "insertions")){
				logger.error(c("invalid signal region type", signal))
			}
			if (type==signal){
				logger.error(paste0("Region type ('", type, "') is not allowed to be the same as signal ('", signal, "')"))
			}
		}
		if (is.element(type, getRegionTypes(.object))){
			logger.warning(c("Overwriting aggregated region type:", type))
		}
		#sort the regions
		coordOrder <- order(as.integer(seqnames(regGr)), start(regGr), end(regGr), as.integer(strand(regGr)))
		regGr <- regGr[coordOrder]

		.object@coord[[type]] <- regGr	

		#initialize empty count matrix for new region type
		nSamples <- length(getSamples(.object))
		nRegs    <- length(regGr)
		# emptyVec <- rep(as.integer(NA), nRegs)
		# .object@counts[[type]] <- data.table(emptyVec)
		# for (i in 1:nSamples){
		# 	.object@counts[[type]][[i]] <- emptyVec
		# }
		if (.object@sparseCounts){
			.object@counts[[type]] <- sparseMatrix(i=c(), j=c(), x=1, dims=c(nRegs,nSamples))
		} else {
			.object@counts[[type]] <- matrix(as.integer(NA), nrow=nRegs, ncol=nSamples)
		}
		
		if (.object@diskDump) .object@counts[[type]] <- as(.object@counts[[type]], "HDF5Array")
		colnames(.object@counts[[type]]) <- getSamples(.object)
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
			for (sid in getSamples(.object)){
				.object@counts[[type]][,sid] <- as.matrix(countOverlaps(regGr, getInsertionSites(.object, samples=sid)[[1]], ignore.strand=TRUE))
			}
		}
		if (doAggr){
			logger.info(c("Aggregated signal counts across", nrow(.object@counts[[type]]), "regions"))
			# rows2keep <- rowAnys(!is.na(.object@counts[[type]]))
			naMat <- !is.na(.object@counts[[type]])
			if (.object@sparseCounts) naMat <- naMat && .object@counts[[type]] != 0
			rows2keep <- rowSums(naMat) > 0
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
#' @param mergeGroups  factor or character vector or column name in sample annotation table
#' @param countAggrFun aggregation function for signal counts.
#'                Currently \code{sum}, \code{mean} and \code{median} (default) are supported.
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
		countAggrFun="median"
	) {
		if (!is.element(countAggrFun, c("sum", "mean", "median"))){
			logger.error(c("Unknown signal count aggregation function:", countAggrFun))
		}
		if (.object@sparseCounts && is.element(countAggrFun, c("mean", "median"))){
			logger.warning("Mean and median merging of samples can be slow due to conversion of sparse matrices")
		}
		sampleNames <- getSamples(.object)
		nSamples <- length(sampleNames)
		ph <- getSampleAnnot(.object)
		if (is.character(mergeGroups) && length(mergeGroups) == 1 && is.element(mergeGroups, colnames(ph))){
			mergeGroups <- ph[,mergeGroups]
		}
		if ((!is.factor(mergeGroups) && !is.character(mergeGroups)) || length(mergeGroups) != nrow(ph)){
			logger.error("Invalid merge groups")
		}
		if (is.factor(mergeGroups)) mergeGroups <- as.character(mergeGroups)

		phm <- muRtools::aggregateDf(ph, mergeGroups)
		mgL <- lapply(rownames(phm), FUN=function(mg){which(mergeGroups==mg)})
		names(mgL) <- rownames(phm)

		.object@sampleAnnot <- phm

		#count data
		regTypes <- getRegionTypes(.object)
		for (rt in regTypes){
			cm <- ChrAccR::getCounts(.object, rt, asMatrix=FALSE) #design question: should we set the parameter naIsZero==FALSE?
			cmm <- do.call("cbind", lapply(mgL, FUN=function(iis){
				if (.object@sparseCounts && is.element(countAggrFun, c("mean", "median"))){
					aggMat <- ChrAccR::getCounts(.object, rt, j=iis, asMatrix=TRUE) #design question: should we set the parameter naIsZero==FALSE?
				} else {
					aggMat <- cm[,iis,drop=FALSE]
				}
				if(countAggrFun=="sum"){
					return(rowSums(aggMat, na.rm=TRUE))
				} else if(countAggrFun=="mean"){
					return(rowMeans(aggMat, na.rm=TRUE))
				} else if(countAggrFun=="median"){
					return(rowMedians(aggMat, na.rm=TRUE))
				}
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
			insL <- .object@fragments
			.object@fragments <- lapply(mgL, FUN=function(iis){
				rr <- insL[iis]
				rr <- lapply(seq_along(iis), FUN=function(i){
					x <- rr[[i]]
					if (is.character(x)){
						if (file.exists(x)) {
							x <- readRDS(x)
						} else {
							logger.error(c("Could not load fragment data from file:", x))
						}
					}
					elementMetadata(x)[,".sample"] <- sampleNames[iis[i]]
					return(x)
				})
				catRes <- do.call("c", rr)
				if (.object@diskDump.fragments) {
					fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext=".rds")
					saveRDS(catRes, fn)
					catRes <- fn
				}
				return(catRes)
			})
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
		require(GenomicAlignments)
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
		grl
	) {
		sids <- names(grl)
		if (length(sids)!=length(grl)){
			logger.error("The list of GRanges must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}
		if (class(grl)!="GRangesList") grl <- GRangesList(grl)
		gr.c <- unlist(grl, use.names=FALSE)
		if (length(gr.c) < 1) logger.error("[addCountDataFromGRL] invalid GRL: Must be of length 1 or more")
		sampleIds <- rep(names(grl), times=elementNROWS(grl))
		sampleIds.cm <- getSamples(.object)
		rts <- getRegionTypes(.object)

		for (rt in rts){
			logger.status(c("Counting reads in region set:", rt))
			gr.ds <- getCoord(.object, rt)
			oo <- findOverlaps(gr.ds, gr.c, ignore.strand=TRUE)
			idxDt <- as.data.table(cbind(
				queryHits(oo), #row indices (regions) in count matrix
				match(sampleIds[subjectHits(oo)], sampleIds.cm) #column indices (samples) in count matrix
			))
			# count the number of occurrences between each index pair
			idxDt <- idxDt[,.N, by=names(idxDt)]
			idxM <- as.matrix(idxDt[,c(1,2)])
			.object@counts[[rt]][idxM] <- idxDt$N
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
#' Add insertion data to DsATAC object based on bam files
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param fns     a named vector of bam file locations. The names must correspond to sample identifiers in the object
#' @param pairedEnd flag indicating whether the bam files are from paired-end sequencing
#' @param .diskDump for internal use only (should fragment data be stored as RDS instead of GRanges)
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
		.diskDump=.object@diskDump.fragments
	) {
		require(GenomicAlignments)
		sids <- names(fns)
		if (length(sids)!=length(fns)){
			logger.error("The vector specifying filenames (fns) must be named")
		}
		if (!all(sids %in% getSamples(.object))){
			logger.error(c("DsATAC dataset does not contain samples:", paste(setdiff(sids, getSamples(.object)), collapse=", ")))
		}
		for (sid in sids){
			logger.start(c("Reading insertion data for sample:", sid))
				if (!is.null(.object@fragments[[sid]])){
					logger.warning(c("Overwriting insertion data for sample:", sid))
				}
				ga <- NULL
				if (pairedEnd){
					ga <- readGAlignmentPairs(fns[sid], use.names=FALSE, param=ScanBamParam(flag=scanBamFlag(isProperPair=TRUE)))
				} else {
					ga <- readGAlignments(fns[sid], use.names=FALSE)
				}
				ga <- setGenomeProps(ga, .object@genome, onlyMainChrs=TRUE)
				fragGr <- getATACfragments(ga, offsetTn=TRUE)
				if (.diskDump){
					fn <- tempfile(pattern="fragments_", tmpdir=tempdir(), fileext = ".rds")
					saveRDS(fragGr, fn)
					fragGr <- fn
				}
				.object@fragments[[sid]] <- fragGr
			logger.completed()
		}

		return(.object)
	}
)
################################################################################
# Manipulating DsATAC objects
################################################################################

#TODO: not tested yet
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
#' @aliases removeRegions
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
		.object@counts[[type]]  <- .object@counts[[type]][inds2keep,]

		# rts <- setdiff(getRegionTypes(.object), type)
		# if (reaggregate && type == "signal" && length(rts)>0) {
		# 	logger.start("Recomputing region aggregation")
		# 	rtGrl <- .object@coord

		# 	.object@coord <- .object@coord["signal"]
		# 	.object@counts  <- .object@counts["signal"]

		# 	for (rt in rts){
		# 		logger.status(c(rt, "..."))
		# 		.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
		# 	}
		# 	logger.completed()	
		# }
		return(.object)
	}
)
#-------------------------------------------------------------------------------
#TODO: not tested yet
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
		}
		
		return(.object)
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
#' @param method  transformation method to be applied. Currently only 'log2', 'quantile' (quantile normalization), 'vst' (DESeq2 Variance Stabilizing Transformation) and 'RPKM' (RPKM normalization) are supported
#' @param regionTypes character vector specifying a name for the region type in which count data should be normalized(default: all region types)
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
		method="quantileNorm",
		regionTypes=getRegionTypes(.object)
	) {
		if (!all(regionTypes %in% getRegionTypes(.object))){
			logger.error(c("Unsupported region type:", paste(setdiff(regionTypes, getRegionTypes(.object)), collapse=", ")))
		}
		if (!is.element(method, c("quantile", "log2", "RPKM", "vst", "tf-idf"))) logger.error(c("Unsupported normalization method type:", method))

		if (method == "quantile"){
			logger.start(c("Performing quantile normalization"))
				require(preprocessCore)
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					# .object@counts[[rt]] <- data.table(normalize.quantiles(as.matrix(.object@counts[[rt]])))
					.object@counts[[rt]] <- normalize.quantiles(ChrAccR::getCounts(.object, rt, asMatrix=TRUE)) #design question: should we set the parameter naIsZero==FALSE?
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
		} else if (method == "RPKM"){
			logger.start(c("Performing RPKM normalization"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					# cm <- as.matrix(.object@counts[[rt]])
					cm <- .object@counts[[rt]]
					regLen <- width(getCoord(.object, rt))
					sizeFac <- matrix(colSums(cm, na.rm=TRUE), ncol=ncol(cm), nrow=nrow(cm), byrow=TRUE)
					# .object@counts[[rt]] <- data.table(cm/(regLen * sizeFac) * 1e3 * 1e6)
					.object@counts[[rt]] <- cm/(regLen * sizeFac) * 1e3 * 1e6
					if (.object@diskDump) .object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("RPKM", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "log2"){
			c0 <- 1
			logger.start(c("log2 transforming counts"))
				if (.object@sparseCounts) logger.error("Generating sparse matrix for matrix with possible true zero entries (log2)")
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					idx <- .object@counts[[rt]]!=0
					.object@counts[[rt]][idx] <- log2(.object@counts[[rt]][idx] + c0)
					.object@countTransform[[rt]] <- c("log2", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "vst"){
			logger.start(c("Applying DESeq2 VST"))
				if (.object@sparseCounts) logger.warning("Generating sparse matrix for matrix with possible true zero entries (VST)")
				require(DESeq2)
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					dds <- DESeqDataSet(getCountsSE(.object, rt), design=~1)
					# .object@counts[[rt]] <- data.table(assay(vst(dds, blind=TRUE)))
					.object@counts[[rt]] <- assay(vst(dds, blind=TRUE))
					if (!.object@diskDump && .object@sparseCounts){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "sparseMatrix")
					}
					if (.object@diskDump){
						.object@counts[[rt]] <- as(.object@counts[[rt]], "HDF5Array")
					}
					.object@countTransform[[rt]] <- c("deseq.vst", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "tf-idf"){
			# TF-IDF transformation as applied in LSI (Shendure lab) (Cusanovich, et al. (2018). A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. Cell, 1-35)
			# Recycled some code from: https://github.com/shendurelab/mouse-atac
			logger.start(c("Applying TF-IDF transformation"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cm <- !is.na(.object@counts[[rt]]) & .object@counts[[rt]] > 0 #indicator matrix: are there any counts in that region
					cnames <- colnames(cm)
					if (class(cm)=="lgCMatrix"){
						tf <- Matrix::t(Matrix::t(cm) / Matrix::colSums(cm)) #term frequency
						idf <- tf * log(1 + ncol(cm) / Matrix::rowSums(cm)) # inverse document frequency
					} else {
						tf <- t(t(cm) / colSums(cm)) #term frequency
						idf <- tf * log(1 + ncol(cm) / rowSums(cm)) # inverse document frequency
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
			rem <- rowSums(getCounts(.object, rt, naIsZero=TRUE)) < numAllowed
			nRem <- sum(rem)
			nRegs <- getNRegions(.object, rt)
			if (nRem > 0){
				.object <- removeRegions(.object, rem, rt)
			}
			logger.status(c("Removed", nRem, "regions", paste0("(", round(nRem/nRegs, 4)*100, "%)"), "of type", rt))
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
		if (class(rsl)!="GRangesList") rsl <- GRangesList(rsl)

		res <- matrix(as.numeric(NA), nrow=length(rsl), ncol=length(getSamples(.object)))
		if (bySample){
			rslGr <- unlist(rsl, use.names=FALSE)
			idx.rsl <- rep(1:length(rsl), times=elementNROWS(rsl))

			res <- do.call("cbind", lapply(getSamples(.object), FUN=function(sid){
				insGr <- getInsertionSites(.object, sid)
				ov <- countOverlaps(rslGr, insGr, ignore.strand=TRUE)
				rr <- tapply(ov, idx.rsl, sum)
				rr <- rr[as.character(1:length(rsl))] # make sure the counts are returned in the same order as in rsl (not sure, if really necessary)
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
		require(qvalue)
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
		mm <- matchMotifs(mmObjs[["motifs"]], se, genome=mmObjs[["genome"]])
		regionMotifMatch <- as.matrix(motifMatches(mm))

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
				 qvalue(res[,"pVal"])$qvalue,
				 error = function(ee) {
				 	if (ee$message == "missing or infinite values in inputs are not allowed"){
				 		return(qvalue(res[,"pVal"], lambda=0)$qvalue)
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
#' @return a \code{data.frame} summarizing Fisher's Exact Test enrichment statistics for each motif
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
		require(chromVAR)
		require(motifmatchr)
		res <- NULL

		countSe <- getCountsSE(.object, type, naIsZero=TRUE)
		genomeObj <- getGenomeObject(.object@genome)
		genome(countSe) <- providerVersion(genomeObj) # hack to override inconsistent naming of genome versions (e.g. hg38 and GRCh38)

		countSe <- addGCBias(countSe, genome=genomeObj)

		# for motifmatchr
		mmInput <- prepareMotifmatchr(genomeObj, motifs)
		mmObj <- matchMotifs(mmInput[["motifs"]], countSe, genome=genomeObj)

		res <- computeDeviations(object=countSe, annotations=mmObj)

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

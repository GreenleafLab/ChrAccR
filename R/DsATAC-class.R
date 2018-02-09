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
		counts = "list"
	),
	contains = "DsAcc",
	package = "ChrAccR"
)
setMethod("initialize","DsATAC",
	function(
		.Object,
		coord,
		counts,
		sampleAnnot,
		genome
	) {
		.Object@coord       <- coord
		.Object@counts      <- counts
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object
	}
)

#' @param sampleAnnot \code{data.frame} object containing sample annotation
#' @param genome    character string containing genome assembly
#' @noRd
DsATAC <- function(sampleAnnot, genome){
	obj <- new("DsATAC",
		list(),
		list(),
		sampleAnnot,
		genome
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
#' @param asMatrix return a matrix instead of a \code{data.table}
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
		asMatrix=FALSE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		res <- .object@counts[[type]]
		if (asMatrix) res <- as.matrix(res)
		return(res)
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

		cat("DsATAC chromatin accessibility dataset \n")
		cat("contains:\n")
		cat(" * ", length(ss), " samples: ", str.ss, " \n")
		cat(" * ", str.rts, " \n")
	}
)

################################################################################
# Summary functions
################################################################################

#-------------------------------------------------------------------------------
if (!isGeneric("getRegionMapping")) {
	setGeneric(
		"getRegionMapping",
		function(.object, ...) standardGeneric("getRegionMapping"),
		signature=c(".object")
	)
}
#' getRegionMapping-methods
#'
#' Retrieve a mapping from regions to GC indices in the dataset
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param type    character string specifying a name for the region type
#' @return list containing vectors of indices of GCs for each region of the
#'         specified type
#' 
#' @rdname getRegionMapping-DsATAC-method
#' @docType methods
#' @aliases getRegionMapping
#' @aliases getRegionMapping,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("getRegionMapping",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		type
	) {
		baseGr <- getCoord(.object, type="counts")
		regGr  <- getCoord(.object, type=type)
		oo <- findOverlaps(baseGr, regGr)
		sl <- tapply(queryHits(oo), subjectHits(oo), c)

		res <- rep(list(integer(0)), getNRegions(.object, type=type))
		res[as.integer(names(sl))] <- sl
		return(res)
	}
)
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
			if (!is.element(signal, getRegionTypes(.object))){
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
		emptyVec <- rep(as.integer(NA), nRegs)
		.object@counts[[type]] <- data.table(emptyVec)
		for (i in 1:nSamples){
			.object@counts[[type]][[i]] <- emptyVec
		}
		colnames(.object@counts[[type]]) <- getSamples(.object)

		#Aggregate count signal
		if (!is.null(signal) && is.element(aggrFun, c("sum", "mean", "median"))){
			signalGr <- getCoord(.object, "signal")
			if (genome(regGr)[1]!=genome(signalGr)[1]){
				logger.warning(c("Potentially incompatible genome assemblies (object:", genome(signalGr)[1], ", regGr:", genome(regGr)[1], ")"))
			}
			oo <- findOverlaps(signalGr, regGr)
			if (any(duplicated(queryHits(oo)))) logger.info("Some signals map to multiple regions")
			dtC <- data.table(getCounts(.object)[queryHits(oo),], mergedIndex=subjectHits(oo))

			if (aggrFun=="sum") {
				rr <- dtC[, lapply(.SD, sum, na.rm=TRUE), by=.(mergedIndex)]
			} else if (aggrFun=="mean") {
				rr <- dtC[, lapply(.SD, mean, na.rm=TRUE), by=.(mergedIndex)]
			} else if (aggrFun=="median") {
				rr <- dtC[, lapply(.SD, median, na.rm=TRUE), by=.(mergedIndex)]
			} else {
				logger.error(c("Unknown signal count aggregation function:", aggrFun))
			}
			.object@counts[[type]][rr[["mergedIndex"]],] <- rr[,!"mergedIndex"]
			rm(dtC, rr); cleanMem() #clean-up

			logger.info(c("Aggregated signal counts across", nrow(.object@counts[[type]]), "regions"))
			rows2keep <- rowAnys(!is.na(.object@counts[[type]]))
			logger.info(c("  of which", sum(rows2keep), "regions contained signal counts"))
			#discard regions where all signal counts are unobserved
			if (dropEmpty){
				.object@coord[[type]]   <- .object@coord[[type]][rows2keep]
				.object@counts[[type]]  <- .object@counts[[type]][rows2keep,]
			}
		}

		return(.object)
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
setMethod("addCountDataFromBam",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		fns
	) {
		require(GenomicAlignments)
		
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
			.object@counts[[rt]][,sids] <- assays(ov.rse)$count 
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
		rts <- getRegionTypes(.object)
		for (rt in rts){
			logger.status(c("Counting reads in region set:", rt))
			gr.ds <- getCoord(.object, rt)
			for (sid in sids){
				gr.c <- grl[[sid]]
				.object@counts[[rt]][,sid] <- countOverlaps(gr.ds, gr.c, ignore.strand=TRUE)
			}
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
#' @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
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
		type="signal",
		reaggregate=TRUE
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
		.object@counts[[type]]  <- .object@meth[[type]][inds2keep,]

		rts <- setdiff(getRegionTypes(.object), type)
		if (reaggregate && type == "signal" && length(rts)>0) {
			logger.start("Recomputing region aggregation")
			rtGrl <- .object@coord

			.object@coord <- .object@coord["signal"]
			.object@counts  <- .object@counts["signal"]

			for (rt in rts){
				logger.status(c(rt, "..."))
				.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
			}
			logger.completed()	
		}
		return(.object)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("normalizeCounts")) {
	setGeneric(
		"normalizeCounts",
		function(.object, ...) standardGeneric("normalizeCounts"),
		signature=c(".object")
	)
}
#' normalizeCounts-methods
#'
#' Normalize methylation levels
#'
#' @param .object \code{\linkS4class{DsATAC}} object
#' @param method  normalization method to be applied. Currently only 'quantile' is supported
#' @param regionTypes character vector specifying a name for the region type in which count data should be normalized(default: all region types)
#' @return a new \code{\linkS4class{DsATAC}} object with normalized count data
#' 
#' @rdname normalizeCounts-DsATAC-method
#' @docType methods
#' @aliases normalizeCounts
#' @aliases normalizeCounts,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("normalizeCounts",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		method="quantile",
		regionTypes=getRegionTypes(.object)
	) {
		if (!all(regionTypes %in% getRegionTypes(.object))){
			logger.error(c("Unsupported region type:", paste(setdiff(regionTypes, getRegionTypes(.object)), collapse=", ")))
		}
		if (!is.element(method, c("quantile"))) logger.error(c("Unsupported normalization method type:", method))

		if (method == "quantile"){
			logger.start(c("Performing quantile normalization"))
				require(preprocessCore)
				for (rt in regionTypes){
					logger.status(c("Regiomn type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					.object@counts[[rt]] <- data.table(normalize.quantiles(as.matrix(.object@counts[[rt]])))
					colnames(.object@counts[[rt]]) <- cnames
				}
			logger.completed()
		}
		return(.object)
	}
)

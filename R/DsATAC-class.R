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
		insertions = "list",
		counts = "list",
		countTransform = "list"
	),
	contains = "DsAcc",
	package = "ChrAccR"
)
setMethod("initialize","DsATAC",
	function(
		.Object,
		insertions,
		coord,
		counts,
		sampleAnnot,
		genome
	) {
		.Object@insertions  <- insertions
		.Object@coord       <- coord
		.Object@counts      <- counts
		.Object@countTransform <- rep(list(character(0)), length(.Object@counts))
		names(.Object@countTransform) <- names(.Object@counts)
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object@pkgVersion  <- packageVersion("ChrAccR")
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
		asMatrix=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object))) logger.error(c("Unsupported region type:", type))
		res <- .object@counts[[type]]
		if (asMatrix) res <- as.matrix(res)
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
		if (!all(samples %in% names(.object@insertions))) logger.error(c("Object does not contain insertion information for samples:", paste(setdiff(samples, names(.object@insertions)), collapse=", ")))
		res <- list()
		for (sid in samples){
			isW <- width(.object@insertions[[sid]])>1 # the insertion site is already width=1 --> single end. For paired end-data all of these should be TRUE
			grins <- NULL
			if (any(!isW)){
				grins <- .object@insertions[[sid]][isW]
			}
			if (all(isW)){
				# paired-end data - default case
				grins <- c(
					resize(.object@insertions[[sid]], width=1, fix="start"),
					resize(.object@insertions[[sid]], width=1, fix="end")
				)
			} else if (any(isW)){
				# mixed paired-end and single-end data
				logger.warning(c("mixed paired-end and single-end data detected for sample", sid))
				grins <- c(
					grins,
					resize(.object@insertions[[sid]][isW], width=1, fix="start"),
					resize(.object@insertions[[sid]][isW], width=1, fix="end")
				)
			}
			#sort the result
			res[[sid]] <- grins[order(as.integer(seqnames(grins)), start(grins), end(grins), as.integer(strand(grins)))]
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

		cat("DsATAC chromatin accessibility dataset \n")
		cat("contains:\n")
		cat(" * ", length(ss), " samples: ", str.ss, " \n")
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
		emptyVec <- rep(as.integer(NA), nRegs)
		.object@counts[[type]] <- data.table(emptyVec)
		for (i in 1:nSamples){
			.object@counts[[type]][[i]] <- emptyVec
		}
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
		}
		if (!is.null(signal) && signal=="insertions"){
			doAggr <- TRUE
			for (sid in getSamples(.object)){
				.object@counts[[type]][,sid] <- countOverlaps(regGr, getInsertionSites(.object, samples=sid)[[1]], ignore.strand=TRUE)			
			}
		}
		if (doAggr){
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
		# Buenrostro, et al. (2013). Nature Methods, 10(12), 1213–1218. 
		# For peak-calling and footprinting, we adjusted the read start sites to represent the center of the transposon binding event. Previous descriptions of the Tn5 transposase show that the trans- poson binds as a dimer and inserts two adaptors separated by 9 bp (ref. 11). Therefore, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the – strand were offset −5 bp
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
#' @return a new \code{\linkS4class{DsATAC}} object with insertions for each sample
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
		pairedEnd=TRUE
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
				if (!is.null(.object@insertions[[sid]])){
					logger.warning(c("Overwriting insertion data for sample:", sid))
				}
				ga <- NULL
				if (pairedEnd){
					ga <- readGAlignments(fns[sid], use.names=FALSE)
				} else {
					ga <- readGAlignmentPairs(fns[sid], use.names=TRUE)
				}
				.object@insertions[[sid]] <- getATACinsertion(ga, offsetTn=TRUE)
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
#' @param method  transformation method to be applied. Currently only 'log2', 'quantile' (quantile normalization) and 'RPKM' (RPKM normalization) are supported
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
		if (!is.element(method, c("quantile"))) logger.error(c("Unsupported normalization method type:", method))

		if (method == "quantile"){
			logger.start(c("Performing quantile normalization"))
				require(preprocessCore)
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					.object@counts[[rt]] <- data.table(normalize.quantiles(as.matrix(.object@counts[[rt]])))
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("quantileNorm", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "RPKM"){
			logger.start(c("Performing RPKM normalization"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					cnames <- colnames(.object@counts[[rt]])
					cm <- as.matrix(.object@counts[[rt]])
					regLen <- width(getCoord(.object, rt))
					sizeFac <- matrix(colSums(cm, na.rm=TRUE), ncol=ncol(cm), nrow=nrow(cm), byrow=TRUE)
					.object@counts[[rt]] <- data.table(cm/(regLen * sizeFac) * 1e3 * 1e6)
					colnames(.object@counts[[rt]]) <- cnames
					.object@countTransform[[rt]] <- c("RPKM", .object@countTransform[[rt]])
				}
			logger.completed()
		} else if (method == "log2"){
			logger.start(c("log2 transforming counts"))
				for (rt in regionTypes){
					logger.status(c("Region type:", rt))
					.object@counts[[rt]] <- log2(.object@counts[[rt]])
					.object@countTransform[[rt]] <- c("log2", .object@countTransform[[rt]])
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
			rem <- rowSums(getCounts(.object, rt) >= thresh) < numAllowed
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
		cm <- ChrAccR::getCounts(.object, type, asMatrix=TRUE)
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
		se <- SummarizedExperiment(assays=list(counts=cm), rowRanges=coords)
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
			fr <- fisher.test(x=idx, y=regionMotifMatch[,mo], alternative="greater")
			df <- data.frame(
				pVal=fr$p.value,
				oddsRatio=fr$estimate
			)
			return(df)
		}))
		rownames(res) <- motifNames
		res[,"motif"] <- motifNames
		res[,"qVal"] <- qvalue(res[,"pVal"])$qvalue

		return(res)
	}
)

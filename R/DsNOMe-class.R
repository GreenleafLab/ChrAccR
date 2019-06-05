#' DsNOMe
#'
#' A class for storing NOMe accessibility data
#'
#' @section Slots:
#' \describe{
#'   \item{\code{meth}}{
#'		List of GC methylation for sites and summarized regions.
#'   }
#'   \item{\code{covg}}{
#'		List of GC read coverage for sites and summarized regions.
#'   }
#' }
#'
#' @name DsNOMe-class
#' @rdname DsNOMe-class
#' @author Fabian Mueller
#' @exportClass DsNOMe
setClass("DsNOMe",
	slots = list(
		meth        = "list",
		covg        = "list"
	),
	contains = "DsAcc",
	package = "ChrAccR"
)
setMethod("initialize","DsNOMe",
	function(
		.Object,
		coord,
		meth,
		covg,
		sampleAnnot,
		genome
	) {
		.Object@coord       <- coord
		.Object@meth        <- meth
		.Object@covg        <- covg
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object@pkgVersion  <- packageVersion("ChrAccR")
		.Object
	}
)

#' @param siteCoord \code{GRanges} object containing coordinates of GC dinucleotides
#' @param siteMeth  \code{data.table} object containing methylation values for each
#'                  GC dinucleotide and each sample 
#' @param siteCovg  \code{data.table} object containing read coverage values for each
#'                  GC dinucleotide and each sample 
#' @param sampleAnnot \code{data.frame} object containing sample annotation
#' @param genome    character string containing genome assembly
#' @noRd
DsNOMe <- function(siteCoord, siteMeth, siteCovg, sampleAnnot, genome){
	obj <- new("DsNOMe",
		list(sites=siteCoord),
		list(sites=siteMeth),
		list(sites=siteCovg),
		sampleAnnot,
		genome
	)
	return(obj)
}

################################################################################
# Getters
################################################################################
#generic is already defined in minfi package --> redefine
# if (!isGeneric("getMeth")) {
	setGeneric(
		"getMeth",
		function(.object, ...) standardGeneric("getMeth"),
		signature=c(".object")
	)
# }
#' getMeth-methods
#'
#' Return table of methylation values
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @param asMatrix return a matrix instead of a \code{data.table}
#' @return \code{data.table} or \code{matrix} containing methylation levels for
#'         each site/region and sample
#'
#' @rdname getMeth-DsNOMe-method
#' @docType methods
#' @aliases getMeth
#' @aliases getMeth,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getMeth",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type="sites",
		asMatrix=FALSE
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		res <- .object@meth[[type]]
		if (asMatrix) res <- as.matrix(res)
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getCovg")) {
	setGeneric(
		"getCovg",
		function(.object, ...) standardGeneric("getCovg"),
		signature=c(".object")
	)
}
#' getCovg-methods
#'
#' Return table of read coverage values
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @param asMatrix return a matrix instead of a \code{data.table}
#' @return \code{data.table} or \code{matrix} containing read coverage for
#'         each site/region and sample
#'
#' @rdname getCovg-DsNOMe-method
#' @docType methods
#' @aliases getCovg
#' @aliases getCovg,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getCovg",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type="sites",
		asMatrix=FALSE
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		res <- .object@covg[[type]]
		if (asMatrix) res <- as.matrix(res)
		return(res)
	}
)

################################################################################
# Display
################################################################################
setMethod("show","DsNOMe",
	function(object) {
		ss <- getSamples(object)
		str.ss <- paste(getSamples(object), collapse=", ")
		if (length(ss) > 5) str.ss <- paste(c(getSamples(object)[1:5], "..."), collapse=", ")
		rts <- getRegionTypes(object)
		str.rts <- "no region types"
		if (length(rts) > 0) str.rts <- paste0(length(rts), " region types: ", paste(rts, collapse=", "))

		cat("DsNOMe chromatin accessibility dataset \n")
		cat("contains:\n")
		cat(" * ", length(ss), " samples: ", str.ss, " \n")
		cat(" * ", getNRegions(object, "sites"), " GC methylation measurements", "\n")
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
if (!isGeneric("mergeStrands")) {
	setGeneric(
		"mergeStrands",
		function(.object, ...) standardGeneric("mergeStrands"),
		signature=c(".object")
	)
}
#' mergeStrands-methods
#'
#' Merge + and - strands of the dataset by adding read coverage and recomputing
#' Methylation levels
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param reaggregate redo region aggregation (only has an effect if there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsNOMe}} object with the strands merged
#'
#' @rdname mergeStrands-DsNOMe-method
#' @docType methods
#' @aliases mergeStrands
#' @aliases mergeStrands,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("mergeStrands",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		reaggregate=TRUE
	) {
		gr <- getCoord(.object)
		doMerge <- length(unique(seqlevels(gr))) > 1
		if (!doMerge) {
			logger.warning("Nothing to do: data is not from more than one strand")
		} else {
			gr.merged <- gr
			strand(gr.merged) <- "*"
			gr.merged <- unique(gr.merged)
			newN <- length(gr.merged)

			oo <- findOverlaps(gr, gr.merged, type="equal", ignore.strand=TRUE)
			if (queryLength(oo)!=length(gr)) logger.error("Something went wrong in matching coords to merged coords [queryLength]")
			if (length(oo)!=length(gr)) logger.error("Something went wrong in matching coords to merged coords [non-unique]")
			mergedIndex <- subjectHits(oo) # target index in merged GRanges object
			mergedIndex.char <- as.character(mergedIndex)

			#assign new, merged coordinates
			.object@coord[["sites"]] <- gr.merged
			rm(gr.merged, mergedIndex, oo); cleanMem() #clean-up

			#join coverage
			dtC <- getCovg(.object)
			dtC[,mergedIndex:=mergedIndex.char]
			.object@covg[["sites"]] <- dtC[, lapply(.SD, sum, na.rm=TRUE), by=.(mergedIndex)][,mergedIndex:=NULL] #! NAs are now 0s
			dtC[,mergedIndex:=NULL] #reverse concurrent programming effect (added mergedIndex column)

			dtM <- getMeth(.object)
			#compute the mean methylation, weighted by coverage
			emptyVec <- rep(as.numeric(NA), newN)
			rr <- data.table(emptyVec)
			for (i in 1:ncol(dtM)){
				rr[[i]] <- emptyVec
			}
			colnames(rr) <- colnames(dtM)
			for (cn in colnames(dtM)){
				logger.status(c("column:", cn))
				dt <- data.table(meth=dtM[[cn]], covg=dtC[[cn]], mergedIndex=mergedIndex.char)
				dt[,wm:=meth*covg]
				dtr <- dt[,.(sum(wm, na.rm=TRUE), sum(covg, na.rm=TRUE)), by=.(mergedIndex)]
				dtr[,weightedSum:=V1/V2]
				# dtr <- dt[, lapply(.SD, weighted.mean, w=covg), by=.(mergedIndex)] #takes forever
				# dtr <- dt[, lapply(.SD, function(x, w){sum(x*w, na.rm=TRUE)/sum(w, na.rm=TRUE)}, w=covg), by=.(mergedIndex)] #takes forever
				rr[[cn]] <- dtr[["weightedSum"]]

				rm(dt, dtr); cleanMem() #clean-up
			}
			.object@meth[["sites"]] <- rr[,colnames(.object@covg[["sites"]]), with=FALSE]
			rm(dtC, dtM, rr); cleanMem() #clean-up
			
			# if (length(getRegionTypes(.object))>0){
			# 	logger.warning("Detected multiple region types. These will be dropped and not be recalculated")
			# 	.object@coord <- .object@coord["sites"]
			# 	.object@meth  <- .object@meth["sites"]
			# 	.object@covg  <- .object@covg["sites"]
			# }
			rts <- getRegionTypes(.object)
			if (reaggregate && length(rts)>0) {
				logger.start("Recomputing region aggregation")
				rtGrl <- .object@coord

				.object@coord <- .object@coord["sites"]
				.object@meth  <- .object@meth["sites"]
				.object@covg  <- .object@covg["sites"]

				for (rt in rts){
					logger.status(c(rt, "..."))
					.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
				}
				logger.completed()	
			}
		}
		return(.object)
	}
)

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
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param type    character string specifying a name for the region type
#' @return list containing vectors of indices of GCs for each region of the
#'         specified type
#' 
#' @rdname getRegionMapping-DsNOMe-method
#' @docType methods
#' @aliases getRegionMapping
#' @aliases getRegionMapping,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getRegionMapping",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type
	) {
		siteGr <- getCoord(.object)
		regGr  <- getCoord(.object, type=type)
		oo <- findOverlaps(siteGr, regGr)
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
#' Aggregate methylation levels and coverage values accross a set of regions
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param regGr   \code{GRanges} object containing regions to summarize
#' @param type    character string specifying a name for the region type
#' @param methAggrFun aggregation function for methylation levels.
#'                Currently \code{mean}, \code{median}
#'                and \code{weightedMean} (default) are supported.
#' @param dropEmpty discard all regions with no observed methylation levels
#' @return a new \code{\linkS4class{DsNOMe}} object with aggregated regions
#'
#' @details
#' Coverage values are aggregated by summing up coverage values for individual GCs
#' while the aggregation function for methylation levels is specified by the
#' \code{methAggrFun} parameter.
#' 
#' @rdname regionAggregation-DsNOMe-method
#' @docType methods
#' @aliases regionAggregation
#' @aliases regionAggregation,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("regionAggregation",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		regGr,
		type,
		methAggrFun="weightedMean",
		dropEmpty=TRUE
	) {
		if (!is.element(methAggrFun, c("mean", "median", "weightedMean"))){
			logger.error(c("Unknown methylation aggregation function:", methAggrFun))
		}
		if (type=="sites"){
			logger.error("Region type is not allowed to be named 'sites'")
		}
		if (is.element(type, getRegionTypes(.object))){
			logger.warning(c("Overwriting aggregated region type:", type))
		}
		#sort the regions
		coordOrder <- order(as.integer(seqnames(regGr)), start(regGr), end(regGr), as.integer(strand(regGr)))
		regGr <- regGr[coordOrder]

		siteGr <- getCoord(.object)
		if (genome(regGr)[1]!=genome(siteGr)[1]){
			logger.warning(c("Potentially incompatible genome assemblies (object:", genome(siteGr)[1], ", regGr:", genome(regGr)[1], ")"))
		}
		oo <- findOverlaps(siteGr, regGr)
		if (any(duplicated(queryHits(oo)))) logger.info("Some GCs map to multiple regions")

		.object@coord[[type]] <- regGr

		# aggregate coverage (sum)
		dtC <- data.table(getCovg(.object)[queryHits(oo),], mergedIndex=subjectHits(oo))
		rr <- dtC[, lapply(.SD, sum, na.rm=TRUE), by=.(mergedIndex)]

		#initialize empty coverage matrix
		nSamples <- length(getSamples(.object))
		nRegs    <- length(regGr)
		emptyVec <- rep(as.integer(NA), nRegs)
		.object@covg[[type]] <- data.table(emptyVec)
		for (i in 1:nSamples){
			.object@covg[[type]][[i]] <- emptyVec
		}
		colnames(.object@covg[[type]]) <- colnames(.object@covg[["sites"]])

		.object@covg[[type]][rr[["mergedIndex"]],] <- rr[,!"mergedIndex"]
		rm(rr); cleanMem() #clean-up

		# aggregate methylation
		dtM <- data.table(getMeth(.object)[queryHits(oo),], mergedIndex=subjectHits(oo))

		#initialize empty methylation matrix
		emptyVec <- rep(as.numeric(NA), nRegs)
		.object@meth[[type]] <- data.table(emptyVec)
		for (i in 1:nSamples){
			.object@meth[[type]][[i]] <- emptyVec
		}
		colnames(.object@meth[[type]]) <- colnames(.object@meth[["sites"]])

		if (is.element(methAggrFun, c("mean", "median"))){
			if (methAggrFun=="mean") {
				rr <- dtM[, lapply(.SD, mean, na.rm=TRUE), by=.(mergedIndex)]
			} else if (methAggrFun=="median") {
				rr <- dtM[, lapply(.SD, median, na.rm=TRUE), by=.(mergedIndex)]
			} else {
				logger.error(c("Unknown methylation aggregation function:", methAggrFun))
			}
			.object@meth[[type]][rr[["mergedIndex"]],] <- rr[,!"mergedIndex"]
			rm(dtM, dtC, rr); cleanMem() #clean-up
		} else if (is.element(methAggrFun, c("weightedMean"))){
			for (cn in setdiff(colnames(dtM), "mergedIndex")){
				# logger.status(c("sample:", cn))
				dt <- data.table(meth=dtM[[cn]], covg=dtC[[cn]], mergedIndex=dtM[["mergedIndex"]])
				dt[,wm:=meth*covg]
				dtr <- dt[,.(sum(wm, na.rm=TRUE), sum(covg, na.rm=TRUE)), by=.(mergedIndex)]
				dtr[,weightedSum:=V1/V2]
				.object@meth[[type]][dtr[["mergedIndex"]], cn] <- dtr[["weightedSum"]]
				rm(dt, dtr); cleanMem() #clean-up
			}
			rm(dtM, dtC); cleanMem() #clean-up
		} else {
			logger.error(c("Unknown methylation aggregation function:", methAggrFun))
		}

		logger.info(c("Aggregated measurements across", nrow(.object@meth[[type]]), "regions"))
		rows2keep <- rowAnys(!is.na(.object@meth[[type]]))
		logger.info(c("  of which", sum(rows2keep), "regions contained GCs with measurements"))
		#discard regions where all methylation levels are unobserved
		if (dropEmpty){
			.object@coord[[type]] <- .object@coord[[type]][rows2keep]
			.object@meth[[type]]  <- .object@meth[[type]][rows2keep,]
			.object@covg[[type]]  <- .object@covg[[type]][rows2keep,]
		}

		return(.object)
	}
)

################################################################################
# Maniputlating DsNOMe objects
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
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param indices a vector of indices of sites/regions to be removed. Can be numeric, integer or logical.
#' @param type    character string specifying a name for the region type (sefault: sites)
#' @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsNOMe}} object with sites/regions removed
#' 
#' @rdname removeRegions-DsNOMe-method
#' @docType methods
## @aliases removeRegions
#' @aliases removeRegions,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("removeRegions",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		indices,
		type="sites",
		reaggregate=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))

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
		.object@meth[[type]]  <- .object@meth[[type]][inds2keep,]
		.object@covg[[type]]  <- .object@covg[[type]][inds2keep,]

		rts <- setdiff(getRegionTypes(.object), type)
		if (reaggregate && type == "sites" && length(rts)>0) {
			logger.start("Recomputing region aggregation")
			rtGrl <- .object@coord

			.object@coord <- .object@coord["sites"]
			.object@meth  <- .object@meth["sites"]
			.object@covg  <- .object@covg["sites"]

			for (rt in rts){
				logger.status(c(rt, "..."))
				.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
			}
			logger.completed()	
		}
		return(.object)
	}
)

if (!isGeneric("maskMethNA")) {
	setGeneric(
		"maskMethNA",
		function(.object, ...) standardGeneric("maskMethNA"),
		signature=c(".object")
	)
}
#' maskMethNA-methods
#'
#' Set the indices specified in a mask to NA
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param mask    a mask, i.e. a logical matrix of indices to set to NA
#' @param type    character string specifying a name for the region type (default: sites)
#' @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsNOMe}} object with sites/regions masked
#' 
#' @rdname maskMethNA-DsNOMe-method
#' @docType methods
#' @aliases maskMethNA
#' @aliases maskMethNA,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("maskMethNA",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		mask,
		type="sites",
		reaggregate=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))

		if (!is.logical(mask) || nrow(mask)!=nrow(.object@meth[[type]]) || ncol(mask)!=ncol(.object@meth[[type]])){
			logger.error(c("Unsupported type for mask matrix"))
		}

		.object@meth[[type]][mask]  <- NA

		rts <- setdiff(getRegionTypes(.object), type)
		if (reaggregate && type == "sites" && length(rts)>0) {
			logger.start("Recomputing region aggregation")
			rtGrl <- .object@coord

			.object@coord <- .object@coord["sites"]
			.object@meth  <- .object@meth["sites"]
			.object@covg  <- .object@covg["sites"]

			for (rt in rts){
				logger.status(c(rt, "..."))
				.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
			}
			logger.completed()
		}
		return(.object)
	}
)

if (!isGeneric("normalizeMeth")) {
	setGeneric(
		"normalizeMeth",
		function(.object, ...) standardGeneric("normalizeMeth"),
		signature=c(".object")
	)
}
#' normalizeMeth-methods
#'
#' Normalize methylation levels
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param method  normalization method to be applied. Currently only 'quantile' is supported
#' @param type    character string specifying a name for the region type (default: sites)
#' @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsNOMe}} object with normalized methylation levels
#' 
#' @rdname normalizeMeth-DsNOMe-method
#' @docType methods
#' @aliases normalizeMeth
#' @aliases normalizeMeth,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("normalizeMeth",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type="sites",
		method="quantile",
		reaggregate=TRUE
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		if (!is.element(method, c("quantile"))) logger.error(c("Unsupported normalization method type:", method))

		if (method == "quantile"){
			logger.start(c("Performing quantile normalization"))
				require(preprocessCore)
				cnames <- colnames(.object@meth[[type]])
				.object@meth[[type]] <- data.table(normalize.quantiles(as.matrix(.object@meth[[type]])))
				colnames(.object@meth[[type]]) <- cnames
			logger.completed()
		}

		rts <- setdiff(getRegionTypes(.object), type)
		if (reaggregate && type == "sites" && length(rts)>0) {
			logger.start("Recomputing region aggregation")
			rtGrl <- .object@coord

			.object@coord <- .object@coord["sites"]
			.object@meth  <- .object@meth["sites"]
			.object@covg  <- .object@covg["sites"]

			for (rt in rts){
				logger.status(c(rt, "..."))
				.object <- regionAggregation(.object, rtGrl[[rt]], rt, dropEmpty=TRUE)
			}
			logger.completed()
		}
		return(.object)
	}
)

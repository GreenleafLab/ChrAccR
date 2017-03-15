#' DsNOMe
#'
#' A class for storing information on pipeline jobs.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{coord}}{
#'		List of coordinates (GRanges objects) for GC dinucleotides (sites) and summarized regions.
#'   }
#'   \item{\code{meth}}{
#'		List of GC methylation for sites and summarized regions.
#'   }
#'   \item{\code{covg}}{
#'		List of GC read coverage for sites and summarized regions.
#'   }
#'   \item{\code{sampleAnnot}}{
#'		Sample annotation Table
#'   }
#'   \item{\code{genome}}{
#'		Genome assembly
#'   }
#' }
#'
#' @section Methods:
#' \describe{
#'    \item{\code{\link{samples,DsNOMe-method}}}{
#'      Retrieve a vector of sample IDs for the dataset
#'    }
#' }
#'
#' @name DsNOMe-class
#' @rdname DsNOMe-class
#' @author Fabian Mueller
#' @exportClass DsNOMe
setClass("DsNOMe",
	slots = list(
		coord       = "list",
		meth        = "list",
		covg        = "list",
		sampleAnnot = "data.frame",
		genome      = "character"
	),
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
if (!isGeneric("getSamples")) {
	setGeneric(
		"getSamples",
		function(.object) standardGeneric("getSamples"),
		signature=c(".object")
	)
}
#' getSamples-methods
#'
#' Return sample IDs in a dataset
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @return Character vector of sample IDs in the dataset
#'
#' @rdname getSamples-DsNOMe-method
#' @docType methods
#' @aliases getSamples
#' @aliases getSamples,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getSamples",
	signature(
		.object="DsNOMe"
	),
	function(
		.object
	) {
		return(rownames(.object@sampleAnnot))
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getSampleAnnot")) {
	setGeneric(
		"getSampleAnnot",
		function(.object) standardGeneric("getSampleAnnot"),
		signature=c(".object")
	)
}
#' getSampleAnnot-methods
#'
#' Return sample annotation table of a dataset
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @return \code{data.frame} containing sample annotation
#'
#' @rdname getSampleAnnot-DsNOMe-method
#' @docType methods
#' @aliases getSampleAnnot
#' @aliases getSampleAnnot,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getSampleAnnot",
	signature(
		.object="DsNOMe"
	),
	function(
		.object
	) {
		return(.object@sampleAnnot)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getRegionTypes")) {
	setGeneric(
		"getRegionTypes",
		function(.object, ...) standardGeneric("getRegionTypes"),
		signature=c(".object")
	)
}
#' getRegionTypes-methods
#'
#' Return sample IDs in a dataset
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param inclSites include \code{"sites"} in the result
#' @return Character vector of sample IDs in the dataset
#'
#' @rdname getRegionTypes-DsNOMe-method
#' @docType methods
#' @aliases getRegionTypes
#' @aliases getRegionTypes,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getRegionTypes",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		inclSites=FALSE
	) {
		res <- names(.object@coord)
		if (!inclSites) res <- setdiff(res, "sites")
		return(res)
	}
)
#-------------------------------------------------------------------------------
if (!isGeneric("getCoord")) {
	setGeneric(
		"getCoord",
		function(.object, ...) standardGeneric("getCoord"),
		signature=c(".object")
	)
}
#' getCoord-methods
#'
#' Return coordinates of sites/regions in a dataset
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @param asMatrix return a matrix instead of a \code{data.table}
#' @return \code{GRanges} object containing coordinates for covered
#'         sites/regions
#'
#' @rdname getCoord-DsNOMe-method
#' @docType methods
#' @aliases getCoord
#' @aliases getCoord,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getCoord",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type="sites"
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		res <- .object@coord[[type]]
		return(res)
	}
)
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
if (!isGeneric("getNRegions")) {
	setGeneric(
		"getNRegions",
		function(.object, ...) standardGeneric("getNRegions"),
		signature=c(".object")
	)
}
#' getNRegions-methods
#'
#' Return the number of regions of a given type
#'
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @return the number of regions of that type
#'
#' @rdname getNRegions-DsNOMe-method
#' @docType methods
#' @aliases getNRegions
#' @aliases getNRegions,DsNOMe-method
#' @author Fabian Mueller
#' @export
setMethod("getNRegions",
	signature(
		.object="DsNOMe"
	),
	function(
		.object,
		type="sites"
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		return(length(.object@coord[[type]]))
	}
)
#-------------------------------------------------------------------------------

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
		.object
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
			
			if (length(getRegionTypes(.object))>0){
				logger.warning("Detected multiple region types. These will be dropped and not be recalculated")
				.object@coord <- .object@coord["sites"]
				.object@meth  <- .object@meth["sites"]
				.object@covg  <- .object@covg["sites"]
			}
		}
		return(.object)
	}
)
#-------------------------------------------------------------------------------

################################################################################
# Saving and loading DsNOMe objects
################################################################################
#' saveDsNOMe
#' 
#' Save a DsNOMe dataset to disk for later loading
#' @param .object \code{\linkS4class{DsNOMe}} object
#' @param path    destination to save the object to
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @export
saveDsNOMe <- function(.object, path){
	if (dir.exists(path)){
		logger.error("could not save object. Path already exists")
	}
	dir.create(path, recursive=FALSE)
	dsFn <- file.path(path, "ds.rds")
	saveRDS(.object, dsFn)
	invisible(NULL)
}

#' loadDsNOMe
#' 
#' Load a DsNOMe dataset from disk
#' @param path    Location of saved \code{\linkS4class{DsNOMe}} object
#' @return \code{\linkS4class{DsNOMe}} object
#' @author Fabian Mueller
#' @export
loadDsNOMe <- function(path){
	if (!dir.exists(path)){
		logger.error(c("Could not load object. Path does not exist:", path))
	}
	dsFn <- file.path(path, "ds.rds")
	.object <- readRDS(dsFn)
	return(.object)
}

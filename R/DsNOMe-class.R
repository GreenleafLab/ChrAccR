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
		if (!is.element(getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
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
		if (!is.element(getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
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
		if (!is.element(getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		res <- .object@covg[[type]]
		if (asMatrix) res <- as.matrix(res)
		return(res)
	}
)
#-------------------------------------------------------------------------------

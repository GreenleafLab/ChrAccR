#' DsAcc
#'
#' A class for accessibility datasets
#'
#' @section Slots:
#' \describe{
#'   \item{\code{coord}}{
#'		List of coordinates (GRanges objects) for accessibility summarized regions.
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
#'    \item{\code{\link{samples,DsAcc-method}}}{
#'      Retrieve a vector of sample IDs for the dataset
#'    }
#' }
#'
#' @name DsAcc-class
#' @rdname DsAcc-class
#' @author Fabian Mueller
#' @exportClass DsAcc
setClass("DsAcc",
	slots = list(
		coord       = "list",
		sampleAnnot = "data.frame",
		genome      = "character",
		diskDump    = "logical",
		pkgVersion  = "ANY"
	),
	package = "ChrAccR"
)
setMethod("initialize","DsAcc",
	function(
		.Object,
		coord,
		sampleAnnot,
		genome,
		diskDump
	) {
		if (diskDump){
			require(DelayedArray)
			require(HDF5Array)
		}
		.Object@coord       <- coord
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object@diskDump    <- diskDump
		.Object@pkgVersion  <- packageVersion("ChrAccR")
		.Object
	}
)

#' @param siteCoord \code{GRanges} object containing coordinates of GC dinucleotides
#' @param sampleAnnot \code{data.frame} object containing sample annotation
#' @param genome    character string containing genome assembly
#' @param diskDump  should large matrices be stored on disk rather than in main memory
#' @noRd
DsAcc <- function(siteCoord, sampleAnnot, genome, diskDump=FALSE){
	obj <- new("DsAcc",
		siteCoord,
		sampleAnnot,
		genome,
		diskDump
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @return Character vector of sample IDs in the dataset
#'
#' @rdname getSamples-DsAcc-method
#' @docType methods
#' @aliases getSamples
#' @aliases getSamples,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getSamples",
	signature(
		.object="DsAcc"
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @return \code{data.frame} containing sample annotation
#'
#' @rdname getSampleAnnot-DsAcc-method
#' @docType methods
#' @aliases getSampleAnnot
#' @aliases getSampleAnnot,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getSampleAnnot",
	signature(
		.object="DsAcc"
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param inclSites include \code{"sites"} in the result
#' @return Character vector of sample IDs in the dataset
#'
#' @rdname getRegionTypes-DsAcc-method
#' @docType methods
#' @aliases getRegionTypes
#' @aliases getRegionTypes,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getRegionTypes",
	signature(
		.object="DsAcc"
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @param asMatrix return a matrix instead of a \code{data.table}
#' @return \code{GRanges} object containing coordinates for covered
#'         sites/regions
#'
#' @rdname getCoord-DsAcc-method
#' @docType methods
#' @aliases getCoord
#' @aliases getCoord,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getCoord",
	signature(
		.object="DsAcc"
	),
	function(
		.object,
		type
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		res <- .object@coord[[type]]
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param type    character string specifying the rgion type or \code{"sites"} (default)
#' @return the number of regions of that type
#'
#' @rdname getNRegions-DsAcc-method
#' @docType methods
#' @aliases getNRegions
#' @aliases getNRegions,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getNRegions",
	signature(
		.object="DsAcc"
	),
	function(
		.object,
		type="sites"
	) {
		if (!is.element(type, getRegionTypes(.object, inclSites=TRUE))) logger.error(c("Unsupported region type:", type))
		return(length(.object@coord[[type]]))
	}
)

################################################################################
# Display
################################################################################
setMethod("show","DsAcc",
	function(object) {
		ss <- getSamples(object)
		str.ss <- paste(getSamples(object), collapse=", ")
		if (length(ss) > 5) str.ss <- paste(c(getSamples(object)[1:5], "..."), collapse=", ")
		rts <- getRegionTypes(object)
		str.rts <- "no region types"
		if (length(rts) > 0) str.rts <- paste0(length(rts), " region types: ", paste(rts, collapse=", "))

		cat("DsAcc chromatin accessibility dataset \n")
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
# Maniputlating DsAcc objects
################################################################################
if (!isGeneric("addSampleAnnotCol")) {
	setGeneric(
		"addSampleAnnotCol",
		function(.object, ...) standardGeneric("addSampleAnnotCol"),
		signature=c(".object")
	)
}
#' addSampleAnnotCol-methods
#'
#' add a sample annotation column to the sample annotation table
#'
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param name    a name for the new column
#' @param vals    vector of values
#' @return a new \code{\linkS4class{DsAcc}} object with added sample annotation
#' 
#' @rdname addSampleAnnotCol-DsAcc-method
#' @docType methods
#' @aliases addSampleAnnotCol
#' @aliases addSampleAnnotCol,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("addSampleAnnotCol",
	signature(
		.object="DsAcc"
	),
	function(
		.object,
		name,
		vals
	) {
		ph <- .object@sampleAnnot
		if (length(vals)!=ncol(ph)){
			logger.error(c("vals must contain exactly one value for each sample"))
		}
		if (is.element(name, colnames(ph))){
			logger.warning(c("Replacing sample annotation column:", name))
		}
		
		ph[,name] <- vals

		.object@sampleAnnot <- ph
		return(.object)
	}
)
#-------------------------------------------------------------------------------
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
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param indices a vector of indices of sites/regions to be removed. Can be numeric, integer or logical.
#' @param type    character string specifying a name for the region type (sefault: sites)
#' @param reaggregate redo region aggregation (only has an effect if type is sites and there are aggregated regions in the dataset)
#' @return a new \code{\linkS4class{DsAcc}} object with sites/regions removed
#' 
#' @rdname removeRegions-DsAcc-method
#' @docType methods
#' @aliases removeRegions
#' @aliases removeRegions,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("removeRegions",
	signature(
		.object="DsAcc"
	),
	function(
		.object,
		indices,
		type
	) {
		inclSites <- FALSE
		if (type == "sites") inclSites <- TRUE
		if (!is.element(type, getRegionTypes(.object, inclSites=inclSites))) logger.error(c("Unsupported region type:", type))

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
		return(.object)
	}
)
################################################################################
# Saving and loading DsAcc objects
################################################################################
#' saveDsAcc
#' 
#' Save a DsAcc dataset to disk for later loading
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param path    destination to save the object to
#' @param forceDiskDump force large matrices (counts) to be stored as HDF5 (even when the object was not created using \code{diskDump=TRUE})
#' @param updateDiskRef update disk dumped (HDF5) references (e.g. for count data)
#' @return (invisibly) The object (with potentially updated disk dumped references)
#' @author Fabian Mueller
#' @export
saveDsAcc <- function(.object, path, forceDiskDump=FALSE, updateDiskRef=TRUE){
	if (dir.exists(path)){
		logger.error("could not save object. Path already exists")
	}
	dir.create(path, recursive=FALSE)

	# save region count data as HDF5
	if (.hasSlot(.object, "counts") && !is.null(.object@counts) && length(.object@counts) > 0){
		if (forceDiskDump || (.hasSlot(.object, "diskDump") && .object@diskDump)){
			require(DelayedArray)
			require(HDF5Array)
			logger.start("Saving region count data to HDF5")
				countDir <- file.path(path, "countData")
				dir.create(countDir)
				for (i in 1:length(.object@counts)) {
					rt <- names(.object@counts)[i]
					logger.status(c("Region type:", rt))
					# just an assertion to make sure the colnames correspond to the sample ids
					if (!all(colnames(.object@counts[[rt]])==getSamples(.object))) logger.error("Assertion failed: column names do not correspond to sample names (counts)")
					if (updateDiskRef){
						.object@counts[[rt]] <- writeHDF5Array(.object@counts[[rt]], filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
					} else {
						dummy <- writeHDF5Array(.object@counts[[rt]], filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
					}
				}
			logger.completed()
			.object@diskDump <- TRUE
		}
	}

	dsFn <- file.path(path, "ds.rds")
	saveRDS(.object, dsFn)

	invisible(.object)
}

#' loadDsAcc
#' 
#' Load a DsAcc dataset from disk
#' @param path    Location of saved \code{\linkS4class{DsAcc}} object
#' @return \code{\linkS4class{DsAcc}} object
#' @author Fabian Mueller
#' @export
loadDsAcc <- function(path){
	if (!dir.exists(path)){
		logger.error(c("Could not load object. Path does not exist:", path))
	}
	dsFn <- file.path(path, "ds.rds")
	.object <- readRDS(dsFn)

	# load region count data from HDF5
	if (.hasSlot(.object, "diskDump") && .object@diskDump && .hasSlot(.object, "counts") && !is.null(.object@counts) && length(.object@counts) > 0){
		require(DelayedArray)
		require(HDF5Array)
		logger.start("Loading region count data from HDF5")
			countDir <- file.path(path, "countData")
			for (i in 1:length(.object@counts)) {
				rt <- names(.object@counts)[i]
				logger.status(c("Region type:", rt))
				.object@counts[[rt]] <- HDF5Array(filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
				colnames(.object@counts[[rt]]) <- getSamples(.object) # reset the column names (workaround for the issue that writeHDF5Array does not write dimnames to HDF5)
			}
		logger.completed()
	}
	return(.object)
}

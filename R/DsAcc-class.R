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
#'   \item{\code{diskDump}}{
#'		Flag indicating whether large matrices and objects will be kept on disk rather than in main memory.
#'   }
#'   \item{\code{pkgVersion}}{
#'		Version number of the ChrAccR package that created the object
#'   }
#' }
#'
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
			if (!requireNamespace("DelayedArray") || !requireNamespace("HDF5Array")) logger.error(c("Could not load dependency: DelayedArray, HDF5Array"))
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
#' Retrieve the number of samples contained in a DsAcc object
#'
#' @param x DsAcc object
setMethod("length", "DsAcc",
	function(x){
		return(length(getSamples(x)))
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
if (!isGeneric("getGenome")) {
	setGeneric(
		"getGenome",
		function(.object) standardGeneric("getGenome"),
		signature=c(".object")
	)
}
#' getGenome-methods
#'
#' Return the genome assembly
#'
#' @param .object \code{\linkS4class{DsAcc}} object
#' @return Character string containing the genome assembly
#'
#' @rdname getGenome-DsAcc-method
#' @docType methods
#' @aliases getGenome
#' @aliases getGenome,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getGenome",
	signature(
		.object="DsAcc"
	),
	function(
		.object
	) {
		return(.object@genome)
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
		if (length(vals)!=nrow(ph)){
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
			logger.start("Saving region count data to HDF5")
				countDir <- file.path(path, "countData")
				dir.create(countDir)
				for (i in 1:length(.object@counts)) {
					rt <- names(.object@counts)[i]
					logger.status(c("Region type:", rt))
					sampleNames <- getSamples(.object)
					# just an assertion to make sure the colnames correspond to the sample ids
					if (!all(colnames(.object@counts[[rt]])==sampleNames)) logger.error("Assertion failed: column names do not correspond to sample names (counts)")
					if (updateDiskRef){
						.object@counts[[rt]] <- HDF5Array::writeHDF5Array(.object@counts[[rt]], filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
						colnames(.object@counts[[rt]]) <- sampleNames # reset the column names (workaround for the issue that writeHDF5Array does not write dimnames to HDF5)
					} else {
						dummy <- HDF5Array::writeHDF5Array(.object@counts[[rt]], filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
					}
				}
			logger.completed()
			.object@diskDump <- TRUE
		}
	}
	# save fragment data as RDS
	if (.hasSlot(.object, "fragments") && !is.null(.object@fragments) && length(.object@fragments) > 0){
		if (forceDiskDump || (.hasSlot(.object, "diskDump.fragments") && .object@diskDump.fragments)){
			logger.start("Saving fragment data to RDS")
				fragDir <- file.path(path, "fragments")
				dir.create(fragDir)
				nSamples <- length(.object@fragments)
				chunkL <- lapply(1:nSamples, identity)
				chunkedFragmentFiles <- .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
				if (chunkedFragmentFiles){
					createNewChunkL <- TRUE
					isDd <- sapply(.object@fragments, is.character)
					# all disk-dumped?
					if (all(isDd)){
						ddFns <- unlist(.object@fragments)
						chunkL <- tapply(1:length(ddFns), ddFns, identity)
						# check if already chunked
						createNewChunkL <- all(elementNROWS(chunkL)==1)
					}
					if (createNewChunkL){
						chunkL <- split(1:nSamples, rep(1:ceiling(nSamples/.object@diskDump.fragments.nSamplesPerFile), each=.object@diskDump.fragments.nSamplesPerFile)[1:nSamples])
						names(chunkL) <- NULL
					}
				}
				for (k in 1:length(chunkL)) {
					fn <- file.path(fragDir, paste0("fragmentGr_", k, ".rds"))
					iis <- chunkL[[k]]
					fragGrl <- .object@fragments[chunkL[[k]]]
					fn_source <- ""
					isFile <- sapply(fragGrl, is.character)
					isMix <- length(unique(isFile)) != 1
					if (isMix) logger.error("Mix of disk dumped and in-memory fragment objects not supported yet")
					isFileReady <- all(isFile)
					if (isFileReady) {
						# check if all samples in the current chunk are in the same file
						uu <- unique(unlist(fragGrl))
						if (length(uu)==1) {
							fn_source <- uu
						} else {
							logger.error("Repackaging of already packaged fragments not supported yet")
						}
					}
					if (isFileReady){
						if (file.exists(fn_source)) {
							file.copy(fn_source, fn)
						} else {
							logger.error(c("Could not find fragment data file:", fn_source))
						}
					} else {
						saveRDS(fragGrl, fn, compress=TRUE)
					}
					if (updateDiskRef){
						for (i in chunkL[[k]]){
							.object@fragments[[i]] <- fn
						}
					}
				}
			logger.completed()
			.object@diskDump.fragments <- TRUE
		}
	}

	dsFn <- file.path(path, "ds.rds")
	saveRDS(.object, dsFn, compress=TRUE)

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
		logger.start("Loading region count data from HDF5")
			countDir <- file.path(path, "countData")
			for (i in 1:length(.object@counts)) {
				rt <- names(.object@counts)[i]
				logger.status(c("Region type:", rt))
				.object@counts[[rt]] <- HDF5Array::HDF5Array(filepath=file.path(countDir, paste0("regionCounts_", i, ".h5")), name=paste0("count_hdf5_", rt))
				colnames(.object@counts[[rt]]) <- getSamples(.object) # reset the column names (workaround for the issue that writeHDF5Array does not write dimnames to HDF5)
			}
		logger.completed()
	}
	# load fragment data from RDS
	if (.hasSlot(.object, "diskDump.fragments") && .object@diskDump.fragments && .hasSlot(.object, "fragments") && !is.null(.object@fragments) && length(.object@fragments) > 0){
		logger.start("Updating fragment RDS file references")
			nSamples <- length(.object@fragments)
			fragDir <- file.path(path, "fragments")
			
			chunkedFragmentFiles <- .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
			if (!chunkedFragmentFiles){
				for (i in 1:nSamples) {
					fn <- file.path(fragDir, paste0("fragmentGr_", i, ".rds"))
					if (!file.exists(fn)) logger.error(paste0("Invalid save: Could not find fragment file:", fn))
					.object@fragments[[i]] <- fn
				}
			} else {
				fragFns <- list.files(fragDir, pattern="fragmentGr_")
				for (fn in fragFns){
					fn_full <- file.path(fragDir, fn)
					fGrl <- readRDS(fn_full)
					if (!all(names(fGrl) %in% names(.object@fragments))) logger.error(c("Incompatible fragment file: could not find sample names from file:", fn))
					.object@fragments[names(fGrl)] <- rep(list(fn_full), length(fGrl))
				}
			}
		logger.completed()
	}
	return(.object)
}

################################################################################
# Retrieving differential comparison info
################################################################################
if (!isGeneric("getComparisonTable")) {
	setGeneric(
		"getComparisonTable",
		function(.object, ...) standardGeneric("getComparisonTable"),
		signature=c(".object")
	)
}
#' getComparisonTable-methods
#'
#' Retrieve a table describing pairwise comparisons on a \code{\linkS4class{DsAcc}} object
#'
#' @param .object \code{\linkS4class{DsAcc}} object
#' @param cols    column names in the sample annotation table to consider for pairwise comparisons
#' @param cols1vAll  column names in the sample annotation table to consider for 1-vs-all comparisons
#' @param compNames  vector of character strings specifying a fixed comparison names to be parsed (format "$GRP1_NAME vs $GRP1_NAME [$ANNOTATION_COLUMN]")
#' @param minGroupSize Minimum size of a group to be used in comparison. Affects the annotation columns that will be used for comparisons.
#' @return a \code{data.frame} with comparison inforamtion containing columns for the comparison name (\code{compName}), 
#'         column in the annotation table (\code{compCol})
#'         and group names for the two groups in the comparison (\code{grp1Name, grp2Name}),
#' 
#' @rdname getComparisonTable-DsAcc-method
#' @docType methods
#' @aliases getComparisonTable
#' @aliases getComparisonTable,DsAcc-method
#' @author Fabian Mueller
#' @export
setMethod("getComparisonTable",
	signature(
		.object="DsAcc"
	),
	function(
		.object,
		cols=NULL,
		cols1vAll=NULL,
		compNames=NULL,
		minGroupSize=2L
	) {
		colsAdd <- NULL
		fixedCompInfo <- NULL
		# parse fixed comparison names
		if (!is.null(compNames)){
			re <- "^(.+) vs (.+) \\[(.+)\\]$"
			isMatch <- grepl(re, compNames)
			if (any(!isMatch)){
				logger.error(c("The following comparison names could not be parsed:", paste(compNames[!isMatch], collapse=", ")))
			}
			fixedCompInfo <- data.frame(
				compName=compNames,
				compCol=gsub(re, "\\3", compNames),
				# grp1Name=gsub(re, "\\1", compNames),
				# grp2Name=gsub(re, "\\2", compNames),
				stringsAsFactors=FALSE
			)
			colsAdd <- unique(fixedCompInfo[,"compCol"])
			# check for 1-vs-all comparisons
			is1vsAll <- grepl("^\\.ALL", compNames) | grepl(" vs \\.ALL", compNames)
			if (any(is1vsAll)){
				cols1vAll <- union(cols1vAll, unique(fixedCompInfo[is1vsAll,"compCol"]))
			}
		}

		# get comparison info
		sannot <- getSampleAnnot(.object)
		sampleGrps <- getGroupsFromTable(sannot, cols=unique(c(colsAdd, cols, cols1vAll)), minGrpSize=minGroupSize)
		if (length(sampleGrps) < 1) logger.error("No valid comparisons found (to begin with)")
		compTab <- do.call("rbind", lapply(1:length(sampleGrps), FUN=function(i){
			tt <- NULL
			grpNs <- sapply(sampleGrps[[i]], length)
			names(grpNs) <- names(sampleGrps[[i]])
			if (length(sampleGrps[[i]]) == 2) {
				tt <- data.frame(
					compName=paste0(names(sampleGrps[[i]])[1], " vs ", names(sampleGrps[[i]])[2],  " [", names(sampleGrps)[i], "]"),
					compCol=names(sampleGrps)[i],
					grp1Name=names(sampleGrps[[i]])[1],
					grp2Name=names(sampleGrps[[i]])[2],
					stringsAsFactors=FALSE
				)
			} else if (length(sampleGrps[[i]]) > 2) {
				if (is.element(names(sampleGrps)[i], cols1vAll)){
					gns <- names(sampleGrps[[i]])
					tt <- data.frame(
						compName=paste0(gns, " vs ", ".ALL",  " [", names(sampleGrps)[i], "]"),
						compCol=names(sampleGrps)[i],
						grp1Name=gns,
						grp2Name=".ALL",
						stringsAsFactors=FALSE
					)
				} else {
					grpNames <- t(combn(names(sampleGrps[[i]]), 2))
					tt <- data.frame(
						compName=paste0(grpNames[,1], " vs ", grpNames[,2],  " [", names(sampleGrps)[i], "]"),
						compCol=names(sampleGrps)[i],
						grp1Name=grpNames[,1],
						grp2Name=grpNames[,2],
						stringsAsFactors=FALSE
					)
				}
			}
			tt[,"nGrp1"] <- grpNs[tt[,"grp1Name"]]
			tt[,"nGrp2"] <- sapply(1:nrow(tt), FUN=function(i){
				if (tt[i,"grp2Name"] == ".ALL"){
					return(sum(grpNs[names(grpNs)!=tt[i,"grp2Name"]]))
				} else {
					return(grpNs[tt[i,"grp2Name"]])
				}
			})
			
			return(tt)
		}))
		if (is.null(compTab)) logger.error("No valid comparisons found")

		# add comparison info for fixed comparisons
		if (!is.null(compNames)){
			compTabMirror <- do.call("rbind", lapply(1:nrow(compTab), FUN=function(i){
				tt <- data.frame(
					compName=paste0(compTab[i,"grp2Name"], " vs ", compTab[i,"grp1Name"],  " [", compTab[i,"compCol"], "]"),
					compCol=compTab[i,"compCol"],
					grp1Name=compTab[i,"grp2Name"],
					grp2Name=compTab[i,"grp1Name"],
					nGrp1=compTab[i,"nGrp2"],
					nGrp2=compTab[i,"nGrp1"],
					stringsAsFactors=FALSE
				)
			}))
			matchIdx <- match(fixedCompInfo[,"compName"], compTab[,"compName"])
			matchIdxMirror <- match(fixedCompInfo[,"compName"], compTabMirror[,"compName"])
			noMatch <- is.na(matchIdx) & is.na(matchIdxMirror)
			if (any(noMatch)){
				logger.error(c("The following comparison names could not be matched:", paste(fixedCompInfo[noMatch], collapse=", ")))
			}
			fixedCompTab <- rbind(compTab[na.omit(matchIdx),], compTabMirror[na.omit(matchIdxMirror),])

			compTabIdx.dupRem <- !(compTab[,"compCol"] %in% colsAdd) # don't duplicate the columns for which fixed comparisons have been specified
			compTab <- rbind(
				fixedCompTab,
				compTab[compTabIdx.dupRem,]
			)
		}
		return(compTab)
	}
)

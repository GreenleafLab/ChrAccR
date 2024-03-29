################################################################################
# Configuration and global package options
################################################################################
.config <- new.env()
.config$tmpDir <- tempdir()
.config$cleanMem <- TRUE
.config$colorSchemes <- list(
	".default" = c("#009FE3", "#DE7E00", "#8EC041", "#FFCC00", "#951B81", "#BE1716", "#7C83B3", "#671719", "#E0CDA6", "#775725", "#000000")
)
.config$colorSchemesCont <- list(
	".default" = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF"),
	".default.div" = c("#01665E", "#35978F", "#80CDC1", "#C7EAE5", "#F5F5F5", "#F6E8C3", "#DFC27D", "#BF812D", "#8C510A"),
	# ".default.geneactivity" = c("#e0f3db", "#a8ddb5", "#4eb3d3", "#08589e")
	".default.geneactivity" = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D') # solarextra
)
.config$geneModelVersions <- c(
	"hg38"="gencode.v27",
	"hg19"="gencode.v19",
	"mm10"="gencode.vM16",
	"mm9" ="gencode.vM1"
)
.config$analysisName <- "ChrAccR analysis"
.config$regionTypes <- NULL
.config$chromVarMotifs <- c("jaspar_vert")
.config$chromVarRegionTypes <- NULL
.config$chromVarMotifNamesForDimRed <- NULL
.config$genesOfInterest <- NULL
.config$annotationColumns <- NULL
.config$annotationMinGroupSize <- 2L
.config$annotationMaxGroupCount <- NULL
.config$doPeakCalling <- FALSE
.config$peakCallingProfile <- NULL
.config$addCellRangerPeaks <- TRUE
.config$annotationPeakGroupColumn <- NULL
.config$annotationPeakGroupAgreePerc <- 1
.config$filteringCovgCount <- 1L
.config$filteringCovgReqSamples <- 0.75
.config$filteringSexChroms <- FALSE
.config$filteringScMinFragmentsPerCell <- 1000L
.config$filteringScMaxFragmentsPerCell <- Inf
.config$filteringScMinTssEnrichment <- 6
.config$normalizationMethod <- "quantile"
.config$exploratoryLogNormCounts <- TRUE
.config$exploratoryNSubsample <- 2e6
.config$differentialColumns <- NULL
.config$differentialColumns1vsAll <- NULL
.config$differentialCompNames <- NULL
.config$differentialAdjColumns <- NULL
.config$differentialCutoffL2FC <- 2
.config$lolaDbPaths <- NULL
.config$scIterativeLsiRegType <- NULL
.config$scIterativeLsiParams <- list(
	it0clusterResolution = 0.4,
	it1clusterResolution = 0.4,
	it2clusterResolution = 0.4,
	umapParams = list(distMethod="euclidean", min_dist=0.5, n_neighbors=25)
)
.config$scGeneActivity <- FALSE
.config$muPipeR_cmdr <- NULL

#' setConfigElement
#'
#' Set a configuration item to a given value
#'
#' @param name	name of the config item
#' @param value	value of the config item
#' @return nothing of particular interest.
#' @section Options used by the package:
#' \describe{
#'   \item{\bold{\code{tmpDir}}\code{ = temdir()}}{
#'        Directory for temporary files. Must be existing.
#'   }
#'   \item{\bold{\code{cleanMem}}\code{ = TRUE}}{
#'        During runtime, regularly clean-out the memory in order to reduce memory overuse
#'   }
#'   \item{\bold{\code{colorSchemes}}}{
#'       named \code{list} of DISCRETE color schemes to be used for plotting. Each element should be a named vector specifying colors for groups/annotations.
#'   }
#'   \item{\bold{\code{colorSchemesCont}}}{
#'       named \code{list} of CONTINOUS color schemes to be used for plotting. Each element should be a vector specifying a range of colors.
#'   }
#'   \item{\bold{\code{geneModelVersions}}}{
#'       Gene model versions to be used for various genomes
#'   }
#'   \item{\bold{\code{analysisName}}\code{ = "ChrAccR analysis"}}{
#'       A title for the analysis (a string).
#'   }
#'   \item{\bold{\code{regionTypes}}}{
#'       Region types to be used in the analysis
#'   }
#'   \item{\bold{\code{chromVarRegionTypes}}\code{ = NULL}}{
#'       Region types to be used for chromVar analysis. If \code{NULL} (default), ChrAccR will automatically look for region types with the keyword \code{"peak"} in their name.
#'   }
#' 	 \item{\bold{\code{chromVarMotifs}}\code{ = "jaspar_vert"}}{
#'       Character vector of names of TF motif sets to be used in ChromVAR analyses. By default the vertebrate set of the JASPAR database will be used.
#'   }
#'   \item{\bold{\code{chromVarMotifNamesForDimRed}}}{
#'       Names of motifs to be used for dimension reduction plots in the reports. [only relevant for single-cell data]
#'   }
#'   \item{\bold{\code{genesOfInterest}}}{
#'       Names of genes of interest to be highlighted in the reports (e.g. dimension reduction) in the reports.
#'       [currently only relevant for single-cell data and only when \code{scGeneActivity} is activated]
#'   }
#'   \item{\bold{\code{annotationColumns}}}{
#'       Sample annotation columns to be used for reporting
#'   }
#'   \item{\bold{\code{annotationMinGroupSize}}\code{ = 2}}{
#'       Minimum size of a group to be used in the reports. Influences which annotation columns are automatically selected for reporting.
#'   }
#'   \item{\bold{\code{annotationMaxGroupCount}}\code{ = NULL}}{
#'       Maximum number of groups to be used in the reports. Influences which annotation columns are
#'       automatically selected for reporting. If \code{NULL} (default) it will effectively be the
#'       number of samples - 1.
#'   }
#'   \item{\bold{\code{doPeakCalling}}\code{ = FALSE}}{
#'       Perform per-sample peak calling and retrieve consensus peak set. Requires that \code{macs2} is installed and can be called from the command line. [for bulk data analysis only]
#'   }
#'   \item{\bold{\code{peakCallingProfile}}\code{ = NULL}}{
#'       If set to a string describing a valid profile, will apply a special profile for \code{macs2} peak calling.
#'       [only valid in combination with the \code{doPeakCalling} option]
#'   }
#'   \item{\bold{\code{annotationPeakGroupColumn}}}{
#'       Annotation column to base the consensus peak set replication filtering on.
#'   }
#'   \item{\bold{\code{annotationPeakGroupAgreePerc}}\code{ = 1.0}}{
#'       Percent of samples that have to agree to identify consensus peaks. See \code{\link{getConsensusPeakSet}} for details.
#'   }
#'   \item{\bold{\code{filteringCovgCount}}\code{ = 1L}}{
#'       Minimum insertion count to filter count matrices by. See \code{\link{filterLowCovg,DsATAC-method}} for details. [for bulk data analysis only]
#'   }
#'   \item{\bold{\code{filteringCovgReqSamples}}\code{ = 0.75}}{
#'       Minimum required samples to apply low coverage filtering to. See \code{\link{filterLowCovg,DsATAC-method}} for details. [for bulk data analysis only]
#'   }
#'   \item{\bold{\code{filteringSexChroms}}\code{ = FALSE}}{
#'       Flag indicating whether to remove sex chromosomes.
#'   }
#'   \item{\bold{\code{filteringScMinFragmentsPerCell}}\code{ = 1000L}}{
#'       Minimum number of fragments per cell to retain a cell in the analysis. [for single-cell data analysis only]
#'   }
#'   \item{\bold{\code{filteringScMaxFragmentsPerCell}}\code{ = Inf}}{
#'       Maximum number of fragments allowed per cell to retain a cell in the analysis. [for single-cell data analysis only]
#'   }
#'   \item{\bold{\code{filteringScMinTssEnrichment}}\code{ = 6}}{
#'       Minimum TSS enrichment score per cell to retain a cell in the analysis. [for single-cell data analysis only]
#'   }
#'   \item{\bold{\code{normalizationMethod}}\code{ = "quantile"}}{
#'       Normalization method to use for count normalization. Allowed methods include the ones listed in \code{\link{transformCounts,DsATAC-method}}. [for bulk data analysis only]
#'   }
#'   \item{\bold{\code{exploratoryLogNormCounts}}\code{ = TRUE}}{
#'       Should a log-normalization be applied in the exploratory plot sections of the reports (dimension reduction, heatmaps)
#'   }
#'   \item{\bold{\code{exploratoryNSubsample}}\code{ = 2e6}}{
#'       Number of regions to subsample in exploratory analysis in order to increase computational performance.
#'   }
#'   \item{\bold{\code{differentialColumns}}}{
#'       Sample annotation columns to be used for differential testing and reporting
#'   }
#'   \item{\bold{\code{differentialColumns1vsAll}}}{
#'       Sample annotation columns to be used for differential testing and reporting in a 1-vs-all group setting. Should be a subset of \code{differentialColumns}.
#'   }
#'   \item{\bold{\code{differentialCompNames}}}{
#'       Comparison names from which comparison information is derived. Must be in the format of "$GRP1_NAME vs $GRP2_NAME [$ANNOTATION_COLUMN]".
#'   }
#'   \item{\bold{\code{differentialAdjColumns}}}{
#'       Sample annotation columns to be adjusted for in differential testing
#'   }
#'   \item{\bold{\code{differentialCutoffL2FC}}}{
#'       Cutoff on log2 fold-change to be used for reporting differential accessibility.
#'   }
#'   \item{\bold{\code{lolaDbPaths}}}{
#'       Precomputed LOLA databases to be used for enrichment analysis. If \code{NULL} (default), ChrAccR will download an apropriate core database.
#'   }
#'   \item{\bold{\code{scIterativeLsiRegType}}}{
#'       For single-cell analysis only: region type to be used for clustering and dimension reduction using iterative LSI. By default (\code{NULL}),
#'       ChrAccR will look for a region type named \code{"tiling"}.
#'   }
#'   \item{\bold{\code{scIterativeLsiParams}}}{
#'       Parameters to use for iterative LSI. See \code{\link{iterativeLSI,DsATACsc-method}} for details.
#'   }
#'   \item{\bold{\code{scGeneActivity}}\code{ = FALSE}}{
#'       For single-cell analysis only: Compute gene activity from accessibility.
#'       Possible options are \code{"RBF"} for radial-basis-function-weighted count aggregation (default when set to \code{TRUE}) or 
#'       \code{"Cicero"} for Cicero correlation-based aggregation
#'   }
#' }
#' @author Fabian Mueller
#' @export
setConfigElement <- function(name, value){
	if (!exists(name, .config)){
		logger.error(c("No such configuration element:", name))
	}
	# TODO: implement option checker (especially for report-relevant options)
	.config[[name]] <- value
}

#' getConfigElement
#'
#' Get the value for a configuration item
#'
#' @param name	name of the config item
#' @return the value of the config item
#' @author Fabian Mueller
#' @export
getConfigElement <- function(name){
	if (!exists(name, .config)){
		logger.warning(c("No such configuration element:", name, "--> NULL returned"))
	}
	.config[[name]]
}

# convert a named vector to list. If x is already a list, apply it recursively to all elements of x
v2l <- function(x){
	if (is.list(x)){
		return(lapply(x, v2l))
	} else {
		if (length(x) > 1 && !is.null(names(x))){
			return(as.list(x))
		} else {
			return(x)
		}
	}
}
# convert a list with one-element entries into a vector. If not all elements of x have length 1, apply it recursively to all of them
l2v <- function(x, ...){
	if (!is.list(x)) return(x)
	lls <- sapply(x, length)
	if (all(lls<2)){
		return(unlist(x, recursive=FALSE))
	} else {
		idx <- lls > 1
		x[idx] <- lapply(x[idx], l2v)
		return(x)
	}
}
#' saveConfig
#'
#' Save the current configuration to a configuration file (JSON)
#'
#' @param dest	Filename for the config file in JSON format.
#' @return nothing of particular interest.
#'
#' @author Fabian Mueller
#' @export
saveConfig <- function(dest){
	cfgL <- as.list(.config)
	# toJSON does not allow for named vectors (https://github.com/jeroen/jsonlite/issues/76). Here's a workaround
	namedVectors <- intersect(names(cfgL), c("geneModelVersions", "colorSchemes"))
	for (nn in namedVectors){
		cfgL[[nn]] <- v2l(cfgL[[nn]])
	}
	cat(jsonlite::toJSON(cfgL, pretty=TRUE, null="null",na="string"), file=dest)
}

#' loadConfig
#'
#' Sets the configuration from a configuration file (JSON)
#'
#' @param cfgFile	Config file in JSON format. As output by \link{saveConfig}
#' @return nothing of particular interest. The configuration is set for the current environment
#'
#' @author Fabian Mueller
#' @export
loadConfig <- function(cfgFile){
	# toJSON does not allow for named vectors (https://github.com/jeroen/jsonlite/issues/76). Here's a workaround
	namedVectors <- c("geneModelVersions", "colorSchemes")

	cfgList <- jsonlite::fromJSON(cfgFile)
	for (nn in names(cfgList)){
		if (is.element(nn, namedVectors)){
			cfgList[[nn]] <- l2v(cfgList[[nn]])
		}
		if (is.element(nn,ls(.config))){
			.config[[nn]] <- cfgList[[nn]]
		} else {
			logger.warning(c("Ignoring unknown item '",nn,"' in when loading configuration from file",cfgFile))
		}
	}
}

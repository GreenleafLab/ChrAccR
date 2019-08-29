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
	".default.div" = c("#01665E", "#35978F", "#80CDC1", "#C7EAE5", "#F5F5F5", "#F6E8C3", "#DFC27D", "#BF812D", "#8C510A")
)
.config$geneModelVersions <- c(
	"hg38"="gencode.v27",
	"hg19"="gencode.v19",
	"mm10"="gencode.vM16",
	"mm9" ="gencode.vM1"
)
.config$regionTypes <- NULL
.config$chromVarMotifs <- c("jaspar_vert")
.config$chromVarRegionTypes <- NULL
.config$chromVarMotifsForDimRed <- NULL
.config$annotationColumns <- NULL
.config$annotationMinGroupSize <- 2L
.config$exploratoryLogNormCounts <- TRUE
.config$differentialColumns <- NULL
.config$differentialCompNames <- NULL
.config$differentialAdjColumns <- NULL
.config$lolaDbPaths <- NULL
.config$scIterativeLsiRegType <- NULL

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
#'   \item{\bold{\code{regionTypes}}}{
#'       Region types to be used in the analysis
#'   }
#'   \item{\bold{\code{chromVarRegionTypes}}}{
#'       Region types to be used for chromVar analysis. If \code{NULL} (default), ChrAccR will automatically look for region types with the keyword \code{"peak"} in their name.
#'   }
#' 	 \item{\bold{\code{chromVarMotifs}}}{
#'       Character vector of names of TF motif sets to be used in ChromVAR analyses
#'   }
#'   \item{\bold{\code{annotationColumns}}}{
#'       Sample annotation columns to be used for reporting
#'   }
#'   \item{\bold{\code{annotationMinGroupSize}}}{
#'       Minimum size of a group to be used in the reports. Influences which annotation columns are selected for reporting.
#'   }
#'   \item{\bold{\code{exploratoryLogNormCounts}}}{
#'       Should a log-normalization be applied in the exploratory plot sections of the reports (dimension reduction, heatmaps)
#'   }
#'   \item{\bold{\code{differentialColumns}}}{
#'       Sample annotation columns to be used for differential testing and reporting
#'   }
#'   \item{\bold{\code{differentialCompNames}}}{
#'       Comparison names from which comparison information is derived. Must be in the format of "$GRP1_NAME vs $GRP1_NAME [$ANNOTATION_COLUMN]".
#'   }
#'   \item{\bold{\code{differentialAdjColumns}}}{
#'       Sample annotation columns to be adjusted for in differential testing
#'   }
#'   \item{\bold{\code{lolaDbPaths}}}{
#'       Precomputed LOLA databases to be used for enrichment analysis. If \code{NULL} (default), ChrAccR will download an apropriate core database.
#'   }
#'   \item{\bold{\code{scIterativeLsiRegType}}}{
#'       For single-cell analysis only: region type to be used for clustering and dimension reduction using iterative LSI. By default (\code{NULL}),
#'       ChrAccR will look for a region type named \code{"tiling"}.
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
	cat(jsonlite::toJSON(as.list(.config), pretty=TRUE), file=dest)
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
	cfgList <- jsonlite::fromJSON(cfgFile)
	for (nn in names(cfgList)){
		if (is.element(nn,ls(.config))){
			.config[[nn]] <- cfgList[[nn]]
		} else {
			logger.warning(c("Ignoring unknown item '",nn,"' in when loading configuration from file",cfgFile))
		}
	}
}

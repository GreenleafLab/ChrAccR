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
	".default" = c("#01665E", "#35978F", "#80CDC1", "#C7EAE5", "#F5F5F5", "#F6E8C3", "#DFC27D", "#BF812D", "#8C510A")
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
.config$annotationColumns <- NULL

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
#' }
#' @author Fabian Mueller
#' @export
setConfigElement <- function(name, value){
	if (!exists(name, .config)){
		logger.error(c("No such configuration element:", name))
	}
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
	cat(toJSON(as.list(.config), pretty=TRUE), file=dest)
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
	cfgList <- fromJSON(cfgFile)
	for (nn in names(cfgList)){
		if (is.element(nn,ls(.config))){
			.config[[nn]] <- cfgList[[nn]]
		} else {
			logger.warning(c("Ignoring unknown item '",nn,"' in when loading configuration from file",cfgFile))
		}
	}
}

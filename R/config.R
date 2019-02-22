################################################################################
# Configuration and global package options
################################################################################
.config <- new.env()
.config$tmpDir <- tempdir()
.config$cleanMem <- TRUE
.config$geneModelVersions <- c(
	"hg38"="gencode.v27",
	"hg19"="gencode.v19",
	"mm10"="gencode.vM16",
	"mm9" ="gencode.vM1"
)

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
#'        \code{Directory for temporary files. Must be existing.}
#'   }
#'   \item{\bold{\code{cleanMem}}\code{ = TRUE}}{
#'        \code{During runtime, regularly clean-out the memory in order to reduce memory overuse}
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

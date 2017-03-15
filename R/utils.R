#' cleanMem
#' 
#' clean the system mory by invoking the garbage collector
#' @param iter.gc number of times to invoke the garbage collector
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @export
cleanMem <- function(iter.gc=1L){
	if (getConfigElement("cleanMem")){
		for (i in 1:iter.gc){gc()}
	}
	invisible(NULL)
}

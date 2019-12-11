################################################################################
# Various annotation functionality
################################################################################

#' getChrAccRAnnotationPackage
#' 
#' retrieve the corresponding ChrAccRAnnotation package for a given genome
#' @param genome character string specifying the genome
#' @return name of the annotation package, if installed. \code{NULL} and a warning if the package is not installed
#' @author Fabian Mueller
getChrAccRAnnotationPackage <- function(genome){
	# genome <- "hg38"
	pkgName <- paste0("ChrAccRAnnotation", gsub("^([[:alnum:]])","\\U\\1", genome, perl=TRUE)) # upper case first letter
	if (!require(pkgName, quietly=TRUE, character.only=TRUE)) {
		logger.warning(c("Could not find annotation package:", pkgName))
		pkgName <- NULL
	}
	return(pkgName)
	# # to get a function specifically from a package:
	# get("getGenomeParams", asNamespace(pkgName))
}

################################################################################
# TF annotation
################################################################################

#' getTfAnnot
#' 
#' retrieve TF annotation data
#' @param type annotation type. Currently only \code{"humantfs"} (pulls info from humantfs.ccbr.utoronto.ca) is supported
#' @return a data frame of TF annotation
#' @author Fabian Mueller
#' @export
getTfAnnot <- function(type="humantfs"){
	tfAnnot <- NULL
	if (!is.element(type, c("humantfs"))) logger.error(c("Unknown TF annotation type:", type))
	fn <- system.file(file.path("extdata", paste0("tfAnnot_", type, ".rds")), package="ChrAccR")
	if (file.exists(fn)){
		tfAnnot <- readRDS(fn)
	} else {
		if (type=="humantfs"){
			fUrl <- "http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt"
			tt <- readTab(fUrl)
			# add JASPAR motif names
			mml <- prepareMotifmatchr("hg38", "jaspar")$motifs
			tfNames <- TFBSTools::name(mml)
			tfSymbols <- getJasparSymbols(tfNames)
			tfMotifIds <- TFBSTools::ID(mml)
			motifIdL <- rep(list(NULL), nrow(tt))
			names(motifIdL) <- tt[, "HGNC.symbol"]
			for (nn in names(tfSymbols)){
				idx <- which(toupper(names(motifIdL)) %in% toupper(tfSymbols[[nn]])) # case insensitive matching
				for (i in idx){
					motifIdL[[i]] <- c(motifIdL[[i]], nn)
				}
			}
			tt[,"motifIds_jaspar"] <- sapply(motifIdL, FUN=function(x){paste(x, collapse=";")})
			tt[tt[,"motifIds_jaspar"]=="","motifIds_jaspar"] <- NA
			tfAnnot <- tt
		}
	}
	return(tfAnnot)
}

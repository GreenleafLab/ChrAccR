################################################################################
# getting insertions from alignments
################################################################################
#' getATACinsertion
#' 
#' Given a \code{GAlignmentPairs} or \code{GAlignments} object, return a \code{GRanges} object containing the insertion
#' @param ga           \code{GAlignmentPairs} (or \code{GAlignments} for single-end sequencing) object
#' @param offsetTn     apply offsets for Tn5 dimer cut site (+4bp on genomic + strand; -5bp on genomic - strand)
#' @return \code{GRanges} object containing derived insertions. For paired-end data (recommended), the width of the resulting ranges corresponds to the insert size
#'         for single-end data, the width is set to 1bp
#' @author Fabian Mueller
#' @export
getATACinsertion <- function(ga, offsetTn=TRUE){
	isPaired <- FALSE
	if (is.element("GAlignments", class(ga))){
		logger.info("detected single-end data")
	} else if (is.element("GAlignmentPairs", class(ga))){
		logger.info("detected paired-end data")
		isPaired <- TRUE
	} else {
		logger.error(c("Invalid input. Expected GAlignments object"))
	}

	res <- NULL
	if (isPaired){
		r1 <- GenomicAlignments::first(ga)
		r2 <- GenomicAlignments::second(ga)
		coordMat <- cbind(start(r1), end(r1), start(r2), end(r2))
		res <- granges(r1)
		start(res) <- rowMins(coordMat)
		end(res) <- rowMaxs(coordMat)
		if (offsetTn){
			# shift inserts inward due to the Tn5 dimer offset:
			# --> +4bp
			# -------------------------
			# -------------------------
			#                  -5bp <--
			res <- GenomicRanges::resize(res, width(res)-4, fix="end", ignore.strand=TRUE)
			res <- GenomicRanges::resize(res, width(res)-5, fix="start", ignore.strand=TRUE)
		}
	} else {
		#single end insertions just have width 1
		res <- granges(ga)
		res <- promoters(res, upstream=0, downstream=1)
		if (offsetTn){
			# shift due to Tn5 dimer offset
			res <- GenomicRanges::shift(res, 4)
			isNeg <- strand(res)=="-"
			res[isNeg] <- GenomicRanges::shift(res[isNeg], -9)
		}
	}

	return(res)
}

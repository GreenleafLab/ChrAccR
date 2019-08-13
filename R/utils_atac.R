################################################################################
# getting fragments from alignments
################################################################################
#' getATACfragments
#' 
#' Given a \code{GAlignmentPairs} or \code{GAlignments} object, return a \code{GRanges} object containing the fragment (or insertion site for single-end data)
#' @param ga           \code{GAlignmentPairs} (or \code{GAlignments} for single-end sequencing) object
#' @param offsetTn     apply offsets for Tn5 dimer cut site (+4 bp on genomic + strand; -4 bp on genomic - strand)
#' @return \code{GRanges} object containing derived insertions. For paired-end data (recommended), the width of the resulting ranges corresponds to the insert size
#'         for single-end data, the width is set to 1bp
#' @author Fabian Mueller
#' @export
getATACfragments <- function(ga, offsetTn=TRUE){
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
		res <- granges(r1)
		# counts on properly paired reads in immune atlas data
		# --> we have stranded data
		# sum(start(r1) < start(r2) & end(r2) > end(r1) & strand(res)=="+")/length(r1) # 30%
		# sum(start(r1) <= start(r2) & end(r2) >= end(r1) & strand(res)=="+")/length(r1) # 50%
		# sum(start(r2) < start(r1) & end(r2) < end(r1) & strand(res)=="-")/length(r1) # 30%
		# sum(start(r1) < start(r2) & end(r2) > end(r1) & strand(res)=="-")/length(r1) # 0%
		# sum(start(r2) < start(r1) & end(r2) < end(r1) & strand(res)=="+")/length(r1) # 0%
		# sum(start(r1) == start(r2) & end(r1) == end(r2))/length(r1) # 38%
		# sum(start(r1) < end(r2) & strand(res)=="+")/length(r1) # 50%
		# sum(start(r1) < end(r2) & strand(res)=="-")/length(r1) # 28%
		# sum(start(r2) < end(r1) & strand(res)=="-")/length(r1) # 50%
		# sum(start(r2) < end(r1) & strand(res)=="+")/length(r1) # 28%
		start(res) <- ifelse(strand(res)=="-", start(r2), start(r1))
		end(res)   <- ifelse(strand(res)=="-", end(r1), end(r2))
		# coordMat <- cbind(start(r1), end(r1), start(r2), end(r2))
		# start(res) <- rowMins(coordMat)
		# end(res) <- rowMaxs(coordMat)
		if (offsetTn){
			# shift inserts inward due to the Tn5 dimer offset:
			# --> +4 bp
			# -------------------------
			# -------------------------
			#                 -4 bp <--
			res <- GenomicRanges::resize(res, width(res)-4, fix="end", ignore.strand=TRUE)
			# res <- GenomicRanges::resize(res, width(res)-5, fix="start", ignore.strand=TRUE)
			res <- GenomicRanges::resize(res, width(res)-4, fix="start", ignore.strand=TRUE)
		}
	} else {
		#single end: get insertions instead of fragments (just have width 1)
		res <- granges(ga)
		res <- promoters(res, upstream=0, downstream=1)
		if (offsetTn){
			# shift due to Tn5 dimer offset
			res <- GenomicRanges::shift(res, 4)
			isNeg <- strand(res)=="-"
			# res[isNeg] <- GenomicRanges::shift(res[isNeg], -9)
			res[isNeg] <- GenomicRanges::shift(res[isNeg], -8)
		}
	}

	return(res)
}

#' getInsertionSitesFromFragmentGr
#' 
#' Given a \code{GAlignmentPairs} or \code{GAlignments} object, return a \code{GRanges} object containing the fragment (or insertion site for single-end data)
#' @param fragGr       GRanges object containing fragments
#' @return \code{GRanges} object containing derived insertions.
#' @author Fabian Mueller
#' @noRd
getInsertionSitesFromFragmentGr <- function(fragGr){
	isW <- width(fragGr)>1 # the insertion site is already width=1 --> single end. For paired end-data all of these should be TRUE
	grins <- GRanges()
	if (any(!isW)){
		# width==1 --> single-end data
		grins <- fragGr[!isW]
	}
	if (any(isW)){
		peStarts <- GRanges()
		peEnds   <- GRanges()
		if (all(isW)){
			# paired-end data - default case
			peStarts <- resize(fragGr, width=1, fix="start")
			peEnds   <- resize(fragGr, width=1, fix="end")
		} else {
			# mixed paired-end and single-end data
			logger.warning(c("mixed paired-end and single-end data detected for sample", sid))
			peStarts <- resize(fragGr[isW], width=1, fix="start")
			peEnds   <- resize(fragGr[isW], width=1, fix="end")
		}

		# # THIS IS NOT VALID: The Tn5 dimer is NOT always loaded with one read1 and one read2 adapter
		# # avoid double counting on neighboring fragments:
		# # only count those insertion sites once that originate from the fragment on the left and another time from the fragment on the right
		# # These incidences should only be taken into account if the fragments have different orientation (+/- strand) since the Tn5 is loaded
		# # with both read1 and read2 adapters (!WRONG ASSUMPTION!)
		# peEnds.inv <- peEnds
		# strand(peEnds.inv) <- ifelse(strand(peEnds)=="+", "-", ifelse(strand(peEnds)=="-", "+", "*"))
		# peEnds <- peEnds[!overlapsAny(peEnds.inv, peStarts, ignore.strand=FALSE)] # remove insertion sites from fragment end points that can be found as start points of fragments on the opposite strand
		# strand(peStarts)[overlapsAny(peStarts, peEnds.inv, ignore.strand=FALSE)] <- "*" #set the strand to both if it is supported by a forward and a reverse fragment

		grins <- c(
			grins,
			peStarts,
			peEnds
		)
	}
	grins <- grins[order(as.integer(seqnames(grins)), start(grins), end(grins), as.integer(strand(grins)))]
	return(grins)
}

#' readMACS2peakFile
#' 
#' Reads the MACS2 ouput as GRanges
#' @param fn	Filename for MACS2 narrow peak file
#' @return \code{GRanges} object containing peak information
#' @author Fabian Mueller
#' @export
readMACS2peakFile <- function(fn){
	if (!is.character(fn) || length(fn)!=1 || !file.exists(fn)){
		logger.error(c(fn, "is not a valid peak file"))
	}
	if (!grepl("_peaks\\.narrowPeak$", fn)){
		logger.warning(c("MACS2 peak file with unrecognized suffix/file extension:", fn))
	}
	df <- readTab(fn, header=FALSE)
	colnames(df) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "foldChange", "negLog10pval", "negLog10qval", "summitOffset")
	df[df[,"strand"]==".","strand"] <- "*"

	gr <- muRtools::df2granges(df, ids=df[,"name"], chrom.col=1L, start.col=2L, end.col=3L, strand.col=6L, coord.format="B0RE", assembly=NULL, doSort=TRUE, adjNumChromNames=TRUE)
	elementMetadata(gr)[,"summit"] <- start(gr)+elementMetadata(gr)[,"summitOffset"]

	return(gr)
}

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

#-------------------------------------------------------------------------------
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

################################################################################
# Peak helper functions
################################################################################
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

#-------------------------------------------------------------------------------
#' getNonOverlappingByScore
#' 
#' Retrieve the set of non-verlapping regions by iteratively picking the region with maximum score for
#' each set of consecutively overlapping regions
#' @param gr           \code{GRanges} object
#' @param scoreCol     name of the column to be used as score in the \code{elementMetadata} of the \code{gr} object
#' @return \code{GRanges} object containing non-overlapping regions
#' @author Fabian Mueller
#' @export
getNonOverlappingByScore <- function(gr, scoreCol="score"){
	gr.rem <- gr

	res <- GRanges()
	seqlevels(res) <- seqlevels(gr)
	seqlengths(res) <- seqlengths(gr)
	genome(res) <- genome(gr)
	i <- 0
	while (length(gr.rem) > 0){
		i <- i + 1
		# logger.status(c("iteration", i)) #DEBUG
		scs <- elementMetadata(gr.rem)[,scoreCol]
		gr.merged <- reduce(gr.rem, min.gapwidth=0L, with.revmap=TRUE, ignore.strand=TRUE)
		# maxScoreIdx <- sapply(gr.merged, FUN=function(x){
		# 	idx <- elementMetadata(x)[,"revmap"][[1]]
		# 	return(idx[which.max(scs[idx])])
		# }) #too slow
		maxScoreIdx <- sapply(elementMetadata(gr.merged)[,"revmap"], FUN=function(idx){
			return(idx[which.max(scs[idx])])
		})
		bait <- gr.rem[maxScoreIdx]
		res <- c(res, bait)
		gr.rem <- gr.rem[!overlapsAny(gr.rem, bait, ignore.strand=TRUE)]
		# logger.info(c(length(gr.rem), "regions left")) #DEBUG
	}
	return(res)
}

#-------------------------------------------------------------------------------
#' getConsensusPeakSet
#' 
#' Retrieve a consensus peak set from a set of peak lists
#' @param grl		   list or \code{GRangesList} object containing the peak sets for each sample
#' @param mode         consensus mode. Currently only "no_by_score" (non-overlapping; i.e. select the peak with the highest score from each set of
#'                     overlapping peaks) is supported.
#' @param grouping     vector of group memberships (numeric, character or factor). must be of the same length as \code{grl}
#' @param groupAgreePerc percentile of members in a group required to contain a peak in order to keep it.
#'                     E.g. a value of 1 (default) means that all replicates in a group are required to contain that peak in order
#'                     to keep it.
#' @param groupConsSelect if set, the peak set will also be checked for consistency, i.e. in order to retain a peak
#'                     it has to be consistently be present or absent in each group (as specified in \code{groupAgreePerc} percent of samples)
#' @param scoreCol     name of the column to be used as score in the \code{elementMetadata} of the peak sets. This will determine which peak is selected
#'                     if multiple peaks overlap
#' @param keepOvInfo   keep annotation columns in the elementMetadata of the results specifying whether a consensus peak overlaps with a
#'                     peak in each sample
#' @return \code{GRanges} object the containing consensus peak set
#' @author Fabian Mueller
#' @export
getConsensusPeakSet <- function(grl, mode="no_by_score", grouping=NULL, groupAgreePerc=1.0, groupConsSelect=FALSE, scoreCol="score", keepOvInfo=FALSE){
	supportedModes <- c("no_by_score")
	if (!is.element(mode, supportedModes)) logger.error(c("unsupported mode:", mode))
	if (!is.list(grl) && !is.element(class(grl), c("GRangesList", "CompressedGRangesList"))){
		logger.error("Invalid grl argument. expected list or GRangesList")
	}
	nSamples <- length(grl)
	if (length(grouping) > 0 && length(grouping)!=nSamples) logger.error("Grouping info must match the number of samples")
	if (is.null(names(grl))){
		logger.warning("No sample names specified for the peak list. Assigning arbitrary names")
		names(grl) <- paste0("sample", 1:nSamples)
	}
	if (groupAgreePerc > 1 || groupAgreePerc < 0){
		logger.error(c("Invalid value for groupAgreePerc. Must be in [0,1]"))
	}

	sampleIds <- names(grl)

	res <- NULL
	if (mode=="no_by_score"){
		# Testing ALTERNATIVES A and B: investigate whether merging all peaks and getting non-overlapping set (instead of iterating over samples)
		# is i) faster and ii) yields more consistent results
		# ==>
		# i) B is much faster
		# ii) results should be very similar
		# 
		# # ALTERNATIVE (A) iterate over samples
		# i <- 0
		# for (sid in sampleIds){
		# 	i <- i + 1
		# 	logger.status(c("Processing peaks for sample:", sid, paste0("(",i, " of ", length(sampleIds), ")")))
		# 	peakSet.cur <- grl[[sid]]
		# 	if (!is.element(class(peakSet.cur), c("GRanges"))) logger.error("Not a GRanges object")
			
		# 	#add coverage info for all samples
		# 	elementMetadata(peakSet.cur)[,paste0(".ov.", sampleIds)] <- FALSE # as.logical(NA)

		# 	if (is.null(res)){
		# 		#remove overlapping peaks in initial sample based on normalized scores
		# 		peakSet.cur <- getNonOverlappingByScore(peakSet.cur, scoreCol=scoreCol)
		# 		#initialize peak set with all peaks from the first sample
		# 		elementMetadata(peakSet.cur)[,paste0(".ov.", sid)] <- TRUE
		# 		res <- peakSet.cur
		# 	} else {
		# 		# add new peaks and remove the overlapping ones by taking the peaks with the best score
		# 		res <- getNonOverlappingByScore(c(res, peakSet.cur), scoreCol=scoreCol)
		# 		elementMetadata(res)[,paste0(".ov.", sid)] <- overlapsAny(res, peakSet.cur, ignore.strand=TRUE)
		# 	}
		# }
		# ALTERNATIVE (B) joined set
		res <- getNonOverlappingByScore(unlist(grl), scoreCol=scoreCol)
		for (sid in sampleIds){
			elementMetadata(res)[,paste0(".ov.", sid)] <- overlapsAny(res, grl[[sid]], ignore.strand=TRUE)
		}
	}
	nTotal <- length(res)
	ovMat <- as.matrix(elementMetadata(res)[,paste0(".ov.", sampleIds)])
	colnames(ovMat) <- sampleIds

	if (!keepOvInfo){
		#remove the helper columns
		for (sid in sampleIds){
			elementMetadata(res)[,paste0(".ov.", sid)] <- NULL
		}
	}

	if (!is.null(grouping) && groupAgreePerc > 0){
		logger.start("Accounting for peak reproducibility across replicates")
			groupF <- factor(grouping)
			# matrix indicating whether a peak is represented by sufficiently many samples in a group
			gRepMat <- matrix(as.logical(NA), nrow=nTotal, ncol=nlevels(groupF))
			colnames(gRepMat) <- levels(groupF)
			# matrix indicating whether a peak is absent in sufficiently many samples in a group
			gAbsentMat <- matrix(as.logical(NA), nrow=nTotal, ncol=nlevels(groupF))
			colnames(gAbsentMat) <- levels(groupF)
			for (gg in levels(groupF)){
				sidsRepl <- sampleIds[groupF==gg]
				ovMat_cur <- ovMat[, sidsRepl, drop=FALSE]
				nReq <- as.integer(ceiling(groupAgreePerc * length(sidsRepl)))
				gRepMat[,gg] <- rowSums(ovMat_cur) >= nReq
				gAbsentMat[,gg] <- rowSums(!ovMat_cur) >= nReq
			}
			keep <- matrixStats::rowAnys(gRepMat)

			nRem <- nTotal - sum(keep)
			logger.info(c("Removed", nRem, "of", nTotal, "peaks", paste0("(",round(nRem/nTotal*100, 2),"%)"), "because they were not reproduced by more than", round(groupAgreePerc*100, 2), "% of samples in any group"))

			if (groupConsSelect){
				keep2 <- matrixStats::rowAlls(gRepMat | gAbsentMat)

				nRem <- nTotal - sum(keep2)
				logger.info(c("Removed", nRem, "of", nTotal, "peaks", paste0("(",round(nRem/nTotal*100, 2),"%)"), "because they were not consistent for all groups (", round(groupAgreePerc*100, 2), "% of samples in a group consistently contained or did not contain a peak)"))
				keep <- keep & keep2
			}
			res <- res[keep]
		logger.completed()
	}

	#sort
	res <- sortSeqlevels(res)
	res <- sort(res)
	return(res)
}

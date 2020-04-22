#' prepAnaDir
#' 
#' Prepare analysis directory for ChrAccR workflow
#' @param anaDir	analysis directory
#' @param initConfig save config JSON to config directory
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @noRd
prepAnaDir <- function(anaDir, initConfig=TRUE){
	res <- dir.exists(anaDir)
	if (res){
		logger.info(c("Analysis directory already exists"))
	} else {
		dir.create(anaDir)
	}
	dds <- c("config", "log", "data", "reports")
	for (dd in dds){
		if (!dir.exists(file.path(anaDir, dd))) dir.create(file.path(anaDir, dd))
	}
	if (initConfig) saveConfig(file.path(anaDir, "config", "config.json"))
	invisible(res)
}


#' getWfState
#' 
#' get the current state of the ChrAccR workflow from the analysis directory
#' @param anaDir	analysis directory
#' @return An S3 object describing the current state of the workflow
#' @author Fabian Mueller
#' @noRd
getWfState <- function(anaDir){
	res <- list(
		anaDir = anaDir,
		configPath = as.character(NA),
		logDir = as.character(NA),
		dataDir = as.character(NA),
		dsAtacPaths = c(
			raw = as.character(NA),
			filtered = as.character(NA),
			processed = as.character(NA)
		),
		reportDir = as.character(NA),
		existingReports = c(
			summary = FALSE,
			filtering = FALSE,
			normalization = FALSE,
			exploratory = FALSE,
			differential = FALSE,
			index = FALSE
		)
	)
	if (!dir.exists(anaDir)) logger.error(c("Analysis directory does not exist:", anaDir))

	dd <- "log"
	if (dir.exists(file.path(anaDir, dd))) res[["logDir"]] <- dd
	dd <- "data"
	if (dir.exists(file.path(anaDir, dd))) res[["dataDir"]] <- dd
	dd <- "reports"
	if (dir.exists(file.path(anaDir, dd))) res[["reportDir"]] <- dd
	fn <- file.path("config", "config.json")
	if (file.exists(file.path(anaDir, fn))) res[["configPath"]] <- fn

	if (!is.na(res[["reportDir"]])){
		if (file.exists(file.path(anaDir, res[["reportDir"]], "index.html"))) res[["existingReports"]]["index"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "summary.html"))) res[["existingReports"]]["summary"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "filtering.html"))) res[["existingReports"]]["filtering"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "normalization.html"))) res[["existingReports"]]["normalization"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "exploratory.html"))) res[["existingReports"]]["exploratory"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "differential.html"))) res[["existingReports"]]["differential"] <- TRUE
	}

	if (!is.na(res[["dataDir"]])){
		stepNames <- c("raw", "filtered", "processed")
		for (sn in stepNames){
			fp <- file.path(res[["dataDir"]], paste0("dsATAC_", sn))
			if (dir.exists(file.path(anaDir, fp))) res[["dsAtacPaths"]][sn] <- fp
		}
	}
	
	class(res) <- "workflowState"
	return(res)
}

#' resetWfToStage
#' 
#' reset the current stage of the ChrAccR workflow for the analysis directory. !!! will delete content !!!
#' @param anaDir	analysis directory
#' @return Invisible: An S3 object describing the current state of the workflow
#' @author Fabian Mueller
#' @noRd
resetWfToStage <- function(anaDir, resetTo){
	validStages <- c("raw", "filtered", "processed")
	if (!is.element(resetTo, validStages)) logger.error(c("Invalid stage to reset to:", resetTo))
	wfState <- getWfState(anaDir)

	delPaths <- character(0)
	if (is.element(resetTo, c("processed", "filtered", "raw"))){
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "exploratory*"))
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "differential*"))
	}
	if (is.element(resetTo, c("filtered", "raw"))){
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "normalization*"))
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$dsAtacPaths["processed"]))
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$dataDir, "iterativeLSI*"))
	}
	if (is.element(resetTo, c("raw"))){
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$dsAtacPaths["filtered"]))
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "filtering*"))
		delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "summary*"))
	}

	delPaths <- c(delPaths, file.path(wfState$anaDir, wfState$reportDir, "index*"))
	
	unlink(delPaths, recursive=TRUE)

	wfState <- getWfState(anaDir)
	invisible(wfState)
}

#' makeReportIndex
#' 
#' Make a report index based on a report directory
#' @param reportDir	report directory
#' @param anaName	Title of the Index report
#' @param reportIds Identifiers of reports to look for
#' @return (invisible) \code{muReportR::Report} object containing the report
#' @author Fabian Mueller
#' @noRd
makeReportIndex <- function(reportDir, anaName="ChrAccR analysis", reportIds=character(0)){
	reportDesc <- data.frame(
		id = c("summary", "filtering", "normalization", "exploratory", "differential"),
		htmlName = c("summary", "filtering", "normalization", "exploratory", "differential"),
		title = c("Summary", "Filtering", "Normalization", "Exploratory analysis", "Differential analysis"),
		desc = c(
			"A summary and quality control (QC) of the loaded dataset",
			"A summary of the filtering steps applied to the dataset",
			"A summary of the normalization procedure applied to count data",
			"Exploratory analyses including unsupervised methods for clustering and dimension reduction and transcription factor motif activities",
			"Analysis of differential accessibility between groups of samples and characterization of identified differences"
		),
		stringsAsFactors = FALSE
	)
	rownames(reportDesc) <- reportDesc[,"id"]

	if (length(reportIds) < 1 || !all(reportIds %in% rownames(reportDesc))){
		logger.warning("Could not create report summary: invalid report ids")
		invisible(NULL)
	}

	htmlFns <- file.path(reportDir, paste0(reportDesc[reportIds,"htmlName"], ".html"))
	if (!all(file.exists(htmlFns))){
		logger.warning("Could not create report summary: Not all requested reports exist")
		invisible(NULL)
	}

	if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
	initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
	rr <- muReportR::createReport(file.path(reportDir, paste0("index", ".html")), anaName, page.title=anaName, init.configuration=initConfigDir, theme="stanford")

	imgDir <- muReportR::getReportDir(rr, dir="images", absolute=TRUE)
	file.copy(system.file("extdata/dna.png", package="ChrAccR", mustWork=TRUE), imgDir)

	txt <- c(
		"Here are the resultes of your ChrAccR ATAC-seq analysis. The following reports have been generated:"
	)
	rr <- muReportR::addReportSection(rr, "Table of contents", txt, collapsed="never")

	## Build table of contents
	stext <- c("<table class=\"pipeline\" style=\"background-image:url(index_images/dna.png);\">", "<colgroup>",
		"\t<col width=\"650\" />", "</colgroup>")
	for (rid in reportIds) {
		rfile <- file.path(paste0(reportDesc[rid,"htmlName"], ".html"))
		mod.open <- paste("<a href=", rfile, "><div class=", "completed", ">", sep = "\"")
		mod.close <- "</div></a>"
		stext <- c(stext, "<tr>")
		stext <- c(stext, paste0("\t<td>", mod.open))
		stext <- c(stext, paste0("\t\t<h3>", reportDesc[rid, "title"], "</h3>"))
		stext <- c(stext, paste0("\t\t<p>", reportDesc[rid, "desc"], "</p>"))
		stext <- c(stext, paste0("\t", mod.close, "</td>"))
		stext <- c(stext, "</tr>")
	}
	stext <- c(stext, "</table>")

	cat(paste(stext, collapse = "\n"), "\n", file=rr@fname, sep="", append=TRUE)

	muReportR::off(rr)
	invisible(rr)
}

#-------------------------------------------------------------------------------
#' run_atac_qc
#' 
#' Run the summary QC analysis for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing QC statistics and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_qc <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	doReport <- !wfState$existingReports["summary"]

	report <- NULL
	qcStats <- NULL

	if (doReport){
		logger.start("Creating summary report")
			report <- createReport_summary(dsa, file.path(wfState$anaDir, wfState$reportDir))
		logger.completed()
	}
	res <- list(
		qcStats=qcStats,
		report=report
	)
	class(res) <- "ChrAccR_runRes_qc"
	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_peakcalling
#' 
#' Run peak calling for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing the annotated \code{\linkS4class{DsATAC}} object, per-sample peak calls, a consensus peak set and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_peakcalling <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	report <- NULL
	peakGr <- NULL
	peakGrl <- NULL

	isSingleCell <- class(dsa)=="DsATACsc"
	dsn <- dsa

	doPeakCall <- !isSingleCell
	if (doPeakCall){
		logger.start("Peak calling")
			
			logger.start("Per-sample peak sets")
				peakGrl <- callPeaks(dsn)
			logger.completed()

			logger.start("Consensus peak set")
				ggs <- NULL
				cn <- getConfigElement("annotationPeakGroupColumn")
				if (!is.null(cn)){
					aa <- getSampleAnnot(dsn)
					if (is.element(cn, colnames(aa))){
						logger.info(c("Using annotation column '", cn, "' to check for peak consistency in groups"))
						ggs <- aa[,cn]
					} else {
						logger.warning(c("Could not find peak group annotation column:", cn, ". Continuing without group consensus peaks"))
					}
				}
				peakGr <- getConsensusPeakSet(peakGrl, mode="no_by_score", grouping=ggs, groupAgreePerc=getConfigElement("annotationPeakGroupAgreePerc"), groupConsSelect=FALSE, scoreCol="score", keepOvInfo=TRUE)
				logger.info(c("Identified", length(peakGr), "consensus peaks"))
				logger.start("Adding peak set to DsATAC")
					dsn <- regionAggregation(dsn, peakGr, ".peaks.cons", signal="insertions", dropEmpty=FALSE)
				logger.completed()
			logger.completed()
		logger.completed()
	} else {
		logger.warning("Peaks will not be called")
	}

	res <- list(
		ds_anno=dsn,
		peakGr=peakGr,
		peakGrl=peakGrl,
		report=report
	)
	class(res) <- "ChrAccR_runRes_callPeaks"

	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_filtering
#' 
#' Run the filtering for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing the filtered \code{\linkS4class{DsATAC}} object, filtering statistics and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_filtering <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	report <- NULL
	filterStats <- list(
		regionStats=NULL,
		params=list(),
		steps=character(0)
	)

	regTypes <- getRegionTypes(dsa)
	filterStats$regionStats <- data.frame(
		before=sapply(regTypes, FUN=function(rt){getNRegions(dsa, rt)})
	)
	rownames(filterStats$regionStats) <- regTypes

	isSingleCell <- class(dsa)=="DsATACsc"
	dsf <- dsa

	# low coverage filtering
	reqSamples <- getConfigElement("filteringCovgReqSamples")
	thresh <- getConfigElement("filteringCovgCount")
	if (reqSamples > 0 && thresh > 0 && !isSingleCell){
		filterStats$params[["covgReqSamples"]] <- reqSamples
		filterStats$params[["covgCount"]] <- thresh
		filterStats$steps <- c(filterStats$steps, "region_covg")
		logger.start("Coverage filtering")
			logger.info(c("filteringCovgReqSamples:", reqSamples))
			logger.info(c("filteringCovgCount:", thresh))
			dsf <- filterLowCovg(dsa, thresh=thresh, reqSamples=reqSamples)
			filterStats$regionStats[,"after_region_covg"] <- sapply(regTypes, FUN=function(rt){getNRegions(dsf, rt)})
		logger.completed()
	}

	# low coverage filtering
	exclChroms <- c("chrM")
	if (getConfigElement("filteringSexChroms")){
		exclChroms <- union(exclChroms, c("chrX", "chrY"))
	}
	if (length(exclChroms) > 0){
		filterStats$steps <- c(filterStats$steps, "chromosomes")
		filterStats$params[["exclChroms"]] <- exclChroms
		logger.start("Chromosome filtering")
			logger.info(c("Excluding chromosomes:", paste(exclChroms, collapse=",")))
			dsf <- filterChroms(dsf, exclChrom=exclChroms)
			filterStats$regionStats[,"after_chromosomes"] <- sapply(regTypes, FUN=function(rt){getNRegions(dsf, rt)})
		logger.completed()
	}
	
	if (isSingleCell){
		sampleIds <- unique(getSampleAnnot(dsa)[,".sampleId"])
		getNcellsPerSample <- function(dsa){
			sids <- getSampleAnnot(dsa)[,".sampleId"]
			res <- table(sids)[sampleIds]
			return(res)
		}

		filterStats[["cellStats"]] <- data.frame(
			before = getNcellsPerSample(dsa)
		)
		rownames(filterStats[["cellStats"]]) <- sampleIds

		cellQcTab <- getScQcStatsTab(dsa)
		
		keepCell <- rep(TRUE, length(getSamples(dsa)))
		fragTmin <- getConfigElement("filteringScMinFragmentsPerCell")
		if (fragTmin > 0){
			if (any(keepCell & cellQcTab[,"nPass"] < fragTmin)){
				filterStats$steps <- c(filterStats$steps, "sc_minFrags")
				filterStats$params[["scMinFragmentsPerCell"]] <- fragTmin
			}
			keepCell <- keepCell & cellQcTab[,"nPass"] >= fragTmin
		}
		fragTmax <- getConfigElement("filteringScMaxFragmentsPerCell")
		if (fragTmax < Inf){
			if (any(keepCell & cellQcTab[,"nPass"] > fragTmax)){
				filterStats$steps <- c(filterStats$steps, "sc_maxFrags")
				filterStats$params[["scMaxFragmentsPerCell"]] <- fragTmax
			}
			keepCell <- keepCell & cellQcTab[,"nPass"] <= fragTmax
		}
		
		tssT <- getConfigElement("filteringScMinTssEnrichment")
		if (is.element("tssEnrichment", colnames(cellQcTab)) && tssT > 0){
			if (any(keepCell & cellQcTab[,"tssEnrichment"] < tssT)){
				filterStats$steps <- c(filterStats$steps, "sc_tssEnrichment")
				filterStats$params[["scMinTssEnrichment"]] <- tssT
			}
			keepCell <- keepCell & cellQcTab[,"tssEnrichment"] >= tssT
		}
		if (!all(keepCell)){
			logger.start("Single-cell filtering")
				if (fragTmin > 0) logger.info(c("filteringScMinFragmentsPerCell:", fragTmin))
				if (fragTmax < Inf) logger.info(c("filteringScMaxFragmentsPerCell:", fragTmax))
				if (tssT > 0) logger.info(c("filteringScMinTssEnrichment:", tssT))

				dsf <- dsf[keepCell]
				# workaround: currently saving single-cell DsATAC datasets does not support re-chunking of disk-dumped fragment data
				chunkedFragmentFiles <- dsf@diskDump.fragments && .hasSlot(dsf, "diskDump.fragments.nSamplesPerFile") && dsf@diskDump.fragments.nSamplesPerFile > 1
				if (chunkedFragmentFiles){
					logger.start("Undisking ...")
						dsf <- undiskFragmentData(dsf)
					logger.completed()
				}
			logger.completed()
			filterStats[["cellStats"]][,"after"] <- getNcellsPerSample(dsf)
		}
	}

	filterStats$regionStats[,"after"] <- sapply(regTypes, FUN=function(rt){getNRegions(dsf, rt)})

	doReport <- !wfState$existingReports["filtering"]
	if (doReport){
		logger.start("Creating filtering report")
			report <- createReport_filtering(dsf, file.path(wfState$anaDir, wfState$reportDir), unfilteredObj=dsa, filterStats=filterStats)
		logger.completed()
	}

	res <- list(
		ds_filtered=dsf,
		filterStats=filterStats,
		report=report
	)
	class(res) <- "ChrAccR_runRes_filtering"
	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_normalization
#' 
#' Run count normalization for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing the normalized \code{\linkS4class{DsATAC}} object and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_normalization <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	report <- NULL

	isSingleCell <- class(dsa)=="DsATACsc"
	dsn <- dsa

	validNormMethods <- c("quantile", "percentile", "rankPerc", "log2", "RPKM", "CPM", "vst", "tf-idf", "logCPM", "logRPKM", "none")
	nm <- getConfigElement("normalizationMethod")
	if (!is.element(nm, validNormMethods)){
		logger.error(c("Invalid normalization method:", nm))
	}
	doNorm <- !isSingleCell && !is.element(nm, c("none"))
	if (doNorm){
		logger.start("Normalizing count data")
			logger.info(c("normalizationMethod:", nm))
			doLog <- grepl("^log", nm)
			nm <- gsub("^log", "", nm)
			dsn <- transformCounts(dsn, method=nm)
			if (doLog){
				dsn <- transformCounts(dsn, method="log10")
			}
		logger.completed()
	} else {
		logger.warning("Count data will not be normalized")
	}

	doReport <- doNorm && !wfState$existingReports["normalization"]
	if (doReport){
		logger.start("Creating normalization report")
			report <- createReport_normalization(dsn, file.path(wfState$anaDir, wfState$reportDir), unnormObj=dsa)
		logger.completed()
	}

	res <- list(
		ds_norm=dsn,
		report=report
	)
	class(res) <- "ChrAccR_runRes_norm"

	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_sc_unsupervised
#' 
#' Run unsupervised analysis for single-cell ATAC-seq data (i.e. iterative LSI, clustering and cluster peak detection)
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing the annoted \code{\linkS4class{DsATAC}} object, the results of running iterative LSI and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_sc_unsupervised <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	report <- NULL

	isSingleCell <- class(dsa)=="DsATACsc"
	dsan <- dsa
	geneAct <- NULL
	doUnsupervised <- isSingleCell
	if (doUnsupervised){
		findItLsiRt <- function(ds){
			return(findOrderedNames(getRegionTypes(ds), c("tiling", "^t[1-9]"), exact=FALSE, ignore.case=TRUE))
		}
		itLsiRt <- getConfigElement("scIterativeLsiRegType")
		if (length(itLsiRt) < 1) {
			itLsiRt <- findItLsiRt(dsan)
			if (is.na(itLsiRt)) {
				logger.error("Could not determine region type for iterative LSI")
			}
		}
		logger.start("Iterative LSI")
			logger.info(c("Using region type:", itLsiRt))
			argL <- getConfigElement("scIterativeLsiParams")
			argL[[".object"]] <- dsan
			argL[["it0regionType"]] <- itLsiRt
			# logger.info(c("Using cluster resolution:", clustRes))
			# itLsi <- iterativeLSI(
			# 	dsan,
			# 	it0regionType=itLsiRt,
			# 	it0clusterResolution=clustRes,
			# 	it1clusterResolution=clustRes,
			# 	it2clusterResolution=clustRes,
			# 	rmDepthCor=0.5,
			# 	normPcs=FALSE,
			# 	umapParams=umapParams
			# )
			itLsi <- do.call(iterativeLSI, argL)
		logger.completed()
		logger.start("Aggregating counts across initial cluster peaks")
			dsan <- regionAggregation(dsan, itLsi$clusterPeaks_unfiltered, ".peaks.itlsi0", signal="insertions", dropEmpty=FALSE, bySample=FALSE)
		logger.completed()
		logger.start("Aggregating counts across final iterative LSI features (peaks)")
			dsan <- regionAggregation(dsan, itLsi$regionGr, ".itlsi.features", signal="insertions", dropEmpty=FALSE, bySample=FALSE)
		logger.completed()
		dsan <- addSampleAnnotCol(dsan, ".itlsi.clustering", itLsi$clustAss[getSamples(dsan)])
		logger.info("Added clustering info to DsATAC object")

		# gene activity scores
		gaMethod <- getConfigElement("scGeneActivity")
		doGeneAct <- !is.null(gaMethod) && !is.na(gaMethod)
		if (doGeneAct){
			if (is.logical(gaMethod)) {
				doGeneAct <- gaMethod
				if (doGeneAct) gaMethod <- "RBF"
			} else if (is.character(gaMethod)) {
				doGeneAct <- nchar(gaMethod) > 0
			}
		}
		if (doGeneAct){
			logger.start("Computing gene activity scores")
				if (is.element(tolower(gaMethod), c("cicero"))){
					logger.info("Using Cicero")
					geneAct <- ChrAccR::getCiceroGeneActivities(dsan, ".peaks.itlsi0", promoterGr=NULL, dimRedCoord=itLsi$pcaCoord[,itLsi$pcs])
				} else if (is.element(tolower(gaMethod), c("rbf"))){
					logger.info("Using RBF-weighted count aggregation")
					geneAct <- ChrAccR::getRBFGeneActivities(dsan, ".peaks.itlsi0", tssGr=NULL)
				} else {
					logger.error(c("Unknown method for gene activity score computation:", gaMethod))
				}
			logger.completed()
		}
	} else {
		logger.warning("Unsupervised analysis will not be performed")
	}

	res <- list(
		ds_anno=dsan,
		itLsi=itLsi,
		geneActivity=geneAct,
		report=report
	)
	class(res) <- "ChrAccR_runRes_sc_unsupervised"

	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_exploratory
#' 
#' Run exploratory analyses for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @param itLsiObj  [for single-cell only; optional] pre-computed result of a call to \code{iterativeLSI(.object, ...)}
#' @param geneActSe [for single-cell only; optional] pre-computed result of a call to \code{getCiceroGeneActivities(.object, ...)}
#' @return \code{S3} object containing exploratory metrics and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_exploratory <- function(dsa, anaDir, itLsiObj=NULL, geneActSe=NULL){
	wfState <- getWfState(anaDir)
	doReport <- !wfState$existingReports["exploratory"]

	report <- NULL

	if (doReport){
		logger.start("Creating exploratory report")
			report <- createReport_exploratory(dsa, file.path(wfState$anaDir, wfState$reportDir), itLsiObj=itLsiObj, geneActSe=geneActSe)
		logger.completed()
	}
	res <- list(
		metrics=NULL,
		report=report
	)
	class(res) <- "ChrAccR_runRes_exploratory"
	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac_differential
#' 
#' Run differential analyses for ATAC-seq data
#' @param dsa       \code{\linkS4class{DsATAC}} object
#' @param anaDir	analysis directory
#' @return \code{S3} object containing differential analysis results and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_differential <- function(dsa, anaDir){
	wfState <- getWfState(anaDir)
	doReport <- !wfState$existingReports["differential"]

	report <- NULL

	if (doReport){
		logger.start("Creating differential report")
			report <- createReport_differential(dsa, file.path(wfState$anaDir, wfState$reportDir))
		logger.completed()
	}
	res <- list(
		metrics=NULL,
		report=report
	)
	class(res) <- "ChrAccR_runRes_differential"
	return(res)
}

#-------------------------------------------------------------------------------
#' run_atac
#' 
#' Run the complete ChrAccR analysis for ATAC-seq data
#' @param anaDir	analysis directory
#' @param input		Input object. Can be either \code{NULL}, a character string, a \code{DsATAC}. Set to \code{NULL} when you want to continue a previous analysis
#' @param sampleAnnot sample annotation table (\code{data.frame}) or \code{NULL} if continuing existing analysis or input is a \code{DsATAC} object
#' @param genome    genome assembly. Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param sampleIdCol column name in the sample annotation table containing unique sample Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param regionSets a list of GRanges objects which contain region sets over which count data will be aggregated. Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param startStage stage where to start the analysis from. can be one of \code{"raw"}, \code{"filtered"}, \code{"processed"}. Only relevant if not continuing existing analysis.
#' @param resetStage flag indicating whether to reset the analysis directory (i.e. deleting previously generated reports and datasets), when continuing previous analyses (\code{input} argument is \code{NULL}).
#' @return \code{\linkS4class{DsATAC}} object (invisible)
#' @author Fabian Mueller
#' @export
run_atac <- function(anaDir, input=NULL, sampleAnnot=NULL, genome=NULL, sampleIdCol=NULL, regionSets=NULL, startStage="raw", resetStage=NULL){
	doContinue <- prepAnaDir(anaDir, initConfig=TRUE)
	wfState <- getWfState(anaDir)
	
	loadConfig(file.path(wfState$anaDir, wfState$configPath))

	logFn <- file.path(wfState$anaDir, wfState$logDir, paste0(muRtools::getHashString("chraccr"), ".log"))
	logger.start(fname=c(NA, logFn))
	logger.start("ChrAccR analysis")
	############################################################################
	# Prepare DsATAC dataset
	############################################################################
	saveDs <- TRUE
	dsa <- NULL
	dsa_unnorm <- NULL
	dsStages <- c("raw", "filtered", "processed")
	if (!is.element(startStage, dsStages)) logger.error(c("Invalid dataset start stage:", startStage))

	if (is.null(input)){
		if (doContinue) {
			logger.info("Continuing ChrAccR analysis")
			logger.start("Loading DsATAC dataset")
				dsfPath <- file.path(wfState$anaDir, wfState$dsAtacPaths["filtered"])
				if (!is.na(wfState$dsAtacPaths["processed"])){
					logger.info("Continuing from processed dataset")
					dsa <- loadDsAcc(file.path(wfState$anaDir, wfState$dsAtacPaths["processed"]))
					startStage <- "processed"
					saveDs <- FALSE
				} else if (!is.na(wfState$dsAtacPaths["filtered"])){
					logger.info("Continuing from filtered dataset")
					dsa <- dsa_unnorm <- loadDsAcc(dsfPath)
					startStage <- "filtered"
				} else if (!is.na(wfState$dsAtacPaths["raw"])){
					logger.info("Continuing from raw dataset")
					dsa <- loadDsAcc(file.path(wfState$anaDir, wfState$dsAtacPaths["raw"]))
					startStage <- "raw"
				} else {
					logger.error("Could not find existing dataset in analysis directory")
				}
				isSingleCell <- class(dsa)=="DsATACsc"

				if (is.element(startStage, c("processed")) && !is.na(wfState$dsAtacPaths["filtered"]) && !isSingleCell){
					logger.status("Loading filtered, unnormalized dataset")
					dsa_unnorm <- loadDsAcc(dsfPath)
				}
			logger.completed()
		} else {
			logger.error("Cannot continue existing analysis: analysis does not exist")
		}
	} else {
		if (doContinue) {
			logger.error("Starting new analysis although the analysis directory already exists. For continuing an existing analysis 'input' argument has to be NULL")
		}
		if (is.element(class(input), c("DsATAC", "DsATACsc"))){
			sampleAnnot <- getSampleAnnot(input)
		}
		validSampleAnnot <- is.data.frame(sampleAnnot)
		if (is.element(class(input), c("DsATAC", "DsATACsc"))){
			logger.info(c("Detected input type:", class(input)))
			logger.info(c("Start stage:", startStage))
			dsa <- input
			# if (startStage=="processed"){
			# 	saveDs <- FALSE
			# }
		} else if (validSampleAnnot){
			startStage <- "raw"
			inputType <- as.character(NA)
			inputFns <- character(0)
			if (is.character(input)){
				inputFns <- input
				if (length(input) == 1 && is.element(input, colnames(sampleAnnot))){
					inputFns <- sampleAnnot[,input]
				}

				# check if input is cellranger output
				validCellRanger <- sapply(inputFns, FUN=function(x){
					dir.exists(x) && dir.exists(file.path(x, "outs"))
				})

				if (all(validCellRanger)){
					inputType <- "sc_cellranger"
				} else if (all(grepl("\\.fragments\\.tsv(\\.gz)?$", inputFns))){
					inputType <- "sc_fragments"
				} else if (length(inputFns) > 1){
					# assuming bulk
					if (all(grepl("\\.bam$", inputFns))){
						inputType <- "bulk_bam"
					} else if (all(grepl("\\.bed$", inputFns))){
						inputType <- "bulk_fragBed"
					}
				}
			}
			if (is.na(inputType)) logger.error("Could not detect input type")
			logger.info(c("Detected input type:", inputType))
			# detect sample id column if not specified
			if (is.null(sampleIdCol)){
				sampleIdCol <- grep("sample", colnames(sampleAnnot), value=TRUE)
				if (length(sampleIdCol) < 1) {
					logger.error("Could not uniquely determine sample id column from sample annotation column names.")
				} else if (length(sampleIdCol) > 1) {
					logger.warning("Could not uniquely determine sample id column from sample annotation column names. --> Picking the first one")
					sampleIdCol <- sampleIdCol[1]
				}
				logger.info(c("Automatically determined column '", sampleIdCol, "' as sample identifier column"))
			}
			doRepMerge <- anyDuplicated(sampleAnnot[,sampleIdCol]) && grepl("^bulk", inputType)
			repGrps <- NULL
			if (doRepMerge){
				logger.info("Rows in the sample annotation table contain duplicate ids. --> treating these as replicates")
				repGrps <- sampleAnnot[,sampleIdCol]
				sampleAnnot[,".sampleId_made_unique"] <- make.unique(sampleAnnot[,sampleIdCol], sep="_")
				sampleIdCol <- ".sampleId_made_unique"
			}
			logger.start("Preparing DsATAC dataset from input files")
				if (inputType == "bulk_bam"){
					dsa <- DsATAC.bam(sampleAnnot, inputFns, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, diskDump=TRUE, keepInsertionInfo=TRUE, pairedEnd=TRUE)
				} else if (inputType == "bulk_fragBed"){
					dsa <- DsATAC.fragmentBed(sampleAnnot, inputFns, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, diskDump=TRUE, keepInsertionInfo=TRUE)
				} else if (inputType == "sc_fragments"){
					minFrags <- getConfigElement("filteringScMinFragmentsPerCell")
					maxFrags <- getConfigElement("filteringScMaxFragmentsPerCell")
					dsa <- DsATACsc.fragments(sampleAnnot, inputFns, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, minFragsPerBarcode=minFrags, maxFragsPerBarcode=maxFrags, cellAnnot=NULL, keepInsertionInfo=TRUE, diskDump.fragments=FALSE, cellQcStats=TRUE)
				} else if (inputType == "sc_cellranger"){
					dsa <- DsATAC.cellranger(sampleAnnot, input, genome, dataDir="", regionSets=regionSets, addPeakRegions=TRUE, sampleIdCol=sampleIdCol, keepInsertionInfo=TRUE, diskDump.fragments=FALSE)
				}
				if (doRepMerge){
					logger.start("Merging replicates")
						dsa <- mergeSamples(dsa, repGrps, countAggrFun="sum")
					logger.completed()
					dsa@sampleAnnot[,".sampleId_made_unique"] <- NULL # delete artificially unique sample id column
				}
			logger.completed()
		} else {
			logger.error("Invalid input parameter")
		}
	}
	isSingleCell <- class(dsa)=="DsATACsc"

	if (is.null(resetStage)){
		resetStage <- FALSE
		if (is.null(input) && doContinue){
			resetStage <- TRUE
		}
	}
	if (resetStage){
		logger.warning(c("Resetting analysis directory to", startStage, "state"))
		wfState <- resetWfToStage(anaDir, startStage)
	}

	# identifying consensus peak set
	doPeakCalling <- !isSingleCell && getConfigElement("doPeakCalling") && !is.element(".peaks.cons", getRegionTypes(dsa)) && length(dsa@fragments) > 0
	if (doPeakCalling){
		logger.start("Running peak calling analysis")
			res <- run_atac_peakcalling(dsa, anaDir)
			dsa <- res$ds_anno
		logger.completed()
		logger.start("Saving peak data")
			saveRDS(res$peakGr, file.path(wfState$anaDir, wfState$dataDir, paste0("peakGr_consensus", ".rds")))
			saveRDS(res$peakGrl, file.path(wfState$anaDir, wfState$dataDir, paste0("peakGrl_perSample", ".rds")))
		logger.completed()
	}

	saveDs_raw <- saveDs && is.na(wfState$dsAtacPaths["raw"]) && is.element(startStage, c("raw"))
	if (saveDs_raw){
		wfState$dsAtacPaths["raw"] <- file.path(wfState$dataDir, "dsATAC_raw")
		logger.start("Saving raw DsATAC dataset")
			saveDsAcc(dsa, file.path(wfState$anaDir, wfState$dsAtacPaths["raw"]))
		logger.completed()
	}
	############################################################################
	# Run analysis
	############################################################################
	doQc <- TRUE
	if (doQc){
		logger.start("Running summary QC analysis")
			res <- run_atac_qc(dsa, anaDir)
		logger.completed()
	}
	#---------------------------------------------------------------------------
	doFilter <- is.element(startStage, c("raw"))
	if (doFilter){
		logger.start("Running filtering analysis")
			res <- run_atac_filtering(dsa, anaDir)
			dsa <- res$ds_filtered
		logger.completed()
	}
	saveDs_filtered <- saveDs && (doFilter) && is.na(wfState$dsAtacPaths["filtered"])
	if (saveDs_filtered){
		wfState$dsAtacPaths["filtered"] <- file.path(wfState$dataDir, "dsATAC_filtered")
		logger.start("Saving filtered DsATAC dataset")
			saveDsAcc(dsa, file.path(wfState$anaDir, wfState$dsAtacPaths["filtered"]))
		logger.completed()
	}
	#---------------------------------------------------------------------------
	# Bulk: normalize data
	diffWarn <- FALSE
	if (is.null(dsa_unnorm)) {
		diffWarn <- TRUE
		dsa_unnorm <- dsa
	}
	doNorm <- !isSingleCell && is.element(startStage, c("raw", "filtered"))
	if (doNorm){
		logger.start("Running normalization analysis")
			res <- run_atac_normalization(dsa, anaDir)
			dsa <- res$ds_norm
		logger.completed()
	}

	# Single-cell: Perform unsupervised analysis
	itLsi <- NULL
	geneActSe <- NULL
	itLsiFp <- file.path(wfState$anaDir, wfState$dataDir, "iterativeLSI")
	itLsiFn <- paste0(itLsiFp, ".rds")
	gaFn <- file.path(wfState$anaDir, wfState$dataDir, "geneActivitySE.rds")
	doScUnsupervised <- isSingleCell && is.element(startStage, c("raw", "filtered")) && !file.exists(itLsiFn)
	if (doScUnsupervised){
		logger.start("Running unsupervised single-cell analysis")
			res <- run_atac_sc_unsupervised(dsa, anaDir)
			dsa <- res$ds_anno
			itLsi <- res$itLsi
		logger.completed()
		logger.start("Saving unsupervised results")
			saveRDS(itLsi, itLsiFn)
			uwot::save_uwot(itLsi$umapRes, paste0(itLsiFp, "_uwot"))
		logger.completed()
		if (!is.null(res$geneActivity)){
			geneActSe <- res$geneActivity
			logger.start("Saving gene activity scores")
				saveRDS(geneActSe, gaFn)
			logger.completed()
		}
	} else if (isSingleCell && file.exists(itLsiFn)){
		logger.start("Loading iterative LSI result from disk")
			itLsi <- readRDS(itLsiFn)
			itLsi$umapRes <- uwot::load_uwot(paste0(itLsiFp, "_uwot"))
		logger.completed()
		if (file.exists(gaFn)){
			logger.start("Loading gene activity result from disk")
				geneActSe <- readRDS(gaFn)
			logger.completed()
		}
	}

	if (FALSE){
		logger.start("Removing fragment data (to save space)")
			dsa <- removeFragmentData(dsa)
		logger.completed()
	}

	saveDs_processed <- saveDs && (doNorm || doScUnsupervised) && is.na(wfState$dsAtacPaths["processed"])
	if (saveDs_processed){
		wfState$dsAtacPaths["processed"] <- file.path(wfState$dataDir, "dsATAC_processed")
		logger.start("Saving processed DsATAC dataset")
			saveDsAcc(dsa, file.path(wfState$anaDir, wfState$dsAtacPaths["processed"]))
		logger.completed()
	}
	#---------------------------------------------------------------------------
	doExploratory <- TRUE
	if (doExploratory){
		logger.start("Running exploratory analysis")
			res <- run_atac_exploratory(dsa, anaDir, itLsiObj=itLsi, geneActSe=geneActSe)
		logger.completed()
	}
	#---------------------------------------------------------------------------
	doDifferential <- !isSingleCell && (!is.null(getConfigElement("differentialColumns")) || !is.null(getConfigElement("differentialCompNames")))
	if (doDifferential){
		if (startStage=="processed" && diffWarn){
			logger.warning(c("You are starting the analysis from a processed dataset. It is not recommended to compute differential accessibility from normalized counts"))
		}
		logger.start("Running differential analysis")
			res <- run_atac_differential(dsa_unnorm, anaDir)
		logger.completed()
	}

	# update workflow state
	wfState <- getWfState(anaDir)
	
	# Make report index
	# delete existing index
	reportDir <- file.path(wfState$anaDir, wfState$reportDir)
	if (wfState$existingReports[["index"]]){
		unlink(file.path(reportDir, "index*"), recursive=TRUE)
	}
	repIds <- c("summary", "filtering", "normalization", "exploratory", "differential")
	repIds <- intersect(repIds, names(wfState$existingReports)[wfState$existingReports])
	makeReportIndex(reportDir, anaName=getConfigElement("analysisName"), reportIds=repIds)
	
	logger.completed()
	invisible(dsa)
}

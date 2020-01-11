#' prepAnaDir
#' 
#' Prepare analysis directory for ChrAccR workflow
#' @param anaDir	analysis directory
#' @return nothing of particular interest
#' @author Fabian Mueller
#' @noRd
prepAnaDir <- function(anaDir){
	res <- dir.exists(anaDir)
	if (res){
		logger.info(c("Analysis directory already exists"))
	} else {
		dir.create(anaDir)
	}
	dds <- c("config", "log", "data")
	for (dd in dds){
		if (!dir.exists(file.path(anaDir, dd))) dir.create(file.path(anaDir, dd))
	}
	saveConfig(file.path(anaDir, "config", "config.json"))
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
			exploratory = FALSE,
			differential = FALSE
		)
	)
	if (!dir.exists(anaDir)) logger.error(c("Analysis directory does not exist:", anaDir))

	dd <- "log"
	if (dir.exists(file.path(anaDir, dd))) res[["logDir"]] <- dd
	dd <- "data"
	if (dir.exists(file.path(anaDir, dd))) res[["dataDir"]] <- dd
	dd <- "reports"
	if (dir.exists(file.path(anaDir, dd))) res[["reportDir"]] <- dd
	fn <- file.path(anaDir, "config", "config.json")
	if (file.exists(fn)) res[["configPath"]] <- fn

	if (!is.na(res[["reportDir"]])){
		if (file.exists(file.path(anaDir, res[["reportDir"]], "summary.html"))) res[["existingReports"]]["summary"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "exploratory.html"))) res[["existingReports"]]["exploratory"] <- TRUE
		if (file.exists(file.path(anaDir, res[["reportDir"]], "differential.html"))) res[["existingReports"]]["differential"] <- TRUE
	}

	if (!is.na(res[["dataDir"]])){
		stepNames <- c("raw", "filtered", "processed")
		for (sn in stepNames){
			fp <- paste0("dsATAC_", sn)
			if (dir.exists(file.path(anaDir, res[["dataDir"]], fp))) res[["dsAtacPaths"]][sn] <- fp
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
		delPaths <- c(delPaths, file.path(wfState$dataDir, wfState$reportDir, "exploratory*"))
		delPaths <- c(delPaths, file.path(wfState$dataDir, wfState$reportDir, "differential*"))
	}
	if (is.element(resetTo, c("filtered", "raw"))){
		delPaths <- c(delPaths, file.path(wfState$dataDir, wfState$dsAtacPaths["processed"]))
	}
	if (is.element(resetTo, c("raw"))){
		delPaths <- c(delPaths, file.path(wfState$dataDir, wfState$dsAtacPaths["filtered"]))
		delPaths <- c(delPaths, file.path(wfState$dataDir, wfState$reportDir, "summary*"))
	}
	
	unlink(delPaths, recursive=TRUE)

	wfState <- getWfState(anaDir)
	invisible(wfState)
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
		regionStats=NULL
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
		logger.start("Coverage filtering")
			logger.info(c("filteringCovgReqSamples:", reqSamples))
			logger.info(c("filteringCovgCount:", thresh))
			dsf <- filterLowCovg(dsa, thresh=thresh, reqSamples=reqSamples)
			filterStats$regionStats[,"after_covg"] <- sapply(regTypes, FUN=function(rt){getNRegions(dsf, rt)})
		logger.completed()
	}

	# low coverage filtering
	exclChroms <- c("chrM")
	if (getConfigElement("filteringSexChroms")){
		exclChroms <- union(exclChroms, c("chrX", "chrY"))
	}
	if (length(exclChroms) > 0){
		logger.start("Chromosome filtering")
			dsf <- filterChroms(dsf, exclChrom=exclChroms)
			filterStats$regionStats[,"after_chrom"] <- sapply(regTypes, FUN=function(rt){getNRegions(dsf, rt)})
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
			keepCell <- keepCell & cellQcTab[,"nPass"] >= fragTmin
		}
		fragTmax <- getConfigElement("filteringScMinFragmentsPerCell")
		if (fragTmax < Inf){
			keepCell <- keepCell & cellQcTab[,"nPass"] <= fragTmax
		}
		
		tssT <- getConfigElement("filteringScMinTssEnrichment")
		if (is.element("tssEnrichment", colnames(cellQcTab)) && tssT > 0){
			keepCell <- keepCell & cellQcTab[,"tssEnrichment"] > tssT
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
			clustRes <- getConfigElement("scIterativeLsiClusterResolution")
			logger.info(c("Using cluster resolution:", clustRes))
			itLsi <- iterativeLSI(dsan, it0regionType=itLsiRt, it0clusterResolution=clustRes, it1clusterResolution=clustRes, it2clusterResolution=clustRes)
		logger.completed()
		logger.start("Aggregating counts across cluster peaks")
			dsan <- regionAggregation(dsan, itLsi$clusterPeaks_unfiltered, ".peaks.itlsi", signal="insertions", dropEmpty=FALSE, bySample=FALSE)
		logger.completed()
	} else {
		logger.warning("Unsupervised analysis will not be performed")
	}

	res <- list(
		ds_anno=dsan,
		itLsi=itLsi,
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
#' @return \code{S3} object containing exploratory metrics and an analysis report object
#' @author Fabian Mueller
#' @export
run_atac_exploratory <- function(dsa, anaDir, itLsiObj=NULL){
	wfState <- getWfState(anaDir)
	doReport <- !wfState$existingReports["exploratory"]

	report <- NULL

	if (doReport){
		logger.start("Creating exploratory report")
			report <- createReport_exploratory(dsa, file.path(wfState$anaDir, wfState$reportDir), itLsiObj=itLsiObj)
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
#' @param input		Input object. Can be either \code{NULL}, a character sting, a \code{DsATAC}
#' @param sampleAnnot sample annotation table (\code{data.frame}) or \code{NULL} if continuing existing analysis or input is a \code{DsATAC} object
#' @param genome    genome assembly. Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param sampleIdCol column name in the sample annotation table containing unique sample Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param regionSets a list of GRanges objects which contain region sets over which count data will be aggregated. Only relevant if not continuing existing analysis and input is not a \code{DsATAC} object
#' @param startStage stage where to start the analysis from. can be one of \code{"raw"}, \code{"filtered"}, \code{"processed"}. Only relevant if not continuing existing analysis.
#' @return \code{\linkS4class{DsATAC}} object (invisible)
#' @author Fabian Mueller
#' @export
run_atac <- function(anaDir, input=NULL, sampleAnnot=NULL, genome=NULL, sampleIdCol=NULL, regionSets=NULL, startStage="raw"){
	doContinue <- prepAnaDir(anaDir)
	wfState <- getWfState(anaDir)
	
	loadConfig(wfState$configPath)

	############################################################################
	# Prepare DsATAC dataset
	############################################################################
	saveDs <- FALSE
	dsa <- NULL
	dsStages <- c("raw", "filtered", "processed")
	if (!is.element(startStage, dsStages)) logger.error(c("Invalid dataset start stage:", startStage))
	if (is.null(input)){
		if (doContinue) {
			logger.info("Continuing ChrAccR analysis")
			logger.start("Loading DsATAC dataset")
				if (!is.na(wfState$dsAtacPaths["processed"])){
					logger.info("Continuing from processed dataset")
					dsa <- loadDsAcc(wfState$dsAtacPaths["processed"])
					startStage <- "processed"
				} else if (!is.na(wfState$dsAtacPaths["filtered"])){
					logger.info("Continuing from filtered dataset")
					dsa <- loadDsAcc(wfState$dsAtacPaths["filtered"])
					startStage <- "filtered"
				} else if (!is.na(wfState$dsAtacPaths["raw"])){
					logger.info("Continuing from raw dataset")
					dsa <- loadDsAcc(wfState$dsAtacPaths["raw"])
					startStage <- "raw"
				} else {
					logger.error("Could not find existing dataset in analysis directory")
				}
			logger.completed()
		} else {
			logger.error("Cannot continue existing analysis: analysis does not exist")
		}
	} else {
		if (is.element(class(input), c("DsATAC", "DsATACsc"))){
			sampleAnnot <- getSampleAnnot(dsa)
		}
		validSampleAnnot <- is.data.frame(sampleAnnot)
		if (is.element(class(input), c("DsATAC", "DsATACsc"))){
			dsa <- input
		} else if (validSampleAnnot){
			saveDs <- TRUE
			inputType <- as.character(NA)
			inputFns <- character(0)
			if (is.character(input)){
				inputFns <- input
				if (length(input) == 1 && is.element(input, colnames(sampleAnnot))){
					inputFns <- sampleAnnot[,input]
					# check if input is cellranger output
					validCellRanger <- sapply(inputFns, FUN=function(x){
						dir.exists(x) && dir.exists(file.path(x, "outs"))
					})
					if (all(validCellRanger)){
						inputType <- "sc_cellranger"
					}
				}
				if (length(inputFns) > 1){
					if (all(grepl("\\.bam$", inputFns))){
						inputType <- "bulk_bam"
					} else if (all(grepl("\\.fragments\\.tsv(\\.gz)?$", inputFns))){
						inputType <- "sc_fragments"
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
					dsa <- DsATAC.bam(sampleAnnot, inputFns, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, diskDump=FALSE, keepInsertionInfo=TRUE, pairedEnd=TRUE)
				} else if (inputType == "sc_fragments"){
					minFrags <- getConfigElement("filteringScMinFragmentsPerCell")
					maxFrags <- getConfigElement("filteringScMaxFragmentsPerCell")
					dsa <- DsATACsc.fragments(sampleAnnot, inputFns, genome, regionSets=regionSets, sampleIdCol=sampleIdCol, minFragsPerBarcode=minFrags, maxFragsPerBarcode=maxFrags, cellAnnot=NULL, keepInsertionInfo=TRUE, cellQcStats=TRUE)
				} else if (inputType == "sc_cellranger"){
					dsa <- DsATAC.cellranger(sampleAnnot, input, genome, dataDir="", regionSets=regionSets, addPeakRegions=TRUE, sampleIdCol=sampleIdCol, keepInsertionInfo=TRUE)
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

	saveDs_raw <- saveDs && is.na(wfState$dsAtacPaths["raw"])
	if (saveDs_raw){
		wfState$dsAtacPaths["raw"] <- "dsATAC_raw"
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
		saveDs_filtered <- saveDs && (doFilter) && is.na(wfState$dsAtacPaths["filtered"])
		if (saveDs_filtered){
			wfState$dsAtacPaths["filtered"] <- "dsATAC_filtered"
			logger.start("Saving filtered DsATAC dataset")
				saveDsAcc(dsa, file.path(wfState$anaDir, wfState$dsAtacPaths["filtered"]))
			logger.completed()
		}
	}
	#---------------------------------------------------------------------------
	# TODO: identify consensus peak set?
	doNorm <- !isSingleCell && is.element(startStage, c("raw", "filtered"))
	if (doNorm){
		logger.start("Running normalization analysis")
			res <- run_atac_normalization(dsa, anaDir)
			dsa <- res$ds_norm
		logger.completed()
	}

	itLsi <- NULL
	doScUnsupervised <- isSingleCell && is.element(startStage, c("raw", "filtered"))
	if (doScUnsupervised){
		logger.start("Running unsupervised single-cell analysis")
			res <- run_atac_sc_unsupervised(dsa, anaDir)
			dsa <- res$ds_anno
			itLsi <- res$itLsi
		logger.completed()
		logger.start("Saving unsupervised results")
			fp <- file.path(wfState$anaDir, wfState$dataDir, "iterativeLSI")
			saveRDS(itLsi, paste0(fp, ".rds"))
			uwot::save_uwot(itLsi$umapRes, paste0(fp, "_uwot"))
		logger.completed()
	}

	saveDs_processed <- saveDs && (doNorm || doScUnsupervised) && is.na(wfState$dsAtacPaths["processed"])
	if (saveDs_processed){
		wfState$dsAtacPaths["processed"] <- "dsATAC_processed"
		logger.start("Saving processed DsATAC dataset")
			saveDsAcc(dsa, file.path(wfState$anaDir, wfState$dsAtacPaths["processed"]))
		logger.completed()
	}
	#---------------------------------------------------------------------------
	doExploratory <- TRUE
	if (doExploratory){
		logger.start("Running exploratory analysis")
			res <- run_atac_exploratory(dsa, anaDir, itLsiObj=itLsi)
		logger.completed()
	}
	#---------------------------------------------------------------------------
	doDifferential <- !isSingleCell && (!is.null(getConfigElement("differentialColumns")) || !is.null(getConfigElement("differentialCompNames")))
	if (doDifferential){
		logger.start("Running differential analysis")
			res <- run_atac_differential(dsa, anaDir)
		logger.completed()
	}
			
	invisible(dsa)
}

if (!isGeneric("createReport_filtering")) {
	setGeneric(
		"createReport_filtering",
		function(.object, ...) standardGeneric("createReport_filtering"),
		signature=c(".object")
	)
}
#' createReport_filtering-methods
#'
#' Create a report summarizing steps and statistics
#'
#' @param .object    filtered \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @param unfilteredObj unfiltered \code{\linkS4class{DsATAC}} object
#' @param filterStats filtering statistics as output by \code{\link{run_atac_filtering}}
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_filtering-DsATAC-method
#' @docType methods
#' @aliases createReport_filtering
#' @aliases createReport_filtering,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("createReport_filtering",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir,
		unfilteredObj,
		filterStats=NULL
	) {
		if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- muReportR::createReport(file.path(reportDir, paste0("filtering", ".html")), "Filtering summary", page.title = "Filtering", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- muReportR::getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- muReportR::getReportDir(rr, dir="data", absolute=TRUE)

		hasFragments <- length(.object@fragments) > 0
		isSingleCell <- class(.object)=="DsATACsc"
		validFilterStats <- !is.null(filterStats)

		regTypes <- intersect(getRegionTypes(.object), getRegionTypes(unfilteredObj))
		if (length(regTypes) < 1) logger.error("Not enough region sets contained in both the filtered and unfiltered object")

		logger.start("Dataset overview section")
			txt <- c(
				"The filtered ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), ifelse(isSingleCell, " cells.", " samples.")
			)
			rr <- muReportR::addReportSection(rr, "Dataset summary", txt, level=1L, collapsed=FALSE)

			txt <- c("Signal has been summarized for the following region sets:")
			rr <- muReportR::addReportParagraph(rr, txt)

			regCountTab <- data.frame(
				"#regions" = sapply(getRegionTypes(.object), FUN=function(rt){getNRegions(.object, rt)}),
				"transformations" =  sapply(getRegionTypes(.object), FUN=function(rt){paste(rev(.object@countTransform[[rt]]), collapse=" -> ")}),
				check.names=FALSE
			)
			rownames(regCountTab) <- getRegionTypes(.object)

			rr <- muReportR::addReportTable(rr, regCountTab, row.names=TRUE, first.col.header=FALSE)

			hasFragments <- length(.object@fragments) > 0
			txt <- c("Fragment data IS NOT available.")
			if (hasFragments) txt <- c("Fragment data IS available.")
			rr <- muReportR::addReportParagraph(rr, txt)
		logger.completed()

		stepDesc <- list()
		if (validFilterStats){
			getParam <- function(name){
				if (is.element(name, names(filterStats$params))){
					return(filterStats$params[[name]])
				} else {
					return(NA)
				}
			}
			stepDesc <- list(
				"region_covg" = paste("Regions with a coverage of less than", getParam("covgCount"), "in more than", round(1-getParam("covgReqSamples")*100,2), "% of samples have been removed."),
				"chromosomes" = paste("Regions and fragments on chromosomes", paste(getParam("exclChroms"), collapse=", "), "have been removed."),
				"sc_minFrags" = paste("Cells with less than", getParam("scMinFragmentsPerCell"), "unique fragments have been removed"),
				"sc_maxFrags" = paste("Cells with more than", getParam("scMaxFragmentsPerCell"), "unique fragments have been removed (they might be doublets)"),
				"sc_tssEnrichment" = paste("Cells with a low TSS enrichment score (lower than", getParam("scMinTssEnrichment"), ") have been removed")
			)

			txt <- c(
				"The following filtering steps were applied to the dataset:"
			)
			rr <- muReportR::addReportSection(rr, "Filtering steps", txt, level=1L, collapsed=FALSE)
			rr <- muReportR::addReportList(rr, stepDesc[filterStats$steps], type="u")
		}
		

		ft <- do.call("rbind", lapply(regTypes, FUN=function(rt){
			data.frame(
				group=rt,
				nBefore=getNRegions(unfilteredObj, rt),
				nAfter=getNRegions(.object, rt),
				stringsAsFactors=FALSE
			)
		}))
		
		filtTabL <- list(
			overall=ft
		)

		if (validFilterStats){
			stepNames <- setdiff(colnames(filterStats$regionStats), c("before", "after"))
			for (sn in stepNames){
				i <- which(colnames(filterStats$regionStats)==sn)
				sn <- gsub("^after_", "", sn)
				filtTabL[[sn]] <- data.frame(
					group=regTypes,
					nBefore=filterStats$regionStats[,i-1],
					nAfter=filterStats$regionStats[,i],
					stringsAsFactors=FALSE
				)
			}
		}

		addStatFun <- function(ft){
			ft[,"filtered"] <- ft[,"nBefore"] - ft[,"nAfter"]
			ft[,"filtered_perc"] <- ft[,"filtered"]/ft[,"nBefore"] * 100
			ft[,"retained"] <- ft[,"nAfter"]
			ft[,"retained_perc"] <- ft[,"retained"]/ft[,"nBefore"] * 100
			return(ft)
		}
		filtTabL <- lapply(filtTabL, addStatFun)
		
		plotPies <- function(ft){
			grps <- ft[,"group"]
			rownames(ft) <- grps
			ft[,"label_filtered"] <- paste0(ft[,"filtered"], "\nfiltered\n(", round(ft[,"filtered_perc"], 2), "%)")
			ft[,"label_retained"] <- paste0(ft[,"retained"], "\nretained\n(", round(ft[,"retained_perc"], 2), "%)")
			nTypes <- length(grps)
			par(mfrow=c(ceiling(sqrt(nTypes)), ceiling(sqrt(nTypes))))
			for (gg in grps){
				pie(c(ft[gg,"filtered"], ft[gg,"retained"]), labels=c(ft[gg,"label_filtered"], ft[gg,"label_retained"]), col=c("#8d3c1e", "#175e54"), clockwise=TRUE, main=gg)
			}
		}

		txt <- character(0)
		if (validFilterStats) txt <- c(txt, "Regions in the datasets were filtered according to the above criteria. ")
		c(txt, "The following plot shows the number of retained and removed regions per region type in the dataset.")
		rr <- muReportR::addReportSection(rr, "Region filtering", txt, level=1L, collapsed=FALSE)

		plotFn <- paste0("regionFilterPieOverall")
		repPlot <- muReportR::createReportPlot(plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=300L)
			plotPies(filtTabL[["overall"]])
		repPlot <- muReportR::off(repPlot)
		rr <- muReportR::addReportFigure(rr, "Overall retained and removed regions per region type", repPlot)

		stepNames <- setdiff(names(filtTabL), c("overall"))
		if (validFilterStats && length(stepNames) > 0){
			txt <- c(
				"The plot below summarizes the number of retained and removed regions per filtering step."
			)
			rr <- muReportR::addReportParagraph(rr, txt)

			repPlotL <- lapply(stepNames, FUN=function(sn){
				plotFn <- paste0("regionFilterPie_", normalize.str(sn, return.camel=TRUE))
				repPlot <- muReportR::createReportPlot(plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=300L)
					plotPies(filtTabL[[sn]])
				repPlot <- muReportR::off(repPlot)
				return(repPlot)
			})

			figSettings.step <- stepNames
			names(figSettings.step) <- normalize.str(stepNames, return.camel=TRUE)
			figSettings <- list(
				"Filtering step" = figSettings.step
			)
			rr <- muReportR::addReportFigure(rr, "Retained and removed regions per region type and filtering step", repPlotL, figSettings)
		}

		if (isSingleCell){
			getNcellsPerSample <- function(dsa){
				sampleIds <- unique(getSampleAnnot(dsa)[,".sampleId"])
				sids <- getSampleAnnot(dsa)[,".sampleId"]
				res <- table(sids)[sampleIds]
				return(res)
			}

			cps_before <- getNcellsPerSample(unfilteredObj)
			cps_after <- getNcellsPerSample(.object)
			sampleIds <- intersect(rownames(cps_after), rownames(cps_before))

			if (length(sampleIds) > 0){
				ft <- data.frame(
					group=sampleIds,
					nBefore=cps_before[sampleIds],
					nAfter=cps_after[sampleIds],
					stringsAsFactors=FALSE
				)
				if (length(sampleIds) > 1){
					ft <- rbind(ft, data.frame(
						group=".ALL",
						nBefore=sum(cps_before),
						nAfter=sum(cps_after),
						stringsAsFactors=FALSE
					))
				}
				ft <- addStatFun(ft)

				txt <- character(0)
				if (validFilterStats) txt <- c(txt, "Cells in the datasets were filtered according to the above criteria. ")
				c(txt, "The following plot shows the number of retained and removed cells per sample in the dataset.")
				rr <- muReportR::addReportSection(rr, "Cell filtering", txt, level=1L, collapsed=FALSE)

				plotFn <- paste0("cellFilterPieOverall")
				repPlot <- muReportR::createReportPlot(plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=300L)
					plotPies(ft)
				repPlot <- muReportR::off(repPlot)
				rr <- muReportR::addReportFigure(rr, "Retained and removed cells", repPlot)
			}
		}

		muReportR::off(rr)
		invisible(rr)
	}
)

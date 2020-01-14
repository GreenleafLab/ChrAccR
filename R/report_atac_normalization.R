if (!isGeneric("createReport_normalization")) {
	setGeneric(
		"createReport_normalization",
		function(.object, ...) standardGeneric("createReport_normalization"),
		signature=c(".object")
	)
}
#' createReport_normalization-methods
#'
#' Create a report summarizing normalization
#'
#' @param .object    normalized \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @param unnormObj  unnormalized \code{\linkS4class{DsATAC}} object
#' @param filterStats filtering statistics as output by \code{\link{run_atac_filtering}}
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_normalization-DsATAC-method
#' @docType methods
#' @aliases createReport_normalization
#' @aliases createReport_normalization,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("createReport_normalization",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir,
		unnormObj
	) {
		if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- muReportR::createReport(file.path(reportDir, paste0("normalization", ".html")), "Normalization summary", page.title = "Normalization", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- muReportR::getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- muReportR::getReportDir(rr, dir="data", absolute=TRUE)

		hasFragments <- length(.object@fragments) > 0
		isSingleCell <- class(.object)=="DsATACsc"
		if (isSingleCell){
			logger.error("Normalization report not supported for single-cell datasets")
		}


		regTypes <- intersect(getRegionTypes(.object), getRegionTypes(unnormObj))
		if (length(regTypes) < 1) logger.error("Not enough region sets contained in both the normalized and unnormalized objects")
		if (!all(getSamples(.object)==getSamples(unnormObj))) logger.error("The normalized and unnormalized objects should have the same samples")

		logger.start("Dataset overview section")
			txt <- c(
				"The normalized ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), " samples."
			)
			rr <- muReportR::addReportSection(rr, "Normalization summary", txt, level=1L, collapsed=FALSE)

			txt <- c("Signal has been normalized for the following region sets:")
			rr <- muReportR::addReportParagraph(rr, txt)

			regCountTab <- data.frame(
				"#regions" = sapply(getRegionTypes(.object), FUN=function(rt){getNRegions(.object, rt)}),
				"normalization" =  sapply(getRegionTypes(.object), FUN=function(rt){paste(rev(.object@countTransform[[rt]]), collapse=" -> ")}),
				check.names=FALSE
			)
			rownames(regCountTab) <- getRegionTypes(.object)

			rr <- muReportR::addReportTable(rr, regCountTab, row.names=TRUE, first.col.header=FALSE)

			# hasFragments <- length(.object@fragments) > 0
			# txt <- c("Fragment data IS NOT available.")
			# if (hasFragments) txt <- c("Fragment data IS available.")
			# rr <- muReportR::addReportParagraph(rr, txt)
		logger.completed()

		
		quantProbs <- seq(0, 1, length.out=101L)
		doLog <- getConfigElement("exploratoryLogNormCounts")

		sampleIds <- getSamples(.object)
		sampleIds.ext <- c("[ALL]", sampleIds)

		plotL.qq <- list()
		plotL.hm <- list()
		for (rt in regTypes){
			rtString <- muRtools::normalize.str(rt, return.camel=TRUE)
			cmU <- getCounts(unnormObj, rt, asMatrix=TRUE)
			cmN <- getCounts(.object, rt, asMatrix=TRUE)
			countLabel <- "count"
			if (doLog) {
				cmU <- log2(cmU + 1)
				cmN <- log2(cmN + 1)
				countLabel <- paste("log2", countLabel,sep="_")
			}

			qmU <- t(matrixStats::colQuantiles(cmU, probs=quantProbs, na.rm=TRUE))
			qmN <- t(matrixStats::colQuantiles(cmN, probs=quantProbs, na.rm=TRUE))

			qU <- quantile(cmU, probs=quantProbs, na.rm=TRUE)
			qN <- quantile(cmN, probs=quantProbs, na.rm=TRUE)
			maxCount <- max(c(max(cmU, na.rm=TRUE), max(cmN, na.rm=TRUE)))

			# QQ-plot of normalized and unnormalized counts (overall and for each sample)
			df2p <- data.frame(reshape2::melt(qmU), value2=reshape2::melt(qmN)[,"value"], stringsAsFactors=FALSE)
			cns <- c("quantile", "sampleId", paste(countLabel, c("unnormalized", "normalized"), sep="_"))
			colnames(df2p) <- cns

			# add overall quantiles
			df_add <- data.frame(
				quantile=names(qN),
				sampleId="[ALL]",
				count_unnormalized=qU,
				count_normalized=qN
			)
			colnames(df_add) <- cns
			df2p <- rbind(df2p, df_add)

			i <- 0
			for (sid in sampleIds.ext){
				pp <- ggplot(df2p[df2p[,"sampleId"]==sid,]) + aes_string(x=paste(countLabel, "unnormalized", sep="_"), y=paste(countLabel, "normalized", sep="_")) +
					  coord_fixed() +
					  geom_abline(intercept=0, slope=1, color="#A0A0A0") + geom_path(lineend="round", size=2) +
					  theme(aspect.ratio=1)
				repPlot <- muReportR::createReportGgPlot(pp, paste0("qqPlot_", rtString, "_s", i), rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
				repPlot <- muReportR::off(repPlot)
				plotL.qq <- c(plotL.qq, list(repPlot))
				i <- i + 1
			}

			print("blubb2")

			# Side-by-side quantile heatmaps
			cs <- getConfigElement("colorSchemesCont")[[".default"]]
			colScheme <- circlize::colorRamp2(seq(0, maxCount, length.out=length(cs)), cs)

			hmu <- ComplexHeatmap::Heatmap(qmU, name=paste("Unnormalized", countLabel),
				col = colScheme,
				column_title = "Unnormalized",
				row_title = "quantile",
				cluster_rows=FALSE,
				cluster_columns=FALSE,
				row_names_gp = grid::gpar(fontsize=6),
				column_names_gp = grid::gpar(fontsize=6),
				show_row_names = TRUE, show_column_names = TRUE
			)
			hmn <- ComplexHeatmap::Heatmap(qmN, name=paste("Normalized", countLabel),
				col = colScheme,
				column_title = "Normalized",
				row_title = "quantile",
				cluster_rows=FALSE,
				cluster_columns=FALSE,
				row_names_gp = grid::gpar(fontsize=6),
				column_names_gp = grid::gpar(fontsize=6),
				show_row_names = TRUE, show_column_names = TRUE
			)
			repPlot <- muReportR::createReportPlot(paste0("quantileHeatmaps_", rtString), rr, width=10, height=7, create.pdf=TRUE, high.png=300L)
			# pdf(fn, width=10, height=10, onefile=FALSE)
				ComplexHeatmap::draw(hmu + hmn)
			repPlot <- muReportR::off(repPlot)
			plotL.hm <- c(plotL.hm, list(repPlot))
		}

		figSettings.region <- regTypes
		names(figSettings.region) <- normalize.str(regTypes, return.camel=TRUE)
		figSettings.sampleId <- sampleIds.ext
		names(figSettings.sampleId) <- paste0("s", 0:length(sampleIds))

		txt <- c(
			"The plot below shows the quantile-quantile relationship of unnormalized to normalized insertion counts, overall and for each sample."
		)
		rr <- muReportR::addReportSection(rr, "QQ plot", txt, level=1L, collapsed=FALSE)

		figSettings <- list(
			"Region type" = figSettings.region,
			"Sample" = figSettings.sampleId
		)
		rr <- muReportR::addReportFigure(rr, "QQ-plot of unnormalized vs normalized counts", plotL.qq, figSettings)

		txt <- c(
			"The heatmaps below show the unnormalized and normalized insertion counts at each quantile for each sample in the dataset."
		)
		rr <- muReportR::addReportSection(rr, "Quantile heatmaps", txt, level=1L, collapsed=FALSE)

		figSettings <- list(
			"Region type" = figSettings.region
		)
		rr <- muReportR::addReportFigure(rr, "Heatmaps of unnormalized vs normalized counts for each quantile and sample", plotL.hm, figSettings)

		muReportR::off(rr)
		invisible(rr)
	}
)

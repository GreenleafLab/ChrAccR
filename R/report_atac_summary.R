if (!isGeneric("createReport_summary")) {
	setGeneric(
		"createReport_summary",
		function(.object, ...) standardGeneric("createReport_summary"),
		signature=c(".object")
	)
}
#' createReport_summary-methods
#'
#' Create a report summarizing an accessibility dataset
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_summary-DsATAC-method
#' @docType methods
#' @aliases createReport_summary
#' @aliases createReport_summary,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("createReport_summary",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		require(muReportR)
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- createReport(file.path(reportDir, "summary.html"), "Accessibility Summary", page.title = "Summary", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- getReportDir(rr, dir="data", absolute=TRUE)

		hasFragments <- length(.object@fragments) > 0

		sampleIds <- getSamples(.object)
		sannot <- getSampleAnnot(.object)
		fn.sannot <- file.path(rDir.data.abs, "sampleAnnot.tsv")
		writeTab(sannot, fn.sannot)

		txt <- c(
			"This ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), " samples. ",
			"The table below contains an overview of the provided sample annotation. "
		)
		rr <- addReportSection(rr, "Overview", txt, level=1L, collapsed=FALSE)

		rr <- addReportTable(rr, sannot, row.names=TRUE, first.col.header=FALSE)

		txt <- c("Sample annotation table as ",  paste(c("<a href=\"", rDir.data, "/", "sampleAnnot.tsv", "\">","TSV file","</a>"),collapse=""))
		rr <- addReportParagraph(rr, txt)

		txt <- c("Signal has been summarized for the following region sets:")
		rr <- addReportParagraph(rr, txt)

		regCountTab <- data.frame(
			"#regions" = sapply(getRegionTypes(.object), FUN=function(rt){getNRegions(.object, rt)}),
			check.names=FALSE
		)
		rownames(regCountTab) <- getRegionTypes(.object)

		rr <- addReportTable(rr, regCountTab, row.names=TRUE, first.col.header=FALSE)
		# ll <- lapply(getRegionTypes(.object), FUN=function(rt){
		# 	paste0("<b>", rt, ":</b> ", getNRegions(.object, rt), " regions")
		# })
		# rr <- addReportList(rr, ll, type="u")


		txt <- c("Fragment data IS NOT available.")
		if (hasFragments) txt <- c("Fragment data IS available.")
		rr <- addReportParagraph(rr, txt)

		if (hasFragments){
			txt <- c(
				"This section contains per-sample quality control plots and data. The table below contains fragment counts for each sample."
			)
			rr <- addReportSection(rr, "Sample QC", txt, level=1L, collapsed=FALSE)

			logger.start("Summarizing fragment counts")
				countTab <- data.frame(
					"#fragments" = getFragmentNum(.object),
					check.names=FALSE
				)
				rownames(countTab) <- sampleIds
				rr <- addReportTable(rr, countTab, row.names=TRUE, first.col.header=FALSE)
			logger.completed()

			# plot fragment numbers
			df2p <- data.frame(
				sample = rownames(countTab),
				nFragments = countTab[,"#fragments"]
			)
			pp <- ggplot(df2p) + aes(sample, nFragments) + geom_col() + 
				  scale_x_discrete(name="") + 
				  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

			figFn <- paste0("fragNumBar")
			repPlot <- createReportGgPlot(pp, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
			repPlot <- off(repPlot, handle.errors=TRUE)
			rr <- addReportFigure(rr, "Number of fragments per sample", repPlot)

			logger.start("Plotting fragment size distribution")
				txt <- c(
					"The following plot illustrates the fragment size distribution for each sample. ",
					"Periodicity indicating nucleosome-bound DNA is generally an indicator of good data quality."
				)
				rr <- addReportSection(rr, "Fragment size distribution", txt, level=2L, collapsed=FALSE)

				plotL <- lapply(1:length(sampleIds), FUN=function(i){
					logger.status(c("Sample", i, "of", length(sampleIds), "..."))
					pp <- plotInsertSizeDistribution(.object, sampleIds[i])
					figFn <- paste0("fragSizeDistr_s", i)
					repPlot <- createReportGgPlot(pp, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
					repPlot <- off(repPlot, handle.errors=TRUE)
					return(repPlot)
				})
				figSettings.sampleId <- sampleIds
				names(figSettings.sampleId) <- paste0("s", 1:length(sampleIds))
				figSettings <- list(
					"Sample" = figSettings.sampleId
				)
				rr <- addReportFigure(rr, "Fragment size distribution", plotL, figSettings)
			logger.completed()

			logger.start("Plotting TSS enrichment")
				# retrieve corresponding TSS coordinates corresponding to the specified gene model
				tssGr <- NULL
				geneModelVer <- getConfigElement("geneModelVersions")[.object@genome]
				if (grepl("^gencode", geneModelVer)){
					tssGr <- muRtools::getAnnotGrl.gencode(geneModelVer)[["gene"]]
					tssGr <- tssGr[elementMetadata(tssGr)[,"gene_type"]=="protein_coding"]
					if (length(tssGr) < 1) logger.error("Could not retrieve TSS coordinates")
					tssGr <- promoters(tssGr, upstream=0, downstream=1)
				} else {
					logger.error("Incompatible gene model")
				}
				txt <- c(
					"The following plot illustrates genome-wide enrichment of ATAC signal around transcription start sites (TSS). ",
					"It is indicative of the signal-to-noise ratio in the dataset. High enrichment of signal in at TSSs compared to ",
					"background indicates good sample quality. Ideally, there is a dip in the TSS profile corresponding the +1 nucleosome."
				)
				rr <- addReportSection(rr, "TSS profile", txt, level=2L, collapsed=FALSE)

				qcTab <- data.frame(
					sample = rownames(countTab),
					nFragments = countTab[,"#fragments"],
					tssEnrich = as.numeric(NA),
					stringsAsFactors=FALSE
				)
				rownames(qcTab) <- qcTab[,"sample"]

				plotL <- list()
				for (i in 1:length(sampleIds)){
					logger.status(c("Sample", i, "of", length(sampleIds), "..."))
					tsse <- getTssEnrichment(.object, sampleIds[i], tssGr)
					figFn <- paste0("tssProfile_s", i)
					repPlot <- createReportGgPlot(tsse$plot, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
					repPlot <- off(repPlot, handle.errors=TRUE)
					plotL <- c(plotL, list(repPlot))

					qcTab[sampleIds[i], "tssEnrich"] <- tsse$tssEnrichment
				})
				figSettings.sampleId <- sampleIds
				names(figSettings.sampleId) <- paste0("s", 1:length(sampleIds))
				figSettings <- list(
					"Sample" = figSettings.sampleId
				)
				desc <- c("TSS profile. Genome-wide aggregate of TSS +/- 2kb. Signal was normalized to the mean signal in the 100-bp-wide windows in the tails of the plot. A smoothing window of width 25bp has been applied to draw the profile curve.")
				rr <- addReportFigure(rr, desc, plotL, figSettings)


				pp <- ggplot(qcTab) + aes(nFragments, tssEnrich) + geom_point() + geom_text(aes(label=sample), size=1)
				figFn <- paste0("qcScatter")
				repPlot <- createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
				repPlot <- off(repPlot, handle.errors=TRUE)
				rr <- addReportFigure(rr, "Scatterplot of QC metrics", repPlot)
			logger.completed()
		}

		off(rr)
		invisible(rr)
	}
)
# report <- createReport_summary(dsr, file.path(config$.anaDir, "reports"))

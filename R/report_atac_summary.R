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
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' reportDir <- file.path(".", "ChrAccR_reports")
#' createReport_summary(dsa, reportDir)
#' }
setMethod("createReport_summary",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- muReportR::createReport(file.path(reportDir, "summary.html"), "Accessibility Summary", page.title = "Summary", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- muReportR::getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- muReportR::getReportDir(rr, dir="data", absolute=TRUE)

		hasFragments <- length(.object@fragments) > 0
		isSingleCell <- class(.object)=="DsATACsc"

		fn.sannot <- file.path(rDir.data.abs, "sampleAnnot.tsv")
		if (isSingleCell) fn.sannot <- file.path(rDir.data.abs, "cellAnnot.tsv")

		sampleIds <- getSamples(.object)
		sannot <- getSampleAnnot(.object)
		writeTab(sannot, fn.sannot)


		if (isSingleCell){
			txt <- c(
				"This single-cell ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), " cells. "
			)
		} else {
			txt <- c(
				"This ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), " samples. ",
				"The table below contains an overview of the provided sample annotation. "
			)
		}
		rr <- muReportR::addReportSection(rr, "Overview", txt, level=1L, collapsed=FALSE)

		if (isSingleCell){
			txt <- c("Cell annotation table as ",  paste(c("<a href=\"", rDir.data, "/", "cellAnnot.tsv", "\">","TSV file","</a>"),collapse=""))
		} else {
			rr <- muReportR::addReportTable(rr, sannot, row.names=TRUE, first.col.header=FALSE)
			txt <- c("Sample annotation table as ",  paste(c("<a href=\"", rDir.data, "/", "sampleAnnot.tsv", "\">","TSV file","</a>"),collapse=""))
		}
		rr <- muReportR::addReportParagraph(rr, txt)
		
		txt <- c("Signal has been summarized for the following region sets:")
		rr <- muReportR::addReportParagraph(rr, txt)

		regCountTab <- data.frame(
			"#regions" = sapply(getRegionTypes(.object), FUN=function(rt){getNRegions(.object, rt)}),
			"transformations" =  sapply(getRegionTypes(.object), FUN=function(rt){paste(rev(.object@countTransform[[rt]]), collapse=" -> ")}),
			check.names=FALSE
		)
		rownames(regCountTab) <- getRegionTypes(.object)

		rr <- muReportR::addReportTable(rr, regCountTab, row.names=TRUE, first.col.header=FALSE)
		# ll <- lapply(getRegionTypes(.object), FUN=function(rt){
		# 	paste0("<b>", rt, ":</b> ", getNRegions(.object, rt), " regions")
		# })
		# rr <- addReportList(rr, ll, type="u")

		txt <- c("Fragment data IS NOT available.")
		if (hasFragments) txt <- c("Fragment data IS available.")
		rr <- muReportR::addReportParagraph(rr, txt)

		if (isSingleCell){
			txt <- c(
				"This section provides an overview on single-cell QC statistics. Summary plots show the distribution of statistics per sample."
			)
			rr <- muReportR::addReportSection(rr, "Cell QC", txt, level=1L, collapsed=FALSE)

			summaryDf <- getScQcStatsTab(.object)
			cns <- setdiff(colnames(summaryDf), c("cell", "sample"))

			sampleStatsTab <- data.frame(
				sample=unique(summaryDf[,"sample"]),
				stringsAsFactors=FALSE
			)
			sampleStatsTab[,"nCells"] <- table(summaryDf[,"sample"])[sampleStatsTab[,"sample"]]
			for (cn in cns){
				infix <- ""
				if (!is.element(cn, c("tssEnrichment"))){
					infix <- "fragment_"
				}
				sampleSummaryCn <- paste0("median_", infix, cn)
				sampleStatsTab[,sampleSummaryCn] <- sapply(sampleStatsTab[,"sample"], FUN=function(ss){
					median(summaryDf[summaryDf[,"sample"]==ss, cn], na.rm=TRUE)
				})
			}

			sampleAddCns <- c(
				tssEnrichment=findOrderedNames(colnames(sannot), ".CR.sampleQC.tss_enrichment_score"),
				totalFragments=findOrderedNames(colnames(sannot), ".CR.sampleQC.total_usable_fragments")
			)
			for (cn in names(sampleAddCns)){
				if (!is.na(sampleAddCns[cn])){
					sampleAggr <- tapply(sannot[,sampleAddCns[cn]], summaryDf[,"sample"], FUN=function(x){
						uu <- unique(x)
						if (length(uu)==1){
							return(uu)
						} else {
							return(NA)
						}
					})
					if (!all(is.na(sampleAggr))) sampleStatsTab[, paste0("sample_", cn)] <- sampleAggr[sampleStatsTab[,"sample"]]
				}
			}
			rr <- muReportR::addReportTable(rr, sampleStatsTab, row.names=FALSE, first.col.header=TRUE)
			plotL <- lapply(cns, FUN=function(cn){
				pp <- ggplot(summaryDf) + aes_string(x="sample", y=cn) +
				      geom_violin(adjust=1, fill="#4d4f53") +
				      geom_boxplot(aes(fill=NULL), outlier.shape=NA, width=0.2) +
				      ylim(quantile(summaryDf[,cn], 0.98)) + coord_flip() +
				      theme(axis.title.y=element_blank()) + guides(fill=FALSE) 
				figFn <- paste0("statDistr_", cn)
				repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
				repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
				return(repPlot)
			})
			figSettings.stat <- cns
			names(figSettings.stat) <- cns
			figSettings <- list(
				"Statistic" = figSettings.stat
			)
			rr <- muReportR::addReportFigure(rr, "Distribution of cell statistics", plotL, figSettings)

			if (all(c("tssEnrichment") %in% cns)){

				txt <- c("The plot below shows two frequently used statistics that can be used to asses per-cell quality. ",
					     "The normalized enrichment of insertions at transcription start sites is plotted against the absolute number of fragments."
				)
				rr <- muReportR::addReportParagraph(rr, txt)

				cut_x <- getConfigElement("filteringScMinFragmentsPerCell")
				cut_y <- getConfigElement("filteringScMinTssEnrichment")

				sampleIds <- sampleStatsTab[,"sample"]
				plotL <- lapply(1:length(sampleIds), FUN=function(i){
					subDf <- summaryDf[summaryDf[,"sample"]==sampleIds[i],]
					df2p <- data.frame(
						log_nPass = log10(subDf[,"nPass"]),
						tssEnrichment = subDf[,"tssEnrichment"]
					)
					df2p[, "pointDens"] <- muRtools::getPointDensity(df2p[, "log_nPass"], df2p[, "tssEnrichment"], n=100)

					pp <- ggplot(df2p) + aes_string(x="log_nPass", y="tssEnrichment", color="pointDens")
					if (!is.null(cut_x)) pp <- pp + geom_vline(xintercept=cut_x, color="#4d4f53")
					if (!is.null(cut_y)) pp <- pp + geom_hline(yintercept=cut_y, color="#4d4f53")
					pp <- pp + geom_point(size=0.5) + muRtools::ggAutoColorScale(df2p[, "pointDens"], method="color")

					# pp <- muRtools::create.densityScatter(df2p, is.special=NULL, sparse.points=0.01)
					figFn <- paste0("sampleTssEnrich_s", i)
					repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
					repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
					return(repPlot)
				})
				figSettings.sampleId <- sampleIds
				names(figSettings.sampleId) <- paste0("s", 1:length(sampleIds))
				figSettings <- list(
					"Sample" = figSettings.sampleId
				)
				rr <- muReportR::addReportFigure(rr, "Cell QC statistics: number of fragments vs. TSS enrichment", plotL, figSettings)
			}
		} else {
			if (hasFragments){
				txt <- c(
					"This section contains per-sample quality control plots and data. The table below contains fragment counts for each sample."
				)
				rr <- muReportR::addReportSection(rr, "Sample QC", txt, level=1L, collapsed=FALSE)

				logger.start("Summarizing fragment counts")
					countTab <- data.frame(
						"#fragments" = getFragmentNum(.object),
						check.names=FALSE
					)
					rownames(countTab) <- sampleIds
					rr <- muReportR::addReportTable(rr, countTab, row.names=TRUE, first.col.header=FALSE)
				logger.completed()

				# plot fragment numbers
				df2p <- data.frame(
					sample = factor(rownames(countTab), levels=rownames(countTab)[order(countTab[,"#fragments"], decreasing=TRUE)]),
					nFragments = countTab[,"#fragments"]
				)
				pp <- ggplot(df2p) + aes(sample, nFragments) + geom_col() + 
					  scale_x_discrete(name="") + 
					  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

				figFn <- paste0("fragNumBar")
				repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
				repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
				rr <- muReportR::addReportFigure(rr, "Number of fragments per sample", repPlot)

				logger.start("Plotting fragment size distribution")
					txt <- c(
						"The following plot illustrates the fragment size distribution for each sample. ",
						"Periodicity indicating nucleosome-bound DNA is generally an indicator of good data quality."
					)
					rr <- muReportR::addReportSection(rr, "Fragment size distribution", txt, level=2L, collapsed=FALSE)

					plotL <- lapply(1:length(sampleIds), FUN=function(i){
						logger.status(c("Sample", i, "of", length(sampleIds), "..."))
						pp <- plotInsertSizeDistribution(.object, sampleIds[i])
						figFn <- paste0("fragSizeDistr_s", i)
						repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
						repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
						return(repPlot)
					})
					figSettings.sampleId <- sampleIds
					names(figSettings.sampleId) <- paste0("s", 1:length(sampleIds))
					figSettings <- list(
						"Sample" = figSettings.sampleId
					)
					rr <- muReportR::addReportFigure(rr, "Fragment size distribution", plotL, figSettings)
				logger.completed()

				logger.start("Plotting TSS enrichment")
					# retrieve corresponding TSS coordinates corresponding to the specified gene model
					tssGr <- NULL
					geneModelVer <- getConfigElement("geneModelVersions")[.object@genome]
					annoPkg <- getChrAccRAnnotationPackage(.object@genome)
					if (grepl("^gencode", geneModelVer)){
						downloadGencode <- TRUE
						# gencode version already annotated?
						if (!is.null(annoPkg)){
							geneAnnoName <- "gencode_coding"
							annoParams <- get("getGenomeParams", asNamespace(annoPkg))()
							useAnno <- is.element(geneAnnoName, names(annoParams$geneAnno)) && annoParams$geneAnno[[geneAnnoName]]$version==geneModelVer
							if (useAnno) {
								logger.info(c("Using annotation package:", annoPkg))
								tssGr <- get("getGeneAnnotation", asNamespace(annoPkg))(anno=geneAnnoName, type="tssGr")
								downloadGencode <- FALSE
							}
						}
						# gencode version not annotated --> download
						if (downloadGencode){
							tssGr <- muRtools::getAnnotGrl.gencode(geneModelVer)[["gene"]]
							tssGr <- tssGr[elementMetadata(tssGr)[,"gene_type"]=="protein_coding"]
							if (length(tssGr) < 1) logger.error("Could not retrieve TSS coordinates")
							tssGr <- promoters(tssGr, upstream=0, downstream=1)
						}
					} else {
						logger.error("Incompatible gene model")
					}
					txt <- c(
						"The following plot illustrates genome-wide enrichment of ATAC signal around transcription start sites (TSS). ",
						"It is indicative of the signal-to-noise ratio in the dataset. High enrichment of signal in at TSSs compared to ",
						"background indicates good sample quality. Ideally, there is a dip in the TSS profile corresponding the +1 nucleosome."
					)
					rr <- muReportR::addReportSection(rr, "TSS profile", txt, level=2L, collapsed=FALSE)

					qcTab <- data.frame(
						sample = rownames(countTab),
						nFragments = countTab[,"#fragments"],
						tssEnrichment = as.numeric(NA),
						stringsAsFactors=FALSE
					)
					rownames(qcTab) <- qcTab[,"sample"]

					plotL <- list()
					for (i in 1:length(sampleIds)){
						logger.status(c("Sample", i, "of", length(sampleIds), "..."))
						tsse <- getTssEnrichment(.object, sampleIds[i], tssGr)
						figFn <- paste0("tssProfile_s", i)
						repPlot <- muReportR::createReportGgPlot(tsse$plot, figFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
						repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
						plotL <- c(plotL, list(repPlot))

						qcTab[sampleIds[i], "tssEnrichment"] <- tsse$tssEnrichment
					}
					figSettings.sampleId <- sampleIds
					names(figSettings.sampleId) <- paste0("s", 1:length(sampleIds))
					figSettings <- list(
						"Sample" = figSettings.sampleId
					)
					desc <- c("TSS profile. Genome-wide aggregate of TSS +/- 2kb. Signal was normalized to the mean signal in the 100-bp-wide windows in the tails of the plot. A smoothing window of width 25bp has been applied to draw the profile curve.")
					rr <- muReportR::addReportFigure(rr, desc, plotL, figSettings)


					pp <- ggplot(qcTab) + aes(nFragments, tssEnrichment) + geom_point() + geom_text(aes(label=sample), size=1)
					figFn <- paste0("qcScatter")
					repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
					repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
					rr <- muReportR::addReportFigure(rr, "Scatterplot of QC metrics", repPlot)
				logger.completed()
			}
		}

		muReportR::off(rr)
		invisible(rr)
	}
)

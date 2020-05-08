if (!isGeneric("createReport_differential")) {
	setGeneric(
		"createReport_differential",
		function(.object, ...) standardGeneric("createReport_differential"),
		signature=c(".object")
	)
}
#' createReport_differential-methods
#'
#' Create a report summarizing differential accessibility analysis
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_differential-DsATAC-method
#' @docType methods
#' @aliases createReport_differential
#' @aliases createReport_differential,DsATAC-method
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' reportDir <- file.path(".", "ChrAccR_reports")
#' setConfigElement("regionTypes", setdiff(getRegionTypes(dsa), c("promoters_gc_protein_coding", "t10k")))
#' setConfigElement("differentialColumns", c("stimulus", "cellType"))
#' # adjust for the donor annotation in the differential test
#' setConfigElement("differentialAdjColumns", c("donor"))
#' # create the report
#' createReport_differential(dsa, reportDir)
#' }
setMethod("createReport_differential",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- muReportR::createReport(file.path(reportDir, "differential.html"), "Differential Accessibility", page.title = "Differential", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- muReportR::getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- muReportR::getReportDir(rr, dir="data", absolute=TRUE)
		
		compTab <- getComparisonTable(.object,
			cols=getConfigElement("differentialColumns"),
			compNames=getConfigElement("differentialCompNames"),
			cols1vAll=getConfigElement("differentialColumns1vsAll"),
			minGroupSize=getConfigElement("annotationMinGroupSize"),
			maxGroupCount=getConfigElement("annotationMaxGroupCount")
		)

		if (is.null(compTab)) logger.error("No valid comparisons found")
		if (nrow(compTab) > 10) logger.warning("An extensive amount of comparisons will be performed. Consider being more specific in the differentialColumns and differentialCompNames options.")

		regionTypes <- getRegionTypes(.object)
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		wasTransformed <- sapply(regionTypes, FUN=function(rt){length(.object@countTransform[[rt]])>0})
		if (any(wasTransformed)){
			logger.warning(c("Detected tranformed count data. It is not recommended to compute differential accessibility from normalized counts"))
		}

		logger.start("Computing differential accessibility")
			adjCols <- getConfigElement("differentialAdjColumns")
			logger.start("Differential accessibility objects")
				diffObjL <- lapply(regionTypes, FUN=function(rt){
					logger.status(c("Region type:", rt))
					# do <- ChrAccR:::getDESeq2Dataset(.object, rt, designCols)
					do <- getDESeq2Dataset(.object, rt, adjCols, compTab=compTab)
					fn <- file.path(rDir.data.abs, paste0("diffObj_", normalize.str(rt, return.camel=TRUE), ".rds"))
					saveRDS(do, fn)
					return(do)
				})
				names(diffObjL) <- regionTypes
			logger.completed()
			logger.start("Differential accessibility tables")
				diffTabL <- lapply(regionTypes, FUN=function(rt){
					logger.status(c("Region type:", rt))
					gr <- ChrAccR::getCoord(.object, rt)
					coordDf <- data.frame(
						chrom=seqnames(gr),
						chromStart=start(gr)-1,
						chromEnd=end(gr),
						strand=strand(gr)
					)
					ll <- lapply(1:nrow(compTab), FUN=function(i){
						dat <- getDiffAcc(
							.object, rt, compTab[i,"compCol"], grp1Name=compTab[i,"grp1Name"], grp2Name=compTab[i,"grp2Name"],
							adjustCols=setdiff(adjCols, compTab[i,"compCol"]), # does not matter, the adjustment columns are already in the diffObj
							method='DESeq2',
							diffObj=diffObjL[[rt]]
						)
						dm.coord <- data.frame(
							coordDf,
							dat
						)
						fn <- file.path(rDir.data.abs, paste0("diffTab_", i, "_", normalize.str(rt, return.camel=TRUE), ".tsv"))
						write.table(dm.coord, fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
						return(dat)
					})
					names(ll) <- compTab[,"compName"]
					return(ll)
				})
				names(diffTabL) <- regionTypes
			logger.completed()
		logger.completed()

		isDiffFuns <- list(
			cutL2fc2Padj05 = function(dm){
				abs(dm[,"log2FoldChange"]) > 2 & dm[,"padj"] < 0.05
			},
			cutL2fc2Padj05gain = function(dm){
				dm[,"log2FoldChange"] > 2 & dm[,"padj"] < 0.05
			},
			cutL2fc2Padj05loss = function(dm){
				dm[,"log2FoldChange"] < -2 & dm[,"padj"] < 0.05
			},
			cRankTopPerc1 = function(dm){
				dm[,"cRank_rerank"] < quantile(dm[,"cRank_rerank"], prob=0.01, na.rm=TRUE)
			},
			cRankTopPerc5 = function(dm){
				dm[,"cRank_rerank"] < quantile(dm[,"cRank_rerank"], prob=0.05, na.rm=TRUE)
			}
		)
		diffFunDesc <- c(
			cutL2fc2Padj05     = "Differential [|log2(fold change)| > 2; adj. p-value < 0.05]",
			cutL2fc2Padj05gain = "Gain [log2(fold change) > 2; adj. p-value < 0.05]",
			cutL2fc2Padj05loss = "Loss [log2(fold change) < -2; adj. p-value < 0.05]",
			cRankTopPerc1      = "Differential [combined rank in top 1%]",
			cRankTopPerc5      = "Differential [combined rank in top 5%]"
		)

		txt <- c(
			"Differential Accessibility was quantified for the following comparisons:"
		)
		rr <- muReportR::addReportSection(rr, "Comparisons", txt, level=1L, collapsed=FALSE)

		cTab <- compTab[,c("compName", "compCol", "grp1Name", "nGrp1", "grp2Name", "nGrp2")]
		rownames(cTab) <- as.character(1:nrow(cTab))
		colnames(cTab) <- c("Comparison name", "Annotation column", "Group 1", "N1", "Group 2", "N2")
		rr <- muReportR::addReportTable(rr, cTab)

		txt <- c("Tables containing differential accessibility results can be found below.")
		rr <- muReportR::addReportParagraph(rr, txt)

		tabFileTab <- do.call("rbind",lapply(1:nrow(compTab),FUN=function(i){
			sapply(regionTypes,FUN=function(rt){
				fn <- file.path(rDir.data, paste0("diffTab_", i, "_", normalize.str(rt, return.camel=TRUE), ".tsv"))
				txt <- paste(c("<a href=\"", fn, "\">","TSV","</a>"),collapse="")
				return(txt)
			})
		}))
		rownames(tabFileTab) <- compTab[,"compName"]
		colnames(tabFileTab) <- regionTypes
		rr <- muReportR::addReportTable(rr, tabFileTab)

		logger.start("Plotting")
			plotL.ma <- list()
			for (rt in regionTypes){
				logger.start(c("Region type:", rt))
					for (i in 1:nrow(compTab)){
						logger.start(c("comparison", i, ":", compTab[i,"compName"]))
							dm <- diffTabL[[rt]][[i]]
							df2p.ma <- dm[,c("log2BaseMean", "log2FoldChange")]
							for (funName in names(isDiffFuns)){
								logger.start(c("Comparing using method", diffFunDesc[funName]))
									isDiff <- isDiffFuns[[funName]](dm)
									isDiff[is.na(isDiff)] <- FALSE
									pp <- muRtools::create.densityScatter(df2p.ma, is.special=isDiff, sparse.points=0.001)
									plotFn <- paste0("maPlot_", i, "_", normalize.str(rt, return.camel=TRUE), "_", funName)
									repPlot <- muReportR::createReportGgPlot(pp, plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
									repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
									plotL.ma <- c(plotL.ma, list(repPlot))
								logger.completed()
							}
						logger.completed()
					}
				logger.completed()
			}

			figSettings.comp <- compTab[,"compName"]
			names(figSettings.comp) <- as.character(1:nrow(compTab))
			figSettings.region <- regionTypes
			names(figSettings.region) <- normalize.str(regionTypes, return.camel=TRUE)	
			figSettings.diffFun <- diffFunDesc[names(isDiffFuns)]
			names(figSettings.diffFun) <- names(isDiffFuns)
			figSettings <- list(
				"Comparison" = figSettings.comp,
				"Region type" = figSettings.region,
				"Differential method" = figSettings.diffFun
			)
			rr <- muReportR::addReportFigure(rr, "MA plot", plotL.ma, figSettings)

			lolaDbPaths <- getConfigElement("lolaDbPaths")
			doLola <- !is.null(lolaDbPaths) && all(dir.exists(lolaDbPaths))
			if (doLola){
				logger.start("LOLA analysis")
					logger.start("Preparing LOLA database")
						lolaDb <- loadLolaDbs(lolaDbPaths)
					logger.completed()


					lolaRefTxt <- c("Sheffield, & Bock (2016). LOLA: enrichment analysis for genomic region sets and regulatory elements in R and Bioconductor. <i>Bioinformatics</i>, <b>32</b>(4), 587-589.")
					rr <- muReportR::addReportReference(rr, lolaRefTxt)
					txt <- c(
						"LOLA enrichment analysis ", muReportR::getReportReference(rr, lolaRefTxt), " for differentially accessible regions."
					)
					rr <- muReportR::addReportSection(rr, "LOLA enrichment analysis", txt, level=2L, collapsed=FALSE)

					logger.start("Running LOLA and plotting")
						plotL.lb <- list()
						for (rt in regionTypes){
							logger.start(c("Region type:", rt))
								gr <- ChrAccR::getCoord(.object, rt)
								for (i in 1:nrow(compTab)){
									logger.start(c("comparison", i, ":", compTab[i,"compName"]))
										dm <- diffTabL[[rt]][[i]]
										for (funName in names(isDiffFuns)){
											#[TODO] slow: currently takes ~6min per method on a large database -->parallelize (maybe one analysis with different user sets per region type)
											logger.start(c("Comparing using method", diffFunDesc[funName]))
												isDiff <- isDiffFuns[[funName]](dm)
												isDiff[is.na(isDiff)] <- FALSE

												curSuffix <- paste0(i, "_", normalize.str(rt, return.camel=TRUE), "_", funName)
												lolaResFn <- file.path(rDir.data.abs, paste0("lolaRes_", curSuffix, ".rds"))
												
												lolaRes <- LOLA::runLOLA(gr[isDiff], gr, lolaDb, cores=8)
												saveRDS(lolaRes, lolaResFn)

												pp <- lolaBarPlot(lolaDb, lolaRes, scoreCol="log2OR", orderCol="maxRnk", pvalCut=0.01, maxTerms=200)
												plotFn <- paste0("lolaBar_", curSuffix)
												repPlot <- muReportR::createReportGgPlot(pp, plotFn, rr, width=20, height=5, create.pdf=TRUE, high.png=0L)
												repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
												plotL.lb <- c(plotL.lb, list(repPlot))
											logger.completed()
										}
									logger.completed()
								}
							logger.completed()
						}						
					logger.completed()
					figDesc <- c(
						"LOLA bar plot for enrichment of region annotations. ",
						"Bar height denotes log2(odds-ratio). Bars are ordered according to the combined ranking in the LOLA result. ",
						"Up to 200 most enriched categories (q-value < 0.01; Fisher's Exact Test) are shown for each comparison."
					)
					rr <- muReportR::addReportFigure(rr, figDesc, plotL.lb, figSettings)

				logger.completed()
			}
		logger.completed()


		muReportR::off(rr)
		invisible(rr)
	}
)


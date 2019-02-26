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
setMethod("createReport_differential",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		require(muReportR)
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- createReport(file.path(reportDir, "differential.html"), "Differential Accessibility", page.title = "Differential", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- getReportDir(rr, dir="data", absolute=TRUE)

		sampleIds <- getSamples(.object)
		sannot <- getSampleAnnot(.object)
		sampleGrps <- getGroupsFromTable(sannot, cols=getConfigElement("differentialColumns"))
		if (length(sampleGrps) < 1) logger.error("No valid comparisons found (to begin with)")
		compTab <- do.call("rbind", lapply(1:length(sampleGrps), FUN=function(i){
			tt <- NULL
			if (length(sampleGrps[[i]]) == 2) {
				tt <- data.frame(
					compName=paste0(names(sampleGrps[[i]])[1], " vs ", names(sampleGrps[[i]])[2],  " [", names(sampleGrps)[i], "]"),
					compCol=names(sampleGrps)[i],
					grp1Name=names(sampleGrps[[i]])[1],
					grp2Name=names(sampleGrps[[i]])[2],
					stringsAsFactors=FALSE
				)
			} else if (length(sampleGrps[[i]]) > 2) {
				grpNames <- t(combn(names(sampleGrps[[i]]), 2))
				tt <- data.frame(
					compName=paste0(grpNames[,1], " vs ", grpNames[,2],  " [", names(sampleGrps)[i], "]"),
					compCol=names(sampleGrps)[i],
					grp1Name=grpNames[,1],
					grp2Name=grpNames[,2],
					stringsAsFactors=FALSE
				)
			}
			return(tt)
		}))
		if (is.null(compTab)) logger.error("No valid comparisons found")
		if (nrow(compTab) > 10) logger.warning("An extensive amount of comparisons will be performed. Consider being more specific in the differentialColumns option.")

		regionTypes <- getRegionTypes(.object)
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		logger.start("Computing differential accessibility")
			designCols <- c(getConfigElement("differentialAdjColumns"), unique(compTab[,"compCol"]))
			logger.start("Differential accessibility objects")
				diffObjL <- lapply(regionTypes, FUN=function(rt){
					logger.status(c("Region type:", rt))
					do <- getDESeq2Dataset(.object, rt, designCols)
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
							adjustCols=setdiff(designCols, compTab[i,"compCol"]), # does not matter, the adjustment columns are already in the diffObj
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
				dm[,"cRank_rerank"] < quantile(dm[,"cRank_rerank"], prob=0.01)
			},
			cRankTopPerc5 = function(dm){
				dm[,"cRank_rerank"] < quantile(dm[,"cRank_rerank"], prob=0.05)
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
			"Differential Accessibility was quantified for the following comparisons. ",
			"Tables containing differential accessibility results can be found below."
		)
		rr <- addReportSection(rr, "Comparisons", txt, level=1L, collapsed=FALSE)

		tabFileTab <- do.call("rbind",lapply(1:nrow(compTab),FUN=function(i){
			sapply(regionTypes,FUN=function(rt){
				fn <- file.path(rDir.data, paste0("diffTab_", i, "_", normalize.str(rt, return.camel=TRUE), ".tsv"))
				txt <- paste(c("<a href=\"", fn, "\">","TSV","</a>"),collapse="")
				return(txt)
			})
		}))
		rownames(tabFileTab) <- compTab[,"compName"]
		colnames(tabFileTab) <- regionTypes
		rr <- addReportTable(rr, tabFileTab)

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
									pp <- create.densityScatter(df2p.ma, is.special=isDiff, sparse.points=0.001)
									plotFn <- paste0("maPlot_", normalize.str(rt, return.camel=TRUE), "_", i, "_", funName)
									repPlot <- createReportGgPlot(pp, plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
									repPlot <- off(repPlot, handle.errors=TRUE)
									plotL.ma <- c(plotL.ma, list(repPlot))
								logger.completed()
							}
						logger.completed()
					}
				logger.completed()
			}

			figSettings.region <- regionTypes
			names(figSettings.region) <- normalize.str(regionTypes, return.camel=TRUE)
			figSettings.comp <- compTab[,"compName"]
			names(figSettings.comp) <- as.character(1:nrow(compTab))
			figSettings.diffFun <- diffFunDesc[names(isDiffFuns)]
			names(figSettings.diffFun) <- names(isDiffFuns)
			figSettings <- list(
				"Region type" = figSettings.region,
				"Comparison" = figSettings.comp,
				"Differential method" = figSettings.diffFun
			)
			rr <- addReportFigure(rr, "MA plot", plotL.ma, figSettings)

			lolaDbPaths <- getConfigElement("lolaDbPaths")
			doLola <- !is.null(lolaDbPaths) && all(dir.exists(lolaDbPaths))
			if (doLola){
				logger.start("LOLA analysis")
					logger.start("Preparing LOLA database")
						lolaDb <- loadLolaDbs(lolaDbPaths)
					logger.completed()
					
				logger.completed()
			}
		logger.completed()


		off(rr)
		invisible(rr)
	}
)


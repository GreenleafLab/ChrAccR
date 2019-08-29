if (!isGeneric("createReport_exploratory")) {
	setGeneric(
		"createReport_exploratory",
		function(.object, ...) standardGeneric("createReport_exploratory"),
		signature=c(".object")
	)
}
#' createReport_exploratory-methods
#'
#' Create a report summarizing exploratory analyses of an accessibility dataset
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_exploratory-DsATAC-method
#' @docType methods
#' @aliases createReport_exploratory
#' @aliases createReport_exploratory,DsATAC-method
#' @author Fabian Mueller
#' @export
#' 
#' @examples
#' \dontrun{
#' dsa <- ChrAccRex::loadExample("dsAtac_ia_example")
#' dsa_qnorm <- transformCounts(dsa, method="quantile")
#' setConfigElement("annotationColumns", c("cellType", "donor", "stimulus"))
#' setConfigElement("regionTypes", setdiff(getRegionTypes(dsa), c("promoters_gc_protein_coding", "t10k")))
#' reportDir <- file.path(".", "ChrAccR_reports")
#' createReport_exploratory(dsa_qnorm, reportDir)
#' }
setMethod("createReport_exploratory",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		# reportDir <- file.path("~/myscratch/temp", getHashString("report"))
		if (!requireNamespace("muReportR")) logger.error(c("Could not load dependency: muReportR"))
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- muReportR::createReport(file.path(reportDir, "exploratory.html"), "Exploratory Analysis", page.title="Exploratory", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- muReportR::getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- muReportR::getReportDir(rr, dir="data", absolute=TRUE)

		isSingleCell <- class(.object)=="DsATACsc"

		regionTypes <- getRegionTypes(.object)
		specAnnotCols <- getConfigElement("annotationColumns")
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		sannot <- getSampleAnnot(.object)
		sannot[,".ALL"] <- "all"

		cellIds <- getSamples(.object)
		if (isSingleCell){
			findItLsiRt <- function(ds){
				return(findOrderedNames(getRegionTypes(ds), c("tiling", "^t[1-9]"), exact=FALSE, ignore.case=TRUE))
			}
			itLsiRt <- getConfigElement("scIterativeLsiRegType")
			if (length(itLsiRt) < 1) {
				itLsiRt <- findItLsiRt(.object)
				if (is.na(itLsiRt)) {
					logger.error("Could not determine region type for iterative LSI")
				}
			}
			logger.start("Dimension reduction using iterative LSI")
				logger.info(c("Using region type:", itLsiRt))
				dre <- iterativeLSI(.object, it0regionType=itLsiRt, it0clusterResolution=0.6, it1clusterResolution=0.6, it2clusterResolution=0.6)
				logger.status("Saving ...")
				saveRDS(dre, file.path(rDir.data.abs, "dimRed_iterativeLSI_res.rds"))
				uwot::save_uwot(dre$umapRes, file.path(rDir.data.abs, "dimRed_iterativeLSI_res_uwot"))
			logger.completed()
			
			# annotation
			qcDf <- getScQcStatsTab(.object)
			cellIds <- qcDf[,"cell"]
			sannot[,"clusterAssignment"] <- dre$clustAss[cellIds]
			pcns <- paste0("PC", 1:3)
			for (pcn in pcns){
				sannot[,pcn] <- dre$pcaCoord[cellIds, pcn]
			}
			qcAnnotCols <- c(setdiff(colnames(summaryDf), c("cell", "sample")), pcns)
			specAnnotCols <- c(specAnnotCols, qcAnnotCols, "clusterAssignment")
		}

		sampleGrps <- getGroupsFromTable(sannot, cols=specAnnotCols, minGrpSize=getConfigElement("annotationMinGroupSize"))
		sampleGrps <- c(list(".ALL"=c(1:nrow(sannot))), sampleGrps)
		grpNames <- names(sampleGrps)

		numericCols <- sapply(colnames(sannot), FUN=function(cn){is.numeric(sannot[,cn])})
		numericCols <- names(numericCols)[numericCols]
		numericCols <- setdiff(numericCols, grpNames)
		if (length(specAnnotCols) > 0) numericCols <- intersect(numericCols, specAnnotCols)
		plotAnnotCols <- c(grpNames, numericCols)


		colSchemes <- getConfigElement("colorSchemes")
		colSchemesNum <- getConfigElement("colorSchemesCont")
		grpColors <- lapply(plotAnnotCols, FUN=function(cn){
			cs <- c()
			x <- NULL
			if (is.element(cn, names(sampleGrps))) x <- sampleGrps[[cn]]
			isNum <- is.numeric(sannot[,cn])
			useDefault <- (!isNum && !is.element(cn, names(colSchemes))) || (isNum && !is.element(cn, names(colSchemesNum)))
			if (!useDefault && !isNum) {
				useDefault <- !all(names(x) %in% names(colSchemes[[cn]]))
			}
			if (useDefault) {
				if (isNum){
					cs <- colSchemesNum[[".default"]]
				} else {
					cs <- rep(colSchemes[[".default"]], length.out=length(x))
					names(cs) <- names(x)
				}
			} else {
				if (isNum){
					cs <- colSchemesNum[[cn]]
				} else {
					cs <- colSchemes[[cn]][names(x)]
				}
			}
			return(cs)
		})
		names(grpColors) <- plotAnnotCols

		logger.start("Dataset overview section")
			txt <- c(
				"This ATAC-seq dataset contains accessibility profiles for ", length(getSamples(.object)), ifelse(isSingleCell, " cells.", " samples.")
			)
			rr <- muReportR::addReportSection(rr, "Overview", txt, level=1L, collapsed=FALSE)

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

		if (isSingleCell){
			logger.start("Dimension reduction plots")
				txt <- c(
					"Iterative LSI was used for clustering and for dimension reduction."
				)
				rr <- muReportR::addReportSection(rr, "Dimension reduction", txt, level=1L, collapsed=FALSE)

				mnames <- "umap"
				plotL <- list()
				for (gn in plotAnnotCols){
					pp <- getDimRedPlot(dre$umapCoord[cellIds,], annot=sannot, colorCol=cn, shapeCol=FALSE, colScheme=grpColors[[gn]], ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + coord_fixed()
					plotFn <- paste("dimRed", "umap", normalize.str(gn, return.camel=TRUE), sep="_")
					repPlot <- muReportR::createReportGgPlot(pp, plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
					repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
					plotL <- c(plotL, list(repPlot))
				}
				figSettings.annot <- plotAnnotCols
				names(figSettings.annot) <- normalize.str(plotAnnotCols, return.camel=TRUE)
				figSettings.method <- mnames
				names(figSettings.method) <- mnames

				figSettings <- list(
					"Method" = figSettings.method,
					# "Region type" = figSettings.region,
					"Annotation" = figSettings.annot
				)
				rr <- muReportR::addReportFigure(rr, "Dimension reduction", plotL, figSettings)
			logger.completed()
		} else {
			logger.start("Dimension reduction and heatmap generation")
				doLogNorm <- getConfigElement("exploratoryLogNormCounts")

				txt <- c(
					"Read counts are summarized for various region types and the corresponding ",
					"aggregate count matrices are used for dimension reduction."
				)
				if (doLogNorm) {
					txt <- c(txt,
						" Counts have been log-normalized (log2(count+1))."
					)
				}
				rr <- muReportR::addReportSection(rr, "Dimension reduction", txt, level=1L, collapsed=FALSE)

				mnames <- c("pca", "umap")

				linkMethod <- "ward.D"
				corMethod <- "pearson"
				varRankCut <- 1000L
				colors.hm <- getConfigElement("colorSchemesCont")
				if (is.element("heatmap", names(colors.hm))) {
					colors.hm <- colors.hm[["heatmap"]]
				} else {
					colors.hm <- colors.hm[[".default"]]
				}
				sannot.sub <- sannot[,setdiff(plotAnnotCols, ".ALL"), drop=FALSE]

				plotL <- list()
				hmPlotL <- list()
				for (rt in regionTypes){
					logger.start(c("Region type:", rt))
						rtString <- normalize.str(rt, return.camel=TRUE)
						cm <- getCounts(.object, rt, asMatrix=TRUE)
						if (doLogNorm) cm <- log2(cm + 1)
						tcm <- t(cm)
						coords <- list(
							"pca"  = muRtools::getDimRedCoords.pca(tcm),
							"umap" = muRtools::getDimRedCoords.umap(tcm)
						)
						for (gn in plotAnnotCols){
							logger.status(c("Annotation:", gn))
							for (mn in mnames){
								pp <- getDimRedPlot(coords[[mn]], annot=sannot, colorCol=gn, shapeCol=FALSE, colScheme=grpColors[[gn]], addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + theme(aspect.ratio=1)
								plotFn <- paste("dimRed", mn, rtString, normalize.str(gn, return.camel=TRUE), sep="_")
								repPlot <- muReportR::createReportGgPlot(pp, plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
								repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
								plotL <- c(plotL, list(repPlot))
							}
						}

						logger.status(c("Clustered heatmap"))
						mostVarIdx <- which(rank(-matrixStats::rowVars(cm, na.rm=TRUE), na.last="keep",ties.method="min") <= varRankCut)
						cres.col <- as.hclust(getClusteringDendrogram(cm[mostVarIdx,], distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))
						cres.row <- as.hclust(getClusteringDendrogram(tcm[,mostVarIdx], distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))
						plotFn <- paste0("varRegionHeatmap_", rtString)
						repPlot <- muReportR::createReportPlot(plotFn, rr, width=10, height=10, create.pdf=TRUE, high.png=300L)
							pheatmap::pheatmap(
								cm[mostVarIdx,],
								color=grDevices::colorRampPalette(colors.hm)(100),
								border_color=NA,
								cluster_rows=cres.row, cluster_cols=cres.col,
								annotation_col=sannot.sub,
								annotation_colors=grpColors,
								fontsize_row=8, fontsize_col=3
							)
						repPlot <- muReportR::off(repPlot)
						hmPlotL <- c(hmPlotL, list(repPlot))
					logger.completed()
				}
				figSettings.method <- mnames
				names(figSettings.method) <- mnames
				figSettings.region <- regionTypes
				names(figSettings.region) <- normalize.str(regionTypes, return.camel=TRUE)
				figSettings.annot <- plotAnnotCols
				names(figSettings.annot) <- normalize.str(plotAnnotCols, return.camel=TRUE)
				figSettings <- list(
					"Method" = figSettings.method,
					"Region type" = figSettings.region,
					"Annotation" = figSettings.annot
				)
				rr <- muReportR::addReportFigure(rr, "Dimension reduction", plotL, figSettings)

				txt <- c(
					"Samples have been clustered according to the ", varRankCut, " most variable regions for each region type."
				)
				if (doLogNorm) {
					txt <- c(txt,
						" Counts have been log-normalized (log2(count+1))."
					)
				}
				rr <- muReportR::addReportSection(rr, "Clustered heatmaps", txt, level=1L, collapsed=FALSE)
				figSettings <- list(
					"Region type" = figSettings.region
				)
				legTxt <- paste("Clustered heatmap. The", varRankCut, "most variable regions are shown.")
				if (doLogNorm) legTxt <- paste(legTxt, "Counts have been log-normalized (log2(count+1)).")
				rr <- muReportR::addReportFigure(rr, legTxt, hmPlotL, figSettings)
			logger.completed()
		}

		doChromVar <- FALSE
		regionTypes.cv <- getConfigElement("chromVarRegionTypes")
		if (is.null(regionTypes.cv)) regionTypes.cv <- regionTypes[grepl("peak", regionTypes, ignore.case=TRUE)]
		doChromVar <- length(regionTypes.cv) > 0 && all(regionTypes.cv %in% regionTypes)
		if (doChromVar){
			logger.start("Computing chromVAR scores")
				cvResL <- lapply(regionTypes.cv, FUN=function(rt){
					logger.status(c("Region type:", rt))
					cvd <- getChromVarDev(.object, rt, motifs=getConfigElement("chromVarMotifs"))
					fn <- file.path(rDir.data.abs, paste0("chromVarDev_", normalize.str(rt, return.camel=TRUE), ".rds"))
					saveRDS(cvd, fn)
					return(cvd)
				})
				names(cvResL) <- regionTypes.cv
			logger.completed()

			cromVarRefTxt <- c("Schep, Wu, Buenrostro, & Greenleaf (2017). chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. <i>Nature Methods</i>, <b>14</b>(10), 975-978")
			rr <- muReportR::addReportReference(rr, cromVarRefTxt)
			txt <- c(
				"chromVAR ", muReportR::getReportReference(rr, cromVarRefTxt), " analysis. ",
				"The following motif set(s) were used for the analysis: ", paste(getConfigElement("chromVarMotifs"), collapse=", "), ". ",
				"R data files of chromVAR deviation scores have been attached to this report:"
			)
			rr <- muReportR::addReportSection(rr, "chromVAR", txt, level=1L, collapsed=FALSE)

			ll <- lapply(regionTypes.cv, FUN=function(rt){
				paste0("<b>", rt, ":</b> ", paste(c("<a href=\"", rDir.data, "/", paste0("chromVarDev_", normalize.str(rt, return.camel=TRUE), ".rds"), "\">","RDS file","</a>"),collapse=""))
			})
			rr <- muReportR::addReportList(rr, ll, type="u")

			logger.start("Plotting chromVAR results")
				cvMotifsForDimRed <- getConfigElement("chromVarMotifsForDimRed")
				colors.cv <- getConfigElement("colorSchemesCont")
				if (is.element("chromVAR", names(colors.cv))) {
					colors.cv <- colors.cv[["chromVAR"]]
				} else {
					colors.cv <- colors.cv[[".default.div"]]
				}
				linkMethod <- "ward.D"
				corMethod <- "pearson"
				rankCut <- 100L
				sannot.sub <- sannot[,setdiff(plotAnnotCols, ".ALL"), drop=FALSE]
				plotL.var <- list()
				plotL.dimRed <- list()
				plotL.hm <- list()
				for (rt in regionTypes.cv){
					logger.start(c("Region type:", rt))
						rts <- normalize.str(rt, return.camel=TRUE)
						cvd <- cvResL[[rt]]
						cvv <- chromVAR::computeVariability(cvd)

						pp <- chromVAR::plotVariability(cvv, use_plotly=FALSE)
						plotFn <- paste0("chromVarDevVar_", rts)
						repPlot <- muReportR::createReportGgPlot(pp, plotFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
						repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
						plotL.var <- c(plotL.var, list(repPlot))

						devScores <- chromVAR::deviationScores(cvd)
						cres.col <- as.hclust(getClusteringDendrogram(devScores, distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))

						mostVarIdx <- which(rank(-cvv$variability,na.last="keep",ties.method="min") <= rankCut)
						cres.row <- as.hclust(getClusteringDendrogram(t(devScores[mostVarIdx,]), distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))

						maxDev <- max(abs(devScores[mostVarIdx,]), na.rm=TRUE)

						if (isSingleCell){
							cvMotifs <- cvMotifsForDimRed
							plotFn <- paste0("chromVarDevUmap_", rts)
							plotL.dimRed <- c(plotL.dimRed, list(repPlot))
						} else {
							plotFn <- paste0("chromVarDevHeatmap_", rts)
							repPlot <- muReportR::createReportPlot(plotFn, rr, width=10, height=10, create.pdf=TRUE, high.png=300L)
								pheatmap::pheatmap(
									devScores[mostVarIdx,],
									color=grDevices::colorRampPalette(colors.cv)(100),
									breaks=seq(-maxDev, maxDev, length.out=101), 
									border_color=NA,
									cluster_rows=cres.row, cluster_cols=cres.col,
									annotation_col=sannot.sub,
									annotation_colors=grpColors,
									fontsize_row=8, fontsize_col=3
								)
							repPlot <- muReportR::off(repPlot)
							plotL.hm <- c(plotL.hm, list(repPlot))
						}
					logger.completed()
				}
				figSettings.region <- regionTypes.cv
				names(figSettings.region) <- normalize.str(regionTypes.cv, return.camel=TRUE)
				figSettings <- list(
					"Region type" = figSettings.region
				)
				rr <- muReportR::addReportFigure(rr, "chromVAR variability. TF motifs are shown ordered according to their variability across the dataset.", plotL.var, figSettings)
				if (length(plotL.hm) > 0) {
					rr <- muReportR::addReportFigure(rr, "chromVAR deviation scores. The heatmap shows the scores for 100 most variable TF motifs across the dataset.", plotL.hm, figSettings)
				}
			logger.completed()

		}

		muReportR::off(rr)
		invisible(rr)
	}
)

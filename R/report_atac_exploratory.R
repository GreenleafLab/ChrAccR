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
setMethod("createReport_exploratory",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		require(muReportR)
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- createReport(file.path(reportDir, "exploratory.html"), "Exploratory Analysis", page.title="Exploratory", init.configuration=initConfigDir, theme="stanford")
		rDir.data <- getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- getReportDir(rr, dir="data", absolute=TRUE)

		regionTypes <- getRegionTypes(.object)
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		sannot <- getSampleAnnot(.object)
		sannot[,".ALL"] <- "all"
		sampleGrps <- getGroupsFromTable(sannot, cols=getConfigElement("annotationColumns"))
		sampleGrps <- c(list(".ALL"=c(1:nrow(sannot))), sampleGrps)
		grpNames <- names(sampleGrps)

		colSchemes <- getConfigElement("colorSchemes")
		colSchemesNum <- getConfigElement("colorSchemesCont")
		grpColors <- lapply(names(sampleGrps), FUN=function(cn){
			cs <- c()
			x <- sampleGrps[[cn]]
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
		names(grpColors) <- names(sampleGrps)

		logger.start("Dimension reduction")
			txt <- c(
				"Read counts are summarized for various region types and the corresponding ",
				"aggregate count matrices are used for dimension reduction."
			)
			rr <- addReportSection(rr, "Dimension reduction", txt, level=1L, collapsed=FALSE)

			mnames <- c("pca", "umap")
			plotL <- list()
			for (rt in regionTypes){
				logger.start(c("Region type:", rt))
					tcm <- t(getCounts(.object, rt, asMatrix=TRUE))
					coords <- list(
						"pca"  = muRtools::getDimRedCoords.pca(tcm),
						"umap" = muRtools::getDimRedCoords.umap(tcm)
					)
					for (gn in grpNames){
						logger.status(c("Annotation:", gn))
						for (mn in mnames){
							pp <- getDimRedPlot(coords[[mn]], annot=sannot, colorCol=gn, shapeCol=FALSE, colScheme=grpColors[[gn]], addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + theme(aspect.ratio=1)
							plotFn <- paste("dimRed", mn, normalize.str(rt, return.camel=TRUE), normalize.str(gn, return.camel=TRUE), sep="_")
							repPlot <- createReportGgPlot(pp, plotFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
							repPlot <- off(repPlot, handle.errors=TRUE)
							plotL <- c(plotL, list(repPlot))
						}
					}
				logger.completed()
			}
			figSettings.method <- mnames
			names(figSettings.method) <- mnames
			figSettings.region <- regionTypes
			names(figSettings.region) <- normalize.str(regionTypes, return.camel=TRUE)
			figSettings.annot <- grpNames
			names(figSettings.annot) <- normalize.str(grpNames, return.camel=TRUE)
			figSettings <- list(
				"Method" = figSettings.method,
				"Region type" = figSettings.region,
				"Annotation" = figSettings.annot
			)
			rr <- addReportFigure(rr, "Dimension reduction", plotL, figSettings)
		logger.completed()

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
			rr <- addReportReference(rr, cromVarRefTxt)
			txt <- c(
				"chromVAR ", getReportReference(rr, cromVarRefTxt), " analysis. ",
				"The following motif set(s) were used for the analysis: ", paste(getConfigElement("chromVarMotifs"), collapse=", "), ". ",
				"R data files of chromVAR deviation scores have been attached to this report:"
			)
			rr <- addReportSection(rr, "chromVAR", txt, level=1L, collapsed=FALSE)

			ll <- lapply(regionTypes.cv, FUN=function(rt){
				paste0("<b>", rt, ":</b> ", paste(c("<a href=\"", rDir.data, "/", paste0("chromVarDev_", normalize.str(rt, return.camel=TRUE), ".rds"), "\">","RDS file","</a>"),collapse=""))
			})
			rr <- addReportList(rr, ll, type="u")

			logger.start("Plotting chromVAR results")
				require(pheatmap)
				colors.cv <- getConfigElement("colorSchemesCont")
				if (is.element("chromVAR", names(colors.cv))) {
					colors.cv <- colors.cv[["chromVAR"]]
				} else {
					colors.cv <- colors.cv[[".default.div"]]
				}
				linkMethod <- "ward.D"
				corMethod <- "pearson"
				rankCut <- 100L
				sannot.sub <- sannot[,setdiff(names(sampleGrps), ".ALL"), drop=FALSE]
				plotL.var <- list()
				plotL.hm <- list()
				for (rt in regionTypes.cv){
					logger.start(c("Region type:", rt))
						rts <- normalize.str(rt, return.camel=TRUE)
						cvd <- cvResL[[rt]]
						cvv <- computeVariability(cvd)

						pp <- plotVariability(cvv, use_plotly=FALSE)
						plotFn <- paste0("chromVarDevVar_", rts)
						repPlot <- createReportGgPlot(pp, plotFn, rr, width=10, height=5, create.pdf=TRUE, high.png=0L)
						repPlot <- off(repPlot, handle.errors=TRUE)
						plotL.var <- c(plotL.var, list(repPlot))

						devScores <- deviationScores(cvd)
						cres.col <- as.hclust(getClusteringDendrogram(devScores, distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))

						mostVarIdx <- which(rank(-cvv$variability,na.last="keep",ties.method="min") <= rankCut)
						cres.row <- as.hclust(getClusteringDendrogram(t(devScores[mostVarIdx,]), distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))

						plotFn <- paste0("chromVarDevHeatmap_", rts)
						repPlot <- createReportPlot(plotFn, rr, width=10, height=10, create.pdf=TRUE, high.png=300L)
							pheatmap(
								devScores[mostVarIdx,],
								color=colorRampPalette(colors.cv)(100),
								breaks=seq(-max(abs(devScores)), max(abs(devScores)), length.out=101), 
								border_color=NA,
								cluster_rows=cres.row, cluster_cols=cres.col,
								annotation_col=sannot.sub,
								annotation_colors=grpColors,
								fontsize_row=8, fontsize_col=3
							)
						repPlot <- off(repPlot)
						plotL.hm <- c(plotL.hm, list(repPlot))
					logger.completed()
				}
				figSettings.region <- regionTypes.cv
				names(figSettings.region) <- normalize.str(regionTypes.cv, return.camel=TRUE)
				figSettings <- list(
					"Region type" = figSettings.region
				)
				rr <- addReportFigure(rr, "chromVAR variability. TF motifs are shown ordered according to their variability across the dataset.", plotL.var, figSettings)
				rr <- addReportFigure(rr, "chromVAR deviation scores. The heatmap shows the scores for 100 most variable TF motifs across the dataset.", plotL.hm, figSettings)
			logger.completed()

		}

		off(rr)
		invisible(rr)
	}
)

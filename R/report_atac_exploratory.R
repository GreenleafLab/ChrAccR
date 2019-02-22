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
		rr <- createReport(file.path(reportDir, "exploratory.html"), "Accessibility Summary", page.title = "Summary", init.configuration=initConfigDir)
		rDir.data <- getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- getReportDir(rr, dir="data", absolute=TRUE)

		regionTypes <- getRegionTypes(.object)
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		sannot <- getSampleAnnot(.object)
		sannot[,"[ALL]"] <- "all"
		sampleGrps <- getGroupsFromTable(sannot, cols=getConfigElement("annotationColumns"))
		sampleGrps <- c(list("[ALL]"=c(1:nrow(sannot))), sampleGrps)
		grpNames <- names(sampleGrps)

		colSchemes <- getConfigElement("colorSchemes")
		grpColors <- lapply(names(sampleGrps), FUN=function(cn){
			cs <- c()
			x <- sampleGrps[[cn]]
			useDefault <- !is.element(cn, names(colSchemes))
			if (!useDefault) {
				useDefault <- !all(names(x) %in% names(colSchemes[[cn]]))
			}
			if (useDefault) {
				cs <- rep(colSchemes[[".default"]], length.out=length(x))
				names(cs) <- names(x)
			} else {
				cs <- colSchemes[[cn]][names(x)]
			}
			return(cs)
		})
		names(grpColors) <- names(sampleGrps)

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
					for (mn in mnames){
						pp <- getDimRedPlot(coords[[mn]], annot=sannot, colorCol=gn, shapeCol=FALSE, colScheme=grpColors, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + theme(aspect.ratio=1)
						plotFn <- paste("dimRed", mn, normalize.str(rt), normalize.str(gn), sep="_")
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
		names(figSettings.region) <- normalize.str(regionTypes)
		figSettings.annot <- grpNames
		names(figSettings.annot) <- normalize.str(grpNames)
		figSettings <- list(
			"Method" = figSettings.method,
			"Region type" = figSettings.region,
			"Annotation" = figSettings.annot

		)
		rr <- addReportFigure(rr, "Dimension reduction", plotL, figSettings)

		off(rr)
		invisible(rr)
	}
)
# report <- createReport_exploratory(dsr, file.path(config$.anaDir, "reports"))

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
#' @param itLsiObj   [for single-cell only; optional] pre-computed result of a call to \code{iterativeLSI(.object, ...)}
#' @param geneActSe  [for single-cell only; optional] pre-computed result of a call to \code{getCiceroGeneActivities(.object, ...)}
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
		reportDir,
		itLsiObj=NULL,
		geneActSe=NULL
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
		specAnnotCols0 <- specAnnotCols
		rts <- getConfigElement("regionTypes")
		if (length(rts) > 0) regionTypes <- intersect(rts, regionTypes)
		if (length(regionTypes) < 1) logger.error("Not enough region types specified")

		doGeneAct <- FALSE
		if (isSingleCell && !is.null(geneActSe)){
			if (is.element(class(geneActSe), c("SummarizedExperiment", "RangedSummarizedExperiment"))){
				if (ncol(geneActSe)==length(.object) && all(colnames(geneActSe) == getSamples(.object))){
					doGeneAct <- TRUE
				}
			}
			if (!doGeneAct) {
				logger.warning("Detected invalid gene activity object. --> skipping gene activity section")
			}
		}

		sannot <- getSampleAnnot(.object)
		sannot[,".ALL"] <- "all"

		cellIds <- getSamples(.object)
		if (isSingleCell){
			dre <- itLsiObj
			runItLsi <- is.null(itLsiObj)
			if (!runItLsi){
				isCompatible <- class(itLsiObj) == "iterativeLSIResultSc" && 
				                all(cellIds %in% names(itLsiObj$clustAss)) &&
				                all(cellIds %in% rownames(itLsiObj$pcaCoord)) &&
				                all(cellIds %in% rownames(itLsiObj$umapCoord))
				if (!isCompatible){
					logger.warning("Incompatible itLsiObj object. --> running iterative LSI")
					runItLsi <- TRUE
				}
			}
			if (runItLsi){
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
					clustRes <- getConfigElement("scIterativeLsiClusterResolution")
					logger.info(c("Using cluster resolution:", clustRes))
					umapParams <- getConfigElement("scIterativeLsiUmapParams")
					if (is.null(umapParams)) umapParams <- list(distMethod="euclidean", min_dist=0.5, n_neighbors=25)
					dre <- iterativeLSI(.object, it0regionType=itLsiRt, it0clusterResolution=clustRes, it1clusterResolution=clustRes, it2clusterResolution=clustRes, umapParams=umapParams)
				logger.completed()
			} else {
				logger.info("Using pre-computed iterative LSI result")
			}
			logger.start("Saving iterative LSI result")
				saveRDS(dre, file.path(rDir.data.abs, "dimRed_iterativeLSI_res.rds"))
				uwot::save_uwot(dre$umapRes, file.path(rDir.data.abs, "dimRed_iterativeLSI_res_uwot"))
			logger.completed()
			itLsiPeakRt <- ".peaks.itlsi0"
			regionTypes <- c(regionTypes, itLsiPeakRt)
			doAggr <- TRUE
			if (is.element(itLsiPeakRt, getRegionTypes(.object))){
				doAggr <- length(dre$clusterPeaks_unfiltered) != getNRegions(.object, itLsiPeakRt)
			}
			if (doAggr){
				if (!runItLsi) {
					logger.warning("Cluster peaks annotated in the DsATACsc object seem to be incompatible with specified iterative LSI result. --> reaggregating using itLsi result.")
				}
				logger.start("Aggregating counts across initial cluster peaks")
					.object <- regionAggregation(.object, dre$clusterPeaks_unfiltered, itLsiPeakRt, signal="insertions", dropEmpty=FALSE, bySample=FALSE)
				logger.completed()
				# logger.start("Aggregating counts across regions selected for dimension reduction")
				# 	.object <- regionAggregation(.object, dre$regionGr, "dimRedRegs", signal="insertions", dropEmpty=FALSE, bySample=FALSE)
				# logger.completed()
			}
			
			# annotation
			qcDf <- getScQcStatsTab(.object)
			# cellIds <- qcDf[,"cell"]
			sannot[,"clusterAssignment"] <- dre$clustAss[cellIds]
			sannot <- cbind(sannot, qcDf)
			pcns <- paste0("PC", 1:3)
			for (pcn in pcns){
				sannot[,pcn] <- dre$pcaCoord[cellIds, pcn]
			}
			qcAnnotCols <- c(setdiff(colnames(qcDf), c("cell", "sample")), pcns)
			specAnnotCols <- c(specAnnotCols, qcAnnotCols, "clusterAssignment")
		}

		mgc <- nrow(sannot)
		if (is.null(specAnnotCols0)) mgc <- nrow(sannot)-1 # exclude all-unique columns if no columns are explicitely specified
		mgc <- min(getConfigElement("annotationMaxGroupCount"), mgc)
		if (is.null(specAnnotCols0)) {
			# add automatically found columns to group set
			defaultGrps <- getGroupsFromTable(sannot, cols=NULL, minGrpSize=getConfigElement("annotationMinGroupSize"), maxGrpCount=mgc)
			specAnnotCols <- union(names(defaultGrps), specAnnotCols)
		} else if (!all(specAnnotCols %in% colnames(sannot))){
			logger.warning(c("The following annotation columns could not be found in the sample annotation/QC tables and will be discarded:", paste(setdiff(specAnnotCols, colnames(sannot)), collapse=",")))
			specAnnotCols <- intersect(specAnnotCols, colnames(sannot))
			if (length(specAnnotCols) < 1) {
				logger.warning("All annotation columns have been dropped. Resetting the specified annotation columns to default.")
			}
		}
		sampleGrps <- getGroupsFromTable(sannot, cols=specAnnotCols, minGrpSize=getConfigElement("annotationMinGroupSize"), maxGrpCount=mgc)
		sampleGrps <- c(list(".ALL"=list("all"=c(1:nrow(sannot)))), sampleGrps)
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

		if (isSingleCell){
			logger.start("Dimension reduction plots")
				cusRefTxt <- 'Cusanovich, et al. (2018). A Single-Cell Atlas of In Vivo Mammalian Chromatin Accessibility. <i>Cell</i>  <b>174</b>(5), 1309-1324, <a href="https://dx.doi.org/10.1016/j.cell.2018.06.052">doi:10.1016/j.cell.2018.06.052</a>'
				rr <- muReportR::addReportReference(rr, cusRefTxt)
				granjaRefTxt <- 'Granja, Klemm, McGinnis, et al. (2018). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. <i>Nature Biotechnology</i>  <b>37</b>(12), 1309-1324, <a href="https://dx.doi.org/10.1038/s41587-019-0332-7">doi:10.1038/s41587-019-0332-7</a>'
				rr <- muReportR::addReportReference(rr, granjaRefTxt)

				descItLsi <- c("In order to obtain a low dimensional representation of single-cell ATAC datasets in terms of principal components and UMAP coordinates, an iterative application of the Latent Semantic Indexing approach ", muReportR::getReportReference(rr, cusRefTxt), " described in ", muReportR::getReportReference(rr, granjaRefTxt), " was used. This approach also identifies cell clusters and a peak set that represents a consensus peak set of cluster peaks in a given dataset. In brief, in an initial iteration clusters are identified based on the most accessible regions (e.g. genomic tiling regions). Here, the counts are first normalized using the term frequency - inverse document frequency (TF-IDF) transformation and singular values are computed based on these normalized counts in selected regions (i.e. the most accessible regions in the initial iteration). Clusters are identified based on the singular values using Louvain clustering (as implemented in the Seurat package). Peak calling is then performed on the aggregated insertion sites from all cells of each cluster (using MACS2) and a union/consensus set of peaks uniform-length non-overlapping peaks is selected. In a second iteration, the peak regions whose TF-IDF-normalized counts which exhibit the most variability across the initial clusters provide the basis for a refined clustering using derived singular values. In the final iteration, the most variable peaks across the refined clusters are identified as the final peak set and singular values are computed again. Based on these final singular values UMAP coordinates are computed for low-dimensional projection.")
				rr <- muReportR::addReportSection(rr, "Dimension reduction", descItLsi, level=1L, collapsed=FALSE)

				it0Info <- dre$iterationData$iteration0
				pc1rem <- !is.element("pcs", names(it0Info)) || !is.element(1, it0Info[["pcs"]]) # the first condition is for backwards compatibility only, assuming the default parameters for iterativeLSI have been used (which exclude the first component)
				if (pc1rem){
					txt <- c("The first singular value has been removed in the initial iteration (due to high correlation with read depth)")
					rr <- muReportR::addReportParagraph(rr, txt)
				}

				mnames <- "umap"
				plotL <- list()
				for (gn in plotAnnotCols){
					oor <- sample.int(length(cellIds)) # random order of points
					x <- dre$umapCoord[cellIds,][oor,]
					aa <- sannot[oor, ]
					pp <- getDimRedPlot(x, annot=aa, colorCol=gn, shapeCol=FALSE, colScheme=grpColors[[gn]], ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + coord_fixed()
					isNum <- is.numeric(aa[,gn])
					if (!isNum) pp <- pp + guides(colour=guide_legend(override.aes=list(size=5)))
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
				doSubsample <- rep(FALSE, length(regionTypes))
				names(doSubsample) <- regionTypes
				nSub <- getConfigElement("exploratoryNSubsample")
				if (!is.null(nSub) && nSub > 0 && nSub < Inf){
					for (rt in regionTypes){
						doSubsample[rt] <- nSub < getNRegions(.object, rt)
					}
				}
				txt <- c(
					"Read counts are summarized for various region types and the corresponding ",
					"aggregate count matrices are used for dimension reduction."
				)
				if (doLogNorm) {
					txt <- c(txt,
						" Counts have been log-normalized (log10(count+1))."
					)
				}
				if (any(doSubsample)){
					txt <- c(txt, " The following region types have been subsampled to ", nSub, " features: ",
						paste(names(doSubsample)[doSubsample], collapse=", ")
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
						if (doLogNorm) cm <- log10(cm + 1)
						tcm <- t(cm)
						coords <- list()
						subIdx <- NULL
						if (doSubsample[rt]){
							logger.info(c("Subsampling to", nSub, "(of", ncol(tcm), ") regions"))
							subIdx <- sort(sample.int(ncol(tcm), nSub))
							coords <- list(
								"pca"  = muRtools::getDimRedCoords.pca(tcm[,subIdx]),
								"umap" = muRtools::getDimRedCoords.umap(tcm[,subIdx])
							)
						} else {
							coords <- list(
								"pca"  = muRtools::getDimRedCoords.pca(tcm),
								"umap" = muRtools::getDimRedCoords.umap(tcm)
							)
						}
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
				cvMotifsForDimRed <- getConfigElement("chromVarMotifNamesForDimRed")
				getMostVarMotifsForPlot <- isSingleCell && is.null(cvMotifsForDimRed)
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
				if (!is.null(cvMotifsForDimRed)){
					# normalize motif names to get rid of funky artifacts in ggplot
					cvMotifsForDimRed <- gsub("::","_",cvMotifsForDimRed)
					cvMotifsForDimRed <- gsub("[^[:alnum:]_\\.]","",cvMotifsForDimRed)
				}
				plotL.var <- list()
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
						rownames(devScores) <- gsub("::","_",rownames(devScores))
						rownames(devScores) <- gsub("[^[:alnum:]_\\.]","",rownames(devScores))
						mostVarIdx <- which(rank(-cvv$variability,na.last="keep",ties.method="min") <= rankCut)
						maxDev <- max(abs(devScores[mostVarIdx,]), na.rm=TRUE)

						if (!isSingleCell){
							cres.col <- as.hclust(getClusteringDendrogram(devScores, distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))
							cres.row <- as.hclust(getClusteringDendrogram(t(devScores[mostVarIdx,]), distMethod="euclidean", linkMethod=linkMethod, corMethod=corMethod))
							
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

						if (getMostVarMotifsForPlot){
							selIdx <- which(rank(-cvv$variability,na.last="keep",ties.method="min") <= 10)
							cvMotifsForDimRed <- union(cvMotifsForDimRed, rownames(devScores[selIdx,]))
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


				if (isSingleCell){
					logger.start("chromVAR deviations projected on dimension reduction plots")
						figSettings.motif <- cvMotifsForDimRed
						names(figSettings.motif) <- gsub("_", "", normalize.str(cvMotifsForDimRed, return.camel=TRUE))

						plotL.dimRed <- list()
						for (rt in regionTypes.cv){
							logger.start(c("Region type:", rt))
								rts <- normalize.str(rt, return.camel=TRUE)
								cvM <- t(devScores[cvMotifsForDimRed, cellIds])
								maxDev <- quantile(abs(cvM), 0.99, na.rm=TRUE)
								pL <- lapply(colnames(cvM), FUN=function(mm){
									pp <- getDimRedPlot(dre$umapCoord[cellIds,], annot=cvM, colorCol=mm, shapeCol=FALSE, colScheme=NULL, ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + scale_color_gradientn(colours=colors.cv, limits=c(-maxDev, maxDev), na.value = "#C0C0C0") + coord_fixed()
									figFn <- paste0("chromVarDevUmap_", rts, "_", gsub("_", "", normalize.str(mm, return.camel=TRUE)))
									repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
									repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
									return(repPlot)
								})
								plotL.dimRed <- c(plotL.dimRed, pL)
							logger.completed()
						}
						figSettings <- list(
							"Region type" = figSettings.region,
							"Motif" = figSettings.motif
						)
						rr <- muReportR::addReportFigure(rr, "Dimension reduction annotated with chromVAR deviation scores", plotL.dimRed, figSettings)
					logger.completed()
				}
			logger.completed()

		}

		if (doGeneAct){
			logger.start("Gene activity")
				gaM <- SummarizedExperiment::assay(geneActSe)
				geneNameUniv <- names(geneActSe)
				if (!is.null(SummarizedExperiment::rowRanges(geneActSe))){
					geneAnnot <- elementMetadata(SummarizedExperiment::rowRanges(geneActSe))
					cn <- match("gene_name", colnames(geneAnnot))
					if (!is.na(cn)){
						geneNameUniv <- geneAnnot[,cn]
					}
				}
				# gaM <- t(smoothMagic(t(SummarizedExperiment::assay(geneAct)[,cellIds]), X_knn=dro$pcaCoord[cellIds,dro$pcs], k=15, ka=4, td=3)$Xs)
				gaM <- log2(gaM*1e6+1) # normalize to better dynamic range
				idx <- sample(cellIds) # random cell order for plotting

				geneNames <- getConfigElement("genesOfInterest")
				if (length(geneNames) > 0){
					midx <- match(towlower(geneNames), tolower(geneNameUniv))
					if (any(is.na(midx))){
						logger.warning(c("The following gene names could not be found and will be omitted:", paste(geneNames[is.na(midx)], collapse=",")))
						midx <- na.omit(midx)
					}
					geneNames <- geneNameUniv[midx]
				}

				if (length(geneNames) < 1) {
					logger.info("Picking the 10 most variable genes for gene activity reporting")
					# pick the most variable genes
					vv <- matrixStats::rowVars(as.matrix(gaM), na.rm=TRUE)
					selIdx <- which(rank(-vv, na.last="keep",ties.method="min") <= 10)
					geneNames <- geneNameUniv[selIdx]
				}

				gaAnnot <- t(as.matrix(gaM[geneNames, idx]))
				# truncate quantiles
				gaAnnot <- apply(gaAnnot, 2, FUN=function(x){
					qq <- quantile(x, probs=0.99)
					x[x>qq] <- qq
					return(x)
				})

				md <- S4Vectors::metadata(geneActSe)
				methodTxt <- ""
				if (tolower(md$method) == "rbf"){
					methodTxt <- paste0(
						" Fragment counts in peaks within ", md$params[["maxDist"]], "bp to a TSS have been summed up using RBF-based weighting",
						"correlation cutoff: ", md$params[["corCutOff"]], ")"
					)
				} else if (tolower(md$method) == "cicero"){
					methodTxt <- paste0(
						" Peaks within ", md$params[["maxDist"]], "bp to a TSS have been associated to that TSS using Cicero's correlation-based linking (",
						"sigma: ", md$params[["sigma"]], "; baseline weight: ", md$params[["minWeight"]], ")"
					)
				}
				txt <- c(
					"Gene activities have been computed as the aggregated accessibility of TSS-associated peaks.", methodTxt,
					" The resulting scores for single-cells have been rescaled to one million counts and have been log-normalized."
				)
				rr <- muReportR::addReportSection(rr, "Gene activity", txt, level=1L, collapsed=FALSE)

				umapC <- dre$umapCoord[cellIds,]
				
				figSettings.gene <- geneNames
				names(figSettings.gene) <- gsub("_", "", normalize.str(geneNames, return.camel=TRUE))

				pL <- lapply(seq_along(geneNames), FUN=function(i){
					gn <- colnames(gaAnnot)[i]
					pp <- getDimRedPlot(umapC, annot=gaAnnot, colorCol=gn, shapeCol=FALSE, colScheme=c("#e0f3db", "#a8ddb5", "#4eb3d3", "#08589e"), ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + coord_fixed()
					figFn <- paste0("geneActUmap_", names(figSettings.gene)[i])
					repPlot <- muReportR::createReportGgPlot(pp, figFn, rr, width=7, height=7, create.pdf=TRUE, high.png=0L)
					repPlot <- muReportR::off(repPlot, handle.errors=TRUE)
				})

				figSettings <- list(
					"Gene" = figSettings.gene
				)
				rr <- muReportR::addReportFigure(rr, "Dimension reduction annotated with gene activity scores", pL, figSettings)
			logger.completed()
		}

		muReportR::off(rr)
		invisible(rr)
	}
)

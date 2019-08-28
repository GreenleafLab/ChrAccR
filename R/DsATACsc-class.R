#' DsATACsc
#'
#' A class for storing single-cell ATAC-seq accessibility data
#' inherits from \code{\linkS4class{DsATAC}}. Provides a few additional methods
#' but is otherwise identical to \code{\linkS4class{DsATAC}}.
#' 
#' @name DsATACsc-class
#' @rdname DsATACsc-class
#' @author Fabian Mueller
#' @exportClass DsATACsc
setClass("DsATACsc",
	contains = "DsATAC",
	package = "ChrAccR"
)
setMethod("initialize","DsATACsc",
	function(
		.Object,
		fragments,
		coord,
		counts,
		sampleAnnot,
		genome,
		diskDump,
		diskDump.fragments,
		diskDump.fragments.nSamplesPerFile,
		sparseCounts
	) {
		.Object@fragments  <- fragments
		.Object@coord       <- coord
		.Object@counts      <- counts
		.Object@countTransform <- rep(list(character(0)), length(.Object@counts))
		names(.Object@countTransform) <- names(.Object@counts)
		.Object@sampleAnnot <- sampleAnnot
		.Object@genome      <- genome
		.Object@diskDump    <- diskDump
		.Object@diskDump.fragments <- diskDump.fragments
		.Object@diskDump.fragments.nSamplesPerFile <- diskDump.fragments.nSamplesPerFile
		.Object@sparseCounts <- sparseCounts
		.Object@pkgVersion  <- packageVersion("ChrAccR")
		.Object
	}
)
#' @noRd
DsATACsc <- function(sampleAnnot, genome, diskDump=FALSE, diskDump.fragments=TRUE, sparseCounts=TRUE){
	obj <- new("DsATACsc",
		list(),
		list(),
		list(),
		sampleAnnot,
		genome,
		diskDump,
		diskDump.fragments,
		diskDump.fragments.nSamplesPerFile=500L,
		sparseCounts
	)
	return(obj)
}

################################################################################
# Single-cell methods
################################################################################
#-------------------------------------------------------------------------------

if (!isGeneric("getScQcStatsTab")) {
	setGeneric(
		"getScQcStatsTab",
		function(.object, ...) standardGeneric("getScQcStatsTab"),
		signature=c(".object")
	)
}
#' getScQcStatsTab-methods
#'
#' Retrieve a table of QC statistics for single cells
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @return an \code{data.frame} contain QC statistics for each cell
#' 
#' @rdname getScQcStatsTab-DsATACsc-method
#' @docType methods
#' @aliases getScQcStatsTab
#' @aliases getScQcStatsTab,DsATACsc-method
#' @author Fabian Mueller
#' @export
setMethod("getScQcStatsTab",
	signature(
		.object="DsATACsc"
	),
	function(
		.object
	) {
		cellAnnot <- getSampleAnnot(.object)
		sampleIdCn <- findOrderedNames(colnames(cellAnnot), c(".sampleid", "sampleid", ".CR.cellQC.barcode"), ignore.case=TRUE)
		nFragCns <- c(
			total=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.total"),
			pass=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.passed_filters"),
			tss=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.TSS_fragments"),
			peak=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.peak_region_fragments"),
			duplicate=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.duplicate"),
			mito=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.mitochondrial")
		)
		summaryDf <- data.frame(cell=getSamples(.object), sample=rep("sample", nrow(cellAnnot)))
		if (!is.na(sampleIdCn)) summaryDf[,"sample"] <- cellAnnot[,sampleIdCn]
		for (cn in c("total", "pass")){
			summaryDf[,muRtools::normalize.str(paste("n", cn, sep="_"), return.camel=TRUE)] <- cellAnnot[,nFragCns[cn]]
		}
		# to be divided by total reads
		for (cn in c("mito", "duplicate")){
			if (!is.na(nFragCns[cn])) summaryDf[,muRtools::normalize.str(paste("frac", cn, sep="_"), return.camel=TRUE)] <- cellAnnot[,nFragCns[cn]]/summaryDf[,"nTotal"]
		}
		# to be divided by passing reads
		for (cn in c("tss", "peak")){
			if (!is.na(nFragCns[cn])) summaryDf[,muRtools::normalize.str(paste("frac", cn, sep="_"), return.camel=TRUE)] <- cellAnnot[,nFragCns[cn]]/summaryDf[,"nPass"]
		}
		return(summaryDf)
	}
)

if (!isGeneric("unsupervisedAnalysisSc")) {
	setGeneric(
		"unsupervisedAnalysisSc",
		function(.object, ...) standardGeneric("unsupervisedAnalysisSc"),
		signature=c(".object")
	)
}
#' unsupervisedAnalysisSc-methods
#'
#' Perform unsupervised analysis on single-cell data. Performs dimensionality reduction
#' and clustering.
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param regionType character string specifying the region type
#' @param regionIdx  indices of regions to be used (logical or integer vector). If \code{NULL} (default) all regions of the specified regionType will be used.
#' @param dimRedMethod character string specifying the dimensionality reduction method. Currently on \code{"tf-idf_irlba"} is supported
#' @param usePcs     integer vector specifying the principal components to use for UMAP and clustering
#' @param clusteringMethod character string specifying the clustering method. Currently on \code{"seurat_louvain"} is supported
#' @return an \code{S3} object containing dimensionality reduction results and clustering
#' 
#' @rdname unsupervisedAnalysisSc-DsATACsc-method
#' @docType methods
#' @aliases unsupervisedAnalysisSc
#' @aliases unsupervisedAnalysisSc,DsATACsc-method
#' @author Fabian Mueller
#' @export
setMethod("unsupervisedAnalysisSc",
	signature(
		.object="DsATACsc"
	),
	function(
		.object,
		regionType,
		regionIdx=NULL,
		dimRedMethod="tf-idf_irlba",
		usePcs=1:50,
		clusteringMethod="seurat_louvain"
	) {
		if (!is.element(regionType, getRegionTypes(.object))) logger.error(c("Unsupported region type:", regionType))
		if (!is.element(dimRedMethod, c("tf-idf_irlba"))) logger.error(c("Unsupported dimRedMethod:", dimRedMethod))
		if (!is.integer(usePcs)) logger.error(c("usePcs must be an integer vector"))
		if (!is.element(clusteringMethod, c("seurat_louvain"))) logger.error(c("Unsupported clusteringMethod:", clusteringMethod))
		if (!is.null(regionIdx)){
			if (is.logical(regionIdx)) regionIdx <- which(regionIdx)
			if (!is.integer(regionIdx) || any(regionIdx < 1) || any(regionIdx > getNRegions(.object, regionType))) logger.error("Invalid regionIdx")
		}

		dsn <- .object
		if (!is.null(regionIdx)){
			nRegs <- getNRegions(.object, regionType)
			logger.info(c("Retaining", length(regionIdx), "regions for dimensionality reduction"))
			idx <- rep(TRUE,  nRegs)
			idx[regionIdx] <- FALSE
			dsn <- removeRegions(.object, idx, regionType)
		}
		if (dimRedMethod=="tf-idf_irlba"){
			logger.start(c("Performing dimensionality reduction using", dimRedMethod))			
				if (length(dsn@countTransform[[regionType]]) > 0) logger.warning("Counts have been pre-normalized. dimRedMethod 'tf-idf_irlba' might not be applicable.")
				
				if (!is.element("tf-idf", dsn@countTransform[[regionType]])){
					dsn <- transformCounts(dsn, method="tf-idf", regionTypes=regionType)
				}

				cm <- ChrAccR::getCounts(dsn, regionType, asMatrix=TRUE)
				pcaCoord <- muRtools::getDimRedCoords.pca(t(cm), components=1:max(usePcs), method="irlba_svd")
			logger.completed()
		}
		cellIds <- colnames(cm)
		logger.start(c("Getting UMAP coordinates"))	
			umapCoord <- muRtools::getDimRedCoords.umap(pcaCoord[,usePcs])
			umapRes <- attr(umapCoord, "umapRes")
			attr(umapCoord, "umapRes") <- NULL
		logger.completed()

		if (clusteringMethod=="seurat_louvain"){
			logger.start(c("Performing clustering using", clusteringMethod))	
				if (!requireNamespace("Seurat")) logger.error(c("Could not load dependency: Seurat"))
				# Louvain clustering using Seurat
				dummyMat <- matrix(11.0, ncol=length(cellIds), nrow=11)
				colnames(dummyMat) <- cellIds
				rownames(dummyMat) <- paste0("df", 1:nrow(dummyMat))
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcaCoord, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=usePcs, k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10)
				clustAss <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss) <- names(clustRes@active.ident)
			logger.completed()
		}
		res <- list(
			pcaCoord=pcaCoord,
			umapCoord=umapCoord,
			umapRes=umapRes,
			clustAss=clustAss,
			regionType=regionType,
			regionIdx=regionIdx
		)
		class(res) <- "unsupervisedAnalysisResultSc"
		return(res)
	}
)
#-------------------------------------------------------------------------------

if (!isGeneric("iterativeLSI")) {
	setGeneric(
		"iterativeLSI",
		function(.object, ...) standardGeneric("iterativeLSI"),
		signature=c(".object")
	)
}
#' iterativeLSI-methods
#'
#' EXPERIMENTAL: Perform iterative LSI clustering as described in doi:10.1101/696328
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param it0regionType character string specifying the region type to start with
#' @param it0nMostAcc the number of the most accessible regions to consider in iteration 0
#' @param it0pcs      the principal components to consider in iteration 0
#' @param it0clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in iteration 0
#' @param it0nTopPeaksPerCluster the number of best peaks to be considered for each cluster in the merged peak set (iteration 0)
#' @param it1pcs      the principal components to consider in iteration 0
#' @param it1clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in iteration 1
#' @param it1mostVarPeaks the number of the most variable peaks to consider after iteration 1
#' @param it2pcs      the principal components to consider in the final iteration (2)
#' @param it2clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in the final iteration (2)
#' @return an \code{S3} object containing dimensionality reduction results and clustering
#' 
#' @rdname iterativeLSI-DsATACsc-method
#' @docType methods
#' @aliases iterativeLSI
#' @aliases iterativeLSI,DsATACsc-method
#' @author Fabian Mueller
#' @export
setMethod("iterativeLSI",
	signature(
		.object="DsATACsc"
	),
	function(
		.object,
		it0regionType="t5k",
		it0nMostAcc=20000L,
		it0pcs=2:25,
		it0clusterResolution=0.8,
		it0nTopPeaksPerCluster=2e5,
		it1pcs=1:50,
		it1clusterResolution=0.8,
		it1mostVarPeaks=50000L,
		it2pcs=1:50,
		it2clusterResolution=0.8
	) {
		callParams <- as.list(match.call())
		callParams <- callParams[setdiff(names(callParams), ".object")]
		cellIds <- getSamples(.object)
		if (length(.object@fragments) != length(cellIds)) logger.error("Object does not contain fragment information for all samples")

		ph <- getSampleAnnot(.object)
		depthCol <- colnames(ph) %in% c("numIns", ".CR.cellQC.passed_filters", ".CR.cellQC.total")
		depthV <- NULL
		if (any(depthCol)){
			depthV <- ph[,colnames(ph)[depthCol][1]]
		}

		logger.start("Iteration 0")
			dsr <- .object
			for (rt in setdiff(getRegionTypes(dsr), it0regionType)){
				dsr <- removeRegionType(dsr, rt)
			}
			if (!is.null(it0nMostAcc)){
				regAcc <- safeMatrixStats(ChrAccR::getCounts(dsr, it0regionType, allowSparseMatrix=TRUE), statFun="rowMeans", na.rm=TRUE)
				if (it0nMostAcc < length(regAcc)){
					idx2rem <- rank(-regAcc, na.last="keep", ties.method="min") > it0nMostAcc
					logger.info(c("Retaining the", sum(!idx2rem), "most accessible regions for dimensionality reduction"))
					dsr <- removeRegions(dsr, idx2rem, it0regionType)
				}
			}
			logger.start(c("Performing TF-IDF-based dimension reduction"))
				if (length(dsr@countTransform[[it0regionType]]) > 0) logger.warning("Counts have been pre-normalized. 'tf-idf' might not be applicable.")
				dsn <- transformCounts(dsr, method="tf-idf", regionTypes=it0regionType)

				cm <- ChrAccR::getCounts(dsn, it0regionType, allowSparseMatrix=TRUE)
				pcaCoord_it0 <- muRtools::getDimRedCoords.pca(safeMatrixStats(cm, "t"), components=1:max(it0pcs), method="irlba_svd")
				if (!is.null(depthV)){
					cc <- cor(pcaCoord_it0[,1], depthV, method="spearman")
					logger.info(c("Correlation (Spearman) of PC1 with cell fragment counts:", round(cc, 4)))
				}
				pcaCoord_it0 <- pcaCoord_it0[, it0pcs, drop=FALSE]
			logger.completed()
			logger.start(c("Clustering"))	
				if (!requireNamespace("Seurat")) logger.error(c("Could not load dependency: Seurat"))
				# Louvain clustering using Seurat
				dummyMat <- matrix(11.0, ncol=length(cellIds), nrow=11)
				colnames(dummyMat) <- cellIds
				rownames(dummyMat) <- paste0("df", 1:nrow(dummyMat))
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcaCoord_it0, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcaCoord_it0), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it0clusterResolution)
				clustAss_it0 <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss_it0) <- names(clustRes@active.ident)
			logger.completed()
			logger.start(c("Peak calling"))
				logger.start("Creating cluster pseudo-bulk samples")
					dsr <- addSampleAnnotCol(dsr, "clustAss_it0", as.character(clustAss_it0[cellIds]))
					dsrClust <- mergeSamples(dsr, "clustAss_it0", countAggrFun="sum")
				logger.completed()
				logger.start("Calling peaks")
					clustPeakGrl <- callPeaks(dsrClust)
					if (!is.null(it0nTopPeaksPerCluster)){
						logger.info(paste0("Selecting the", it0nTopPeaksPerCluster, " peaks with highest score for each cluster"))
						clustPeakGrl <- GRangesList(lapply(clustPeakGrl, FUN=function(x){
							idx <- rank(-elementMetadata(x)[,"score_norm"], na.last="keep", ties.method="min") <= it0nTopPeaksPerCluster
							x[idx]
						}))
					}
					
					peakUnionGr <- getNonOverlappingByScore(unlist(clustPeakGrl), scoreCol="score_norm")
					peakUnionGr <- sortGr(peakUnionGr)
				logger.completed()
				logger.start("Aggregating counts for union peak set")
					# dsrClust <- regionAggregation(dsrClust, peakUnionGr, "clusterPeaks", signal="insertions", dropEmpty=FALSE)
					dsr <- regionAggregation(dsr, peakUnionGr, "clusterPeaks", signal="insertions", dropEmpty=FALSE, bySample=FALSE)
				logger.completed()
			logger.completed()
		logger.completed()

		logger.start("Iteration 1")
			it1regionType <- "clusterPeaks"
			logger.start(c("Performing TF-IDF-based dimension reduction"))
				dsr <- removeRegionType(dsr, it0regionType)
				dsn <- transformCounts(dsr, method="tf-idf", regionTypes=it1regionType) #TODO: renormalize based on sequencing depth rather than aggregated counts across peaks only?
				cm <- ChrAccR::getCounts(dsn, it1regionType, allowSparseMatrix=TRUE)
				pcaCoord_it1 <- muRtools::getDimRedCoords.pca(safeMatrixStats(cm, "t"), components=1:max(it1pcs), method="irlba_svd")[, it1pcs, drop=FALSE]
			logger.completed()

			logger.start(c("Clustering"))	
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcaCoord_it1, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcaCoord_it1), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it1clusterResolution)
				clustAss_it1 <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss_it1) <- names(clustRes@active.ident)
			logger.completed()

			if (!is.null(it1mostVarPeaks) && it1mostVarPeaks < nrow(cm)){
				logger.start(c("Identifying cluster-variable peaks"))
					logger.start("Creating cluster pseudo-bulk samples")
						dsr <- addSampleAnnotCol(dsr, "clustAss_it1", as.character(clustAss_it1[cellIds]))
						dsrClust <- mergeSamples(dsr, "clustAss_it1", countAggrFun="sum")
					logger.completed()
					logger.start("Identifying target peaks")
						dsnClust <- transformCounts(dsrClust, method="RPKM", regionTypes=it1regionType)
						l2cpm <- log2(ChrAccR::getCounts(dsnClust, it1regionType) / 1e3 + 1) # compute log2(CPM) from RPKM
						peakVar <- matrixStats::rowVars(l2cpm, na.rm=TRUE)
						if (it1mostVarPeaks < length(peakVar)){
							idx2rem <- rank(-peakVar, na.last="keep", ties.method="min") > it1mostVarPeaks
							logger.info(c("Retaining the", sum(!idx2rem), "most variable peaks"))
							dsr <- removeRegions(dsr, idx2rem, it1regionType)
						}
						peakCoords <- ChrAccR::getCoord(dsr, it1regionType)
					logger.completed()
				logger.completed()
			}
		logger.completed()

		logger.start("Iteration 2")
			it2regionType <- it1regionType
			logger.start(c("Performing TF-IDF-based dimension reduction"))
				bcm_unnorm <- ChrAccR::getCounts(dsr, it2regionType, allowSparseMatrix=TRUE) > 0 # unnormalized binary count matrix
				idfBase <- log(1 + ncol(bcm_unnorm) / safeMatrixStats(bcm_unnorm, "rowSums", na.rm=TRUE))
				dsn <- transformCounts(dsr, method="tf-idf", regionTypes=it2regionType) #TODO: renormalize based on sequencing depth rather than aggregated counts across peaks only?
				cm <- ChrAccR::getCounts(dsn, it2regionType, allowSparseMatrix=TRUE)
				pcaCoord <- muRtools::getDimRedCoords.pca(safeMatrixStats(cm, "t"), components=1:max(it2pcs), method="irlba_svd")
				pcaCoord_sel <- pcaCoord[, it2pcs, drop=FALSE]
			logger.completed()
			logger.start(c("Clustering"))
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcaCoord_sel, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcaCoord_sel), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it2clusterResolution)
				clustAss <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss) <- names(clustRes@active.ident)

				dsr <- addSampleAnnotCol(dsr, "clustAss_it2", as.character(clustAss[cellIds]))
			logger.completed()
			logger.start(c("UMAP coordinates"))	
				umapCoord <- muRtools::getDimRedCoords.umap(pcaCoord_sel)
				umapRes <- attr(umapCoord, "umapRes")
				attr(umapCoord, "umapRes") <- NULL
			logger.completed()
		logger.completed()

		res <- list(
			pcaCoord=pcaCoord,
			pcs = it2pcs,
			idfBase=idfBase,
			umapCoord=umapCoord,
			umapRes=umapRes,
			clustAss=clustAss,
			regionGr=peakCoords,
			clusterPeaks_unfiltered=peakUnionGr,
			iterationData = list(
				iteration0 = list(
					pcaCoord=pcaCoord_it0,
					clustAss=clustAss_it0
				),
				iteration1 = list(
					pcaCoord=pcaCoord_it1,
					clustAss=clustAss_it1
				)
			),
			.params=callParams
		)
		class(res) <- "iterativeLSIResultSc"
		return(res)
	}
)


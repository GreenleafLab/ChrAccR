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
if (!isGeneric("simulateDoublets")) {
	setGeneric(
		"simulateDoublets",
		function(.object, ...) standardGeneric("simulateDoublets"),
		signature=c(".object")
	)
}
#' simulateDoublets-methods
#'
#' EXPERIMENTAL: Simulate doublets by adding counts in matrices for each region set
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param byGroup    sample by group. Must be a column name in the sample annotation table
#' @param n          number of doublets to simulate (per group)
#' @param sampleRatio fraction of non-zero events to subsample. If \code{0 < sampleRatio < 1}, the individual cell counts will be subsampled
#'                   for each region.
#' @return an \code{\linkS4class{DsATACsc}} object containing counts for simulated doublets. Fragment data will be discarded.
#' 
#' @rdname simulateDoublets-DsATACsc-method
#' @docType methods
#' @aliases simulateDoublets
#' @aliases simulateDoublets,DsATACsc-method
#' @author Fabian Mueller
#' @export
#' @noRd
setMethod("simulateDoublets",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		byGroup=NULL,
		n=10000,
		sampleRatio=1.0
	) {
		.sampleSparseMat <- function(mat, sampleRatio=0.5){
			total <- length(mat@x)
			sampleTo <- floor(total * (1-sampleRatio))
			mat@x[sample(seq_len(total), sampleTo)] <- 0
			mat <- drop0(mat)
			mat
		}
		.sampleDensMat <- function(mat, sampleRatio=0.5){
			nonZero <- which(mat > 0)
			sampleTo <- floor(length(nonZero) * (1-sampleRatio))
			mat[sample(nonZero, sampleTo)] <- 0
		}
		
		if (class(.object)!="DsATACsc"){
			logger.warning("Doublet detection is intended for single-cell datasets [DsATACsc] only. Applying it to general DsATAC objects is intended for backwards compatibility only.")
		}

		res <- removeFragmentData(.object)
		nCells <- length(getSamples(.object))
		if (is.null(byGroup)){
			idx <- cbind(sample(seq_len(nCells), n, replace=TRUE), sample(seq_len(nCells), n, replace=TRUE))
		} else {
			ggs <- getSampleAnnot(.object)
			if (!is.element(byGroup, colnames(ggs))){
				logger.error("byGroup must be a valid column name in the sample annotation")
			}
			logger.info(paste0("Simulating doublets by group: '", byGroup, "'..."))
			gIdx <- getGroupsFromTable(ggs[,byGroup, drop=FALSE], minGrpSize=1)[[1]]
			idx <- do.call("rbind", lapply(gIdx, FUN=function(x){
				 cbind(sample(x, n, replace=TRUE), sample(x, n, replace=TRUE))
			}))
		}
		colnames(idx) <- c("idx1", "idx2")
		n <- nrow(idx)

		ph <- data.frame(doubletId=paste0("d", 1:n), idx)
		res@sampleAnnot <- ph
		rownames(res@sampleAnnot) <- ph[,"doubletId"]
		res@diskDump <- FALSE

		regTypes <- getRegionTypes(.object)
		for (rt in regTypes){
			logger.status(paste0("Simulating doublets for region set: '", rt, "'..."))
			cm <- ChrAccR::getCounts(.object, rt, allowSparseMatrix=TRUE)
			pkg <- attr(class(cm), "package")
			isSparse <- is.character(pkg) && pkg=="Matrix"

			cm1 <- cm[, idx[, 1], drop=FALSE]
			cm2 <- cm[, idx[, 2], drop=FALSE]
			if (sampleRatio > 0 && sampleRatio < 1){
				logger.info(c("Subsampling to ratio:", sampleRatio))
				if (isSparse){
					cm1 <- .sampleSparseMat(cm1, sampleRatio)
					cm2 <- .sampleSparseMat(cm2, sampleRatio)
				} else {
					cm1 <- .sampleDensMat(cm1, sampleRatio)
					cm2 <- .sampleDensMat(cm2, sampleRatio)
				}
			}
			cmm <- cm1 + cm2
			colnames(cmm) <- ph[,"doubletId"]

			if (res@sparseCounts & !isSparse) {
				res@counts[[rt]] <- as(cmm, "sparseMatrix")
				res@counts[[rt]] <- drop0(res@counts[[rt]])
			} else {
				res@counts[[rt]] <- cmm
			}
			if (res@diskDump) res@counts[[rt]] <- as(res@counts[[rt]], "HDF5Array")
		}

		return(res)
	}
)
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
			total=findOrderedNames(colnames(cellAnnot), c(".CR.cellQC.total", "nFrags")),
			pass=findOrderedNames(colnames(cellAnnot), c(".CR.cellQC.passed_filters", "nFrags")),
			tss=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.TSS_fragments"),
			peak=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.peak_region_fragments"),
			duplicate=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.duplicate"),
			mito=findOrderedNames(colnames(cellAnnot), ".CR.cellQC.mitochondrial")
		)
		summaryDf <- data.frame(cell=getSamples(.object), sample=rep("sample", nrow(cellAnnot)))
		if (!is.na(sampleIdCn)) summaryDf[,"sample"] <- cellAnnot[,sampleIdCn]
		if (is.na(nFragCns["total"]) && is.na(nFragCns["pass"])){
			logger.info("Number of fragments not annotated. --> trying to count from fragment data")
			hasFragments <- length(.object@fragments) > 0
			if (!hasFragments) logger.error("No fragment data found")
			cellAnnot[,".countedFragments"] <- getFragmentNum(.object)
			nFragCns["total"] <- ".countedFragments"
			nFragCns["pass"] <- ".countedFragments"
		}
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
		# tssEnrichment
		cn <- findOrderedNames(colnames(cellAnnot), c(".tssEnrichment", ".tssEnrichment_smoothed", "tssEnrichment", ".tssEnrichment_unsmoothed"))
		if (!is.na(cn)){
			summaryDf[,"tssEnrichment"] <- cellAnnot[,cn]
		}
		return(summaryDf)
	}
)

#-------------------------------------------------------------------------------
if (!isGeneric("filterCellsTssEnrichment")) {
	setGeneric(
		"filterCellsTssEnrichment",
		function(.object, ...) standardGeneric("filterCellsTssEnrichment"),
		signature=c(".object")
	)
}
#' filterCellsTssEnrichment-methods
#'
#' Filter out cells with low TSS enrichment
#'
#' @param .object  \code{\linkS4class{DsATAC}} object
#' @param cutoff   TSS enrichment cutoff to filter cells
#' @return modified \code{\linkS4class{DsATAC}} object without filtered cells
#' 
#' @rdname filterCellsTssEnrichment-DsATAC-method
#' @docType methods
#' @aliases filterCellsTssEnrichment
#' @aliases filterCellsTssEnrichment,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("filterCellsTssEnrichment",
	signature(
		.object="DsATACsc"
	),
	function(
		.object,
		cutoff=6
	) {
		cellAnnot <- getSampleAnnot(.object)
		cn <- findOrderedNames(colnames(cellAnnot), c(".tssEnrichment", ".tssEnrichment_smoothed", "tssEnrichment", ".tssEnrichment_unsmoothed"))
		if (is.na(cn)){
			logger.info("TSS enrichment not annotated. Computing TSS enrichment ...")
			tsseRes <- getTssEnrichmentBatch(.object, tssGr=NULL)
			.object <- addSampleAnnotCol(.object, ".tssEnrichment_unsmoothed", tsseRes$tssEnrichment)
			.object <- addSampleAnnotCol(.object, ".tssEnrichment", tsseRes$tssEnrichment.smoothed)
			cn <- ".tssEnrichment_unsmoothed"
			cellAnnot <- getSampleAnnot(.object)
		}
		tsse <- cellAnnot[,cn]

		.object <- .object[tsse >= cutoff]

		# workaround: currently saving single-cell DsATAC datasets does not support re-chunking of disk-dumped fragment data
		chunkedFragmentFiles <- .object@diskDump.fragments && .hasSlot(.object, "diskDump.fragments.nSamplesPerFile") && .object@diskDump.fragments.nSamplesPerFile > 1
		if (chunkedFragmentFiles){
			logger.start("Undisking ...")
				.object <- undiskFragmentData(.object)
			logger.completed()
		}

		return(.object)
	}
)

#-------------------------------------------------------------------------------

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
if (!isGeneric("dimRed_UMAP")) {
	setGeneric(
		"dimRed_UMAP",
		function(.object, ...) standardGeneric("dimRed_UMAP"),
		signature=c(".object")
	)
}
#' dimRed_UMAP-methods
#'
#' Retrieve dimension reduction embedding and object using UMAP
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param regions    character string specifying the region type to retrieve the UMAP coordinates from. Alternatively, a \code{GRanges} object specifying coordinates that fragment counts will be aggregated over
#' @param tfidf      normalize the counts using TF-IDF transformation
#' @param pcs        components to use to compute the SVD
#' @param normPcs    flag indicating whether to apply z-score normalization to PCs for each cell
#' @param umapParams  parameters to compute UMAP coordinates (passed on to \code{muRtools::getDimRedCoords.umap} and further to \code{uwot::umap})
#' @return an \code{S3} object containing dimensionality reduction results
#' 
#' @details
#' The output object includes the final singular values/principal components (\code{result$pcaCoord}), the low-dimensional coordinates (\code{result$umapCoord}) as well as region set that provided the basis for the dimension reduction (\code{result$regionGr}).
#' 
#' @rdname dimRed_UMAP-DsATACsc-method
#' @docType methods
#' @aliases dimRed_UMAP
#' @aliases dimRed_UMAP,DsATACsc-method
#' @author Fabian Mueller
#' @export
setMethod("dimRed_UMAP",
	signature(
		.object="DsATACsc"
	),
	function(
		.object,
		regions,
		tfidf=TRUE,
		pcs=1:50,
		normPcs=TRUE,
		umapParams=list(
			distMethod="euclidean",
			min_dist=0.5,
			n_neighbors=25
		)
	) {
		rt <- regions
		gr <- NULL
		if (is.character(rt)){
			if (!is.element(rt, getRegionTypes(.object))){
				logger.error(c("Invalid region type:", rt))
			}
		}
		if (is.element(class(rt), c("GRanges"))){
			logger.start("Aggregating fragment counts")
				gr <- rt
				rt <- ".regionsForDimRed"
				.object <- regionAggregation(.object, gr, rt, signal="insertions", dropEmpty=FALSE, bySample=FALSE)
			logger.completed()
		}

		dsn <- .object
		idfBase <- NULL
		if (tfidf){
			logger.start("Transforming counts using TF-IDF")
				bcm_unnorm <- ChrAccR::getCounts(dsn, rt, allowSparseMatrix=TRUE) > 0 # unnormalized binary count matrix
				idfBase <- log(1 + ncol(bcm_unnorm) / safeMatrixStats(bcm_unnorm, "rowSums", na.rm=TRUE))
				dsn <- transformCounts(dsn, method="tf-idf", regionTypes=rt) #TODO: renormalize based on sequencing depth rather than aggregated counts across peaks only?
			logger.completed()
		}

		gr <- getCoord(dsn, rt)
		cm <- ChrAccR::getCounts(dsn, rt, allowSparseMatrix=TRUE)

		mat <- cm
		pcaCoord <- NULL
		if (length(pcs) > 1){
			logger.start("SVD")
				pcaCoord <- muRtools::getDimRedCoords.pca(safeMatrixStats(cm, "t"), components=1:max(pcs), method="irlba_svd")
				mat <- pcaCoord
				if (normPcs) {
					logger.info("Scaling SVDs")
					mat <- rowZscores(mat, na.rm=TRUE)
				}
				mat <- mat[, pcs, drop=FALSE]
			logger.completed()
		}
		
		logger.start(c("UMAP dimension reduction"))
			paramL <- c(list(X=mat), umapParams)
			umapCoord <- do.call(muRtools::getDimRedCoords.umap, paramL)
			umapRes <- attr(umapCoord, "umapRes")
			attr(umapCoord, "umapRes") <- NULL
		logger.completed()

		res <- list(
			pcaCoord=pcaCoord,
			pcs = pcs,
			idfBase=idfBase,
			umapCoord=umapCoord,
			umapRes=umapRes,
			regionGr=gr,
			.params=list(normPcs=normPcs)
		)
		class(res) <- "DimRed_UMAP_sc"
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
#' Perform iterative LSI clustering and dimension reduction as described in doi:10.1038/s41587-019-0332-7
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param it0regionType character string specifying the region type to start with
#' @param it0nMostAcc the number of the most accessible regions to consider in iteration 0
#' @param it0pcs      the principal components to consider in iteration 0
#' @param it0clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in iteration 0
#' @param it0clusterMinCells the minimum number of cells in a cluster in order for it to be considered in peak calling (iteration 0)
#' @param it0nTopPeaksPerCluster the number of best peaks to be considered for each cluster in the merged peak set (iteration 0)
#' @param it1pcs      the principal components to consider in iteration 0
#' @param it1clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in iteration 1
#' @param it1mostVarPeaks the number of the most variable peaks to consider after iteration 1
#' @param it2pcs      the principal components to consider in the final iteration (2)
#' @param it2clusterResolution resolution paramter for Seurat's  clustering (\code{Seurat::FindClusters}) in the final iteration (2)
#' @param rmDepthCor  coreelation cutoff to be used to discard principal components associated with fragment depth (iteration 0)
#' @param normPcs     flag indicating whether to apply z-score normalization to PCs for each cell (all iterations)
#' @param umapParams  parameters to compute UMAP coordinates (passed on to \code{muRtools::getDimRedCoords.umap} and further to \code{uwot::umap})
#' @return an \code{S3} object containing dimensionality reduction results, peak sets and clustering
#' 
#' @details
#' In order to obtain a low dimensional representation of single-cell ATAC datasets in terms of principal components and UMAP coordinates, we recommend an iterative application of the Latent Semantic Indexing approach [10.1016/j.cell.2018.06.052] described in [doi:10.1038/s41587-019-0332-7]. This approach also identifies cell clusters and a peak set that represents a consensus peak set of cluster peaks in a given dataset. In brief, in an initial iteration clusters are identified based on the most accessible regions (e.g. genomic tiling regions). Here, the counts are first normalized using the term frequency–inverse document frequency (TF-IDF) transformation and singular values are computed based on these normalized counts in selected regions (i.e. the most accessible regions in the initial iteration). Clusters are identified based on the singular values using Louvain clustering (as implemented in the \code{Seurat} package). Peak calling is then performed on the aggregated insertion sites from all cells of each cluster (using MACS2) and a union/consensus set of peaks uniform-length non-overlapping peaks is selected. In a second iteration, the peak regions whose TF-IDF-normalized counts which exhibit the most variability across the initial clusters provide the basis for a refined clustering using derived singular values. In the final iteration, the most variable peaks across the refined clusters are identified as the final peak set and singular values are computed again. Based on these final singular values UMAP coordinates are computed for low-dimensional projection.
#' 
#' The output object includes the final singular values/principal components (\code{result$pcaCoord}), the low-dimensional coordinates (\code{result$umapCoord}), the final cluster assignment of all cells (\code{result$clustAss}), the complete, unfiltered initial cluster peak set (\code{result$clusterPeaks_unfiltered}) as well as the final cluster-variable peak set (\code{result$regionGr}).
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
		it0pcs=1:25,
		it0clusterResolution=0.8,
		it0clusterMinCells=200L,
		it0nTopPeaksPerCluster=2e5,
		it1pcs=1:50,
		it1clusterResolution=0.8,
		it1mostVarPeaks=50000L,
		it2pcs=1:50,
		it2clusterResolution=0.8,
		rmDepthCor=0.5,
		normPcs=TRUE,
		umapParams=list(
			distMethod="euclidean",
			min_dist=0.5,
			n_neighbors=25
		)
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
				pcs <- pcaCoord_it0
				if (normPcs) {
					logger.info("Scaling SVDs")
					pcs <- rowZscores(pcs, na.rm=TRUE) # z-score normalize PCs for each cell
				}
				it0fragCountCor <- NULL
				if (!is.null(depthV) && rmDepthCor > 0 && rmDepthCor < 1){
					it0fragCountCor <- apply(pcs, 2, FUN=function(x){
						cor(x, depthV, method="spearman")
					})
					idx <- which(it0fragCountCor > rmDepthCor)
					if (length(idx) > 0){
						rmStr <- paste(paste0("PC", idx, " (r=", round(it0fragCountCor[idx], 4), ")"), collapse=", ")
						logger.info(c("The following PCs are correlated (Spearman) with cell fragment counts and will be removed:", rmStr))
						it0pcs <- setdiff(it0pcs, idx)
					}
				}
				pcs <- pcs[, it0pcs, drop=FALSE]
				
			logger.completed()
			logger.start(c("Clustering"))	
				if (!requireNamespace("Seurat")) logger.error(c("Could not load dependency: Seurat"))
				# Louvain clustering using Seurat
				dummyMat <- matrix(11.0, ncol=length(cellIds), nrow=11)
				colnames(dummyMat) <- cellIds
				rownames(dummyMat) <- paste0("df", 1:nrow(dummyMat))
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")				
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcs, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcs), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it0clusterResolution)
				clustAss_it0 <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss_it0) <- names(clustRes@active.ident)
				logger.info(c("Number of clusters found:", nlevels(clustAss_it0)))
				ct <- table(clustAss_it0)
				peakCallClusters <- names(ct)[ct >= it0clusterMinCells]
				doExcludeClusters <- !all(levels(clustAss_it0) %in% peakCallClusters)
				if (doExcludeClusters){
					logger.info(c("Considering the following clusters for peak calling:", paste(peakCallClusters, collapse=",")))
				}
			logger.completed()
			logger.start(c("Peak calling"))
				logger.start("Creating cluster pseudo-bulk samples")
					ca <- as.character(clustAss_it0[cellIds])
					dsr <- addSampleAnnotCol(dsr, "clustAss_it0", ca)
					dsm <- dsr
					if (doExcludeClusters){
						dsm <- dsm[ca %in% peakCallClusters]
					}
					dsrClust <- mergeSamples(dsm, "clustAss_it0", countAggrFun="sum")
				logger.completed()
				logger.start("Calling peaks")
					clustPeakGrl <- callPeaks(dsrClust)
					if (!is.null(it0nTopPeaksPerCluster)){
						logger.info(paste0("Selecting the ", it0nTopPeaksPerCluster, " peaks with highest score for each cluster"))
						clustPeakGrl <- GRangesList(lapply(clustPeakGrl, FUN=function(x){
							idx <- rank(-elementMetadata(x)[,"score_norm"], na.last="keep", ties.method="min") <= it0nTopPeaksPerCluster
							x[idx]
						}))
					}
					
					peakUnionGr <- getNonOverlappingByScore(unlist(clustPeakGrl), scoreCol="score_norm")
					peakUnionGr <- sortGr(peakUnionGr)
					names(peakUnionGr) <- NULL
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
				pcaCoord_it1 <- muRtools::getDimRedCoords.pca(safeMatrixStats(cm, "t"), components=1:max(it1pcs), method="irlba_svd")
				pcs <- pcaCoord_it1
				if (normPcs) {
					logger.info("Scaling SVDs")
					pcs <- rowZscores(pcs, na.rm=TRUE)
				}
				pcs <- pcs[, it1pcs, drop=FALSE]
			logger.completed()

			logger.start(c("Clustering"))	
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcs, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcs), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it1clusterResolution)
				clustAss_it1 <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss_it1) <- names(clustRes@active.ident)
				logger.info(c("Number of clusters found:", nlevels(clustAss_it1)))
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
				umapRes <- dimRed_UMAP(dsr, it2regionType, tfidf=TRUE, pcs=it2pcs, normPcs=normPcs, umapParams=umapParams)
				pcaCoord_sel <- umapRes$pcaCoord[, umapRes$pcs, drop=FALSE]
				if (normPcs) pcaCoord_sel <- rowZscores(pcaCoord_sel, na.rm=TRUE)
			logger.completed()

			logger.start(c("Clustering"))
				sObj <- Seurat::CreateSeuratObject(dummyMat, project='scATAC', min.cells=0, min.features=0, assay="ATAC")
				sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings=pcaCoord_sel, key="PC_", assay="ATAC")
				sObj <- Seurat::FindNeighbors(sObj, reduction="pca", assay="ATAC", dims=1:ncol(pcaCoord_sel), k.param=30)
				clustRes <- Seurat::FindClusters(sObj, k.param=30, algorithm=1, n.start=100, n.iter=10, resolution=it2clusterResolution)
				clustAss <- factor(paste0("c", clustRes@active.ident), levels=paste0("c", levels(clustRes@active.ident)))
				names(clustAss) <- names(clustRes@active.ident)
				logger.info(c("Number of clusters found:", nlevels(clustAss)))

				dsr <- addSampleAnnotCol(dsr, "clustAss_it2", as.character(clustAss[cellIds]))
			logger.completed()
		logger.completed()

		res <- list(
			pcaCoord=umapRes$pcaCoord,
			pcs = umapRes$pcs,
			idfBase=umapRes$idfBase,
			umapCoord=umapRes$umapCoord,
			umapRes=umapRes$umapRes,
			clustAss=clustAss,
			regionGr=peakCoords,
			clusterPeaks_unfiltered=peakUnionGr,
			iterationData = list(
				iteration0 = list(
					pcaCoord=pcaCoord_it0,
					clustAss=clustAss_it0,
					pcs=it0pcs,
					pcCorFragmentCount=it0fragCountCor,
					nMostAcc=it0nMostAcc,
					clusterResolution=it0clusterResolution,
					clusterMinCells=it0clusterMinCells,
					nTopPeaksPerCluster=it0nTopPeaksPerCluster
				),
				iteration1 = list(
					pcaCoord=pcaCoord_it1,
					clustAss=clustAss_it1,
					pcs=it1pcs,
					clusterResolution=it1clusterResolution,
					mostVarPeaks=it1mostVarPeaks
				)
			),
			.params=c(list(normPcs=normPcs), callParams)
		)
		class(res) <- "iterativeLSIResultSc"
		return(res)
	}
)

#-------------------------------------------------------------------------------

#' getDiffAcc-methods
#'
#' Compute differential accessibility for single-cell datasets by randomly drawing cells from each group and aggregating them into pseudo-bulk samples
#' which are then compared using bulk differential methods
#'
#' @param .object    \code{\linkS4class{DsATACsc}} object
#' @param regionType character string specifying the region type
#' @param comparisonCol column name in the cell annotation table to base the comparison on. Alternatively, a vector with one element for each
#'                   cell in the dataset that can be coerced to a factor
#' @param grp1Name   name of the first group in the comparison. if not specified, it will be taken as the first factor level specified in the 
#'                   cell annotation (see \code{'comparisonCol'}).
#' @param grp2Name   name of the second group (reference) in the comparison. if not specified, it will be taken as the first factor level specified in the 
#'                   cell annotation (see \code{'comparisonCol'}).
#' @param nCellsPerBulk number of cells to sample to create each pseudo-bulk sample
#' @param nBulkPerGroup number of pseudo-bulk samples to create for each group
#' @param method     Method for determining differential accessibility. Currently only \code{'DESeq2'} is supported
#' @return a \code{data.frame} containing differential accessibility information
#' 
#' @rdname getDiffAcc-DsATACsc-method
#' @docType methods
#' @aliases getDiffAcc
#' @aliases getDiffAcc,DsATAC-method
#' @author Fabian Mueller
#' @export
#' @noRd
setMethod("getDiffAcc",
	signature(
		.object="DsATACsc"
	),
	function(
		.object,
		regionType,
		comparisonCol,
		grp1Name=NULL,
		grp2Name=NULL,
		nCellsPerBulk=100,
		nBulkPerGroup=20,
		method='DESeq2'
	) {
		ph <- getSampleAnnot(.object)
		if (!is.element(method, c("DESeq2"))) logger.error(c("Invalid method for calling differential accessibility:", method))
		if (is.character(comparisonCol) && length(comparisonCol)==1){
			if (!is.element(comparisonCol, colnames(ph))) logger.error(c("Comparison column not found in sample annotation:", comparisonCol))
			contrastF <- factor(ph[,comparisonCol])
		} else if (length(comparisonCol)==nrow(ph)){
			contrastF <- factor(comparisonCol)
		} else {
			logger.error("Invalid value for comparisonCol")
		}
		if (length(levels(contrastF)) < 2)  logger.error(c("Invalid comparison column. There should be at least 2 groups."))

		if (is.null(grp1Name)) grp1Name <- levels(contrastF)[1]
		if (is.null(grp2Name)) grp2Name <- levels(contrastF)[2]
		if (!is.element(grp1Name, c(levels(contrastF), ".ALL"))) logger.error(c("Invalid group name (1). No cells annotated with that group:", grp1Name))
		if (!is.element(grp2Name, c(levels(contrastF), ".ALL"))) logger.error(c("Invalid group name (2). No cells annotated with that group:", grp2Name))
		cidx.grp1 <- which(contrastF==grp1Name)
		if (grp1Name==".ALL") cidx.grp1 <- which(contrastF!=grp2Name)
		cidx.grp2 <- which(contrastF==grp2Name)
		if (grp2Name==".ALL") cidx.grp2 <- which(contrastF!=grp1Name)

		if (method=="DESeq2"){
			logger.info(c("Using method:", method))
			cm <- ChrAccR::getCounts(.object, regionType, allowSparseMatrix=TRUE)

			logger.start("Creating pseudo-bulk samples")
				nCells <- min(c(length(cidx.grp1), length(cidx.grp2)))
				doBoostrap <- FALSE
				if (nCells < nCellsPerBulk){
					logger.warning(c("Few cells detected per group", "--> selecting only", nCells, "cells for sampling"))
					doBoostrap <- TRUE
				} else {
					nCells <- nCellsPerBulk
				}
				logger.info(c("Using", nCells, "cells per sample"))
				logger.info(c("Using", nBulkPerGroup, "samples per group"))

				cidxL.grp1 <- lapply(1:nBulkPerGroup, FUN=function(i){
					sample(cidx.grp1, nCells, replace=doBoostrap)
				})
				cidxL.grp2 <- lapply(1:nBulkPerGroup, FUN=function(i){
					sample(cidx.grp2, nCells, replace=doBoostrap)
				})
				cm.grp1 <- do.call("cbind", lapply(cidxL.grp1, FUN=function(cids){
					safeMatrixStats(cm[,cids,drop=FALSE], statFun="rowSums", na.rm=TRUE)
				}))
				colnames(cm.grp1) <- paste(grp1Name, "sample", 1:nBulkPerGroup, sep="_")
				cm.grp2 <- do.call("cbind", lapply(cidxL.grp2, FUN=function(cids){
					safeMatrixStats(cm[,cids,drop=FALSE], statFun="rowSums", na.rm=TRUE)
				}))
				colnames(cm.grp2) <- paste(grp2Name, "sample", 1:nBulkPerGroup, sep="_")
			logger.completed()

			logger.start("Creating DESeq2 dataset")
				designF <- as.formula(paste0("~", paste("group", collapse="+")))
				sannot <- data.frame(sampleId=c(colnames(cm.grp1), colnames(cm.grp2)), group=rep(c(grp1Name, grp2Name), times=rep(nBulkPerGroup, 2)))
				dds <- DESeq2::DESeqDataSetFromMatrix(
					countData=cbind(cm.grp1, cm.grp2),
					colData=sannot,
					design=designF
				)
				rowRanges(dds) <- getCoord(.object, regionType)
				dds <- DESeq2::DESeq(dds)
			logger.completed()

			logger.start("Differential table")
				diffRes <- DESeq2::results(dds, contrast=c("group", grp1Name, grp2Name))
				dm <- data.frame(diffRes)
				rankMat <- cbind(
					# rank(-dm[,"baseMean"]), na.last="keep", ties.method="min"),
					rank(-abs(dm[,"log2FoldChange"]), na.last="keep", ties.method="min"),
					rank(dm[,"pvalue"], na.last="keep", ties.method="min")
				)
				dm[,"cRank"] <- matrixStats::rowMaxs(rankMat, na.rm=FALSE)
				# dm[,"cRank"] <- rowMaxs(rankMat, na.rm=TRUE)
				dm[!is.finite(dm[,"cRank"]),"cRank"] <- NA
				dm[,"cRank_rerank"] <- rank(dm[,"cRank"], na.last="keep", ties.method="min")

				sidx.grp1 <- which(sannot[,"group"]==grp1Name)
				sidx.grp2 <- which(sannot[,"group"]==grp2Name)
				l10fpkm <- log10(DESeq2::fpkm(dds, robust=TRUE)+1)
				grp1.m.l10fpkm <- rowMeans(l10fpkm[, sidx.grp1, drop=FALSE], na.rm=TRUE)
				grp2.m.l10fpkm <- rowMeans(l10fpkm[, sidx.grp2, drop=FALSE], na.rm=TRUE)
				vstCounts <- assay(DESeq2::vst(dds, blind=FALSE))
				grp1.m.vst <- rowMeans(vstCounts[, sidx.grp1, drop=FALSE], na.rm=TRUE)
				grp2.m.vst <- rowMeans(vstCounts[, sidx.grp2, drop=FALSE], na.rm=TRUE)

				res <- data.frame(
					log2BaseMean=log2(dm[,"baseMean"]),
					meanLog10FpkmGrp1=grp1.m.l10fpkm,
					meanLog10FpkmGrp2=grp2.m.l10fpkm,
					meanVstCountGrp1=grp1.m.vst,
					meanVstCountGrp2=grp2.m.vst,
					dm
				)
				# add group names to column names
				for (cn in c("meanLog10FpkmGrp", "meanVstCountGrp")){
					colnames(res)[colnames(res)==paste0(cn,"1")] <- paste0(cn, "1_", grp1Name)
					colnames(res)[colnames(res)==paste0(cn,"2")] <- paste0(cn, "2_", grp2Name)
				}
			logger.completed()
		}
		return(res)
	}
)

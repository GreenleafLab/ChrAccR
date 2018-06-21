################################################################################
# Utilities and helper functions for TF motif processing
################################################################################

################################################################################
# motifmatchr  helpers
################################################################################
#' prepareMotifmatchr
#' 
#' prepare objects for a \code{motifmatchr} analysis
#' @param genome character string specifying genome assembly
#' @param motifs either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'               that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object)
#' @return a list containing objects to be used as arguments for \code{motifmatchr}
#' @author Fabian Mueller
#' @export
prepareMotifmatchr <- function(genome, motifs){
	require(motifmatchr)
	res <- list()

	# get the species name and the genome sequence object based on the object
	genomeObj <- genome
	if (!is.element("BSgenome", class(genomeObj))){
		genomeObj <- getGenomeObject(genome)
	}
	spec <- organism(genomeObj)

	# get the motif PWMs
	motifL <- TFBSTools::PWMatrixList()
	if (is.character(motifs)){
		if (is.element("jaspar", motifs)){
			# copied code from chromVAR, but updated the JASPAR version
			opts <- list()
			opts["species"] <- spec
			opts["collection"] <- "CORE"
			# gets the non-redundant set by default
			mlCur <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
			if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))){
				names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
			} 
			motifL <- c(motifL, TFBSTools::toPWM(mlCur))
		}
		if (is.element("jaspar2016", motifs)){
			require(chromVAR)
			motifL <- c(motifL, TFBSTools::toPWM(getJasparMotifs(species=spec)))
		}
		if (is.element("homer", motifs)){
			require(chromVARmotifs)
			data("homer_pwms")
			motifL <- c(motifL, homer_pwms)
		}
		if (is.element("encode", motifs)){
			require(chromVARmotifs)
			data("encode_pwms")
			motifL <- c(motifL, encode_pwms)
		}
		if (is.element("cisbp", motifs)){
			require(chromVARmotifs)
			if (spec == "Mus musculus"){
				data("mouse_pwms_v2")
				motifL <- c(motifL, mouse_pwms_v2)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v2")
				motifL <- c(motifL, human_pwms_v2)
			} else {
				logger.warning(c("Could not find cisBP annotation for species", spec))
			}
		} 
		if (length(motifL) < 1) {
			logger.error(c("No motifs were loaded. Unsupported motifs (?) :", motifs))
		}	
	} else if (is.element("PWMatrixList", class(motifs)) || is.element("PFMatrixList", class(motifs))) {
		motifL <- motifs
	} else {
		logger.error(c("unsupported value for motifs:", motifs))
	}	
	res[["genome"]] <- genomeObj
	res[["motifs"]] <- motifL
	return(res)
}
################################################################################
# Motif similarity methods
################################################################################
#' getJasparSymbols
#'
#' Retrieve the TF names (symbols) from a JASPAR identifier
#'
#' @param ss     character vector or JASPAR identifiers
#' @return a list of TF names (symbols) for each identifier
getJasparSymbols <- function(ss){
	ss <- gsub("^MA[0-9]+\\.[1-9]_", "", ss) # jaspar prefix
	ss <- gsub("\\(var\\..*\\)", "", ss) # jaspar variation names
	return(strsplit(ss, "::")) # jaspar multiple TFs (separated by ::)
}

#' getMotifDistMat.jaspar
#'
#' Retrieve motif a comparison table from JASPAR annotation website and construct a dissimilarity matrix for given
#' motif IDs
#'
#' @param motifIds     string vector of motif ids whose dissimilarities are retrieved
#' @param scoreCol     namew of the annotation column in the JASPAR annotation that contains the motif similarity
#' @return a matrix of motif DISsimilarities
getMotifDistMat.jaspar <- function(motifIds=NULL, scoreCol="Ncor"){
	#motif comparison table from JASPAR matrix clustering results
	compFn <- "http://folk.uio.no/azizk/JASPAR_2018_clustering/results/JASPAR_2018_matrix_clustering/vertebrates/JASPAR_2018_matrix_clustering_vertebrates_tables/pairwise_compa.tab"
	# Format:
	# ;mode: matches	thresholds:	cor=-1.000000	ncor=-1.000000	w=0	ncor1=-1.000000	ncor2=-1.000000
	# #id1	id2	name1	name2	cor	Ncor	Ncor1	Ncor2	w1	w2	w	W	Wr	wr1	wr2	strand	offset	uncounted
	# JASPAR_2018_vertebrates_m1_MA0002.2	JASPAR_2018_vertebrates_m1_MA0002.2	RUNX1	RUNX1	1.000000	1.000000	1.000000	1.000000	11	11	11	11	1.000000	1.000000	1.000000	D	0	41
	# JASPAR_2018_vertebrates_m1_MA0002.2	JASPAR_2018_vertebrates_m2_MA0003.3	RUNX1	TFAP2A	0.518599	0.432166	0.471454	0.471454	11	11	10	12	0.833333	0.909091	0.909091	R	-1	41
	# JASPAR_2018_vertebrates_m1_MA0002.2	JASPAR_2018_vertebrates_m3_MA0004.1	RUNX1	Arnt	0.468110	0.255333	0.255333	0.468110	11	6	6	11	0.545455	0.545455	1.000000	D	1	31
	# JASPAR_2018_vertebrates_m1_MA0002.2	JASPAR_2018_vertebrates_m4_MA0006.1	RUNX1	Ahr::Arnt	0.533152	0.290810	0.290810	0.533152	11	6	6	11	0.545455	0.545455	1.000000	D	1	31
	# ...
	compTab <- read.table(compFn, header=TRUE, comment.char=";", "\t", check.names=FALSE, stringsAsFactors=FALSE)
	colnames(compTab) <- gsub("^#(.+)$", "\\1", colnames(compTab))
	extractJasparId <- function(ss){
		sapply(strsplit(ss, "_"), FUN=function(x){x[length(x)]})
	}
	compTab[,"motif1"] <- extractJasparId(compTab[,"id1"])
	compTab[,"motif2"] <- extractJasparId(compTab[,"id2"])
	motifIdsFromTab <- sort(union(compTab[,"motif1"], compTab[,"motif2"]))
	if (is.null(motifIds)){
		motifIds <- motifIdsFromTab
	} else {
		unknownMotifIds <- setdiff(motifIds, motifIdsFromTab)
		if (length(unknownMotifIds) > 0) logger.warning(c("The following motif ids were not found in the JASPAR clustering result table:", paste(unknownMotifIds, collapse=",")))
	}
	compTab <- compTab[compTab[,"motif1"] %in% motifIds & compTab[,"motif2"] %in% motifIds,]

	scoreMat <- matrix(as.numeric(NA), nrow=length(motifIds), ncol=length(motifIds))
	colnames(scoreMat) <- rownames(scoreMat) <- motifIds
	for(i in 1:nrow(compTab)){
		# if (1 %% 100 == 0) print(i)
		scoreMat[compTab[i, "motif1"], compTab[i, "motif2"]] <- compTab[i, scoreCol]
	}
	distMat <- scoreMat
	if (is.element(scoreCol, c("cor", "Ncor"))){
		# for correlation-based similarities the distance is 1-cor
		distMat <- 1-scoreMat
	}
	return(distMat)
}

#' getMotifDistMat
#'
#' Retrieve motif dissimilarity/distance matrix for TF motifs 
#'
#' @param assembly     genome assembly for which the motifs dissimilarity should be retrieved. Only the species information
#'                     of the assembly is really relevant
#' @param mmObj        optional motifmatchr object as returned by \code{ChrAccR::prepareMotifmatchr}
#' @param method       method of dissimilarity quantification. Currently only \code{'jaspar'} (retrieve motif similarities from the annotation of the JASPAR website) is supported.
#' @return a matrix of motif DISsimilarities (\code{dist} object)
#' @author Fabian Mueller
#' @export
getMotifDistMat <- function(assembly="hg38", mmObj=NULL, method="jaspar"){
	require(TFBSTools)
	if (method=="jaspar"){
		spec <- muRtools::normalize.str(organism(getGenomeObject(assembly)))
		if (is.null(mmObj))	mmObj <- prepareMotifmatchr(assembly, "jaspar")
		pwmL <- mmObj$motifs
		fn <- system.file(file.path("extdata", paste0("motifDistMat_jaspar_", spec, ".rds")), package="ChrAccR")
		if (file.exists(fn)){
			mDist <- readRDS(fn)
			lls <- labels(mDist)
			if (!all(lls==names(pwmL))){
				mn_missing <- setdiff(names(pwmL), lls)
				if (length(mn_missing) > 0) logger.error(c("The following motifs could not be found in the JASPAR annotation:", paste(mn_missing, collapse=",")))
				mDist <- as.dist(as.matrix(mDist)[names(pwmL),names(pwmL)])
			}
		} else {
			logger.warning(c("ChrAccR currently does not contain a precomputed JASPAR motif dissimilarity matrix for species", spec, "--> dissimilarities will be retrieved from the JASPAR website"))		
			mIds <- TFBSTools::"ID"(pwmL)
			mNames <- TFBSTools::"name"(pwmL)

			mDist <- getMotifDistMat.jaspar(mIds)
			colnames(mDist) <- rownames(mDist) <- names(pwmL)
			mDist <- as.dist(mDist)
		}
	} else {
		logger.error(c("Unsupported method of motif dissiliarities:", method))
	}
	return(mDist)
}

#' getMotifClustering
#'
#' Retrieve motif clustering of TF motifs
#'
#' @param distM        distance matrix (\code{dist} object) containing motif dissimilarities/distances
#' @param assembly     genome assembly for which the motifs dissimilarity should be retrieved. Only the species information
#'                     of the assembly is really relevant
#' @param motifs either a character string (currently only "jaspar" is supported) or an object containing PWMs
#'               that can be used by \code{motifmatchr::matchMotifs} (\code{PWMatrixList} object)
#' @param distMethod     method of dissimilarity quantification. Currently only \code{'jaspar'} (retrieve motif similarities from the annotation of the JASPAR website) is supported.
#' @param clusterMethod  method to be used for motif clustering (currently only \code{'pam'} (PAM - partitioning around medoids) is supported)
#' @param k            number of clusters. \code{k<1} will result in an automatically selected clustering which is precomputed and stored in \code{ChrAccR}.
#'                     For \code{distMethod=="jaspar"} and \code{clusterMethod=="pam"} this corresponds to the k corresponding to the best silhouette value before a drop (in the silhouette elbow-curve) occurs
#' @return a matrix of motif DISsimilarities (\code{dist} object)
#' @author Fabian Mueller
#' @export
getMotifClustering <- function(distM, assembly="hg38", motifs="jaspar", clusterMethod="pam", k=0){
	if (motifs != "jaspar") logger.error(c("Currently motif clustering is only supported for JASPAR motifs"))
	if (clusterMethod != "pam") logger.error(c("Currently motif clustering is only supported using the PAM clustering method"))
	if (!is.element("dist", class(distM))) logger.error(c("motif distance matrix must be a dist object"))

	spec <- muRtools::normalize.str(organism(getGenomeObject(assembly)))
	motifIds <- labels(distM)

	if (motifs=="jaspar"){
		if (clusterMethod == "pam"){
			if (k<1){
				#auto
				fn <- system.file(file.path("extdata", paste0("motifClustRes_bestSil_jaspar_", spec, ".rds")), package="ChrAccR")
				if (!file.exists(fn)) logger.error(c("ChrAccR currently does not contain a precomputed JASPAR motif clustering for species", spec))
				cr <- readRDS(fn)
			} else {
				require(cluster)
				clustRes.pam <- pam(distM, k=k)
				clustAssign <- clustRes.pam$medoids[clustRes.pam$clustering]
				names(clustAssign) <- motifIds
				clustAssignL <- lapply(clustRes.pam$medoids, FUN=function(mm){names(clustAssign)[clustAssign==mm]})
				names(clustAssignL) <- clustRes.pam$medoids
				clusterNames <- sapply(seq_along(clustAssignL), FUN=function(i){
					paste0(sapply(getJasparSymbols(names(clustAssignL)[i]), paste, collapse="_"), ":", paste(sort(unique(unlist(getJasparSymbols(clustAssignL[[i]])))), collapse="|"))
				})
				names(clusterNames) <- names(clustAssignL)
				cr <- list(k=k, clustAssign=clustAssign, clustAssignL=clustAssignL, clustNames=clusterNames, clustRes=clustRes.pam)
			}
		}
	}
	return(cr)
}

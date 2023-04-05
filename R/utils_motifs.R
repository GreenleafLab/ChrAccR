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
		if (is.element("jaspar_vert", motifs)){
			# JASPER for all vertebrate TFBS
			opts <- list()
			opts["tax_group"] <- "vertebrates"
			opts["collection"] <- "CORE"
			# gets the non-redundant set by default
			mlCur <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
			if (!isTRUE(all.equal(TFBSTools::name(mlCur), names(mlCur)))){
				names(mlCur) <- paste(names(mlCur), TFBSTools::name(mlCur), sep = "_")
			} 
			motifL <- c(motifL, TFBSTools::toPWM(mlCur))
		}
		if (is.element("jaspar2016", motifs)){
			motifL <- c(motifL, TFBSTools::toPWM(chromVAR::getJasparMotifs(species=spec)))
		}
		if (is.element("homer", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			data("homer_pwms")
			motifL <- c(motifL, chromVARmotifs::homer_pwms)
		}
		if (is.element("encode", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			data("encode_pwms")
			motifL <- c(motifL, chromVARmotifs::encode_pwms)
		}
		if (is.element("cisbp", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			if (spec == "Mus musculus"){
				data("mouse_pwms_v1")
				motifL <- c(motifL, chromVARmotifs::mouse_pwms_v1)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v1")
				motifL <- c(motifL, chromVARmotifs::human_pwms_v1)
			} else {
				logger.warning(c("Could not find cisBP annotation for species", spec))
			}
		}
		if (is.element("cisbp_v2", motifs)){
			if (!requireNamespace("chromVARmotifs")) logger.error(c("Could not load dependency: chromVARmotifs"))
			if (spec == "Mus musculus"){
				data("mouse_pwms_v2")
				motifL <- c(motifL, chromVARmotifs::mouse_pwms_v2)
			} else if (spec == "Homo sapiens"){
				data("human_pwms_v2")
				motifL <- c(motifL, chromVARmotifs::human_pwms_v2)
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

#' getMotifOccurrences
#' 
#' Find occurrences of motifs in a given genome
#' @param motifNames character vector of motif names
#' @param motifDb either a character string (currently only "jaspar" and sets contained in \code{chromVARmotifs} ("homer", "encode", "cisbp") are supported) or an object containing PWMs
#'               that can be used by \code{motifmatchr::matchMotifs} (such as an \code{PFMatrixList} or \code{PWMatrixList} object)
#' @param genome character string specifying genome assembly
#' @return a \code{GenomicRangesList} containing motif occurrences
#' @author Fabian Mueller
#' @export
getMotifOccurrences <- function(motifNames=NULL, motifDb="jaspar", genome="hg38"){
	mmArgs <- prepareMotifmatchr(genome, motifDb)
	if (is.null(motifNames)) motifNames <- names(mmArgs$motifs)
	motifNames_in <- motifNames
	motifNames <- intersect(motifNames, names(mmArgs$motifs))
	missingMotifs <- setdiff(motifNames_in, motifNames)
	if (length(missingMotifs) > 0){
		logger.warning(c(length(missingMotifs), "of", length(motifNames_in), "motifs were not found in the annotation database:", paste(missingMotifs, collapse=", ")))
	}

	seqlengths <- muRtools::getSeqlengths4assembly(genome, onlyMainChrs=TRUE)
	validChroms <- names(seqlengths)
	# for each chromosome, get the motif matches
	gr <- do.call("c", lapply(validChroms, FUN=function(chromName){
		# logger.status(c("chr:", chromName))#TODO: remove comment
		mmRes <- motifmatchr::matchMotifs(mmArgs$motifs[motifNames], mmArgs$genome[[chromName]], out="positions")
		# convert the resulting IRangesList to GRanges
		return(do.call("c", lapply(names(mmRes), FUN=function(motifName){
			x <- unlist(mmRes[[motifName]])
			rr <- GRanges()
			if (length(x) > 0){
				elementMetadata(x)[,"motifName"] <- motifName
				strands <- elementMetadata(x)[,"strand"]
				elementMetadata(x)[,"strand"] <- NULL
				rr <- GRanges(seqnames=chromName,ranges=x,strand=strands)
			}
			seqlevels(rr) <- validChroms
			seqlengths(rr) <- seqlengths
			genome(rr) <- genome
			return(rr)
		})))
	}))
	elementMetadata(gr)[,"motifName"] <- factor(elementMetadata(gr)[,"motifName"], levels=motifNames)
	resGrl <- split(gr, elementMetadata(gr)[,"motifName"])
	return(resGrl)
}

#' permutePWMatrix
#'
#' randomly permute the columns of a \code{TFBSTools::PWMatrix} object
#' (strangely the TFBSTools package has a method for \code{PFMatrix} objects, but not for \code{PWMatrix} objects)
#'
#' @param pwm   \code{PFMatrix} or \code{PWMatrix} object whose columns will be permuted
#' @param nperm number of permutations
#' @return a list of (length \code{nperm}) \code{PFMatrix} or \code{PWMatrix} objects representing the permuted versions
#' @author Fabian Mueller
#' @noRd
## @export
permutePWMatrix <- function(pwm, nperm=100){
	res <- list()
	if (class(pwm)=="PFMatrix"){
		res <- lapply(seq_len(nperm), FUN=function(i){TFBSTools::permuteMatrix(pwm, type="intra")})
	} else if (class(pwm)=="PWMatrix"){
		mm <- TFBSTools::Matrix(pwm)
		mLength <- ncol(mm)
		res <- lapply(seq_len(nperm), FUN=function(i){
			x <- pwm
			TFBSTools::Matrix(x) <- mm[,sample.int(mLength)]
			return(x)
		})
	} else {
		logger.error(c("Unable to permute motif matrix of type", class(pwm)))
	}
	return(res)
}
#' permutePWMatrixList
#'
#' randomly permute the columns of each element in a list of motif matrix objects
#'
#' @param pwmL  a list of \code{PFMatrix} or \code{PWMatrix} objects or a \code{PFMatrixList} or \code{PWMatrixList} object with the matrices to permute
#' @param nperm number of permutations applied to each matrix
#' @return a list in which each element is a lists of (length \code{nperm}) \code{PFMatrix} or \code{PWMatrix} objects representing the permuted versions for
#'         each element in \code{pwmL}
#' @author Fabian Mueller
#' @noRd
## @export
permutePWMatrixList <- function(pwmL, nperm=100){
	res <- lapply(pwmL, FUN=function(x){
		permutePWMatrix(x, nperm=nperm)
	})
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
	return(strsplit(ss, "::")) # jaspar multiple TFs (separated by ::); # TODO: treat fusion genes, e.g. "MA0149.1_EWSR1-FLI1"
}

#' getJasparAnnot
#' 
#' retrieve motif annotation data
#' @param ss     character vector or JASPAR identifiers
#' @param type annotation type. Currently only \code{"humantfs"} (pulls info from humantfs.ccbr.utoronto.ca) is supported
#' @return list of data frames of TF annotation (a motif can have multiple annotated TFs) 
#' @author Fabian Mueller
#' @export
getJasparAnnot <- function(ss, type="humantfs"){
	tfa <- getTfAnnot(type)
	res <- lapply(ss, FUN=function(x){
		idx <- grepl(x, tfa[, "motifIds_jaspar"], fixed=TRUE)
		return(tfa[idx,])
	})
	names(res) <- ss
	return(res)
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
	fn <- system.file(file.path("extdata", paste0("motifDistMat_jaspar_", scoreCol, ".rds")), package="ChrAccR")
	if (file.exists(fn)){
		distMat <- readRDS(fn)
		if (!is.null(motifIds)){
			distMat <- distMat[motifIds, motifIds]
		}
	} else {
		logger.info(c("ChrAccR currently does not contain a precomputed JASPAR motif dissimilarity matrix for scoretype", scoreCol, "--> dissimilarities will be retrieved from the JASPAR website"))	
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
	}
	return(distMat)
}

#' getMotifDistMat
#'
#' Retrieve motif dissimilarity/distance matrix for TF motifs 
#'
#' @param assembly     genome assembly for which the motifs dissimilarity should be retrieved. Only the species information
#'                     of the assembly is really relevant. Can be \code{"vert"} for all vertebrate motifs.
#' @param mmObj        optional motifmatchr object as returned by \code{ChrAccR::prepareMotifmatchr}
#' @param method       method of dissimilarity quantification. Currently only \code{'jaspar'} (retrieve motif similarities from the annotation of the JASPAR website) is supported.
#' @return a matrix of motif DISsimilarities (\code{dist} object)
#' @author Fabian Mueller
#' @export
getMotifDistMat <- function(assembly="hg38", mmObj=NULL, method="jaspar"){
	if (method=="jaspar"){
		spec <- assembly
		if (assembly=="vert") {
			if (is.null(mmObj))	mmObj <- prepareMotifmatchr("hg38", "jaspar_vert")
		} else {
			spec <- muRtools::normalize.str(organism(getGenomeObject(assembly)))
			if (is.null(mmObj))	mmObj <- prepareMotifmatchr(assembly, "jaspar")
		}
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
#' @param k            number of clusters. \code{k<1} will result in an automatically selected clustering which is precomputed and stored in \code{ChrAccR}.
#'                     For \code{distMethod=="jaspar"} and \code{clusterMethod=="pam"} this corresponds to the k corresponding to the best silhouette value before a drop (in the silhouette elbow-curve) occurs
#' @param distM        distance matrix (\code{dist} object) containing motif dissimilarities/distances. Only required if \code{k>0}.
#' @param assembly     genome assembly for which the motifs dissimilarity should be retrieved. Only the species information
#'                     of the assembly is really relevant. Can be \code{"vert"} for all vertebrate motifs. Only required if for automatic mode (i.e. \code{k<1}).
#' @param motifs       a character string specifying the motif set (currently only "jaspar" is supported)
#' @param clusterMethod  method to be used for motif clustering (currently only \code{'pam'} (PAM - partitioning around medoids) is supported)

#' @return a list structure containing the clustering result
#' @author Fabian Mueller
#' @export
getMotifClustering <- function(k=0, distM=NULL, assembly="hg38", motifs="jaspar", clusterMethod="pam"){
	if (motifs != "jaspar") logger.error(c("Currently motif clustering is only supported for JASPAR motifs"))
	if (clusterMethod != "pam") logger.error(c("Currently motif clustering is only supported using the PAM clustering method"))
	if (k>0 && !is.element("dist", class(distM))) logger.error(c("motif distance matrix must be a dist object"))

	if (motifs=="jaspar"){
		if (clusterMethod == "pam"){
			if (k<1){
				#auto
				spec <- assembly
				if (assembly != "vert")	spec <- muRtools::normalize.str(organism(getGenomeObject(assembly)))
				fn <- system.file(file.path("extdata", paste0("motifClustRes_bestSil_jaspar_", spec, ".rds")), package="ChrAccR")
				if (!file.exists(fn)) logger.error(c("ChrAccR currently does not contain a precomputed JASPAR motif clustering for species", spec))
				cr <- readRDS(fn)
			} else {
				clustRes.pam <- cluster::pam(distM, k=k)
				clustAssign <- clustRes.pam$medoids[clustRes.pam$clustering]
				names(clustAssign) <- labels(distM)
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

#' collapseMotifMatrix
#'
#' Collapse TF motif matrix of arbitrary values by aggregating values over motif cluster assignment
#'
#' @param X            matrix to be collapsed. Must have the motif names as rownames. E.g. matrix obtained by \code{chromVAR::deviationScores}
#' @param motifClust   optional: motif clustering computed by \code{\link{getMotifClustering}}. If \code{NULL} (default) the default clustering will be retrieved
#' @param assembly     genome assembly for which the motif clustering should be retrieved. Only required if for automatic mode (i.e. \code{motifClust=NULL}).
#' @param motifs       a character string specifying the motif set (currently only "jaspar" is supported)
#' @param aggrFun      function to use to aggregate values
#' @param usePresent  if TRUE and, not all the motifClust motifs present within matrix it'll use only the present ones.
#' @return list containing two elements: \code{X}: Collapsed matrix containing motif cluster aggregated values;
#'         \code{clustering}: clustering result used for aggregation (see\code{\link{getMotifClustering}} for details)
#' @author Fabian Mueller
#' @export
collapseMotifMatrix <- function(X, motifClust=NULL, assembly="hg38", motifs="jaspar", aggrFun=mean,usePresent=FALSE){
	if (motifs != "jaspar") logger.error(c("Currently motif clustering is only supported for JASPAR motifs"))

	if (is.null(motifClust)){
		motifClust <- getMotifClustering(k=0, distM=NULL, assembly=assembly, motifs=motifs, clusterMethod="pam")
	}
	
	if(usePresent){
		names <- rownames(X)[which(rownames(X) %in% names(motifClust$clustAssign))]
		X <- X[which(rownames(X) %in% names),]
		motifClust$clustAssig <- motifClust$clustAssign[which(names(motifClust$clustAssign) %in% names)]
	}
	
	if (is.null(rownames(X)) || !all(rownames(X) %in% names(motifClust$clustAssign))){
		logger.error("X must have valid rownames that are all contained in the clustering")
	}

	df <- reshape2::melt(X)
	colnames(df)[1:2] <- c("motif", "sample")
	df[,"cluster"] <- motifClust$clustAssign[as.character(df[,"motif"])]

	Xc <- reshape2::acast(df, formula=cluster~sample, fun.aggregate=aggrFun, value.var="value")

	if (!is.null(motifClust$clustNames)){
		rownames(Xc) <- motifClust$clustNames[rownames(Xc)]
	}

	res <- list(
		X = Xc,
		clustering = motifClust
	)
	class(res) <- "CollapsedMotifMatrix"
	return(res)
}

################################################################################
# Motif plotting
################################################################################
#' PWMatrixToProbMatrix
#'
#' convert a log2probratio PWM (\code{PWMatrix} from TFBSTools package) to a matrix containing probabilities in [0,1]
#'
#' @param x log2probratio PWM (\code{PWMatrix} from TFBSTools package)
#' @return PWM probability matrix with values in 
#' @author Fabian Mueller
PWMatrixToProbMatrix <- function(x){
	if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
	(2^as(x, "matrix"))*TFBSTools::bg(x)/sum(TFBSTools::bg(x))
}

#' hmSeqLogo
#'
#' Draw a sequence motif logo in a Complex Heatmap using grid.
#' adapted from \code{seqLogo::seqLogo()}
#'
#' @param pwm   PWM (from TFBSTools package)
#' @param x     x center coordinate where the motif should be drawn
#' @param y     y center coordinate where the motif should be drawn
#' @param width drawing width
#' @param height drawing height
#' @param ic.scale \code{logical} If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.
#' @return Draws the motif
#' @author Fabian Mueller
#' @export
#' @examples
#' \dontrun{
#' mm <- prepareMotifmatchr("hg38", "jaspar")$motifs[["MA0137.3_STAT1"]]
#' hmSeqLogo(mm, unit(0.5, "npc"), unit(0.5, "npc"), 0.5, 0.5, ic.scale=TRUE)
#' }
hmSeqLogo <- function(pwm, x=unit(0.5, "npc"), y=unit(0.5, "npc"), width=1, height=1, ic.scale=TRUE){
	if (!requireNamespace("grid")) logger.error(c("Could not load dependency: grid"))
	# # convert units to numbers
	# unitType <- attr(x, "unit")
	unitType <- "npc"
	x <- as.numeric(x)
	y <- as.numeric(y)
	width <- as.numeric(width)
	height <- as.numeric(height)

	# convert the PWM to matrix
	if (class(pwm) == "pwm") {
		pwm <- pwm@pwm
	} else if (class(pwm) == "PWMatrix") {
		pwm <- PWMatrixToProbMatrix(pwm)
	} else if (class(pwm) == "data.frame") {
		pwm <- as.matrix(pwm)
	} else if (class(pwm) != "matrix"){
		stop("pwm must be of class matrix or data.frame")
	}
	if (any(abs(1 - apply(pwm,2,sum)) > 0.01)) stop("Columns of PWM must add up to 1.0")

	chars <- c("A","C","G","T")
	letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
	npos <- ncol(pwm)

	if (ic.scale) {
		facs <- seqLogo:::pwm2ic(pwm)
		facs <- facs/max(facs) # scale columns to max information content
	} else {
		facs <- rep(1, npos)
	}

	wt <- width / npos
	x.pos <- x - width/2
	for (j in 1:npos) {
		column <- pwm[,j]
		hts <- 0.99*column*facs[j]*height
		letterOrder <- order(hts)

		y.pos <- y-height/2
		for (i in 1:length(chars)) {
			letter <- chars[letterOrder[i]]
			ht <- hts[letterOrder[i]]
			if (ht>0) letters <- seqLogo:::addLetter(letters, letter, x.pos, y.pos, ht, wt, fill=c("A"="green","C"="red","G"="blue","T"="orange"))
			y.pos <- y.pos + ht #+ 0.01
		}
		x.pos <- x.pos + wt
	}
	# print(str(letters))
	grid::grid.polygon(x=unit(letters$x, unitType), y=unit(letters$y, unitType), id=letters$id, gp=grid::gpar(fill=letters$fill,col="transparent"))
}

################################################################################
# Footprint plotting
################################################################################
#' plotFootprint
#'
#' Plot motif footprinting results
#'
#' @param footprintDf \code{data.frame} used for plotting (as output by \code{getMotifFootprints(...)$footprintDf})
#' @param pwm   PWM used for plotting the sequence logo (from TFBSTools package or output from \code{prepareMotifmatchr(...)$motifs})
#' @return A \code{ggplot2} object containing the motif plot
#' @author Fabian Mueller
#' @noRd
plotFootprint <- function(footprintDf, pwm=NULL, colorMap=NULL){
	motifLen <- motifUp <- motifDown <- NULL
	if (!is.null(pwm)){
		motifLen <- length(pwm)
		motifUp <- -ceiling(motifLen/2) + 1
		motifDown <- floor(motifLen/2)
	}
	# plot
	if (is.null(colorMap)){
		ggCs <- muRtools::ggAutoColorScale(footprintDf[,"sampleId"], method="color")
	} else {
		ggCs <- scale_color_manual(values=colorMap)
	}
	ppm <- ggplot(footprintDf, aes(x=pos, y=countNormBiasCor, color=sampleId, group=sampleId, fill=sampleId)) + ggCs
	ppb <- ggplot(footprintDf, aes(x=pos, y=Tn5biasNorm, color=sampleId, group=sampleId, fill=sampleId)) + ggCs

	if (!is.null(motifLen)){
		ppm <- ppm + annotate("rect", xmin=motifUp, xmax=motifDown, ymin=-Inf, ymax=Inf, fill="#d9d9d9") 
		ppb <- ppb + annotate("rect", xmin=motifUp, xmax=motifDown, ymin=-Inf, ymax=Inf, fill="#d9d9d9")
	}

	ppm <- ppm + geom_line() #+ geom_smooth(aes(group=cellType),size=2,alpha=0.15,method="gam",formula=y ~ s(x, bs = "cs"))
	ppb <- ppb + geom_line()

	if (!is.null(pwm)){
		# add sequence logo
		logo <- hmSeqLogo(pwm, unit(0.5, "npc"), unit(0.5, "npc"), 1, 1, ic.scale=TRUE)
		logoHeight <- (max(footprintDf[,"countNormBiasCor"], na.rm=TRUE) - min(footprintDf[,"countNormBiasCor"], na.rm=TRUE)) * 0.2
		logoY <- max(footprintDf[,"countNormBiasCor"], na.rm=TRUE) - logoHeight
		ppm <- ppm + cowplot::draw_plot(logo, 100, logoY, 100, logoHeight)
	}

	pp <- cowplot::plot_grid(ppm, ppb + theme(legend.position="none"), ncol=1, rel_heights=c(2, 1), align="v", axis="lr")
	return(pp)
}

################################################################################
# Motif clusters
################################################################################
# chromVAR for Altius institute motif clusters
#' getMotifClusterAnnot_altius
#' 
#' Retrieve non-redundant motif cluster annotation curated by the altius institute (Vierstra et al., Nature, 2020; doi:10.1038/s41586-020-2528-x)
#' @param makeUnique further remove redundancy in motif occurences by retaining only the best matching motif when motif overlaps are present
#' @param genome genome assembly. Currently only \code{"hg38"} is supported
#' @param annotFn Filename for the annotation. Useful to skip automatic download of the large dataset
#' @return Motif cluster annotation
#' @author Fabian Mueller
#' @noRd
getMotifClusterAnnot_altius <- function(makeUnique=TRUE, genome="hg38", annotFn=NULL){
	if (is.null(annotFn)){
		if (genome=="hg38"){
			furl <- "https://muellerf.s3.amazonaws.com/data/ChrAccR/data/annotation/tfMotifClusters_hg38.rds"
			logger.start("Downloading annotation")
				annotFn <- tempfile("annot_", fileext=".rds")
				download.file(furl, annotFn)
			logger.completed()
		} else {
			logger.error("Unsupported genome")
		}
	}
	logger.start("Loading Altius TF motif cluster annotation data")
		mcAnnot <- readRDS(annotFn)
	logger.completed()
	if (makeUnique){
		logger.start("Making overlapping motifs unique")
			mcGr <- unlist(mcAnnot$clusterOcc)
			mcGr_u <- getNonOverlappingByScore(mcGr, scoreCol="score")
			mcGrl <- split(mcGr_u, elementMetadata(mcGr_u)[,"cid"])
			mcGrl <- mcGrl[names(mcAnnot$clusterOcc)]
			rm(mcGr, mcGr_u)
		logger.completed()
		mcAnnot$clusterOcc <- mcGrl
	}
	return(mcAnnot)
}

#' computeDeviations_altius
#' 
#' Compute chromVAR deviations from non-redundant motif clusters curated by the altius institute (Vierstra et al., Nature, 2020; doi:10.1038/s41586-020-2528-x)
#' @param dsa   \code{\linkS4class{DsATAC}} object
#' @param type  haracter string specifying the region type
#' @param mcAnnot motif cluster annotation object. If \code{NULL}, it will be automatically downloaded
#' @return Deviations object as returned by \code{chromVAR::computeDeviations}
#' @author Fabian Mueller
#' @noRd
computeDeviations_altius <- function(dsa, type, mcAnnot=NULL){
	if (is.null(mcAnnot)){
		mcAnnot <- getMotifClusterAnnot_altius(genome=getGenome(dsa))
	}
	logger.start("Preparing chromVAR object")
		gr <- getCoord(dsa, type)
		countSe <- getCountsSE(dsa, type, naIsZero=TRUE)
		motifClusterOccM <- as(do.call(cbind, lapply(mcAnnot$clusterOcc, FUN=function(x){
			overlapsAny(gr, x, ignore.strand=TRUE)
		})), "sparseMatrix")
		mmObj_mc <- SummarizedExperiment::SummarizedExperiment(assays=list(motifMatches=motifClusterOccM), rowRanges=gr, colData=S4Vectors::DataFrame(mcAnnot$clusterAnnot))
		genomeObj <- getGenomeObject(getGenome(dsa))
		genome(countSe) <- S4Vectors::metadata(genomeObj)$genome #override inconsistent naming of genome versions (e.g. hg38 and GRCh38)
		countSe <- chromVAR::addGCBias(countSe, genome=genomeObj)

		ridx <- safeMatrixStats(assay(countSe), "rowSums") > 0 & safeMatrixStats(assay(mmObj_mc), "rowSums") > 0 # only consider regions where there is an actual motif match and counts
		logger.status("Computing deviations")
		cv_mc <- chromVAR::computeDeviations(object=countSe[ridx,], annotations=mmObj_mc[ridx,])
	logger.completed()
	return(cv_mc)
}

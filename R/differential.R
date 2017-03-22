################################################################################
# Methods for differential chromatin accessibility analysis
################################################################################
#' getComparisonInfo
#'
#' retrieve the comparison information for an DsNOMe object. Analogous to \code{RnBeads::get.comparison.info}
#' @param dsn         \code{\linkS4class{DsNOMe}} object
#' @param cmpNames    column names of the sample annotation of the dataset that will be used for comparison
#' @param regionTypes which region types should be processed for differential analysis.
#' @param allPairs    Logical indicating whether all pairwise comparisons should be conducted, when more than 2 groups are present
#' @param adjPairCols argument passed on to \code{rnb.sample.groups}. See its documentation for details.
#' @param minGrpSize  Minimum number of samples required to form a group in comparison
#' @param maxGrpCount maximum number of groups to consider for comparisons
#' @return a list containing one element for each comparison to be conducted. Each element is again a list containing:
#' \describe{
#'   \item{\code{comparison}}{the name of the comparison}
#'   \item{\code{pheno.colname}}{the column name of the sample annotation table the comparison is derived from}
#'   \item{\code{group.names}}{the names of the two groups being compared}
#'   \item{\code{group.inds}}{the sample indices of the samples belonging to the two groups}
#'   \item{\code{paired}}{flag indicating whether paired analysis is conducted}
#'   \item{\code{adj.sva}}{flag indicating whether adjustment for SVA is conducted}
#'   \item{\code{adj.celltype}}{flag indicating whether adjustment for cell type is conducted}
#'   \item{\code{adjustment.table}}{the covariate adjustment table. \code{NULL} if the comparison is not adjusted}
#'   \item{\code{region.types}}{the region types applicable to the analysis}
#' }
#' @author Fabian Mueller
#' @export
getComparisonInfo <- function(dsn, cmpNames=NULL, regionTypes=getRegionTypes(dsn), allPairs=TRUE, adjPairCols=NULL, minGrpSize=1L, maxGrpCount=NULL){
	sannot <- getSampleAnnot(dsn)
	grpInfo <- rnb.sample.groups(sannot, cmpNames, columns.pairs=adjPairCols, min.group.size=minGrpSize, max.group.count=maxGrpCount)

	res <- list()
	for (i in 1:length(grpInfo)){
		groups <- grpInfo[[i]]
		is.paired <- attr(grpInfo,"paired")[i]
		cc <- names(grpInfo)[i]

		if (length(groups)>2){
			if (allPairs) {
				#perform all pairwise comparisons
				pps <- combn(1:length(groups),2)
				for (j in 1:ncol(pps)){
					pp1 <- pps[1,j]
					pp2 <- pps[2,j]
					ll1 <- names(groups)[pp1]
					ll2 <- names(groups)[pp2]
					grps <- c(ll1,ll2)
					comparison <- paste(paste(grps,collapse=" vs. ")," (based on ",cc,")",sep="")

					inds1 <- groups[[pp1]]
					inds2 <- groups[[pp2]]
					adjTab <- NULL

					res.cur <- list(
						comparison=comparison,
						pheno.colname=cc,
						group.names=grps,
						group.inds=list(group1=inds1,group2=inds2),
						paired=is.paired,
						adj.sva=FALSE,
						adj.celltype=FALSE,
						adjustment.table=adjTab,
						region.types=regionTypes
					)
					res.append <- list(res.cur)
					names(res.append) <- res.cur$comparison
					res <- c(res,res.append)
				}
			} else {
				for (j in 1:length(groups)){
					ll <- names(groups)[j]
					grps <- c(ll,paste("non.",ll,sep=""))
					if (is.paired){
						logger.warning(c("Paired analysis is not supported annotations with more than 2 categories and comparing one group vs. all others. ",
									"--> Using unpaired analysis.",
									"Consider enabling the differential.comparison.columns.all.pairwise option or reducing the number of groups ",
									"in this column to 2."))
						is.paired <- FALSE
					}
					comparison <- paste(paste(grps,collapse=" vs. ")," (based on ",cc,")",sep="")
					
					inds1 <- groups[[j]]
					inds2 <- setdiff(unlist(groups),inds1)
					
					adjTab <- NULL
					res.cur <- list(
						comparison=comparison,
						pheno.colname=cc,
						group.names=grps,
						group.inds=list(group1=inds1,group2=inds2),
						paired=is.paired,
						adj.sva=FALSE,
						adj.celltype=FALSE,
						adjustment.table=adjTab,
						region.types=regionTypes
					)
					res.append <- list(res.cur)
					names(res.append) <- res.cur$comparison
					res <- c(res,res.append)
				}
			}
		} else { # length(groups) == 2
			ll1 <- names(groups)[1]
			ll2 <- names(groups)[2]
			grps <- c(ll1,ll2)
			comparison <- paste0(paste(grps,collapse=" vs. ")," (based on ",cc,")")
			inds1 <- groups[[1]]
			inds2 <- groups[[2]]
			
			adjTab <- NULL

			res.cur <- list(
				comparison=comparison,
				pheno.colname=cc,
				group.names=grps,
				group.inds=list(group1=inds1,group2=inds2),
				paired=is.paired,
				adj.sva=FALSE,
				adj.celltype=FALSE,
				adjustment.table=adjTab,
				region.types=regionTypes
			)
			res.append <- list(res.cur)
			names(res.append) <- res.cur$comparison
			res <- c(res,res.append)
		}
	}
	return(res)
}

#-------------------------------------------------------------------------------
### computeDiffAcc.rnb.nome.bin.region
###
### computes a differential accessibitity (methylation in NOMe datasets) in the binary case (2 groups) on the region level.
### @author Fabian Mueller
### @aliases computeDiffAcc.rnb.nome
### @param dsn \code{\linkS4class{DsNOMe}} object
### @param dmtp differential methylation table on the site level (as obtained from \code{\link{RnBeads:::computeDiffMeth.bin.site}})
### @param inds.g1 column indices in \code{b} of group 1 members
### @param inds.g2 column indices in \code{b} of group 2 members
### @details
### Analogous to \code{RnBeads}' \code{computeDiffMeth.bin.region} function
### @return list of differential methylation tables
computeDiffAcc.rnb.nome.bin.region <- function(dsn, dmtp, inds.g1, inds.g2, regionTypes=getRegionTypes(dsn), ...){
	#sanity checks
	if (length(union(inds.g1,inds.g2)) != (length(inds.g1)+length(inds.g2))){
		logger.error("Overlapping sample sets in differential methylation analysis")
	}
	logger.start('Computing Differential Methylation Tables (Region Level)')
	skipSites <- FALSE
	if (is.null(dmtp)){
		logger.info("Computing differential methylation for regions directly (NOT using site-specific differential methylation)")
		skipSites <- TRUE
	}
	diffTabs <- list()
	for (rt in regionTypes){
		if (skipSites){
			covMat <- getCovg(dns, rt, asMatrix=TRUE)
			dmtr <- RnBeads:::computeDiffMeth.bin.site(getMeth(dns,rt), inds.g1, inds.g2, covg=covMat, ...)
		} else {
			inclCov <- !is.null(getCovg(dns, "sites"))
			regions2sites <- getRegionMapping(dns, rt)
			dmtr <- RnBeads:::computeDiffTab.default.region(dmtp, regions2sites, includeCovg=inclCov)
			dmtr4ranks <- RnBeads:::extractRankingCols.region(dmtr)
			combRank <- RnBeads:::combinedRanking.tab(dmtr4ranks, rerank=FALSE)
			dmtr$combinedRank <- combRank
		}
		diffTabs <- c(diffTabs,list(dmtr))
		logger.status(c("Computed table for", rt))
	}
	names(diffTabs) <- regionTypes
	logger.completed()
	return(diffTabs)
}

#' computeDiffAcc.rnb.nome
#'
#' computes differential accessibility for NOMe datasets using \code{RnBeads} functionality
#' @author Fabian Mueller
#' @aliases rnb.execute.computeDiffMeth
#' @param dsn         \code{\linkS4class{DsNOMe}} object
#' @param cmpCols     column names of the sample annotation of the dataset that will be used for comparison
#' @param regionTypes which region types should be processed for differential analysis.
#' @param covgThres   coverage threshold for computing the summary statistics. See \code{\link{RnBeads::computeDiffTab.extended.site}} for details.
#' @param allPairs    Logical indicating whether all pairwise comparisons should be conducted, when more than 2 groups are present
#' @param adjPairCols argument passed on to \code{rnb.sample.groups}. See its documentation for details.
#' @param adjCols     not used yet
#' @param skipSites   flag indicating whether differential methylation in regions should be computed directly and not from sites. This leads to skipping of site-specific differential methylation
#' @param disk.dump Flag indicating whether the resulting differential methylation object should be file backed, ie.e the matrices dumped to disk
#' @param disk.dump.dir disk location for file backing of the resulting differential methylation object. Only meaningful if \code{disk.dump=TRUE}.
#' 						must be a character specifying an NON-EXISTING valid directory.
#' @param ... arguments passed on to binary differential methylation calling. See \code{\link{RnBeads::computeDiffTab.extended.site}} for details.
#' @return an \code{\linkS4class{RnBDiffMeth}} object. See class description for details.
#' @author Fabian Mueller
#' @export
computeDiffAcc.rnb.nome <- function(dsn, cmpCols, regionTypes=getRegionTypes(dsn), covgThres=5L,
		allPairs=TRUE, adjPairCols=NULL,
		adjCols=NULL,
		skipSites=FALSE,
		disk.dump=rnb.getOption("disk.dump.big.matrices"),disk.dump.dir=tempfile(pattern="diffMethTables_"),
		...){

	logger.start("Retrieving comparison info")
	cmpInfo <- getComparisonInfo(dsn, cmpNames=cmpCols, regionTypes=regionTypes, allPairs=allPairs, adjPairCols=adjPairCols, minGrpSize=1L, maxGrpCount=NULL)
	logger.completed()
	if (is.null(cmpInfo)) {
		return(NULL)
	}

	diff.method <- "limma"
	logger.start("Computing differential methylation tables")

	diffmeth <- new("RnBDiffMeth",site.test.method=diff.method,disk.dump=disk.dump,disk.path=disk.dump.dir)
	
	for (i in 1:length(cmpInfo)){
		cmpInfo.cur <- cmpInfo[[i]]
		logger.start(c("Comparing:",cmpInfo.cur$comparison))
		if (cmpInfo.cur$paired){
			logger.status("Conducting PAIRED analysis")
		}

		if (skipSites){
			logger.info("Skipping site-specific differential methylation calling")
			dm <- NULL
		} else {
			dm <- RnBeads:::computeDiffMeth.bin.site(
					getMeth(dsn, asMatrix=TRUE), inds.g1=cmpInfo.cur$group.inds$group1, inds.g2=cmpInfo.cur$group.inds$group2,
					covg=getCovg(dsn, asMatrix=TRUE), covg.thres=covgThres,
					paired=cmpInfo.cur$paired, adjustment.table=cmpInfo.cur$adjustment.table,
					...
			)
			diffmeth <- addDiffMethTable(diffmeth,dm,cmpInfo.cur$comparison,"sites",cmpInfo.cur$group.names)
		}
		cleanMem()
		if (length(cmpInfo.cur$region.types)>0){
			if (skipSites){
				dmr <- computeDiffAcc.rnb.nome.bin.region(dsn, NULL,
					cmpInfo.cur$group.inds$group1, cmpInfo.cur$group.inds$group2,
					regionTypes=cmpInfo.cur$region.types,
					covg.thres=covgThres,
					paired=cmpInfo.cur$paired, adjustment.table=cmpInfo.cur$adjustment.table,
					...
				)
			} else {
				dmr <- computeDiffAcc.rnb.nome.bin.region(dsn,dm,
					cmpInfo.cur$group.inds$group1,cmpInfo.cur$group.inds$group2,
					regionTypes=cmpInfo.cur$region.types
				)	
			}		
			for (rt in cmpInfo.cur$region.types){
				diffmeth <- addDiffMethTable(diffmeth,dmr[[rt]],cmpInfo.cur$comparison, 
					rt, cmpInfo.cur$group.names
				)
			}
		}
		logger.completed()
	}

	diffmeth <- addComparisonInfo(diffmeth,cmpInfo)
	logger.completed()
	return(diffmeth)
}

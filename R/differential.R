################################################################################
# Methods for differential chromatin accessibility analysis - General
################################################################################
#' getComparisonInfo
#'
#' retrieve the comparison information for an DsAcc object. Analogous to \code{RnBeads::get.comparison.info}
#' @param dsa         \code{\linkS4class{DsAcc}} object
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
getComparisonInfo <- function(dsa, cmpNames=NULL, regionTypes=getRegionTypes(dsa), allPairs=TRUE, adjPairCols=NULL, minGrpSize=1L, maxGrpCount=NULL){
	sannot <- getSampleAnnot(dsa)
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

#' fastqDirToTable
#' 
#' scan a directory containing fastq files and create a sample annotation table from the file names
#' @param fqDir   string specifying a directory with fastq files
#' @param tabFn   filename specifying the where the table should be written to. If NULL (default), the table will just be returned
#'                as data frame
#' @return data frame of parsed annotation
#' @author Fabian Mueller
#' @export
fastqDirToTable <- function(fqDir, tabFn=NULL){
	fqs <- list.files(fqDir, pattern=".fastq")
	baseNames <- gsub("_R[12].*\\.fastq(\\.gz)?$", "", fqs)
	baseNamesU <- unique(baseNames)
	readNum <- as.integer(gsub("^.*_R([12]).*\\.fastq(\\.gz)?$", "\\1", fqs))
	fq1names <- sapply(baseNamesU, FUN=function(bn){
		fqs[baseNames==bn & readNum==1][1]
	})
	fq2names <- sapply(baseNamesU, FUN=function(bn){
		fqs[baseNames==bn & readNum==2][1]
	})
	res <- data.frame(
		sampleId=baseNamesU,
		sampleId_org=baseNamesU,
		fileName_fastq1=fq1names,
		fileName_fastq2=fq2names,
		stringsAsFactors=FALSE
	)
	if (!is.null(tabFn)){
		write.table(res, tabFn, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
	}
	return(res)
}


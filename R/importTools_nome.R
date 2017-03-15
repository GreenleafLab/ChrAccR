#' mergeBisSNPMethCalls
#' 
#' Wrapper around command line tools to merge methylation data in bisSNP output
#' format
#' @param fns     a NAMED vector of input file names
#' @param outPrefix prefix for output files
#' @param wd      working directory
#' @param cleanup remove temporary files
#' @details
#' requires avalid installation of \code{bedtools}. Make sure it is included in your \code{PATH} variable
#' @return \code{invisible(NULL)}; The output files with extensions \code{.methlvl.bed} and \code{.covg.bed}
#'         will be generated according to the prefix. They contain genomic coordinates and methylation levels
#'         or read coverage
#' @author Fabian Mueller
#' @noRd
mergeBisSNPMethCalls <- function(fns, outPrefix, wd=tempdir(), cleanup=TRUE){
	require(tools)
	require(R.utils)

	sysWrap <- function(cmd, args){
		logger.info(c("sys:", cmd, args))
		system2(cmd, args, wait=TRUE, stdout="", stderr="")
	}

	# set working directory to the directory where the temp files will be stored
	# (makes things clearer)
	old.wd <- getwd()
	# print(old.wd)
	wd <- tempfile(pattern="mbmc", tmpdir=wd)
	dir.create(wd)
	setwd(wd)

	nns <- names(fns)
	if (is.null(nns)){
		logger.info("no names detcted --> using filenames")
		nns <- gsub(".bed$", "", gsub(".gz$", "", basename(fns[fns.isGz])))
	}
	# extract the input files
	fns.isGz <- file_ext(fns)=="gz"
	if (any(fns.isGz)){
		logger.status("Extracting input files...")
		fns.repl <- gsub(".gz$", "", basename(fns[fns.isGz]))
		for (i in 1:sum(fns.isGz)){
			# sysRes <- sysWrap("zcat", c(fns[fns.isGz][i], "|", "head", "-n", "1000", ">", fns.repl[i]))
			sysRes <- sysWrap("gunzip", c("-c", fns[fns.isGz][i], ">", fns.repl[i]))
		}
		fns[fns.isGz] <- fns.repl
	}

	i <- 0
	for (fn in fns){
		i <- i + 1
		curName <- nns[i]
		logger.start(c("Processing", curName))
			logger.status("Sorting...")
			cmd <- "awk"
			args <- c("'{if(NR>1)print}'", fn, "|", "sort", "-k", "1,1", "-k2,2n", ">", paste0(curName,".sorted.bed"))
			sysRes <- sysWrap(cmd, args)
			logger.status("Split by strand...")
			for (ss in c("+","-")){
				#create empty file in case one of the strands is not present in the input file
				file.create(paste0(curName, ".", ss, ".bed"))
			}
			cmd <- "awk"
			args <- c(paste0("'{print > ", '"', curName, '."$6".bed"}', "'"), paste0(curName,".sorted.bed"))
			sysRes <- sysWrap(cmd, args)
			logger.status("Extracting coordinates, methylation levels and coverage...")
			for (ss in c("+","-")){
				cmd <- "cut"
				args <- c("-f1-3,4", paste0(curName, ".", ss, ".bed"), ">", paste0(curName,".",ss,".methlvl.bg"))
				sysRes <- sysWrap(cmd, args)
				args <- c("-f1-3,5", paste0(curName, ".", ss, ".bed"), ">", paste0(curName,".",ss,".covg.bg"))
				sysRes <- sysWrap(cmd, args)
			}
		logger.completed()
	}
	logger.status("Merging methylation levels and coverage...")
	cmd <- "bedtools"
	for (ss in c("+","-")){
		#create empty files in case one of the strands is not present in all input files
		file.create(paste0("merged.",ss,".methlvl.bed"))
		file.create(paste0("merged.",ss,".covg.bed"))

		args <- c("unionbedg", "-filler", "NA", "-i", paste0(nns, ".",ss,".methlvl.bg"), ">", paste0("merged.",ss,".methlvl.bed"))
		sysRes <- sysWrap(cmd, args)
		args <- c("unionbedg", "-filler", "NA", "-i", paste0(nns, ".",ss,".covg.bg"), ">", paste0("merged.",ss,".covg.bed"))
		sysRes <- sysWrap(cmd, args)
	}
	logger.status("Merging strands...")
	cmd <- "awk"
	args <- c('-F"\\t"', '\'OFS="\\t"{$3=$3"\\t+"}1\'', "merged.+.methlvl.bed", ">", "merged.methlvl.unsorted.bed")
	sysRes <- sysWrap(cmd, args)
	args <- c('-F"\\t"', '\'OFS="\\t"{$3=$3"\\t-"}1\'', "merged.-.methlvl.bed", ">>", "merged.methlvl.unsorted.bed")
	sysRes <- sysWrap(cmd, args)
	args <- c('-F"\\t"', '\'OFS="\\t"{$3=$3"\\t+"}1\'', "merged.+.covg.bed", ">", "merged.covg.unsorted.bed")
	sysRes <- sysWrap(cmd, args)
	args <- c('-F"\\t"', '\'OFS="\\t"{$3=$3"\\t-"}1\'', "merged.-.covg.bed", ">>", "merged.covg.unsorted.bed")
	sysRes <- sysWrap(cmd, args)
	logger.status("Sorting to output files...")
	cmd <- "sort"
	args <- c("-k", "1,1", "-k2,2n", "merged.methlvl.unsorted.bed", ">", paste0(outPrefix, ".methlvl.bed"))
	sysRes <- sysWrap(cmd, args)
	args <- c("-k", "1,1", "-k2,2n", "merged.covg.unsorted.bed", ">", paste0(outPrefix, ".covg.bed"))
	sysRes <- sysWrap(cmd, args)

	#prepend header line #sed -i -e '1iHere is my new top line\' filename
	headLine <- paste(c("chrom", "start", "end", "strand", nns), collapse="\t")
	cmd <- "sed"
	args <- c("-i", "-e", paste0("'1i", headLine, "\\'"), paste0(outPrefix, ".methlvl.bed"))
	sysRes <- sysWrap(cmd, args)
	args <- c("-i", "-e", paste0("'1i", headLine, "\\'"), paste0(outPrefix, ".covg.bed"))
	sysRes <- sysWrap(cmd, args)
	
	if (cleanup){
		unlink(wd, recursive=TRUE)
	}
	setwd(old.wd)
	invisible(NULL)
}

#-------------------------------------------------------------------------------
#' DsNOMe.bisSNP
#' 
#' Create a DsNOMe dataset from multiple input files in bisSNP output format
#' @param inputFns     a NAMED vector of input file names
#' @param sampleAnnot  data.frame spcifying the sample annotation table
#' @param genome       genome assembly
#' @param sampleIds    character vector of sample names
#' @return \code{\linkS4class{Job}} object
#' @author Fabian Mueller
#' @export
DsNOMe.bisSNP <- function(inputFns, sampleAnnot, genome, sampleIds=rownames(sampleAnnot)){
	#TODO:convenience: check if files exists and have the correct format before importing
	if (!all(file.exists(inputFns))){
		logger.error("Not all input files exist")
	}
	if (is.null(names(inputFns))) logger.error("inputFns must be a NAMED vector of file names for input type 'bisSNPfiles'")
	if (is.null(sampleIds)){
		logger.warning("Sample IDs not specified --> Assuming annotation table is in the same order as input file vector")
		rownames(sampleAnnot) <- names(inputFns)
	}
	sampleIds.annot <- intersect(sampleIds, rownames(sampleAnnot))
	sampleIds.not.annot <- setdiff(sampleIds, rownames(sampleAnnot))
	if (length(sampleIds.not.annot)){
		logger.error(c("The following sample IDs have no annotation", paste(sampleIds.not.annot, collapse=",")))
	}
	sampleAnnot <- sampleAnnot[sampleIds,]
	rownames(sampleAnnot) <- sampleIds

	sampleIds.input <- intersect(sampleIds, names(inputFns))
	sampleIds.not.input <- setdiff(sampleIds, names(inputFns))
	if (length(sampleIds.not.input)){
		logger.error(c("The following sample IDs have no input data", paste(sampleIds.not.input, collapse=",")))
	}
	inputFns <- inputFns[sampleIds]


	logger.start("Merging bed files")
		#setConfigElement("tmpDir","/TL/deep/projects/nobackup/fmueller/tmp/readData_nome_20170313_155826_6e78e54c658")
		prefMerged <- file.path(getConfigElement("tmpDir"), "methCalls")
		mergeBisSNPMethCalls(fns, prefMerged, wd=getConfigElement("tmpDir"))
		# mergeBisSNPMethCalls(fns[c(4,11,20)], ffff, wd="/DEEP_fhgfs/projects/fmueller/tmp/blubb", cleanup=FALSE)
		# mergeBisSNPMethCalls(fns, prefMerged, wd="/DEEP_fhgfs/projects/fmueller/tmp/blubb", cleanup=FALSE)
	logger.completed()
	
	logger.start("Reading methylation data")
		logger.status("Methylation levels...")
		# dtMeth <- read.table(paste0(prefMerged, ".methlvl.bed"), sep="\t", header=TRUE, check.names=FALSE)
		# dtMeth <- fread(paste0(prefMerged, ".methlvl.bed"), sep="\t", header=TRUE, nrows=1000)
		dtMeth <- fread(paste0(prefMerged, ".methlvl.bed"), sep="\t", header=TRUE)
		# format(object.size(dtMeth), units = "Gb")
		#adjust coordinates of - strand
		dtMeth[strand=="+", start:=start-1L]
		dtMeth[strand=="+", end:=end-1L]

		if (!all(sampleIds %in% colnames(dtMeth))) {
			logger.error(c("The following sample IDs are missing from the methylation data:", paste(setdiff(sampleIds, colnames(dtMeth)),collapse=",")))
		}
		#convert to GRanges
		df <- as.data.frame(dtMeth[, .(chrom, start, end, strand)])
		rownames(df) <- NULL
		gr <- df2granges(df, chrom.col=1L, start.col=2L, end.col=3L, strand.col=4L, coord.format="B0RI", assembly=genome, doSort=FALSE, adjNumChromNames=grepl("_chr$",genome))
		cleanMem() #clean-up
		
		mm <- dtMeth[,sampleIds, with=FALSE]/100 #matrix of methylation levels
		rm(dtMeth, df); cleanMem() #clean-up
		if (nrow(mm)!=length(gr)) logger.error("Not all methylation levels match a genomic coordinate. Maybe the input file contains invalid chromosome names?")

		logger.status("Coverage...")
		# dtCovg <- read.table(paste0(prefMerged, ".covg.bed"), sep="\t", header=TRUE, check.names=FALSE)
		dtCovg <- fread(paste0(prefMerged, ".covg.bed"), sep="\t", header=TRUE)
		if (!all(sampleIds %in% colnames(dtCovg))) {
			logger.error(c("The following sample IDs are missing from the methylation data:", paste(setdiff(sampleIds, colnames(dtCovg)),collapse=",")))
		}
		cm <- dtCovg[,sampleIds, with=FALSE] #matrix of coverage values
		rm(dtCovg); cleanMem() #clean-up
	logger.completed()

	logger.start("Creating DsNOMe object")
		if (ncol(mm)!=ncol(cm)) logger.error("Methylation and coverage matrices have different numbers of columns")
		if (!all(colnames(mm)==colnames(cm))) logger.error("Methylation and coverage matrices have incompatible column names")
		if (nrow(mm)!=nrow(cm)) logger.error("Methylation and coverage matrices have different numbers of rows")

		obj <- DsNOMe(gr, mm, cm, sampleAnnot, genome)
		rm(gr, mm, cm); cleanMem() #clean-up
	logger.completed()
	return(obj)
}

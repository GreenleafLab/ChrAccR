% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DsATACsc-class.R
\docType{methods}
\name{mergePseudoBulk,DsATACsc-method}
\alias{mergePseudoBulk,DsATACsc-method}
\alias{mergePseudoBulk}
\title{mergePseudoBulk-methods}
\usage{
\S4method{mergePseudoBulk}{DsATACsc}(.object, mergeGroups, cleanSampleAnnot = TRUE)
}
\arguments{
\item{.object}{\code{\linkS4class{DsATACsc}} object}

\item{mergeGroups}{factor or character vector or column name in sample annotation table.
Can alternatively be a (named) list containing sample indices or names
for each group to merge.}

\item{cleanSampleAnnot}{clean up sample annotation table in the new object}
}
\value{
a new \code{\linkS4class{DsATAC}} object with cells merged into pseudobulk samples
}
\description{
Merge cells into pseudobulk samples based on annotation
}
\author{
Fabian Mueller
}

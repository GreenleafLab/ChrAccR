% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_atac.R
\name{getATACfragments}
\alias{getATACfragments}
\title{getATACfragments}
\usage{
getATACfragments(ga, offsetTn = TRUE)
}
\arguments{
\item{ga}{\code{GAlignmentPairs} (or \code{GAlignments} for single-end sequencing) object}

\item{offsetTn}{apply offsets for Tn5 dimer cut site (+4 bp on genomic + strand; -4 bp on genomic - strand)}
}
\value{
\code{GRanges} object containing derived insertions. For paired-end data (recommended), the width of the resulting ranges corresponds to the insert size
        for single-end data, the width is set to 1bp
}
\description{
Given a \code{GAlignmentPairs} or \code{GAlignments} object, return a \code{GRanges} object containing the fragment (or insertion site for single-end data)
}
\author{
Fabian Mueller
}

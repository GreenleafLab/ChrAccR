% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DsATAC-class.R
\docType{methods}
\name{regionAggregation,DsATAC-method}
\alias{regionAggregation,DsATAC-method}
\title{regionAggregation-methods}
\usage{
\S4method{regionAggregation}{DsATAC}(
  .object,
  regGr,
  type,
  signal = NULL,
  aggrFun = "median",
  dropEmpty = TRUE,
  bySample = TRUE,
  chunkSize = 5000L
)
}
\arguments{
\item{.object}{\code{\linkS4class{DsATAC}} object}

\item{regGr}{\code{GRanges} object containing regions to summarize}

\item{type}{character string specifying a name for the region type}

\item{signal}{character string specifying a name for the region type for the signal to be aggregated
If it is \code{NULL} (default), the new region type will be initialized with NA values.
If it is \code{"insertions"} count data will be initialized from insertion sites (if 
fragment data is present in the object).}

\item{aggrFun}{aggregation function for signal counts. Will only be used if \code{signal!="insertions"}
Currently \code{sum}, \code{mean} and \code{median} (default) are supported.}

\item{dropEmpty}{discard all regions with no observed signal counts}

\item{bySample}{[only relevant if \code{signal=="insertions"}]. Process data sample-by-sample to save memory.}

\item{chunkSize}{[only relevant if \code{signal=="insertions" & !bySample}] number of samples to process per chunk (saves memory).
If \code{NULL} or larger than the number of samples, only one chunk will be processed.}
}
\value{
a new \code{\linkS4class{DsATAC}} object with aggregated signal counts per regions
}
\description{
Aggregate signal counts across a set of regions
}
\author{
Fabian Mueller
}

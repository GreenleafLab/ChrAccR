% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DsNOMe-class.R
\docType{methods}
\name{getMeth,DsNOMe-method}
\alias{getMeth,DsNOMe-method}
\alias{getMeth}
\title{getMeth-methods}
\usage{
\S4method{getMeth}{DsNOMe}(.object, type = "sites", asMatrix = FALSE)
}
\arguments{
\item{.object}{\code{\linkS4class{DsNOMe}} object}

\item{type}{character string specifying the rgion type or \code{"sites"} (default)}

\item{asMatrix}{return a matrix instead of a \code{data.table}}
}
\value{
\code{data.table} or \code{matrix} containing methylation levels for
        each site/region and sample
}
\description{
Return table of methylation values
}
\author{
Fabian Mueller
}

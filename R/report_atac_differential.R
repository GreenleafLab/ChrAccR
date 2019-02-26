if (!isGeneric("createReport_differential")) {
	setGeneric(
		"createReport_differential",
		function(.object, ...) standardGeneric("createReport_differential"),
		signature=c(".object")
	)
}
#' createReport_differential-methods
#'
#' Create a report summarizing differential accessibility analysis
#'
#' @param .object    \code{\linkS4class{DsATAC}} object
#' @param reportDir  directory in which the report will be created
#' @return (invisible) \code{muReportR::Report} object containing the report
#' 
#' @rdname createReport_differential-DsATAC-method
#' @docType methods
#' @aliases createReport_differential
#' @aliases createReport_differential,DsATAC-method
#' @author Fabian Mueller
#' @export
setMethod("createReport_differential",
	signature(
		.object="DsATAC"
	),
	function(
		.object,
		reportDir
	) {
		require(muReportR)
		initConfigDir <- !dir.exists(file.path(reportDir, "_config"))
		rr <- createReport(file.path(reportDir, "differential.html"), "Differential Accessibility", page.title = "Differential", init.configuration=initConfigDir)
		rDir.data <- getReportDir(rr, dir="data", absolute=FALSE)
		rDir.data.abs <- getReportDir(rr, dir="data", absolute=TRUE)

		hasFragments <- length(.object@fragments) > 0

		sampleIds <- getSamples(.object)
		sannot <- getSampleAnnot(.object)
		sampleGrps <- getGroupsFromTable(sannot, cols=getConfigElement("differentialColumns"))

		txt <- c(
			"Differential Accessibility was quantified for the following comparisons:"
		)
		rr <- addReportSection(rr, "Comparisons", txt, level=1L, collapsed=FALSE)


		off(rr)
		invisible(rr)
	}
)


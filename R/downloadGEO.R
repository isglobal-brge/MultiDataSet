#' Creation of a MultiDataSet with content of GEO GSE tags
#'
#' Given a vector of GSE tags with a single datasets. The function downloads
#' the dataset. It converts the downloaded \code{ExpressionSet} to the
#' correct set (\code{ExpressionSet} for expression data, \code{SnpSet} for
#' SNP data and \code{MethylationSet} for methylation data) and it
#' creates a MultiDataSet with them.
#'
#' @param gse_list String vector with the GSE tags.
#' @param gse_type String vector with the GES data types (allowed:
#' code{"expression"}, \code{"snp"} and \code{"methylation"}).
#' @examples
#' \dontrun{
#' downloadGEO(
#' 	gse_list = c("GSE40576", "GSE40732"), 
#' 	gse_type = c("methylation", "expression")
#' )
#' }
#' @export
downloadGEO <- function(gse_list, gse_type) {
	if(! "GEOquery" %in% utils::installed.packages()[,"Package"]) {
		stop("Required package 'GEOquery' not installed.")
	} else {
		message("Required package 'GEOquery' will be loaded.")
	}
	## we should check that the gse_list only contains gse tags
	gse_sets <- lapply(gse_list, function(gse) {
	    GEOquery::getGEO(gse, GSEMatrix = TRUE)[[1]]
	})
	md <- createMultiDataSet()
	rr <- lapply(1:length(gse_sets), function(ii) {
		if(gse_type[[ii]] == "expression") {
			md <- add_genexp(md, gse_sets[[ii]])
		} else if(gse_type[[ii]] == "methylation") {
			meth <- methylationSet(
				exprs(gse_sets[[ii]]),
				phenoData(gse_sets[[ii]]),
				featureData(gse_sets[[ii]]),
				annotation(gse_sets[[ii]])
			)
			md <- add_methy(md, meth)
		} else if(gse_type[[ii]] == "snp") {
			snps <- new('SnpSet', 
				phenoData = phenoData(gse_sets[[ii]]),
				featureData = featureData(gse_sets[[ii]]),
				annotation = annotation(gse_sets[[ii]]),
				call = exprs(gse_sets[[ii]])
			)
			md <- add_methy(md, snps)
		} else {
			stop("Invalid data type at location [", ii, "] for '", gse_type[ii], "'.")
		}
	})
	rm(rr)
	return(md)
}

# This SuperSeries is composed of the following SubSeries:
# GSE40576	DNA Methylation Changes and Childhood Asthma in the Inner City [methylation]
# GSE40732	DNA Methylation Changes and Childhood Asthma in the Inner City [gene expression]


# gse.m <- GEOquery::getGEO("GSE40576", GSEMatrix = TRUE)
# gse.g <- GEOquery::getGEO("GSE40732", GSEMatrix = TRUE)
# gse.s <- GEOquery::getGEO("GSE53261", GSEMatrix = TRUE)

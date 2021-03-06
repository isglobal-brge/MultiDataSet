#' MultiDataSet instances
#'
#' The class \code{MultiDataSet} is a superior class to store multiple
#' datasets in form of triplets (assayData-phenoData-featureData). The datasets can be \code{eSet} or
#' \code{SummarizedExperiment} derived or \code{matrices}.
#'
#' The names of the three lists (\code{assayData}, \code{phenoData} and
#' \code{featureData}) must be the same.
#'
#' @name MultiDataSet-class
#' @rdname MultiDataSet-class
#' @exportClass MultiDataSet
#' @slot assayData List of \code{assayData} elements.
#' @slot phenoData List of \code{AnnotatedDataFrame} containing the phenoData
#' of each \code{assayData}.
#' @slot featureData List of \code{AnnotatedDataFrame} containing the featureData
#' of each \code{assayData}.
#' @slot rowRanges List of \code{GenomicRanges} containing the rowRanges
#' of each \code{assayData}.
#' @slot extraData List of other slots of the original object.
#' @slot return_method List of functions used to create the original object.
#' @seealso \code{\link{add_eset}}, \code{\link{add_rse}}
setClass(
  Class = "MultiDataSet",
  representation = representation(
    assayData = "list",
    phenoData = "list",
    featureData = "list",
    rowRanges = "list",
    extraData = "list",
    return_method = "list"
  )
)

#' Method to add an \code{eSet} to \code{MultiDataSet}.
#'
#' This method adds or overwrites a slot of a \code{MultiDataSet} with the content 
#' of the given \code{eSet}.
#'
#' @rdname add_eset
#' @aliases add_eset
#' @param object \code{MultiDataSet} that will be filled.
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.type Character with the type of data of the omic set (e.g. expression, methylation...)
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' are several sets of the same type (e.g. multiple expression assays)
#' @param sample.tables Character with the names of the slots with sample data besides phenoData.
#' @param feature.tables Character with the names of the slots with feature data besides featureData.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten. 
#' @param GRanges \code{GenomicRanges} to be included in rowRanges slot. 
#' @return A new \code{MultiDataSet} with a slot filled.
#' @seealso \code{\link{add_methy}}, \code{\link{add_genexp}}, \code{\link{add_rnaseq}}, \code{\link{add_snps}}
#' @exportMethod add_eset
#' @examples 
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(10), 5))
#' multi <- add_eset(multi, eset, "exampledata", GRanges = NA)
setGeneric("add_eset", function(object, set, dataset.type, dataset.name = NULL, 
                                sample.tables = NULL, feature.tables = NULL,
                                warnings = TRUE, overwrite = FALSE, GRanges) standardGeneric("add_eset"))

#' Method to add an expression microarray dataset to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"expression"} of an
#' \code{MultiDataSet} with the content of the given \code{ExpressionSet}. The fData of 
#' the \code{ExpressionSet} must contain the columns chromosome, start and end.
#' 
#' @aliases add_genexp
#' @rdname add_genexp
#' @exportMethod add_genexp
#' @param object \code{MultiDataSet} that will be filled.
#' @param gexpSet \code{ExpressionSet} to be used to fill the slot.
#' @param ... Arguments to be passed to \code{add_eset}.
#' @return A new \code{MultiDataSet} with the slot \code{"expression"}
#' filled.
#' @examples 
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(4), 2))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), start = c(12414, 1234321),
#'  end = c(121241, 124124114), stringsAsFactors = FALSE)
#' multi <- add_genexp(multi, eset)
setGeneric("add_genexp", function(object, gexpSet, ...) standardGeneric("add_genexp"))

#' Method to add an expression RNA seq dataset to \code{MultiDataSet}.
#' 
#' This method adds or overwrites the slot \code{"rnaseq"} of an
#' \code{MultiDataSet} with the content of the given \code{ExpressionSet}. The fData of 
#' the \code{ExpressionSet} must contain the columns chromosome, start and end.
#' 
#' @aliases add_rnaseq
#' @rdname add_rnaseq-methods
#' @exportMethod add_rnaseq
#' @param object \code{MultiDataSet} that will be filled.
#' @param rnaSet \code{ExpressionSet} to be used to fill the slot.
#' @param ... Arguments to be passed to \code{add_eset}.
#' @return A new \code{MultiDataSet} with the slot \code{"rnaseq"}
#' filled.
#' @examples 
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(4), 2))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), start = c(12414, 1234321),
#'  end = c(121241, 12122414), stringsAsFactors = FALSE)
#' multi <- add_genexp(multi, eset)
setGeneric("add_rnaseq", function(object, rnaSet, ...) standardGeneric("add_rnaseq"))

#' Method to add a slot of methylation to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"methylation"} of an
#' \code{MultiDataSet} with the content of the given \code{MethylationSet} or \code{RatioSet}. 
#' The fData of the input object must contain the columns chromosome and position.
#'
#' @rdname add_methy
#' @aliases add_methy
#' @param object \code{MultiDataSet} that will be filled.
#' @param methySet \code{MethylationSet} or \code{RatioSet} to be used to fill the slot.
#' @param ... Further arguments to be passed to \code{add_eset}.
#' @return A new \code{MultiDataSet} with the slot \code{"methylation"}
#' filled.
#' @exportMethod add_methy
#' @examples 
#' if (require(brgedata)){
#'  multi <- createMultiDataSet()
#'  multi <- add_methy(multi, brge_methy[1:100, ])
#' }
setGeneric("add_methy", function(object, methySet, ...) standardGeneric("add_methy"))

#' Method to add a \code{RangedSummarizedExperiment} to \code{MultiDataSet}.
#'
#' This method adds or overwrites a slot of a \code{MultiDataSet} with the content 
#' of the given  \code{RangedSummarizedExperiment}.
#'
#' @rdname add_rse
#' @aliases add_rse
#' @param object \code{MultiDataSet} that will be filled.
#' @param set Object derived from \code{RangedSummarizedExperiment} to be used to fill the slot.
#' @param dataset.type Character with the type of data of the omic set (e.g. expression, methylation...)
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' are several sets of the same type (e.g. multiple expression assays)
#' @param sample.tables Character with the names of the slots with sample data besides colData.
#' @param feature.tables Character with the names of the slots with feature data besides rowData.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten. 
#' @return A new \code{MultiDataSet} with a slot filled.
#' @exportMethod add_rse
#' @examples 
#' if (require(GenomicRanges) & require(SummarizedExperiment)){
#' multi <- createMultiDataSet()
#' counts <- matrix(runif(200 * 6, 1, 1e4), 200)
#' rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
#'                      IRanges(floor(runif(200, 1e5, 1e6)), width=100),
#'                      strand=sample(c("+", "-"), 200, TRUE),
#'                      feature_id=sprintf("ID%03d", 1:200))
#' colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
#'                     row.names=LETTERS[1:6], id = LETTERS[1:6])
#' names(rowRanges) <- 1:200
#' rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
#'                             rowRanges=rowRanges, colData=colData)
#' multi <- add_rse(multi, rse, "rseEx")
#' }
setGeneric("add_rse", function(object, set, dataset.type, dataset.name = NULL, 
                               sample.tables = NULL, feature.tables = NULL,
                               warnings = TRUE, overwrite = FALSE) standardGeneric("add_rse"))


#' Method to add a \code{SummarizedExperiment} to \code{MultiDataSet}.
#'
#' This method adds or overwrites a slot of a \code{MultiDataSet} with the content 
#' of the given \code{SummarizedExperiment}.
#'
#' @rdname add_se
#' @aliases add_se
#' @param object \code{MultiDataSet} that will be filled.
#' @param set Object derived from \code{SummarizedExperiment} to be used to fill the slot.
#' @param dataset.type Character with the type of data of the omic set (e.g. expression, methylation...)
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' are several sets of the same type (e.g. multiple expression assays)
#' @param sample.tables Character with the names of the slots with sample data besides colData.
#' @param feature.tables Character with the names of the slots with feature data besides rowData.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten. 
#' @param GRanges \code{GenomicRanges} to be included in rowRanges slot. 
#' @return A new \code{MultiDataSet} with a slot filled.
#' @exportMethod add_se
#' @examples 
#' multi <- createMultiDataSet()
#' se <- SummarizedExperiment::SummarizedExperiment(matrix(runif(10), 5))
#' multi <- add_se(multi, se, "exampledata", GRanges = NA)
setGeneric("add_se", function(object, set, dataset.type, dataset.name = NULL, 
                              sample.tables = NULL, feature.tables = NULL,
                              warnings = TRUE, overwrite = FALSE, GRanges) standardGeneric("add_se"))

#' Method to add a slot of SNPs to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"snps"} of an
#' \code{MultiDataSet} with the content of the given \code{SnpSet}. The fData of the
#' \code{SnpSet} must contain the columns chromosome and position.
#'
#' @rdname add_snps
#' @aliases add_snps
#' @param object \code{MultiDataSet} that will be filled.
#' @param snpSet \code{SnpSet} to be used to fill the slot.
#' @param ... Arguments to be passed to \code{add_eset}.
#' @return A new \code{MultiDataSet} with the slot \code{"snps"}
#' filled.
#' @exportMethod add_snps
#' @examples 
#' multi <- createMultiDataSet()
#' geno <- matrix(c(3,1,2,1), ncol = 2)
#' colnames(geno) <- c("VAL0156", "VAL0372")
#' rownames(geno) <- c("rs3115860", "SNP1-1628854")
#' map <- AnnotatedDataFrame(data.frame(chromosome = c("chr1", "chr2"), position = c(12414, 1234321),
#'      stringsAsFactors = FALSE))
#' rownames(map) <- rownames(geno)
#' snpSet <- new("SnpSet", call = geno, featureData = map)
#' pheno <- data.frame(id = c("VAL0156", "VAL0372"))
#' rownames(pheno) <- c("VAL0156", "VAL0372")
#' pData(snpSet) <- pheno
#' multi <- add_snps(multi, snpSet)
setGeneric("add_snps", function(object, snpSet, ...) standardGeneric("add_snps"))

#' Method to add a \code{matrix} to \code{MultiDataSet}.
#'
#' This method adds or overwrites a slot of a \code{MultiDataSet} with the content 
#' of the given \code{matrix}.
#'
#' @rdname add_table
#' @aliases add_table
#' @param object \code{MultiDataSet} that will be filled.
#' @param set \code{matrix} used to fill the slot.
#' @param dataset.type Character with the type of data 
#' @param dataset.name Character with the specific name for this set (NULL by 
#' default). It is useful when there are several sets of the same type.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten. 
#' @return A new \code{MultiDataSet} with a slot filled.
#' @exportMethod add_table
#' @examples 
#' multi <- createMultiDataSet()
#' mat <- matrix(runif(12), nrow = 3)
#' colnames(mat) <- paste0("S", 1:4)
#' rownames(mat) <- paste0("F", 1:3)
#' multi <- add_table(multi, mat, "exampledata")
setGeneric("add_table", function(object, set, dataset.type, dataset.name = NULL,
                                 warnings = TRUE, overwrite = FALSE) standardGeneric("add_table"))


#' Get the name of the ids common to all datasets
#'
#' @rdname commonIds
#' @aliases commonIds
#' @param object \code{MultiDataSet} that will be filtered.
#' @return Character vector with the common ids. 
#' @exportMethod commonIds
#' @examples  
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(9), ncol = 3))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1"), 
#'                           start = c(1, 5, 10),end = c(4, 6, 14), 
#'                           stringsAsFactors = FALSE)
#' sampleNames(eset) <- c("S1", "S2", "S3")
#' pData(eset) <- data.frame(id = c("S1", "S2", "S3"))
#' rownames(pData(eset)) <- c("S1", "S2", "S3")
#' multi <- add_genexp(multi, eset, dataset.name = "g1")
#' eset <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
#'                           start = c(1, 14, 25, 104),end = c(11, 16, 28, 115),
#'                           stringsAsFactors = FALSE)
#' sampleNames(eset) <- c("S1", "G2")
#' pData(eset) <- data.frame(id = c("S1", "G2"))
#' rownames(pData(eset)) <- c("S1", "G2")
#' 
#' multi <- add_genexp(multi, eset, dataset.name="g2")
#' commonIds(multi)
setGeneric("commonIds", function(object) standardGeneric("commonIds"))

#' Method to select samples that are present in all datasets.
#'
#' This method subsets the datasets to only contain the samples that are in all datasets. All sets 
#' will have the samples in the same order, taking into account that there can be duplicates.
#' 
#' If unify.names is TRUE, the sample names of the sets will be unified using the id column of 
#' phenodata. This option is only possible when there are no duplicated ids. 
#'
#' @rdname commonSamples
#' @aliases commonSamples
#' @param object \code{MultiDataSet} that will be filtered.
#' @param unify.names Logical indicating if sample names of the sets should be unified.
#' @return A new \code{MultiDataSet} with only the common samples.
#' @exportMethod commonSamples
#' @examples  
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(9), ncol = 3))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1"), 
#'                           start = c(1, 5, 10),end = c(4, 6, 14), 
#'                           stringsAsFactors = FALSE)
#' sampleNames(eset) <- c("S1", "S2", "S3")
#' pData(eset) <- data.frame(id = c("S1", "S2", "S3"))
#' rownames(pData(eset)) <- c("S1", "S2", "S3")
#' multi <- add_genexp(multi, eset, dataset.name = "g1")
#' eset <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
#'                           start = c(1, 14, 25, 104),end = c(11, 16, 28, 115),
#'                           stringsAsFactors = FALSE)
#' sampleNames(eset) <- c("S1", "G2")
#' pData(eset) <- data.frame(id = c("S1", "G2"))
#' rownames(pData(eset)) <- c("S1", "G2")
#' 
#' multi <- add_genexp(multi, eset, dataset.name="g2")
#' commonSamples(multi)
setGeneric("commonSamples", function(object, unify.names = FALSE) standardGeneric("commonSamples"))

#' Get the name of the datasets that have rowRanges
#' @export
#' 
#' @param object \code{MultiDataSet}
#' @return Character vector with the slots that have rowRanges.
#' @examples 
#' multi <- createMultiDataSet()
#' eset <- new("ExpressionSet", exprs = matrix(runif(10), 5))
#' eset2 <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
#' fData(eset2) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
#'                           start = c(1, 14, 25, 104),end = c(11, 16, 28, 115),
#'                           stringsAsFactors = FALSE)
#' multi <- add_eset(multi, eset, "exampledata", GRanges = NA)
#' multi <- add_genexp(multi, eset2)
#' rowRangesElements(multi)
setGeneric("rowRangesElements", function(object) standardGeneric("rowRangesElements"))

#' Apply iClusterPlus clustering method to a MultiDataSet object
#' 
#' Method \link[iClusterPlus]{iClusterPlus} is applied on a \link{MultiDataSet} object after
#' getting the common samples along all the contained datasets.
#'
#' @rdname w_iclusterplus
#' @aliases w_iclusterplus
#' @param object \code{MultiDataSet}
#' @param commonSamples Logical to indicate if common samples are selected
#' @param ... Arguments passed to function \link[iClusterPlus]{iClusterPlus}
#' @note Argument \code{type} for \link[iClusterPlus]{iClusterPlus} is filled within the method.
#' @return A list of results from \link[iClusterPlus]{iClusterPlus}
#' @export
setGeneric("w_iclusterplus", function(object, commonSamples=TRUE, ...)
    standardGeneric("w_iclusterplus")
)

#' Apply mcia integration method to a MultiDataSet object
#' 
#' Method \link[omicade4]{mcia} is applied on a \link{MultiDataSet} object after
#' getting the common samples along all the contained datasets.
#' 
#' @rdname w_mcia
#' @aliases w_mcia
#' @param object \code{MultiDataSet}
#' @param ... Arguments passed to function \link[omicade4]{mcia}
#' @return A list of results from \link[omicade4]{mcia}
#' @export
setGeneric("w_mcia", function(object, ...)
    standardGeneric("w_mcia")
)
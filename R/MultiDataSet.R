#' MultiDataSet: Implementation of the BRGE's basic classes
#' 
#'  Implementation of the BRGE's (Bioinformatic Research Group in Epidemiology from
#'  Center for Research in Environmental Epidemiology) MultiDataSet and MethylationSet. MultiDataSet
#'  is designed for integrating multi omics data sets and MethylationSet to contain normalized methylation data. MultiDataSet for integrating multi omics data sets
#' 
#' 
#' @docType package
#' @name MultiDataSet
#' 
#' @import BiocGenerics
#' @import Biobase
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment 
#' @import methods
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps 
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom utils installed.packages
#' @importFrom minfi RatioSet GenomicRatioSet
#' @seealso \linkS4class{MultiDataSet}
NULL
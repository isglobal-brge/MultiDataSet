#' Convert a \code{MultiAssayExperiment} to a \code{MultiDataSet}
#' 
#' This function creates a \code{MultiDataSet} using the data of a \code{MultiAssayExperiment}. 
#' 
#' @export mae2mds
#' @param MAE a \code{MultiAssayExperiment} 
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return \code{MultiDataSet} with the of the incoming \code{MultiAssayExperiment}.
mae2mds <- function(MAE, warnings = TRUE){
    
    MDS <- createMultiDataSet()
    
    exps <- MultiAssayExperiment::experiments(MAE)
    map <- MultiAssayExperiment::sampleMap(MAE)
    map <- MultiAssayExperiment::mapToList(map, colnames(map)[1])
    pheno <- as.data.frame(MultiAssayExperiment::pData(MAE))
    
    for (expname in names(exps)){
        set <- exps[[expname]]
        setmap <- map[[expname]]
        rownames(setmap) <- setmap$colname
        if (class(set) == "ExpressionSet"){
            pData(set) <- cbind(pData(set), pheno[setmap[rownames(pData(set)), "primary"], ])
            set$id <- setmap[rownames(pData(set)), "primary"]
            MDS <- add_eset(MDS, set, dataset.type = expname, GRanges = NA)
        } else if (class(set) == "RangedSummarizedExperiment"){
            SummarizedExperiment::colData(set) <- 
                cbind(SummarizedExperiment::colData(set), pheno[setmap[rownames(SummarizedExperiment::colData(set)), "primary"], ])
            set$id <- setmap[rownames(SummarizedExperiment::colData(set)), "primary"]
            MDS <- add_rse(MDS, set, dataset.type = expname)
        } else if (class(set) == "SummarizedExperiment"){
            SummarizedExperiment::colData(set) <- 
                cbind(SummarizedExperiment::colData(set), pheno[setmap[rownames(SummarizedExperiment::colData(set)), "primary"], ])
            set$id <- setmap[rownames(SummarizedExperiment::colData(set)), "primary"]
            MDS <- add_se(MDS, set, dataset.type = expname, GRanges = NA)
        }else if (class(set) == "matrix"){
            pset <- pheno[setmap[colnames(set), "primary"], ]
            pset$id <- setmap[colnames(set), "primary"]
            rownames(pset) <- colnames(set)
            eSet <- ExpressionSet(set)
            pData(eSet) <- pset
            MDS <- add_eset(MDS, eSet, dataset.type = expname, GRanges = NA)
        } else{
            if (warnings){
                warning("MultiDataSet does not support ", class(set), ". ", expname, " will not be added.")
            }
        }
    }
    MDS
}

#' Convert a \code{MultiDataSet} to a \code{MultiAssayExperiment}
#' 
#' This function creates a \code{MultiAssayExperiment} using the data of a \code{MultiDataSet}. 
#' 
#' @export mds2mae
#' @param MDS a \code{MultiDataSet} 
#' @return \code{MultiAssayExperiment} with the of the incoming \code{MultiDataSet}.
mds2mae <- function(MDS){
    
    objlist <- lapply(names(MDS), function(x) MDS[[x]])
    names(objlist) <- names(MDS)
    
    sampleMap <- lapply(pData(MDS), function(x) 
        data.frame(primary = x$id, assay = rownames(x), stringsAsFactors = FALSE))
    dfmap <- MultiAssayExperiment::listToMap(sampleMap)
    
    phenoDataVars <- Reduce(intersect, lapply(pData(MDS), colnames))
    
    pData <- Reduce(rbind, lapply(pData(MDS), function(x) x[, phenoDataVars, drop = FALSE]))
    pData <- pData[!duplicated(pData$id), ]
    rownames(pData) <- pData$id
    pData <- pData[, colnames(pData) != "id"]
    
    MAE <- MultiAssayExperiment::MultiAssayExperiment(objlist, pData, dfmap)
    MAE
}


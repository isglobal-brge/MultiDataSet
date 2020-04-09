context("check MultiDataSet functions")


# Setters Getters ####
test_that("Setters Getters", {
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
                              start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
                              status = c("case", "case", "control"), criteria = c("a", "b", "a"),
                              stringsAsFactors = FALSE)
    pData(eset) <- data.frame(id = letters[1:2])
    
    # Check rowRangesElements ####
    ## Empty dataset
    expect_equal(rowRangesElements(multi), character(0))
    
    ## One dataset with rowRanges
    multi2 <- add_eset(multi, eset, "rnaseq")
    expect_equal(rowRangesElements(multi2), "rnaseq")
    
    ## One dataset without rowRanges
    multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
    expect_equal(rowRangesElements(multi2), character(0))
    
    ## First dataset without rowRanges, second with
    multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
    multi2 <- add_eset(multi2, eset, "rnaseq2")
    expect_equal(rowRangesElements(multi2), "rnaseq2")
    
    ## First dataset with rowRanges, second without
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset, "rnaseq2", GRanges = NA)
    expect_equal(rowRangesElements(multi2), "rnaseq")
    
    ## Both datasets with rowRanges
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset, "rnaseq2")
    expect_equal(rowRangesElements(multi2), c("rnaseq", "rnaseq2"))
    
    
    # Check ncols ####
    ## Empty dataset
    expect_equal(ncol(multi), numeric(0))
    
    ## One dataset 
    multi2 <- add_eset(multi, eset, "rnaseq")
    expect_equal(ncol(multi2), c(rnaseq = 2))
    
    ## Two datasets
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset[, 1], "rnaseq2")
    expect_equal(ncol(multi2), c(rnaseq = 2, rnaseq2 = 1))
    
    # Check nrows ####
    ## Empty dataset
    expect_equal(nrow(multi), numeric(0))
    
    ## One dataset 
    multi2 <- add_eset(multi, eset, "rnaseq")
    expect_equal(nrow(multi2), c(rnaseq = 3))
    
    ## Two datasets
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset[1:2, ], "rnaseq2")
    expect_equal(nrow(multi2), c(rnaseq = 3, rnaseq2 = 2))
    
    ## Check as.list
    multi2 <- add_eset(multi, eset, "rnaseq")
    l <- as.list(multi2)
    expect_equal(names(l), "rnaseq")
    expect_equal(dim(l[[1]]), c(3, 2))
    
    ## Check pData/phenoData
    expect_equal(pData(multi2), list(rnaseq = pData(eset)))
    expect_equal(phenoData(multi2), 
                 list(rnaseq = list(main = phenoData(eset), protocolData = protocolData(eset))))
    
    ## Check fData/featureData
    expect_equal(fData(multi2), list(rnaseq = fData(eset)))
    expect_equal(featureData(multi2), 
                 list(rnaseq = list(main = featureData(eset))))
    
    ## Check sampleNames
    expect_equal(sampleNames(multi), NULL)
    expect_equal(sampleNames(multi2), list(rnaseq = c("a", "b")))
    
    ## Check dims
    expect_equal(dims(multi2), list(rnaseq = c(Features = 3, Samples = 2)))
})

# Show ####
test_that("Show", {
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
                              start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
                              status = c("case", "case", "control"), criteria = c("a", "b", "a"),
                              stringsAsFactors = FALSE)
    
    ## One dataset without rowRanges
    multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
    show(multi2)
    
    ## One dataset with rowRanges
    multi2 <- add_eset(multi, eset, "rnaseq")
    show(multi2)
    
    pData(eset) <- data.frame(id = letters[1:2], sex = c("H", "M"))
    multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
    show(multi2)
})

context("check MethylationSet functions")

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
    expect_equal(ncols(multi), numeric(0))
    
    ## One dataset 
    multi2 <- add_eset(multi, eset, "rnaseq")
    expect_equal(ncols(multi2), c(rnaseq = 2))
    
    ## Two datasets
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset[, 1], "rnaseq2")
    expect_equal(ncols(multi2), c(rnaseq = 2, rnaseq2 = 1))
    
    # Check nrows ####
    ## Empty dataset
    expect_equal(nrows(multi), numeric(0))
    
    ## One dataset 
    multi2 <- add_eset(multi, eset, "rnaseq")
    expect_equal(nrows(multi2), c(rnaseq = 3))
    
    ## Two datasets
    multi2 <- add_eset(multi, eset, "rnaseq")
    multi2 <- add_eset(multi2, eset[1:2, ], "rnaseq2")
    expect_equal(nrows(multi2), c(rnaseq = 3, rnaseq2 = 2))
    
    
})
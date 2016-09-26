context("Subsetting")
library(GenomicRanges)

test_that("Basic subseting", {
    ## Create MULTI
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(9), ncol = 3))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1"), 
                              start = c(1, 5, 10),end = c(4, 6, 14), 
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "S2", "S3")
    pData(eset) <- data.frame(id = c("S1", "S2", "S3"))
    rownames(pData(eset)) <- c("S1", "S2", "S3")
    
    multi <- add_genexp(multi, eset, dataset.name = "g1")
    
    
    eset <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
                              start = c(1, 14, 25, 104),end = c(11, 16, 28, 115),
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "G2")
    pData(eset) <- data.frame(id = c("S1", "G2"))
    rownames(pData(eset)) <- c("S1", "G2")
    
    multi <- add_genexp(multi, eset, dataset.name="g2")
    ## /
    
    ## Check NAMES
    expect_equal(names(multi), c("expression+g1", "expression+g2"))
    ## /
    
    ## Check SAMPLES
    expect_equal(commonIds(multi[c("S1", "S2"), ]), "S1")
    expect_equal(unique(unlist(sampleNames(commonSamples(multi)))), "S1")
    
    ## /
    
    ## Check RANGES
    expect_equivalent(nrow(multi[,, GRanges("chr1:1-7")][["expression+g1"]]), 2L)
    expect_equivalent(nrow(multi[,, GRanges("chr1:1-7")][["expression+g2"]]), 1L)
    ## /
    
})


test_that("Advanced subseting", {
    ## Create MULTI
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(9), ncol = 3))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1"), 
                              start = c(1, 5, 10),end = c(4, 6, 14), gene = letters[1:3],  
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "S2", "S3")
    pData(eset) <- data.frame(id = c("S1", "S2", "S3"), age = c(12, 20, 40))
    rownames(pData(eset)) <- c("S1", "S2", "S3")
    
    multi <- add_genexp(multi, eset, dataset.name = "g1")
    
    
    eset <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
                              start = c(1, 14, 25, 104), end = c(11, 16, 28, 115),
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "G2")
    pData(eset) <- data.frame(id = c("S1", "G2"))
    rownames(pData(eset)) <- c("S1", "G2")
    
    multi <- add_genexp(multi, eset, dataset.name = "g2")
    ## /
    
    ## Check Features
    multi2 <- subset(multi, gene == "a")
    expect_equal(names(multi2), c("expression+g1", "expression+g2"))
    expect_equivalent(nrow(multi2[["expression+g1"]]), 1)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)
    expect_equivalent(nrow(multi2[["expression+g2"]]), 4)
    expect_equivalent(ncol(multi2[["expression+g2"]]), 2)
    
    multi2 <- subset(multi, gene == "a", keep = FALSE)
    expect_equal(names(multi2), "expression+g1")
    expect_equivalent(nrow(multi2[["expression+g1"]]), 1)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)

    
    multi2 <- subset(multi, gene == "A")
    expect_equal(names(multi2), c("expression+g1", "expression+g2"))
    expect_equivalent(nrow(multi2[["expression+g1"]]), 0)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)
    expect_equivalent(nrow(multi2[["expression+g2"]]), 4)
    expect_equivalent(ncol(multi2[["expression+g2"]]), 2)
    
    multi2 <- subset(multi, gene == "A", keep = FALSE)
    expect_equal(names(multi2), "expression+g1")
    expect_equivalent(nrow(multi2[["expression+g1"]]), 0)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)

    expect_error(subset(multi, cot > 10), "feat expression could not be applied to any of the sets.")
    ## /
    
    ## Check PHENOTYPES
    multi2 <- subset(multi, , age > 10)
    expect_equal(names(multi2), c("expression+g1", "expression+g2"))
    expect_equivalent(nrow(multi2[["expression+g1"]]), 3)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)
    expect_equivalent(nrow(multi2[["expression+g2"]]), 4)
    expect_equivalent(ncol(multi2[["expression+g2"]]), 2)
    
    multi2 <- subset(multi, , age > 10 , keep = FALSE)
    expect_equal(names(multi2), "expression+g1")
    expect_equivalent(nrow(multi2[["expression+g1"]]), 3)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 3)
    
    
    multi2 <- subset(multi, , age < 10)
    expect_equal(names(multi2), c("expression+g1", "expression+g2"))
    expect_equivalent(nrow(multi2[["expression+g1"]]), 3)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 0)
    expect_equivalent(nrow(multi2[["expression+g2"]]), 4)
    expect_equivalent(ncol(multi2[["expression+g2"]]), 2)
    
    multi2 <- subset(multi, , age < 10, keep = FALSE)
    expect_equal(names(multi2), "expression+g1")
    expect_equivalent(nrow(multi2[["expression+g1"]]), 3)
    expect_equivalent(ncol(multi2[["expression+g1"]]), 0)
    
    expect_error(subset(multi, , cot > 10), "phe expression could not be applied to any of the sets.")
    ## /
    
})

test_that("wrong variables",{
    
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(9), ncol = 3))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1"), 
                              start = c(1, 5, 10),end = c(4, 6, 14), gene = letters[1:3],  
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "S2", "S3")
    pData(eset) <- data.frame(id = c("S1", "S2", "S3"), age = c(12, 20, 40))
    rownames(pData(eset)) <- c("S1", "S2", "S3")
    
    multi <- add_genexp(multi, eset, dataset.name = "g1")
    
    
    eset <- new("ExpressionSet", exprs = matrix(runif(8), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"), 
                              start = c(1, 14, 25, 104), end = c(11, 16, 28, 115),
                              stringsAsFactors = FALSE)
    sampleNames(eset) <- c("S1", "G2")
    pData(eset) <- data.frame(id = c("S1", "G2"), age = c(15, 30))
    rownames(pData(eset)) <- c("S1", "G2")
    
    multi <- add_genexp(multi, eset, dataset.name = "g2")
    expect_error(subset(multi, 1243), "'feat' must be a logical expression")
    expect_error(subset(multi, , 1243), "'phe' must be a logical expression")
    
})
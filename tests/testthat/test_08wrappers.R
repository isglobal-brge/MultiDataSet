context("check wrappers")

# Wrappers ####
test_that("Multiple omic analysis", {
    multi <- createMultiDataSet()
    
    set.seed(0)
    eset <- new("ExpressionSet", exprs = matrix(runif(20), ncol = 4))
    fData(eset) <- data.frame(chromosome = rep("chr1", 5), 
                              start = 1:5, end = 11:15, stringsAsFactors = FALSE)
    
    eset2 <- new("ExpressionSet", exprs = matrix(runif(20), ncol = 4))
    fData(eset2) <- fData(eset)
    multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
    multi2 <- add_eset(multi2, eset2, "rnaseq2", GRanges = NA)
    
    # iClusterPlus
    ic <- w_iclusterplus(multi2)
    
    expect_equal(class(ic), "list")
    expect_equal(names(ic), c("alpha", "beta", "clusters", "centers", "meanZ", "BIC", "dev.ratio", "dif"))

    
    multi2 <- add_eset(multi2, eset2, "rnaseq3", GRanges = NA)
    multi2 <- add_eset(multi2, eset2, "rnaseq4", GRanges = NA)
    multi2 <- add_eset(multi2, eset2, "rnaseq5", GRanges = NA)
    expect_error(w_iclusterplus(multi2))
    # MCIA
    mc <- w_mcia(multi2)
    
    expect_equal(class(mc), "mcia")
    expect_equal(names(mc), c("call", "mcoa", "coa"))

})

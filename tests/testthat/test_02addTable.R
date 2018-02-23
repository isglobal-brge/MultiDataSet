context("add eSet")

test_that("MultiSet", {
    multi <- createMultiDataSet()
    
    mat <- matrix(runif(12), nrow = 3)
    
    expect_error(add_table(multi, mat, "exampledata"), "Set must contain colnames.")
    
    colnames(mat) <- paste0("S", 1:4)
    expect_error(add_table(multi, mat, "exampledata"), "Set must contain rownames.")
    
    rownames(mat) <- paste0("F", 1:3)
    multi <- add_table(multi, mat, "exampledata")
    expect_equal(class(multi[["exampledata"]]), "matrix")
    
    ## Check overwrite
    expect_error(add_table(multi, mat, "exampledata", overwrite = FALSE))
    multi <- add_table(multi, mat, "exampledata", overwrite = TRUE)
    expect_equal(class(multi[["exampledata"]]), "matrix")
    
    
    multi2 <- multi["S1"]
    expect_equal(class(multi2[["exampledata"]]), "matrix")
    expect_equal(sampleNames(multi2)[["exampledata"]], "S1")
    
    
    mat2 <- matrix(runif(30), nrow = 6)
    colnames(mat2) <- paste0("S", c(4, 5, 6, 1, 2))
    rownames(mat2) <- paste0("F", 1:6)
    multi <- add_table(multi, mat2, "exampledata", "bis")
    multi2 <- commonSamples(multi)
    
    expect_equal(class(multi2[["exampledata"]]), "matrix")
    expect_equal(class(multi2[["exampledata+bis"]]), "matrix")
    
    expect_equal(colnames(multi2[["exampledata"]]), c("S1", "S2", "S4"))
    expect_equal(colnames(multi2[["exampledata+bis"]]), c("S1", "S2", "S4"))
    
    expect_equal(rownames(multi2[["exampledata"]]), c("F1", "F2", "F3"))
    expect_equal(rownames(multi2[["exampledata+bis"]]), c("F1", "F2", "F3", "F4", "F5", "F6"))
    
    ## Check duplicated colnames
    mat2 <- matrix(runif(12), nrow = 3)
    colnames(mat2) <- paste0("S", c(4, 4, 6, 2))
    rownames(mat2) <- paste0("F", 1:3)
    
    expect_error(add_table(multi, mat2, "cot"))
    
    colnames(mat2) <- paste0("S", c(4, 3, 6, 2))
    rownames(mat2) <- paste0("F", rep(1, 3))
    expect_error(add_table(multi, mat2, "cot"))
    
})
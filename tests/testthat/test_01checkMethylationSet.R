context("check MethylationSet functions")
library(S4Vectors)

test_that("MethylationSet creation", {
    beta_matrix <- matrix(runif(4), nrow = 2)
    colnames(beta_matrix) <- c("H", "M")
    rownames(beta_matrix) <- c("cg00050873", "cg00212031")
    phenotypes <- data.frame(age = c(12, 23))
    rownames(phenotypes) <- c("H", "M")
    annot <- data.frame(chr = c("chr1", "chr2"), pos = c(12490, 124109), stringsAsFactors = FALSE) 
    rownames(annot) <- c("cg00050873", "cg00212031")
    expect_match(class(methylationSet(beta_matrix, phenotypes = phenotypes,  
                              annotationDataFrame = annot)), "MethylationSet")
    
    # Check Constructor ####
    expect_match(class(methylationSet(beta_matrix, phenotypes = DataFrame(phenotypes),  
                                      annotationDataFrame = annot)), "MethylationSet")
    
    # Check Initializer
    expect_match(class(new("MethylationSet")), "MethylationSet")
    print(new("MethylationSet"))
})


test_that("Setters Getters", {
    beta_matrix <- matrix(runif(4), nrow = 2)
    colnames(beta_matrix) <- c("H", "M")
    rownames(beta_matrix) <- c("cg00050873", "cg00212031")
    phenotypes <- data.frame(age = c(12, 23))
    rownames(phenotypes) <- c("H", "M")
    annot <- data.frame(chr = c("chr1", "chr2"), pos = c(12490, 124109), stringsAsFactors = FALSE) 
    rownames(annot) <- c("cg00050873", "cg00212031")
    methSet <- methylationSet(beta_matrix, phenotypes = phenotypes,  
                              annotationDataFrame = annot)
    expect_equal(betas(methSet), beta_matrix)
    expect_equal(getMs(methSet), minfi::logit2(beta_matrix))
    
})


test_that("MethylationSet creation errors", {
    beta_matrix <- matrix(runif(4), nrow = 2)
    colnames(beta_matrix) <- c("H", "M")
    rownames(beta_matrix) <- c("cg00050873", "cg00212031")
    phenotypes <- data.frame(age = c(12, 23))
    rownames(phenotypes) <- c("H", "M")
    annot <- data.frame(chr = c("chr1", "chr2"), pos = c(12490, 124109), stringsAsFactors = FALSE) 
    rownames(annot) <- c("cg00050873", "cg00212031")
  
    # Phenotypes
    expect_error(methylationSet(beta_matrix, phenotypes = phenotypes[1, , drop = FALSE],  
                                      annotationDataFrame = annot), "phenotypes must have as rows as columns in betas.")
    
    expect_error(methylationSet(beta_matrix, phenotypes = phenotypes[1, , drop = TRUE],  
                                annotationDataFrame = annot), "phenotypes must be a data.frame or an AnnotatedDataFrame.")
    
    # Feature data
    expect_error(methylationSet(beta_matrix, phenotypes = phenotypes,  
                                annotationDataFrame = annot[1, , drop = FALSE]), "annotationDataFrame must have the same rows than betas.")
    
    expect_error(methylationSet(beta_matrix, phenotypes = phenotypes,  
                                annotationDataFrame = annot[1, , drop = TRUE]), "annotationDataFrame must be a data.frame or an AnnotatedDataFrame.")
    
})

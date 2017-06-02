context("add eSet")
library(GenomicRanges)

# Add eset ####
test_that("Add eSet to MultiSet", {
  multi <- createMultiDataSet()
  eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
  fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
    start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
    status = c("case", "case", "control"), criteria = c("a", "b", "a"),
    stringsAsFactors = FALSE)
  pData(eset) <- data.frame(id = letters[1:2])

  multi2 <- add_eset(multi, eset, "rnaseq")
  expect_s4_class(multi2, "MultiDataSet")
  expect_equal(names(multi2), "rnaseq")
  
  ## User GRanges 
  gr <- makeGRangesFromDataFrame(fData(eset))
  multi2 <- add_eset(multi, eset, "rnaseq", GRanges = gr)
  expect_s4_class(multi2, "MultiDataSet")
  expect_equal(names(multi2), "rnaseq")
  
  ## NA GRanges 
  multi2 <- add_eset(multi, eset, "rnaseq", GRanges = NA)
  expect_match(class(multi2), "MultiDataSet")
  expect_equal(names(multi2), "rnaseq")
  
  ## Check pData
  eset2 <- eset
  pData(eset2) <- data.frame(name = letters[1:2])
  expect_warning(multi2 <- add_eset(multi, eset2, "rnaseq", GRanges = NA), "No id column found in pData. The id will be equal to the sampleNames")
  expect_equal(pData(multi2[["rnaseq"]])$id, c("1", "2"))
})  
  
# Add eSet wrong ####
test_that("add_eset wrong variables",{
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
                              start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
                              status = c("case", "case", "control"), criteria = c("a", "b", "a"),
                              stringsAsFactors = FALSE)
    pData(eset) <- data.frame(id = letters[1:2])
    
    
    ## NA GRanges 
    expect_error(add_eset(multi, eset, "rnaseq", GRanges = "cot"), "GRanges should be a GenomicRanges or NA.")
})

# Specific functions ####
test_that("specific functions based on eSet", {
    
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
                              start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
                              status = c("case", "case", "control"), criteria = c("a", "b", "a"),
                              stringsAsFactors = FALSE)
    pData(eset) <- data.frame(id = letters[1:2])
    multi <- createMultiDataSet()
    
    # Check add_genexp ####
    multi <- add_genexp(multi, eset)
    expect_equal(names(multi), "expression")
    
    # Check add_rnaseq  ####
    multi <- add_rnaseq(multi, eset)
    expect_equal(names(multi), c("expression", "rnaseq"))
  
  
    expect_error(add_genexp(multi, eset), "There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
    expect_warning(add_genexp(multi, eset, overwrite = TRUE), "Slot 'expression' is already set in 'MultiDataSet'. Previous content will be overwritten.")
  
    # Check add_snp  ####
    geno <- matrix(c(3,1,2,1), ncol = 2)
    colnames(geno) <- c("VAL0156", "VAL0372")
    rownames(geno) <- c("rs3115860", "SNP1-1628854")
    map <- AnnotatedDataFrame(data.frame(chromosome = c("chr1", "chr2"), position = c(12414, 1234321),
                                         stringsAsFactors = FALSE))
    rownames(map) <- rownames(geno)
    snpSet <- new("SnpSet", call = geno, featureData = map)
    
    pheno <- data.frame(id = c("VAL0156", "VAL0372"))
    rownames(pheno) <- c("VAL0156", "VAL0372")
    pData(snpSet) <- pheno
    multi2 <- add_snps(multi, snpSet)
    expect_equal(names(multi2), c("expression", "rnaseq", "snps"))
    
    expect_s4_class(multi2[["expression"]], "ExpressionSet")
    expect_s4_class(multi2[["snps"]], "SnpSet")
    
    # Wrong fData ####
    snpSet2 <- snpSet
    fvarLabels(snpSet2)[1] <- "chr"
    expect_error(add_snps(multi, snpSet2), "fData of snpSet must contain columns chromosome and position")
    
    snpSet2 <- snpSet
    fvarLabels(snpSet2)[2] <- "cot"
    expect_error(add_snps(multi, snpSet2), "fData of snpSet must contain columns chromosome and position")
    
    
  # library(methylumi)
  # samps <- read.table(system.file("extdata/samples.txt",
  #                                 package = "methylumi"),sep="\t",header=TRUE)
  # mldat <- methylumi::methylumiR(system.file('extdata/exampledata.samples.txt',package='methylumi'),
  #                     qcfile=system.file('extdata/exampledata.controls.txt',package="methylumi"),
  #                     sampleDescriptions=samps)
  # fvarLabels(mldat) <- tolower(fvarLabels(mldat))
  # fData(mldat)$position <- fData(mldat)$cpg_coordinate
  # pData(mldat)$id <- pData(mldat)$sampleID
  # multi <- add_eset(multi, mldat, "cot", GRanges = NA)
  # expect_equal(names(multi), c("expression", "snps", "cot"))
  # expect_equal(names(multi[, c("cot", "snps")]), c("cot", "snps"))
  # expect_is(multi[, "snps", drop = TRUE], "SnpSet")
  # 
  
  
  beta_matrix <- matrix(runif(4), nrow = 2)
  colnames(beta_matrix) <- c("H", "M")
  rownames(beta_matrix) <- c("cg00050873", "cg00212031")
  phenotypes <- data.frame(age = c(12, 23))
  rownames(phenotypes) <- c("H", "M")
  annot <- data.frame(chromosome = c("chr1", "chr2"), 
             position = c(12414, 1234321))
  mset <- methylationSet(betas = beta_matrix, phenotypes = phenotypes, annotationDataFrame = annot)
  
  multi <- createMultiDataSet()
  multi <- add_methy(multi, mset)
  expect_equal(names(multi), "methylation")
  
  multi <- createMultiDataSet()
  msetbad <- mset
  colnames(fData(msetbad))[1] <- "chr"
  expect_error(multi <- add_methy(multi, msetbad), "fData of methySet must contain columns chromosome and position")
  
  multi <- createMultiDataSet()
  msetbad <- mset
  colnames(fData(msetbad))[2] <- "pos"
  expect_error(multi <- add_methy(multi, msetbad), "fData of methySet must contain columns chromosome and position")
  
  library(minfiData)
  minfiset <- ratioConvert(MsetEx[1:2, ])
  rowData(minfiset) <- fData(mset)
  
  multi <- createMultiDataSet()
  multi <- add_methy(multi, minfiset)
  expect_equal(names(multi), "methylation")
  
  
  multi <- createMultiDataSet()
  msetbad <- minfiset
  colnames(rowData(msetbad))[1] <- "chr"
  expect_error(multi <- add_methy(multi, msetbad), "fData of methySet must contain columns chromosome and position")
  
  multi <- createMultiDataSet()
  msetbad <- minfiset
  colnames(rowData(msetbad))[2] <- "pos"
  expect_error(multi <- add_methy(multi, msetbad), "fData of methySet must contain columns chromosome and position")
  
})

# Specific Wrong ####
test_that("specific functions based on eSet Wrong", {
    multi <- createMultiDataSet()
    eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
    fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
                              start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
                              status = c("case", "case", "control"), criteria = c("a", "b", "a"),
                              stringsAsFactors = FALSE)
    pData(eset) <- data.frame(id = letters[1:2])
    multi <- createMultiDataSet()
    
    # check add_genexp  ####
    eset2 <- eset
    colnames(fData(eset2))[1] <- "chr"
    expect_error(add_genexp(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
    eset2 <- eset
    colnames(fData(eset2))[2] <- "cot"
    expect_error(add_genexp(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
    eset2 <- eset
    colnames(fData(eset2))[3] <- "cot"
    expect_error(add_genexp(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
    
    # check add_rnaseq ####
    eset2 <- eset
    colnames(fData(eset2))[1] <- "chr"
    expect_error(add_rnaseq(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
    eset2 <- eset
    colnames(fData(eset2))[2] <- "cot"
    expect_error(add_rnaseq(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
    eset2 <- eset
    colnames(fData(eset2))[3] <- "cot"
    expect_error(add_rnaseq(multi, eset2), "fData of gexpSet must contain columns chromosome, start and end")
    
})

    
    
#   multiset <- createMultiDataSet()
#   
#   geno <- matrix(c(3,1,2,1), ncol = 2)
#   snpSet <- new("SnpSet", call = geno)
#   expect_warning(multi <- add_snps(multiset, snpSet), "No id column found in pData. The id will be equal to the sampleNames")
#   pheno <- data.frame(id = c("VAL0156", "VAL0372"))
#   pData(snpSet) <- pheno
#   
#   expect_error(add_snps(multiset, snpSet), "Set must contain a fData with column chromosome.")
#   fData(snpSet) <- data.frame(chr = c("chr1", "chr2"), pos = c(12414, 1234321),
#                               stringsAsFactors = FALSE)
#   expect_error(add_snps(multiset, snpSet), "Set must contain a fData with column chromosome.")
#   fData(snpSet) <- data.frame(chromosome = c("chr1", "chr2"), pos = c(12414, 1234321),
#                               stringsAsFactors = FALSE)
#   expect_error(add_snps(multiset, snpSet), "Set must contain a fData with columns position or start.")
#   fData(snpSet) <- data.frame(chr = c("chr1", "chr2"), position = c(12414, 1234321),
#                               stringsAsFactors = FALSE)
#   expect_error(add_snps(multiset, snpSet), "Set must contain a fData with column chromosome.")
#   fData(snpSet) <- data.frame(chr = c("chr1", "chr2"), start = c(12414, 1234321),
#                               stringsAsFactors = FALSE)
#   expect_error(add_snps(multiset, snpSet), "Set must contain a fData with column chromosome.")
#   fData(snpSet) <- data.frame(chromosome = c("chr1", "chr2"), start = c(12414, 1234321), position = c(12414, 1234321),
#                               stringsAsFactors = FALSE)
#   expect_error(add_snps(multiset, snpSet), "Set cannot contain a fData with columns position and start. Only one of them is allowed.")
#   
#   
#   eset <- new("ExpressionSet", exprs = matrix(runif(4), ncol = 2))
#   expect_error(multi <- add_genexp(multiset, eset), "pData of set must contain a column called id.")
#   pData(eset) <- data.frame(id = letters[1:2])
#   
#   expect_error(add_genexp(multiset, eset), "Set must contain a fData with column chromosome.")
#   fData(eset) <- data.frame(chr = c("chr1", "chr2"), pos = c(12414, 1234321),
#                             stringsAsFactors = FALSE)
#   expect_error(add_genexp(multiset, eset), "Set must contain a fData with column chromosome.")
#   fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), pos = c(12414, 1234321),
#                             stringsAsFactors = FALSE)
#   expect_error(add_genexp(multiset, eset), "Set must contain a fData with columns position or start.")
#   fData(eset) <- data.frame(chr = c("chr1", "chr2"), position = c(12414, 1234321),
#                             stringsAsFactors = FALSE)
#   expect_error(add_genexp(multiset, eset), "Set must contain a fData with column chromosome.")

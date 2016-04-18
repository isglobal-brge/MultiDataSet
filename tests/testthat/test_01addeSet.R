context("add eSet")

test_that("MultiSet", {
  multi <- createMultiDataSet()
  eset <- new("ExpressionSet", exprs = matrix(runif(6), ncol = 2))
  fData(eset) <- data.frame(chromosome = c("chr1", "chr2", "chr2"), 
    start = c(12414, 1234321, 1234328),end = c(121241, 1212412414, 1234378), 
    status = c("case", "case", "control"), criteria = c("a", "b", "a"),
    stringsAsFactors = FALSE)
  pData(eset) <- data.frame(id = letters[1:2])

  multi <- add_genexp(multi, eset)
  expect_equal(names(multi), "expression")
  
  expect_error(add_genexp(multi, eset), "There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
  expect_warning(add_genexp(multi, eset, overwrite = TRUE), "Slot 'expression' is already set in 'MultiDataSet'. Previous content will be overwritten.")
  
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
  multi <- add_snps(multi, snpSet)
  expect_equal(names(multi), c("expression", "snps"))
  
  expect_is(multi[["expression"]], "ExpressionSet")
  expect_is(multi[["snps"]], "SnpSet")
  
  
  library(methylumi)
  samps <- read.table(system.file("extdata/samples.txt",
                                  package = "methylumi"),sep="\t",header=TRUE)
  mldat <- methylumiR(system.file('extdata/exampledata.samples.txt',package='methylumi'),
                      qcfile=system.file('extdata/exampledata.controls.txt',package="methylumi"),
                      sampleDescriptions=samps)
  fvarLabels(mldat) <- tolower(fvarLabels(mldat))
  fData(mldat)$position <- fData(mldat)$cpg_coordinate
  pData(mldat)$id <- pData(mldat)$sampleID
  multi <- add_eset(multi, mldat, "cot", GRanges = NA)
  expect_equal(names(multi), c("expression", "snps", "cot"))
  expect_equal(names(multi[, c("cot", "snps")]), c("cot", "snps"))
  expect_is(multi[, "snps", drop = TRUE], "SnpSet")
})

test_that("subseting", {
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

test_that("wrong variables",{
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
})
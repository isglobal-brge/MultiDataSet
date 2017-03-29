context("add rse")
library(GenomicRanges)


test_that("MultiSet", {
  multi <- createMultiDataSet()
  
  counts <- matrix(runif(200 * 6, 1, 1e4), 200)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                       strand=sample(c("+", "-"), 200, TRUE),
                       feature_id=sprintf("ID%03d", 1:200))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6], id = LETTERS[1:6])
  names(rowRanges) <- 1:200
  rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                              rowRanges=rowRanges, colData=colData)
  
  multi <- add_rse(multi, rse, "rseEx")
  expect_equal(names(multi), "rseEx")
  
  expect_error(add_rse(multi, rse, "rseEx"), "There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
  expect_warning(add_rse(multi, rse, "rseEx", overwrite = TRUE), "Slot 'rseEx' is already set in 'MultiDataSet'. Previous content will be overwritten.")
  
  
  
  expect_is(multi[["rseEx"]], "RangedSummarizedExperiment")

  library(minfiData)
  data("MsetEx")
  GRset <- mapToGenome(MsetEx)
  colData(GRset)$id <- letters[1:6]
  multi <- add_rse(multi, GRset, "GRSet")
  expect_is(multi[["GRSet"]], "GenomicMethylSet")
  
  
})

test_that("wrong variables",{
  #TO DO
})
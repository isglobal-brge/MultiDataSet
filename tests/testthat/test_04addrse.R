context("add SE and RSE")
library(GenomicRanges)
library(SummarizedExperiment)


# SummarizedExperiment ####
test_that("add SummarizedExperiment", {
    multi <- createMultiDataSet()
    
    nrows <- 10; ncols <- 6
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                         row.names=LETTERS[1:6])
    rowData <- DataFrame(chr = "chr1", start = 1:10, end = 11:20)
    se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                colData=colData, rowData = rowData)
    
    expect_s4_class(multi2 <- add_se(multi, se0, "seEx"), "MultiDataSet")
    expect_equal(names(multi2), "seEx")
    expect_s4_class(multi2[["seEx"]], "SummarizedExperiment")
    
    
    # Check Overwrite ####
    expect_error(add_se(multi2, se0, "seEx"), "There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
    expect_warning(add_se(multi2, se0, "seEx", overwrite = TRUE), "Slot 'seEx' is already set in 'MultiDataSet'. Previous content will be overwritten.")
    
    # Check GRanges
    gr <- makeGRangesFromDataFrame(rowData)
    expect_s4_class(multi2 <- add_se(multi, se0, "seEx", GRanges = gr), "MultiDataSet")
    
    expect_s4_class(multi2 <- add_se(multi, se0, "seEx", GRanges = NA), "MultiDataSet")
    
    expect_error(multi2 <- add_se(multi, se0, "seEx", GRanges = "cot"), "GRanges should be a GenomicRanges or NA.")
})


# RangedSummarizedExperiment ####
test_that("add RangedSummarizedExperiment", {
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
  
  
  
  expect_s4_class(multi[["rseEx"]], "RangedSummarizedExperiment")

  library(minfiData)
  data("MsetEx")
  GRset <- mapToGenome(MsetEx[1:10, 1:2])
  
  expect_warning(multi2 <- add_rse(multi, GRset, "GRSet"), "No id column found in rowRanges. The id will be equal to the sampleNames")
  expect_equal(sampleNames(multi2[["GRSet"]]), c("5723646052_R02C02", "5723646052_R04C01"))
  expect_equal(multi2[["GRSet"]]$id, c("5723646052_R02C02", "5723646052_R04C01"))
  
  colData(GRset)$id <- letters[1:2]
  multi2 <- add_rse(multi, GRset, "GRSet")
  expect_is(multi2[["GRSet"]], "GenomicMethylSet")
  expect_equal(sampleNames(multi2[["GRSet"]]), c("5723646052_R02C02", "5723646052_R04C01"))
  expect_equal(multi2[["GRSet"]]$id, c("a", "b"))
})

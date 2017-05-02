context("Conversion")

library(GenomicRanges)
library(Biobase)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(MEALData)

# Create a MultiAssayExperiment

test_that("mae2mds", {
    arraydat <- matrix(seq(101, 108), ncol=4,
                       dimnames=list(c("ENST00000294241", "ENST00000355076"),
                                     c("array1", "array2", "array3", "array4")))
    arraypdat <- as(data.frame(slope53=rnorm(4),
                               row.names=c("array1", "array2", "array3",
                                           "array4")), "AnnotatedDataFrame")
    exprdat <- ExpressionSet(assayData=arraydat, phenoData=arraypdat)
    
    patient.data <- data.frame(sex=c("M", "F", "M", "F"),
                               age=38:41,
                               row.names=c("Jack", "Jill", "Bob", "Barbara"))
    exprmap <- data.frame(primary=rownames(patient.data)[c(1, 2, 4, 3)],
                          assay=c("array1", "array2", "array3", "array4"),
                          stringsAsFactors = FALSE)
    methyldat <-
        matrix(1:10, ncol=5,
               dimnames=list(c("ENST00000355076", "ENST00000383706"),
                             c("methyl1", "methyl2", "methyl3",
                               "methyl4", "methyl5")))
    methylmap <- data.frame(primary = c("Jack", "Jack", "Jill", "Barbara", "Bob"),
                            assay = c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5"),
                            stringsAsFactors = FALSE)
    microdat <- matrix(201:212, ncol=3,
                       dimnames=list(c("hsa-miR-21", "hsa-miR-191",
                                       "hsa-miR-148a", "hsa-miR148b"),
                                     c("micro1", "micro2", "micro3")))
    micromap <- data.frame(primary = c("Jack", "Barbara", "Bob"),
                           assay = c("micro1", "micro2", "micro3"),
                           stringsAsFactors = FALSE)
    # gr1 <-
    #     GRanges(seqnames = "chr3", ranges = IRanges(58000000, 59502360),
    #             strand = "+", score = 5L, GC = 0.45)
    # gr2 <-
    #     GRanges(seqnames = c("chr3", "chr3"),
    #             ranges = IRanges(c(58493000, 3), width=9000),
    #             strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
    # gr3 <-
    #     GRanges(seqnames = c("chr1", "chr2"),
    #             ranges = IRanges(c(1, 4), c(3, 9)),
    #             strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
    # grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
    # names(grl) <- c("snparray1", "snparray2", "snparray3")
    # 
    # rangemap <- data.frame(primary = c("Jack", "Jill", "Jill"),
    #                        assay = c("snparray1", "snparray2", "snparray3"),
    #                        stringsAsFactors = FALSE)
    
    nrows <- 5; ncols <- 4
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    rowRanges <- GRanges(rep(c("chr1", "chr2"), c(2, nrows - 2)),
                         IRanges(floor(runif(nrows, 1e5, 1e6)), width=100),
                         strand=sample(c("+", "-"), nrows, TRUE),
                         feature_id=sprintf("ID\\%03d", 1:nrows))
    names(rowRanges) <- letters[1:5]
    colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 2),
                         row.names= c("mysnparray1", "mysnparray2",
                                      "mysnparray3", "mysnparray4"))
    rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                rowRanges=rowRanges, colData=colData)
    rangemap2 <-
        data.frame(primary = c("Jack", "Jill", "Bob", "Barbara"),
                   assay = c("mysnparray1", "mysnparray2", "mysnparray3",
                             "mysnparray4"), stringsAsFactors = FALSE)
    listmap <- list(exprmap, methylmap, micromap, rangemap2)
    names(listmap) <- c("Affy", "Methyl 450k", "Mirna", "CNV gistic2")
    dfmap <- listToMap(listmap)
    
    
    objlist <- list("Affy" = exprdat, "Methyl 450k" = methyldat,
                    "Mirna" = microdat, "CNV gistic2" = rse)
    myMultiAssay <- MultiAssayExperiment(objlist, patient.data, dfmap)

    mds <- mae2mds(myMultiAssay)
    expect_is(mds, "MultiDataSet")
    expect_equal(names(mds), c("Affy", "Methyl 450k", "Mirna", "CNV gistic2"))
    expect_is(mds[["Affy"]], "ExpressionSet")
    expect_is(mds[["Methyl 450k"]], "ExpressionSet")
    expect_is(mds[["Mirna"]], "ExpressionSet")
    expect_is(mds[["CNV gistic2"]], "SummarizedExperiment")
    expect_equal(commonIds(mds), c("Jack", "Barbara", "Bob"))
    
})


test_that("mds2mae", {
    data(eset)
    data(pheno)
    fvarLabels(eset)[1] <- "chromosome"
    
    mset <- prepareMethylationSet(betavals, pheno)
    multi <- createMultiDataSet()
    multi <- add_methy(multi, mset)
    multi <- add_genexp(multi, eset)
    multi <- add_eset(multi, eset, dataset.type = "test", GRanges = NA)
    
    mae <- mds2mae(multi)
    
    expect_is(mae, "MultiAssayExperiment")
    expect_equal(names(mae), c("methylation", "expression", "test"))
    expect_equal(nrow(pData(mae)), length(Reduce(union, sampleNames(multi))))
    
    expect_is(experiments(mae)[[1]], "MethylationSet")
    expect_is(experiments(mae)[[2]], "ExpressionSet")
    expect_is(experiments(mae)[[3]], "ExpressionSet")

})
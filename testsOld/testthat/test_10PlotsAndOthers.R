context("Miscellaneous")

test_that("chrNumToChar", {
    chromosomes <- c(1, 3, 4, 23, 15, 24, 25, 26)
    expect_equal(chrNumToChar(chromosomes), 
                 c("chr1", "chr3", "chr4", "chrX", "chr15", "chrY", "chrXY", "chrMT"))
})

test_that("volcano plot", {
    pvals <- c(1e-3, 1e-4, 1e-5)
    names(pvals) <- c("a", "b", "c")
    fc <- c(-2, 4, -5)
    volcano_plot(pvals, fc)
    a <- volcano_plot(pvals, fc, show.effect = TRUE)
    a <- volcano_plot(pvals, fc, names = letters[1:3])
    a <- volcano_plot(pvals, fc, tFC = NULL, show.labels = TRUE)
    expect_equal(fc, fc)
})
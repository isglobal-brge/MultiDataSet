context("check ResultSet")

test_that("ResultSet basics", {

    library(limma)
    library(GenomicRanges)
    ## Create ResultSet with data.frame
    rset <- create_resultset("test", list(main = list(result = mtcars)), list())
    show(rset)
    expect_equal(getAssociation(rset), mtcars)
    
    ## Check length
    expect_equal(length(rset), 1)

    ## Check names
    expect_equal(names(rset), "main")
    
    ## Check opt
    expect_equal(opt(rset), list(fun_origin = "test"))
    
    ## Test fData
    expect_equal(fData(rset), list())    
    
    ## Check plot
    expect_error(plot(rset))
    
    ## Test show
    rset <- create_resultset("association", list(main = list(result = mtcars)), 
                             list(data.frame(chr = 1, start = 1:10)), lOptions = list(names = c("a", "b")))
    show(rset)
    
    ## Test varLabels
    expect_equal(varLabels(rset), list(main = NULL))
    
    ## Check opt
    expect_equal(opt(rset), list(fun_origin = "association", names = c("a", "b")))
    
    ## Create ResultSet with limma (use lmFit examples)
    sd <- 0.3*sqrt(4/rchisq(10, df = 4))
    y <- matrix(rnorm(10*6, sd = sd), 10, 6)
    rownames(y) <- paste("Gene", 1:10)
    y[1:2, 4:6] <- y[1:2, 4:6] + 2
    design <- cbind(Grp1 = 1, Grp2vs1 = c(0, 0, 0, 1, 1, 1))
    fit <- lmFit(y, design)
    fite <- eBayes(fit)
    df <- topTable(fite, coef = 2, number = Inf)
    
    fdata <- data.frame(chr = "1", start = 1:10, stringsAsFactors = FALSE)
    rownames(fdata) <- paste("Gene", 1:10)
    
    rset <- create_resultset("crossomics", list(pac = list(result = fit, error = NA), 
                                                cot = list(result = mtcars)), 
                             list(pac = fdata, cot = data.frame(chr = 1, start = 1:10)), 
                             list(method = "met", package = "metpac"))
    show(rset)
    
    ## Check getAssociation
    expect_equal(getAssociation(rset), df)
    expect_equal(getAssociation(rset, coef = 1, contrast = cbind(First = c(0,1))), df)
    expect_equal(getAssociation(rset, fNames = c("chr", "start")), cbind(df, fdata[rownames(df), ]))
    expect_error(getAssociation(rset, fNames = "chra"))
    
    ## Check length
    expect_equal(length(rset), 2)
    
    ## Check names
    expect_equal(names(rset), c("pac", "cot"))
    
    ## Check varLabels
    expect_equal(varLabels(rset), list(pac = c("Grp1", "Grp2vs1"), cot = NULL))
    
    ## Test fData
    expect_equal(fData(rset), list(pac = fdata, cot = data.frame(chr = 1, start = 1:10)))    
    
    ## Check plot
    plot(rset, type = "manhattan", fNames = c("chr", "start"))
    plot(rset, type = "manhattan", fNames = c("chr", "start"), highlight = GRanges("1:2-5"))
    plot(rset, type = "manhattan", fNames = c("chr", "start"), subset = GRanges("1:2-5"))
    plot(rset, type = "qq")
    plot(rset, type = "volcano")
})
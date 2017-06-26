## LOAD LIBRARIES
## ------------------------------------------------------------------------- ##
library(stringr)
library(limma)

## CREATE FEATURE DATA
## ------------------------------------------------------------------------- ##
fd <- data.frame(
    feature = c( paste0("ft", seq(1, 100)), paste0("ft", seq(501, 600)) ),
    chromosome = c( rep("chr1", 100), rep("chr2", 100)),
    position = c(1:100, 1:100),
    gene = c(
        paste0(sample(LETTERS[1:5], size=100, replace = T),
            sample(LETTERS[1:21], size=100, replace = T)),
        paste0(sample(LETTERS[15:20], size=100, replace = T),
           sample(LETTERS[1:21], size=100, replace = T))
    ),
    stringsAsFactors = FALSE
)
rownames(fd) <- fd$feature

## CREATE ASSAY DATA
## ------------------------------------------------------------------------- ##
as <- matrix(rnorm(n = 200 * 20), nrow = 200, ncol = 20)
as[2:4, 1:5] <- rnorm(n = 15, mean=1, sd=0.5)
rownames(as) <- fd$feature
colnames(as) <- paste0("s", str_pad(1:20, width = 3, pad="0"))

## PERFORM ASSOCIATION ANALYSIS
## ------------------------------------------------------------------------- ##
cc <- c(rep("case", 5), rep("control", 15))

dm <- model.matrix(~cc)
md <- lmFit(as, dm)

## CREATE RESULT SET
## ------------------------------------------------------------------------- ##
rset <- create_resultset(
    fOrigin = "dummy",
    lResult = list(dummy=list(result=md,error=NA)),
    fData = list(dummy=fd),
    lOptions = list(sva=FALSE)
)

plot(rset, type="manhattan", fNames=c("chromosome", "position"))

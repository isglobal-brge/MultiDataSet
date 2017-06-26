#' @describeIn ResultSet Allows to plot a series of plots (QQ plot, Manhattan
#' plot and Volcano plot) depending on the results stored in the
#' \code{ResultSet}.
#' @param y -
#' @param rid Name or index of the internal result to be used
#' @param coef Coefficient to be returne, usually 2
#' @param type Type of plot to be drawn
#' @param tPV Threshold for P-Value
#' @param tFC Threshold for log FC of effect
#' @param show.effect (default: \code{TRUE}). Used in volcano plot. If \code{TRUE}, 
#' effect is shown as FC instead of logFC. 
#' @param show.lambda (default: \code{TRUE}) If \code{TRUE} shows lambda
#' score for the given model.
#' @export
setMethod(
    f = "plot",
    signature = "ResultSet",
    definition = function(x, y, rid = 1, coef = 2, contrast = NULL, type, 
                          tFC = 2, tPV = -log10(0.001),
                          show.effect = FALSE, show.lambda = TRUE,
                          fNames = c("chromosome", "start"), ...) {
        
        if (class(x@results[[rid]]$result) != "MArrayLM"){
            stop("plot function is only available for results in a MArrayLM object.")
        }
        
        type <- tolower(type)
        type <- match.arg(type, choices = c("qq", "volcano", "manhattan"))
        
        if(type == "qq") {
            dta <- getAssociation(x, rid = rid, coef = coef,
                                  contrast = contrast, fNames = NULL)
            qq_plot(dta$P.Value, show.lambda = show.lambda)
        } else if(type == "manhattan") {
            dta <- getAssociation(x, rid = rid, coef = coef, contrast = contrast, 
                                  fNames = fNames)
            # Select columns P.Value, Chromosome and "Position"
            dta <- dta[ , c("P.Value", fNames)]
            # Add column with rs name per SNP
            dta$SNP <- dta
            # Update colnames of 'dta' by the ones used on qqman::manhattan
            colnames(dta) <- c("P", "CHR", "BP", "SNP")
            # Two steps:
            # 1) If col CHR has more than one chromsome joint by "_", 
            #    select only first
            # 2) If col CHR follows format chr1 ... chr21, remove "chr" term
            dta$CHR <- gsub("chr", "", sapply(strsplit(dta$CHR, "_"), "[[", 1))
            # If col CHR has chromosome X and Y instead of 23 and 24,
            # reformat to 23 and 24
            dta$CHR <- gsub("X", "23", gsub("Y", "24", dta$CHR))
            # Select only the SNPs in chromosome 1 to 24
            dta <- dta[dta$CHR %in% as.character(1:24), ]
            # Convert both pposition and chromosome to numeric
            dta$CHR <- as.numeric(dta$CHR)
            dta$BP <- as.numeric(dta$BP)
            # Plot using qqman
            qqman::manhattan(dta, ylab = "-log10(P.Value)", ...)
        } else if(type == "volcano") {
            dta <- getAssociation(x, rid = rid, coef = coef, contrast = contrast,
                                  fNames = NULL)
            volcano_plot(
                pval = dta$P.Value,
                fc = dta$logFC,
                names = rownames(dta),
                tFC = tFC,
                tPV = tPV,
                show.effect = show.effect
            )
        } else {
            stop("Invalid type of plot ('", type, "').")
        }
    })

#' @describeIn ResultSet Allows to plot a series of plots (QQ plot, Manhattan
#' plot and Volcano plot) depending on the results stored in the
#' \code{ResultSet}.
#' @param rid Name or index of the internal result to be used
#' @param coef Coefficient to be returne, usually 2
#' @param contrast If coefficient to be used was multicategorical, number
#' of the contrast to be returned.
#' @param type Type of plot to be drawn
#' @param tPV Threshold for P-Value
#' @param tFC Threshold for log FC of effect
setMethod(
    f = "plotAssociation",
    signature = "ResultSet",
    definition = function(object, rid = 1, coef = 2, contrast = NULL, type = "volcano", 
                          tPV = NULL, tFC = NULL, show.effect = FALSE) {
        
        type <- tolower(type)
        type <- match.arg(type, choices = c("qq", "volcano", "manhattan"))
        
        dta <- getAssociation(object, rid=rid, coef = coef, contrast = contrast, 
                        fNames = c("chromosome", "position"))
        
        if(type == "qq") {
            qq_plot(dta$P.Value)
        } else if(type == "manhattan") {
            dta$SNP <- rownames(dta)
            colnames(dta) <- c("P", "CHR", "BP", "SNP")
            dta$CHR <- gsub("chr", "", sapply(strsplit(dta$chromosome, "_"), "[[", 1))
            dta$CHR <- gsub("X", "23", gsub("Y", "24", dta$CHR))
            dta <- dta[dta$CHR %in% as.character(1:24), ]
            dta$CHR <- as.numeric(dta$CHR)
            dta$BP <- as.numeric(dta$start)
            
            qqman::manhattan(dta, ylab="-log10(P.Value)", ...)
        } else if(type == "volcano") {
            volcano_plot(
                pval = dta$P.Value,
                fc = dta$logFC,
                names = rownames(dta),
                tFC = tFC,
                tPV = tPV,
                show.effect = show.effect
            )
        }  else {
            stop("Invalid type of plot ('", type, "').")
        }
    
    })

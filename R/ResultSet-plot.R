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
    f = "plot",
    signature = "ResultSet",
    definition = function(x, y, rid = 1, coef = 2, contrast = 1, type, tPV, tFC, show.effect=FALSE) {
        
        if(type == "qq") {
            dta <- topTable(object, rid=rid, coef=coef, contrast=contrast,
                            fData = NULL)
            qq_plot(dta$P.Value)
        } else if(type == "manhattan") {
            dta <- topTable(object, rid=rid, coef=coef, contrast=contrast)
            # Select columns P.Value, Chromosome and "Position"
            dta <- dta[ , c("P.Value", "chromosome", "start")]
            # Add column with rs name per SNP
            dta$SNP <- rownames(object@results[[rid]]$result)
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
            qqman::manhattan(dta, ylab="-log10(P.Value)", ...)
        } else if(type == "volcano") {
            dta <- topTable(object, rid=rid, coef=coef, contrast=contrast,
                            fData = NULL)
            volcano_plot(
                pval=dta$P.Value,
                fc=dta$logFC,
                names=rownames(dta),
                tFC=tFC,
                tPV=tPV,
                show.effect=effect
            )
        } else {
            stop("Invalid type of plot ('", type, "').")
        }
    })

#' Chemical Data Selection within the Floating Percentile Model
#' 
#' Generate a dataset of concentrationso for chemicals with significantly higher concentrations 
#' among toxic samples by comparison to non-toxic samples
#' 
#' @param  data data.frame containing, at a minimum, chemical concentrations as columns and a logical \code{Hit} column classifying toxicity
#' @param  paramList character vector naming columns of \code{data} containing concentrations
#' @param  plot logical value indicating whether or not to generate figures comparing concentrations between subsets where \code{Hit == TRUE} and \code{Hit == FALSE} (default = \code{FALSE})
#' @param  ... additional arguments passed to \code{chemSig}
#' @details \code{chemSigSelect} is used within \code{FPM} to quickly run \code{chemSig} and export concentrations from \code{data} for the
#' significant chemicals to include in the floating percentile model algorithm. Results are exported in list with two data.frames, separated into chemicals that were significant (\code{"sig"}) and non-significant (\code{"nonsig"})
#' For information on \code{data} and \code{paramList}, see details of \code{?FPM}. 
#' 
#' If \code{plot = TRUE}, then a series of plots will be exported comparing the Hit/No-hit data subsets for all chemicals in \code{paramList}. Plots with blue lines are generated from
#' significant chemicals, and plots with green lines are generated from non-significant chemicals. The user currently has only
#' limited control over graphical parameters of plots, and attempts to change graphical parameters are likely to cause errors.
#' @seealso plot.chemSigSelect, FPM, chemSig
#' @return list of two data.frames (\code{"sig"} and \code{"nonsig"}), \code{chemSigSelect} class object
#' @examples
#' paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn")
#' chemSigSelect(h.tristate, paramList)$sig[1:6,]
#' chemSigSelect(h.tristate, paramList, testType = "p")$sig[1:6,]
#' @export
chemSigSelect <- function(data, 
                          paramList, 
                          plot = FALSE, 
                          ...){
    sig <- chemSig(data = data, paramList = paramList, ...)
    outP <- names(sig)[sig]
    outP2 <- paramList[!paramList %in% outP]
    
    out <- data.frame(data[, c(outP, "Hit")])

    attr(out, "nChem") <- length(outP)
    if(attr(out, "nChem") == 0) {
        stop("No significant chemicals were selected")
    }
    
    out2 <- data.frame(data[, c(outP2, "Hit")])
    
    ## reordering paramList to include significant, then non-significant chemicals
    p.paramList <- c(paramList[which(paramList %in% colnames(out))], paramList[which(!paramList %in% colnames(out))])
    
    out3 <- list(sig = out, nonsig = out2)
    class(out3) <- "chemSigSelect"
    
    if(plot) plot.chemSigSelect(out3, ...)
    return(out3)
}## end code

#' Plot method for 'chemSigSelect'
#' 
#' The plot method for 'chemSigSelect' objects
#' 
#' @param  x \code{chemSigSelect} class object
#' @param  paramList character vector of column names of chemical concentration variables in \code{data}; values will be inherited from \code{x} if not specified
#' @param  type character vector indicating the type of plot to generate (see Details; default = \code{c("boxplot", "ecdf", "density")})
#' @param  logAxes logical value; whether axes should be projected in log-space (default = \code{TRUE})
#' @param  ... additional arguments passed to \code{plot}. The user currently has limited capacity to change graphical parameters.
#' @details \code{plot.chemSigSelect} is a class-specific plotting function that generates boxplots, empirical cumulative density functions, or
#' kernel density plots comparing the \code{Hit == TRUE} and \code{Hit == FALSE} subsets for all chemicals in \code{paramList}.
#' By default, all plot types are generated, but \code{type} can be changed to generate specific plot types. The available types are
#' specified by \code{"type = boxplot"}, \code{"type = ecdf"}, and \code{"type = density"}, though character inputs are flexible; 
#' for example, other reasonable versions of these types include \code{"Boxplot"} and \code{"ECDF"}.
#' Note that because \code{boxplot}, \code{plot.ecdf}, and \code{plot.density} functions use plotting arguments in different ways, supplying
#' user-defined plotting arguments may cause unintended changes (or possibly errors or warnings) when plotting multiple types. For that reason, we recommend that the user
#' select a specific plot type (using the \code{type} argument) before adjusting plot arguments (via \code{...}) that apply to that type. 
#' Moreover, many of the plotting arguments are specified by \code{plot.chemSigSelect} and, therefore, cannot be 
#' adjusted by the user. These include but are not limited to parameters like \code{col} and \code{lwd}.
#' @seealso chemSigSelect, boxplot, plot.ecdf, plot.density
#' @return base graphics
#' @importFrom grDevices rgb
#' @importFrom graphics boxplot
#' @importFrom graphics points
#' @importFrom stats ecdf
#' @importFrom stats plot.ecdf
#' @importFrom graphics legend
#' @importFrom stats density 
#' @importFrom graphics lines 
#' @importFrom graphics mtext
#' @examples
#' x <- chemSigSelect(data = h.tristate, paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn"))
#' plot(x, type = "boxplot")
#' @export
plot.chemSigSelect <- function(x, 
                               paramList = NULL,
                               type = c("boxplot", "ecdf", "density"), 
                               logAxes = TRUE,
                               ...) {
    if(is.null(paramList)){
        paramList <- c(colnames(x[[1]]), colnames(x[[2]]))
    }
    paramList <- paramList[paramList != "Hit"]

    
    # Boxplots ------------------------------------------------------------------------
    if("boxplot" %in% type|
       "Boxplot" %in% type|
       "box" %in% type|
       "Box" %in% type|
       "b" %in% type|
       "B" %in% type|
       "box plot" %in% type|
       "Box plot" %in% type|
       "Box Plot" %in% type) {
        for(i.parm in paramList) {
            
            if(i.parm %in% colnames(x[[1]])) {
                col.string <- c("royalblue", "lightskyblue", "lightskyblue3")
                rb.t <- rgb(65, 105, 225, maxColorValue = 255, alpha = 150, names = "royalblue50")
                lsb.t <- rgb(135, 206, 250, maxColorValue = 255, alpha = 150, names = "lightskyblue50")
                over.t <- rgb(100, 155.5, 237.5, maxColorValue = 255, alpha = 150, names = "overlap.col50")
                t.col.string <- c(rb.t, lsb.t)
                tmp <- x[[1]][,c(i.parm, "Hit")]
            } else {
                col.string <- c("springgreen4", "springgreen1", "springgreen2")
                sg4.t <- rgb(0, 139, 69, maxColorValue = 255, alpha = 150, names = "springgreen4.50")
                sg1.t <- rgb(0, 255, 127, maxColorValue = 255, alpha = 150, names = "springgreen1.50")
                over.tg <- rgb(0, 197, 98, maxColorValue = 255, alpha = 150, names = "overlap2.col50")
                t.col.string <- c(sg4.t, sg1.t, over.tg)
                tmp <- x[[2]][,c(i.parm, "Hit")]
            }
           
            tmp$num <- ifelse(tmp$Hit=="TRUE", 1, 2)
            tmp$col <- ifelse(tmp$Hit=="TRUE", t.col.string[1], t.col.string[2])
             
            
            if(logAxes){logAx = "y"} else {logAx = ""}
            boxplot(tmp[,i.parm] ~ tmp$num, outline = FALSE, log = logAx,
                    ylim = c(min(tmp[,i.parm]) * 0.8, max(tmp[, i.parm]) * 1.2),
                    names = FALSE, 
                    xlab = "", ylab = paste0("[ ", i.parm, " ]"), cex.lab = 1.5,
                    main = i.parm, cex.main = 2,
                    col = "white", ...)
            points(y = tmp[, i.parm],
                   x = jitter(tmp[,"num"], 0.75),
                   pch = 21, bg = tmp$col, col = tmp$col)
            mtext("Hit", side = 1, at = 1, line = 3, cex = 1.5) 
            mtext("No Hit", side = 1, at = 2, line = 3, cex = 1.5)
        }
    }
    
    
    # ECDFs ------------------------------------------------------------------------
    if("ecdf" %in% type|
       "ECDF" %in% type|
       "CDF" %in% type|
       "cdf" %in% type) {
        for(i.parm in paramList) {
            
            if(i.parm %in% colnames(x[[1]])) {
                tmp <- x[[1]][,c(i.parm, "Hit")]
            } else {
                tmp <- x[[2]][,c(i.parm, "Hit")]
            }
            
            h.ecdf <- ecdf(tmp[tmp$Hit == TRUE, i.parm])
            nh.ecdf <- ecdf(tmp[tmp$Hit == FALSE, i.parm])
            
            if(i.parm %in% colnames(x[[1]])) {
                tmp$col <- ifelse(tmp$Hit=="TRUE", "royalblue", "lightskyblue")
            } else {
                tmp$col <- ifelse(tmp$Hit=="TRUE", "springgreen4", "springgreen1") 
            }
            
            if(logAxes){logAx = "x"} else {logAx = ""}
            plot.ecdf(h.ecdf, verticals = TRUE, pch = NA, lwd = 2, log = logAx,
                      xlim = c(min(tmp[,i.parm]), max(tmp[,i.parm])),
                      col = ifelse(i.parm %in% colnames(x[[1]]), "royalblue", "springgreen4"),
                      main = i.parm,
                      xlab = paste0("[ ", i.parm, " ]"), ...)
            plot.ecdf(nh.ecdf, add = TRUE, verticals = TRUE, pch = NA, lwd = 2, 
                      col = ifelse(i.parm %in% colnames(x[[1]]), "lightskyblue", "springgreen1"))
            legend("bottomright", c("Hit", "No Hit"), lty = 1, lwd = 2, 
                   col = c(ifelse(i.parm %in% colnames(x[[1]]), "royalblue", "springgreen4"), 
                           ifelse(i.parm %in% colnames(x[[1]]), "lightskyblue", "springgreen1")))
            
        }
    }
    
    
    # Kernel density ------------------------------------------------------------------------
    if("density" %in% type|
       "dens" %in% type|
       "Density" %in% type|
       "d" %in% type|
       "D" %in% type) {
        for(i.parm in paramList) {
            
            if(i.parm %in% colnames(x[[1]])) {
                col.string <- c("royalblue", "lightskyblue", "lightskyblue3")
                rb.t <- rgb(65, 105, 225, maxColorValue = 255, alpha = 150, names = "royalblue50")
                lsb.t <- rgb(135, 206, 250, maxColorValue = 255, alpha = 150, names = "lightskyblue50")
                over.t <- rgb(100, 155.5, 237.5, maxColorValue = 255, alpha = 150, names = "overlap.col50")
                t.col.string <- c(rb.t, lsb.t)
                tmp <- x[[1]][,c(i.parm, "Hit")]
            } else {
                col.string <- c("springgreen4", "springgreen1", "springgreen2")
                sg4.t <- rgb(0, 139, 69, maxColorValue = 255, alpha = 150, names = "springgreen4.50")
                sg1.t <- rgb(0, 255, 127, maxColorValue = 255, alpha = 150, names = "springgreen1.50")
                over.tg <- rgb(0, 197, 98, maxColorValue = 255, alpha = 150, names = "overlap2.col50")
                t.col.string <- c(sg4.t, sg1.t, over.tg)
                tmp <- x[[2]][, c(i.parm, "Hit")]
            }
            
            h.dens <- density(tmp[tmp$Hit == TRUE, i.parm], from = 10^(floor(log10(min(tmp[tmp$Hit == TRUE, i.parm])))))
            nh.dens <- density(tmp[tmp$Hit == FALSE, i.parm], from = 10^(floor(log10(min(tmp[tmp$Hit == TRUE, i.parm])))))
            
            if(logAxes){
                logAx <- "x"
            } else {
                logAx <- ""
            }
            
            plot(h.dens, 
                 xlim = c(min(tmp[,i.parm]), max(tmp[,i.parm])), log = logAx,
                 ylim = c(min(h.dens$y, nh.dens$y), max(h.dens$y, nh.dens$y)),
                 col = col.string[1],
                 lwd = 2,
                 xlab = paste0("[ ", i.parm, " ]"),
                 main = i.parm, ...)
            lines(nh.dens, col = col.string[2])
            legend("topright", c("Hit", "No Hit"),
                   col = col.string, lwd = 2, lty = 1,
                   cex = 1)
        }
    }
}## end code
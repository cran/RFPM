#' Cross-Validation Optimization of the Floating Percentile Model
#' 
#' Use leave-one-out (LOO) or k-folds cross-validation methods to calculate parameter inputs that optimize benchmark performance
#' while attempting to account for out-of-sample uncertainty
#' 
#' @param  data data.frame containing, at a minimum, chemical concentrations as columns and a logical \code{Hit} column classifying toxicity
#' @param  paramList character vector of column names of chemical concentration variables in \code{data}
#' @param  FN_crit numeric vector of values between 0 and 1 indicating false negative threshold for floating percentile model benchmark selection (default = \code{seq(0.1, 0.9, 0.05)})
#' @param  alpha.test numeric vector of values between 0 and 1 indicating type-I error rate for chemical selection (default = \code{seq(0.05, 0.5, by = 0.05)})
#' @param  k integer with length = 1 and value > 1 indicating how many folds to include in k-folds type cross-validation method (default = \code{5})
#' @param  seed integer with length = 1 indicating the random seed to set for assigning k classes for k-folds cross-validation (default = \code{1})
#' @param  plot whether to plot the output of \code{cvFPM} (default = \code{TRUE})
#' @param  which numeric or character indicating which type of plot to generate (see Details; default = \code{c(1, 2, 3, 4)})
#' @param  colors values recognizible as colors to be passed to \code{colorRampPalette} (via \code{colorGradient}) to generate a palette for plotting (default = \code{heat.colors(10)})
#' @param  colsteps integer; number of discrete steps to interpolate colors in \code{colorGradient} (default = \code{100})
#' @param  ... additional arguments passed to \code{chemSig} and \code{FPM}
#' @details 
#' \code{cvFPM} allows users to "tune" the \code{FN_crit} and \code{alpha.test} arguments used by \code{FPM} (via \code{chemSig}). This is achieved by splitting the empirical dataset into "test" and "training" subsets,
#' calculating benchmarks for the training set, and then evaluating the benchmarks' ability to predict Hits in the out-of-sample test set. The output of \code{cvFPM} is similar to \code{optimFPM}: optimal \code{FN_crit} and \code{alpha.test} inputs
#' based on several classification metrics (see \code{?optimFPM} for more details). The key difference between \code{cvFPM} and \code{optimFPM} is that \code{cvFPM} attempts to account for
#' out-of-sample uncertainty, whereas \code{optimFPM} is specific (and potentially overly specific) to the full FPM dataset. Because the primary use of FPM SQBs will be to predict
#' toxicity in sediment samples where toxicity is not measured, the FPM should be parameterized in a way that best accounts for out-of-sample uncertainty. In other words, while FPM
#' generates classification metrics like "overall reliability" for SQBs, they are unlikely to achieve the expected reliability when applied to new samples. This is an inherent limitation of SQBs,
#' which the FPM cannot fully address but that \code{cvFPM} considers.
#' 
#' Two cross-validation methods are available, controlled through the \code{k} argument. If the user specifies \code{k = NULL} or \code{k = nrow(data)}, then leave-one-out (LOO) is used. 
#' LOO is computationally intensive but is better suited to small datasets.
#' Other values of \code{k} (e.g., the default \code{k = 5}) will result in applying a k-folds cross-validation method, which uses larger test subsets (and smaller training sets), 
#' evaluates fewer scenarios, and greatly improves runtime for large datasets. The \code{seed} argument can be used to establish a consistent result; if \code{is.null(seed)}, the result will vary based on randomization.
#' Allowing for randomization may be desireable to understand between-run variability in \code{cvFPM} output caused by re-sampling of training/test sets.
#' 
#' By setting \code{plot = TRUE} (the default), the outcome of cross-validation can be visualized over the range of \code{FN_crit} values considered. Visualizing the results
#' can inform the user about variability in the cross-validation process, ranges of potentially reasonable \code{FN_crit} values, etc. Graphical output depends on
#' whether many \code{FN_crit} and/or many \code{alpha.test} are evaluated, with line plots or heat plots alternately generated.
#' 
#' IMPORTANT: \code{cvFPM} is not in itself optimized for runtime - running \code{cvFPM} can take a long time
#' 
#' The \code{which} argument can be used to specify which of the metric-specific plots should be generated when \code{plot = TRUE}. Inputs
#' to \code{which} are, by default, \code{c(1, 2, 3, 4)}.
#' 
#' @seealso chemSig, FPM, optimFPM
#' @return data.frame of metric output, base R graphical output
#' @importFrom stats median
#' @importFrom graphics par
#' @examples
#' paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn")
#' par(mfrow = c(2,2))
#' cvFPM(h.tristate, paramList, seq(0.1, 0.9, 0.1), 0.05)
#' @export

cvFPM <- function(data,
                  paramList,
                  FN_crit = seq(0.1, 0.9, by = 0.05),
                  alpha.test = seq(0.05, 0.5, by = 0.05),
                  k = 5,
                  seed = 1,
                  plot = TRUE,
                  which = c(1, 2, 3, 4),
                  colors = heat.colors(10),
                  colsteps = 100,
                  ...){
    
    if(any(FN_crit >= 1 | FN_crit <= 0)){
        stop("FN_crit must only include values between 0 and 1")
    }
    if(any(alpha.test >= 1 | alpha.test <= 0)){
        stop("alpha.test must only include values between 0 and 1")
    }
    if(nrow(data[data$Hit,])<=3 | nrow(data[!data$Hit,])<=3){
        stop("data must contain more than 3 Hit and 3 No-Hit samples")
    }
    if(!is.null(k)){
        if(k <= 1) {stop("k must be a single value >1 and <=nrow(data)")}
    }
    
    ## Initilize inputs/outputs
    grid <- expand.grid(FN_crit, alpha.test)
    out <- list()
    if(!is.null(seed)){
        set.seed(seed)
    }
    ## Attempt to speed function by default using subsampling + LOO - blends k-folds and LOO
    if(!is.null(k)) {
        hitPortion <- ceiling(mean(data$Hit) * (nrow(data) / k)) ## per-fold sample size for Hits
        nohitPortion <- floor(mean(!data$Hit) * (nrow(data) / k)) ## per-fold samlpe size for No-Hits
        
        if(hitPortion < 3 | nohitPortion < 3) {stop("k-folds method not possible - Hit or No-hit N < 3; reduce k and try again")}
        
        data2 <- split(data, data$Hit)
        kHit <- rep(1:k, each = hitPortion)
        if(length(kHit) < nrow(data2[["TRUE"]])) {
            kHit <- c(kHit, sample(1:k, nrow(data2[["TRUE"]]) - length(kHit), replace = F))
        }
        data2[["TRUE"]]$k <- NA
        data2[["TRUE"]]$k <- sample(kHit, size = nrow(data2[["TRUE"]]), replace = F)

        kNoHit <- rep(1:k, each = nohitPortion)
        if(length(kNoHit) < nrow(data2[["FALSE"]])) {
            kNoHit <- c(kNoHit, sample(1:k, nrow(data2[["FALSE"]]) - length(kNoHit), replace = F))
        }
        data2[["FALSE"]]$k <- NA
        data2[["FALSE"]]$k <- sample(kNoHit, size = nrow(data2[["FALSE"]]), replace = F)
        
        data2 <- rbind(data2[["TRUE"]], data2[["FALSE"]])
    } else {
        k <- nrow(data2) # k is never null - LOO: k = nrow(data)
        data2 <- data
        data2$k <- 1:nrow(data2)
    }
    
    ## Results for full dataset - for plotting to see loss of info from out-of-sample uncertainty
    fpm.full <- do.call(rbind, apply(grid, 1, function(i){
        FPM(data, paramList = paramList, FN_crit = i[[1]], alpha.test = i[[2]], ...)[["FPM"]][,c("sens", "spec", "OR", "FM", "MCC")]
    }))
    fpm.full$sensSpecRatio <- signif(apply(fpm.full[,c("sens", "spec")], 1, function(x) min(x)/max(x)), 3)
    fpm.full <- fpm.full[, c("sensSpecRatio", "OR", "FM", "MCC")]
    
    out <- apply(grid, 1, function(i){
        out.i <- lapply(1:k, function(j) {
            train <- data2[!data2$k == j,]
            test <- data2[data2$k == j,]
            if(sum(train$Hit)<3 | sum(!train$Hit) <3) stop("Error in subset - Hit or No-hit N<3; reduce k and try again")

            fpm.j <- FPM(train, paramList = paramList, FN_crit = i[[1]], alpha.test = i[[2]], ...)[["FPM"]]
            n.j <- names(fpm.j)[names(fpm.j) %in% paramList]
            
            Hit_result <- data.frame(
                pHit = logical(nrow(test)),
                tHit = logical(nrow(test)),
                compHit = logical(nrow(test))
            )

            # determinig if Hit classification is correct or not
            Hit_result[, "pHit"] <- predict.FPM(fpm.j, newdata = test)# predicted Hit class
            Hit_result[, "tHit"] <- test$Hit # actual Hit class
            Hit_result[, "compHit"] <- Hit_result[, "tHit"] == Hit_result[, "pHit"] # correct Hit class?
            
            TP <- sum(Hit_result[,"pHit"] & Hit_result[,"tHit"])
            TN <- sum(!Hit_result[,"pHit"] & !Hit_result[,"tHit"])
            FP <- sum(Hit_result[,"pHit"] & !Hit_result[,"tHit"])
            FN <- sum(!Hit_result[,"pHit"] & Hit_result[,"tHit"])
            
            sens <- TP/(TP + FN)
            spec <- TN/(TN + FP)
            sensSpecRatio <- min(c(sens, spec))/max(c(sens, spec))
            if(is.nan(sensSpecRatio)) sensSpecRatio <- 0 ## allowing for NAs to be input as zero values - worst possible outcome
            ppr <- TP/(TP + FP)
            OR <- mean(Hit_result[,"compHit"])
            FM <- exp(mean(log(c(sens, ppr))))
            if(is.nan(FM)) FM <- 0 ## allowing for NAs to be input as zero values - worst possible outcome
            MCC <- (TP * TN - FP * FN)/sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN))
            # not allowing this to be zero or -1 to avoid confusion in interpretion
            return(data.frame(FN_crit = i[[1]], alpha.test = i[[2]], sensSpecRatio, OR, FM, MCC))
        })
        k_runs <- do.call(rbind, out.i)
        return(
            with(k_runs,
                data.frame(
                    FN_crit = unique(FN_crit),
                    alpha.test = unique(alpha.test),
                    sensSpecRatio.min = min(sensSpecRatio, na.rm = T),
                    sensSpecRatio.max = max(sensSpecRatio, na.rm = T),
                    sensSpecRatio.med = median(sensSpecRatio, na.rm = T),
                    sensSpecRatio.mean = mean(sensSpecRatio, na.rm = T),
                    OR.min = min(OR, na.rm = T),
                    OR.max = max(OR, na.rm = T),
                    OR.med = median(OR, na.rm = T),
                    OR.mean = mean(OR, na.rm = T),
                    FM.min = min(FM, na.rm = T),
                    FM.max = max(FM, na.rm = T),
                    FM.med = median(FM, na.rm = T),
                    FM.mean = mean(FM, na.rm = T),
                    MCC.min = min(MCC, na.rm = T),
                    MCC.max = max(MCC, na.rm = T),
                    MCC.med = median(MCC, na.rm = T),
                    MCC.mean = mean(MCC, na.rm = T)
                )
            )
        )
    })
    
    out <- do.call(rbind, out); rownames(out) <- NULL
    
    # use for plotting
    metrics <- c("sensSpecRatio", "OR", "FM", "MCC")
    Summary <- do.call(c, lapply(metrics, function(m){
        tmp <- rowSums(out[,grep(x = names(out), pattern = m)], na.rm = T)
        min(which(tmp == max(tmp, na.rm = T)), na.rm = T)
        }))
    names(grid) <- c("FN_crit", "alpha.test")
    pickArgs <- grid[Summary, 1:2]
    rownames(pickArgs) <- c("sensSpecRatio", "OR", "FM", "MCC")
    
    ## ------------------------------------------------------
    ## Plotting
    
    if(plot){
        if(!any(which %in% c(1, 2, 3, 4))){
            stop("which argument must equal or include 1, 2, 3, and/or 4")
        }
        
        if(length(FN_crit) > 1 & length(alpha.test) == 1){
            if(1 %in% which){
                outData <- out[,grep(x = names(out), pattern = "sensSpecRatio")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,1], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "False Negative Limit", ylab = "Sensitivity/Specificity Ratio",
                      ylim = c(0,1))
                 abline(v = pickArgs[1, "FN_crit"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,1], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,1], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = fpm.full[,1], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, 
                        x.intersp = 0.25, xpd = TRUE, lty = c(1, 1, 3, 2), 
                        col = c("dodgerblue3", "darkorange2", "gray20", "black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
    
            if(2 %in% which){
                outData <- out[,grep(x = names(out), pattern = "OR")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,1], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "False Negative Limit", ylab = "Overall Reliability",
                      ylim = c(0,1))
                 abline(v = pickArgs[2, "FN_crit"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,1], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,1], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = fpm.full[,2], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
            
            if(3 %in% which){
                outData <- out[,grep(x = names(out), pattern = "FM")]
                names(outData) <- c("min", "max", "median", "mean")
                    
                 plot(x = grid[,1], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "False Negative Limit", ylab = "Fowlkes-Mallows Index",
                      ylim = c(0,1))
                 abline(v = pickArgs[3, "FN_crit"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,1], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,1], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = fpm.full[,3], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black","black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
            
            
            if(4 %in% which){
                outData <- out[,grep(x = names(out), pattern = "MCC")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,1], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "False Negative Limit", ylab = "Matthew's Correlation Coefficient",
                      ylim = c(-0.5,1))
                 abline(v = pickArgs[4, "FN_crit"], col = "black", lty = 2, lwd = 1)
                 abline(h = 0, col = "lightgray", lty=3)
                 lines(x = grid[,1], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,1], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,1], y = fpm.full[,4], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
        } else if(length(FN_crit) == 1 & length(alpha.test) > 1){
            if(1 %in% which){
                outData <- out[,grep(x = names(out), pattern = "sensSpecRatio")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,2], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "Test Alpha", ylab = "Sensitivity/Specificity Ratio",
                      ylim = c(0,1))
                 abline(v = pickArgs[1, "alpha.test"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,2], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,2], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = fpm.full[,1], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
    
            if(2 %in% which){
                outData <- out[,grep(x = names(out), pattern = "OR")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,2], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "Test Alpha", ylab = "Overall Reliability",
                      ylim = c(0,1))
                 abline(v = pickArgs[2, "alpha.test"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,2], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,2], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = fpm.full[,2], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
            
            if(3 %in% which){
                outData <- out[,grep(x = names(out), pattern = "FM")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = grid[,2], y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "Test Alpha", ylab = "Fowlkes-Mallows Index",
                      ylim = c(0,1))
                 abline(v = pickArgs[3, "alpha.test"], col = "black", lty = 2, lwd = 1)
                 lines(x = grid[,2], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,2], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = fpm.full[,3], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
            
            
            if(4 %in% which){
                outData <- out[,grep(x = names(out), pattern = "MCC")]
                names(outData) <- c("min", "max", "median", "mean")
                
                 plot(x = alpha.test, y = outData$median, type = "l", col = "dodgerblue3", lwd = 1,
                      xlab = "Test Alpha", ylab = "Matthew's Correlation Coefficient",
                      ylim = c(-0.5, 1) )
                 abline(v = pickArgs[4, "alpha.test"], col = "black", lty = 2, lwd = 1)
                 abline(h = 0, col = "lightgray", lty=3)
                 lines(x = grid[,2], y = outData$mean, col = "darkorange2")
                 lines(x = grid[,2], y = outData$min, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = outData$max, col = "gray20", lty = 3)
                 lines(x = grid[,2], y = fpm.full[,4], col = "black", lwd = 1.5)
                 legend("bottom", bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                        lty = c(1,1,3,2), 
                        col = c("dodgerblue3", "darkorange2", "gray20","black", "black"),
                        legend = c("Median", "Mean", "Min/Max", "Optimized value", "All-sample result"))
            }
        } else {
            heads <- c("Sensitivity / Specificity Ratio", "Overall Reliability", 
                       "Fowlkes-Mallows Index", "Matthew's Correlation Coefficient")
            
            for(i in which){
                outData <- out[,grep(x = names(out), pattern = c("sensSpecRatio", "OR", "FM", "MCC")[i])]
                names(outData) <- c("min", "max", "median", "mean")                
                
                cols <- colorGradient(x = 1 - outData$mean, colors = colors, colsteps = colsteps, na.rm = T)
                
                plot(grid[,1], grid[,2], pch = 15, cex = 3, col = cols,
                     xlab = "False Negative Limit", ylab = "Test Alpha", main = heads[i])
                points(pickArgs$FN_crit[i], pickArgs$alpha.test[i], pch = 0, cex = 4, col = "black")
                mtext(side = 3, text = paste0("Mean range: ", 
                                              signif(min(outData$mean), 2), "-", 
                                              signif(max(outData$mean), 2),"; ",
                                              "optimized mean: ", signif(outData$mean[Summary[i]], 2)
                                              ))
            }
        }
    }
    return(pickArgs)
}# end code

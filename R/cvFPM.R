#' Cross-Validation Optimization of the Floating Percentile Model
#' 
#' Use k-folds cross-validation method to calculate parameter inputs that optimize benchmark performance
#' while attempting to account for out-of-sample error
#' 
#' 
#' @param  data data.frame containing, at a minimum, chemical concentrations as columns and a logical \code{Hit} column classifying toxicity
#' @param  paramList character vector of column names of chemical concentration variables in \code{data}
#' @param  FN_crit numeric vector of values between 0 and 1 indicating false negative threshold for floating percentile model benchmark selection (default = \code{seq(0.1, 0.9, 0.05)})
#' @param  k numeric value >1 indicating the number of k-folds to include in cross-validation (default = \code{5})
#' @param  seed random seed to set for reproducible results (default = \code{NULL}, i.e., no seed)
#' @param  plot whether to plot the output of \code{cvFPM} (default = \code{TRUE})
#' @param  tryStop specifies the number of times the cross-validation algorithm will try to run before ending (see Details; default = \code{10})
#' @param  simplify logical; whether to return just the optimized \code{FN_crit} value or more detailed diagnostic information
#' @param  which numeric or character indicating which type of plot to generate (see Details; default = \code{c(1, 2)})
#' @param  ... additional arguments passed to \code{chemSig} and \code{FPM}
#' @details 
#' \code{cvFPM} allows users to "tune" \code{FN_crit} argument in \code{FPM}. This is achieved by splitting the empirical dataset into "test" and "training" subsets,
#' calculating benchmarks for the training set, and then calculating the benchmarks' prediction errors using the test set.
#' This process is repeated several times (the number depending on the size of \code{k} relative to the sample size), and the results are summarized statistically. 
#' Lastly, this process is repeated for each \code{FN_crit} value specified by the user, resulting in comparable
#' statistics for each \code{FN_crit}. The output in the console indicates which \code{FN_crit} value resulted in the consistently optimal benchmarks 
#' (meaning highest overall reliability or most balanced errors). 
#' By setting \code{plot = TRUE} (the default), the outcome of cross-validation can be visualized over the range of \code{FN_crit} values considered. Visualizing the results
#' can inform the user about variability in the cross-validation process, ranges of potentially reasonable \code{FN_crit} values, etc.
#' 
#' \code{cvFPM} does not currently support optimization of the \code{alpha} parameter of \code{FPM}; 
#' \code{optimFPM} allows the user to optimize \code{alpha} but only using the empirical data (not through cross-validation).
#' 
#' Errors may be encountered by setting the value of \code{k} too high or too low, resulting in an inability of \code{cvFPM} 
#' to generate meaningful subsets for testing and floating percentile model calculations. Groups for subsetting are roughly
#' evenly applied within the cross-validation method, so it is reasonable to expect that \code{ceiling(nrow(data)/k)} is the number of
#' samples in any given test subset, with \code{nrow(data) - ceiling(nrow(data)/k)} being the size of the training subset. If a large number
#' of samples still generates an error, consider increasing the \code{tryStop} value and rerunning \code{cvFPM}. The easiest way
#' to avoid this type of error is to keep \code{k} low relative to \code{nrow(data)} (bearing in mind that \code{k} must be >1).
#' 
#' The \code{which} argument can be used to specify which of the two plots should be generated when \code{plot = TRUE}. These plots include
#' the optimization results based on the overall reliability metric or the balanced rate of false positives and false negatives. Inputs
#' to \code{which} are, by default, \code{c(1, 2)}, but flexible character inputs also can be used, for example \code{which = "OR"} or \code{which = "balanced"}.
#' 
#' @seealso chemSig, FPM, seq
#' @return list with 1 or 3 objects (depending on whether or not \code{simplify = TRUE});
#' these include 1) \code{optim_FN}, the optimized \code{FN_crit} value; 2) \code{CV_OR}, the detailed breakdown of overall
#' reliability values for each \code{FN_crit} value; and 3) \code{CV_FPM}, floating percentile model benchark statistics based on
#' all cross-validation runs for each \code{FN_crit}
#' @importFrom stats median
#' @importFrom graphics par
#' @examples
#' paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn")
#' cvFPM(data = h.tristate, paramList = paramList, FN_crit = seq(0.1, 0.9, 0.1), which = "OR")
#' @export
cvFPM <- function(data,
                  paramList,
                  FN_crit = seq(0.1, 0.9, 0.05),
                  k = 5,
                  seed = NULL,
                  plot = TRUE,
                  tryStop = 10,
                  simplify = TRUE,
                  which = c(1, 2),
                   ...){
    if(any(FN_crit > 1)){
        stop("FN_crit must include numbers between 0 and 1")
    }
    if(k <= 0){
        stop("k must be a positive number between 1 and nrow(data)")
    }
    if(!is.null(seed) & !is.numeric(seed)){
        stop("seed must be a positive number")
    }

    if(!is.null(seed)){
        set.seed(seed)
    }
    
    hitRatio <- mean(data$Hit)
    foldN <- ceiling(nrow(data)/k)
    
    foldNHit <- ceiling(hitRatio * foldN)
    foldNNoHit <- ceiling((1 - hitRatio) * foldN)
    

    data$fold <- NA 
    
    tryFold <- 0
    repeat{
        tryFold <- tryFold + 1
        data[data$Hit,]$fold <- sample(rep(1:k, each = foldNHit), 
                                       size = nrow(data[data$Hit,]), replace = FALSE)
        data[!data$Hit,]$fold <- sample(rep(1:k, each = foldNNoHit), 
                                       size = nrow(data[!data$Hit,]), replace = FALSE)
        
        check <- table(data$fold, data$Hit)
        
        if(!any(apply(check, 2, function(x) {
                any(x == 0)})) & 
           nrow(check) == k){
            break
        }
        
        if(tryFold == tryStop){
            stop("Reduce k and try again")
        }
    }
    
    suppressMessages(expr = {
        paramListFix <- paramList[chemSig(data = data, paramList = paramList, ...)]

    })
    
    sims <- list()
    fpms <- list()
    for(FN in 1:length(FN_crit)){
        initial <- FPM(data, paramList = paramList, FN_crit = FN_crit[FN], ...)[["FPM"]]

        OR.k <- numeric()
        Balance.k <- numeric()
        fpms[[FN]] <- list()
        
        for(k in 1:k){
            kTrain <- data[data$fold != k,]
            kTest <- data[data$fold == k,]
            kTestHit <- kTest[kTest$Hit, ]
            kTestNoHit <- kTest[!kTest$Hit, ]
            
            fpm.k <- FPM(kTrain, paramList = paramListFix,
                FN_crit = FN_crit[FN], paramOverride = TRUE, ...)[["FPM"]][paramListFix]

            fpms[[FN]][[k]] <- fpm.k
            
            TP.k <- apply(as.matrix(apply(as.matrix(fpm.k, ncol = length(paramListFix)), 1, function(x) {
                        apply(as.matrix(kTestHit[, paramListFix], ncol = length(paramListFix)), 1, function(y) {
                            any(y > x)})})), 2, sum)
            
            TN.k <- apply(as.matrix(apply(as.matrix(fpm.k, ncol = length(paramListFix)), 1, function(x) {
                        apply(as.matrix(kTestNoHit[, paramListFix]), 1, function(y) {
                            all(y <= x)})})), 2, sum)
            
            FP.k <- apply(as.matrix(apply(as.matrix(fpm.k, ncol = length(paramListFix)), 1, function(x) {
                        apply(as.matrix(kTestNoHit[, paramListFix]), 1, function(y) {
                            any(y > x)})})), 2, sum)
            FN.k <- apply(as.matrix(apply(as.matrix(fpm.k, ncol = length(paramListFix)), 1, function(x) {
                        apply(as.matrix(kTestHit[, paramListFix]), 1, function(y) {
                            all(y <= x)})})), 2, sum)
            
            Balance.k <- c(Balance.k, abs(FP.k - FN.k))/nrow(kTest)
            OR.k <- c(OR.k, (TP.k + TN.k)/nrow(kTest))
        }
        
        fpms[[FN]] <- do.call(rbind, list(mean = apply(do.call(rbind, fpms[[FN]]), 2, mean),
                            median = apply(do.call(rbind, fpms[[FN]]), 2, median),
                            min = apply(do.call(rbind, fpms[[FN]]), 2, min),
                            max = apply(do.call(rbind, fpms[[FN]]), 2, max)))
        names(fpms)[FN] <- paste0("FN_", FN_crit[FN])
        sims[[FN]] <- data.frame(Bal.i = signif(abs(initial[, "FP"] - initial[, "FN"])/sum(initial[,c("TP", "TN", "FP", "FN")]), 2),
                                 Bal.kmin = signif(min(Balance.k), 2),
                                 Bal.kmax = signif(max(Balance.k), 2),
                                 Bal.kmean = signif(mean(Balance.k), 2),
                                 Bal.kmed = signif(median(Balance.k), 2),
                                 OR.i = signif(initial[, "OR"], 2),
                                 OR.kmin = signif(min(OR.k), 2),
                                 OR.kmax = signif(max(OR.k), 2),
                                 OR.kmean = signif(mean(OR.k), 2),
                                 OR.kmed = signif(median(OR.k), 2))
    }
    out <- data.frame(FN_crit = FN_crit, do.call(rbind, sims))

    out_select <- c(
        balance_FN.FP = out$FN_crit[min(which.min(rowSums(out[, c("Bal.kmin", 
                                                            "Bal.kmax", 
                                                            "Bal.kmean", 
                                                            "Bal.kmed")])))],
        max_OR = out$FN_crit[min(which.max(rowSums(out[, c("OR.kmin", 
                                                            "OR.kmax", 
                                                            "OR.kmean", 
                                                            "OR.kmed")])))]
        )
    
    if(plot){
        # prevent cvFPM from overwriting user parameters
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(par(oldpar))
        
        if(!any(which %in% c(1,2)) & 
           !any(which %in% c("OR", "or", "reliability", "Reliability", "rel", "Rel",
                             "Bal", "bal", "balance", "balanced", "Balance", "Balanced"))){
            stop("which argument must equal 1, 2, or c(1,2)"
            )}
        
        if(any(c(1, "OR", "or", "reliability", "Reliability", "rel", "Rel") %in% which)){
            with(out, expr = {
             par(oma = c(0, 0, 3, 0))
             plot(x = FN_crit, y = OR.i, type = "l", col = 1, lwd = 2,
                  xlab = "False Negative Limit", ylab = "Overall Reliability",
                  ylim = c(min(OR.kmin), max(OR.kmax)))
             abline(v = out_select["max_OR"], col = "darkorange", lty = 2, lwd = 2)
             lines(x = FN_crit, y = OR.kmean, col = 2)
             lines(x = FN_crit, y = OR.kmed, col = 4)
             lines(x = FN_crit, y = OR.kmin, col = 4, lty = 3)
             lines(x = FN_crit, y = OR.kmax, col = 4, lty = 3)
             legend("top", inset = -0.15, bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                    lty = c(1,1,1,3,2), col = c(1,2,4,4,"darkorange"), lwd = c(2,1,1,1,2),
                    legend = c("Actual OR", "Mean test OR", "Median test OR", "Min/Max test OR", "Optimal FN"))
            })
        }
        if(any(c(2, "Bal", "bal", "balance", "balanced", "Balance", "Balanced") %in% which)){
            with(out, expr = {
             par(oma = c(0, 0, 3, 0))
             plot(x = FN_crit, y = Bal.i, type = "l", col = 1, lwd = 2,
                  xlab = "False Negative Limit", ylab = "Balance between FP and FN",
                  ylim = c(min(Bal.kmin), max(Bal.kmax)))
             abline(h = 0, lty = 3, col = "gray")
             abline(v = out_select["balance_FN.FP"], col = "darkorange", lty = 2, lwd = 2)
             lines(x = FN_crit, y = Bal.kmean, col = 2)
             lines(x = FN_crit, y = Bal.kmed, col = 4)
             lines(x = FN_crit, y = Bal.kmin, col = 4, lty = 3)
             lines(x = FN_crit, y = Bal.kmax, col = 4, lty = 3)
             legend("top", inset = -0.15, bg = "white", ncol = 2, cex = 0.8, x.intersp = 0.25, xpd = TRUE,
                    lty = c(1,1,1,3,2), col = c(1,2,4,4,"darkorange"), lwd = c(2,1,1,1,2),
                    legend = c("Actual balance", "Mean test balance", "Median test balance", "Min/Max test balance", "Optimized FN_crit"))
            })
        }
    }
    
    if(simplify){
        return(list(optim_FN = out_select))
    } else {
        return(list(CV_Stat = out, CV_FPM = fpms, optim_FN = out_select))
    }
}## end code
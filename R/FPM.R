#' Floating Percentile Model
#' 
#' Generate sediment quality benchmarks using the floating percentile model algorithm
#' 
#' @param  data data.frame containing, at a minimum, chemical concentrations as columns and a logical \code{Hit} column classifying toxicity
#' @param  paramList character vector of column names of chemical concentration variables in \code{data}
#' @param  paramFixed character vector of column names of chemical concentration variables to retain, bypassing testing for specific chemicals (default = \code{NULL}). See Details.
#' @param  paramOverride logical; whether to retain every chemical variable in \code{paramList} (default = \code{FALSE}). See Details.
#' @param  increment numeric value greater than 1; number of increments to evaluate (default = \code{10}). See Details.
#' @param  precision numeric value between 0 and 1 (default = \code{0.1})
#' @param  FN_crit numeric vector of values between 0 and 1 indicating false negative threshold(s) for benchmark selection (default = \code{0.2})
#' @param  seed random seed to set for reproducible results; only for handling edge cases of ranking ties (default = \code{1})
#' @param  lockInfo logical; whether to return the reason for and order in which benchmarks were "locked" within the model algorithm (default = \code{FALSE}). See Details.
#' @param  ... additional argument passed to \code{chemSigSelect} and \code{chemSig}
#' @details \code{FPM} is the main function provided in 'RFPM', which was developed firstly as a redevelopment of the Washington Department of Ecology's Excel-based
#' floating percentile model tool (Avocet 2003; Ecology 2011), and secondly as a means to evaluate uncertainties and sensitivities associated with the model. \code{FPM} generates 
#' sediment quality benchmarks for chemicals with significantly higher concentrations among \code{Hit} samples (meaning they were determined to be categorically toxic).
#' 
#' \code{FPM} is an algorithmic approach to setting sediment quality benchmarks using sediment chemistry data and toxicity test results.
#' Toxicity is treated as a binary classification - either a \code{Hit == TRUE} or \code{Hit == FALSE} (meaning toxic or non-toxic) by some user-defined definition.
#' The most important input to \code{FPM} apart from the empirical data is \code{FN_crit}, which determines an upper limit for false negative errors associated with floating percentile model benchmarks.
#' The default \code{FN_crit} recommended by the Department of Ecology is \code{0.2}; though intended to be protective, the value of \code{0.2} is arbitrary. We recommend
#' that the user run the \code{optimFPM} and/or \code{cvFPM} functions to find the \code{FN_crit} value(s) that optimize benchmark performance within an acceptable error range for the site.
#' \code{optimFPM} can also help users optimize the \code{alpha} parameter (see \code{?chemSig}), which is also somewhat arbitrarily set at a conventional default of \code{0.05}.
#' 
#' There are two arguments that have defaults in \code{FPM} that the user may desire to change in certain circumstances, but that we
#' generally recommend not changing without good reason. These are \code{paramFixed} and \code{paramOverride}, which override the chemical 
#' selection process, resulting in potentially non-toxic chemicals being assigned benchmarks.
#' The \code{paramFixed} argument, which only forces named chemicals into the model algorithm, is looser than \code{paramOverride}, which forces all chemicals in paramList into the model algorithm.
#' See \code{?chemSig} for more information regarding default parameters used within \code{FPM}.
#' Even if chemical names are supplied to \code{paramFixed}, \code{FPM} will still use hypothesis testing methods to consider all other chemicals for inclusion.
#' 
#' \code{increment} determines (inversely) how large or small values should be that are added to percentile values in the model algorithm.
#' A larger \code{increment} results in smaller incremental additions and vice-versa. The WA Department of Ecology recommends a default of
#' \code{increment = 10}. This is a reasonable value, and we recommend not decreasing \code{increment} below \code{10}. 
#' Increasing \code{increment} will increase computation time, and may or may not result in more accurate benchmarks. So, we recommend not
#' increasing \code{increment} much higher than \code{10}.
#' 
#' \code{precision} determines how many iterative loops will be attempted within the model algorithm when trying to increase each benchmark. If increasing the
#' benchmark would increase the false negative rate above \code{FN_crit}, the benchmark would then be decreased, the increment size is divided by \code{increment}, and
#' then the smaller incremental addition is used to increase the benchmark. This process repeats for a fixed number of iterations, which is related to \code{precision}.
#' If the benchmark cannot be increased after the fixed number of iterations, the benchmark is locked in place.
#' The default value for \code{precision} is \code{0.1}, but the value could be lower, if desired. Lowering the value will increase computation time and may or may not
#' result in more accurate benchmarks. In general, we recommend reducing \code{precision} rather than increasing \code{increment} in order to potentially enhance the
#' precision of benchmark calculations.
#' 
#' The \code{lockInfo} argument allows the user to export information about what caused the model algorithm to lock for each
#' chemical. Output options are: \code{"FN"} for exceeding the false negative limit (i.e., \code{FN_crit}), \code{"FP"} if the number of false positives was reduced to zero, 
#' \code{"Max"} if the empirical maximum concentration was exceeded, or \code{Mix} if more than one of the first three options occurred.
#' 
#' The following classification statistics are reported alongside the generated benchmarks:
#' \code{TP}, \code{FN}, \code{TN}, and \code{FP} - the numbers of true positive, false negative, true negative, and false positive predictions
#' \code{pFN} and \code{pFP} - proportions of false predictions (false No-hit and false Hit, respectively)
#' \code{HR} - hit rate or sensitivity; the probability of detecting a Hit
#' \code{NR} - negative rate or specificity; the probability of detecting a No-hit
#' \code{PHR} - positive predictive value or precision; how often were Hit predictions correct?
#' \code{PNR} - negative predictive value; how often were No-hit predictions correct?
#' \code{FHR} - false positive predictive value; how often were Hit predictions incorrect?
#' \code{FNR} - false negative predictive value; how often were No-hit predictions incorrect?
#' \code{OR} - overall reliability; the probability of making a correct prediction (Hit or No-hit)
#' 
#' The second output of \code{FPM} is a metric called \code{chemDensity}. This is a measure of how much
#' the percentile "floated" in the algorithm from the starting position up to the chemical's value at which it was locked in place.
#' Values of \code{chemDensity} closer to 1 floated less and vice-versa. By floating less, this indicates that
#' even small changes in the chemical concentration resulted in one of the acceptance criteria failing (as discussed above with regard to \code{lockInfo}). When comparing
#' the \code{chemDensity} among chemicals, those with lower values might be viewed as having less of an influence on toxicity predictions
#' and vice-versa. For those interested in understanding the relative importance of chemicals among benchmarks, we recommend using \code{chemVI} and considering the \code{MADP} and \code{dOR} outputs.
#' 
#' @seealso optimFPM, cvFPM, chemSig, chemSigSelect, chemVI
#' @return list of 2 or 4 objects (depending on \code{lockInfo}): 
#' 1) Benchmarks and toxicity classification error statistics; 
#' 2) order in which benchmarks were locked in place; 
#' 3) reason for benchmarks being locked in place; and 
#' 4) \code{chemDensity} statistic
#' @importFrom stats quantile
#' @references
#' Avocet. 2003. Development of freshwater sediment quality values for use in Washington State. Phase II report: Development and recommendation of SQVs for freshwater sediments in Washington State. Publication No. 03-09-088. Prepared for Washington Department of Ecology. Avocet Consulting, Kenmore, WA.
#' Ecology. 2011. Development of benthic SQVs for freshwater sediments in Washington, Oregon, and Idaho. Publication no. 11-09-054. Toxics Cleanup Program, Washington State Department of Ecology, Olympia, WA.
#' @examples
#' paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn")
#' FPM(data = h.tristate, paramList = paramList, ExcelMode = TRUE, warn = FALSE)
#' FPM(data = h.tristate, paramList = paramList, FN_crit = c(0.1, 0.2, 0.3))
#' @export
FPM <- function(data, 
                paramList, 
                paramFixed = NULL,
                paramOverride = FALSE,
                increment = 10, 
                precision = 0.1,
                FN_crit = 0.2, 
                seed = 1, 
                lockInfo = FALSE,
                ...) {

    if(nrow(data) == 0){ 
        stop("No data supplied to FPM") 
    }
    if(!is.logical(data$Hit)){
        stop("Hit must be provided as a logical vector (TRUE/FALSE)")
    }
    
    if(increment <= 1) {
        stop("increment must be greater than 1")
    }
    
    if(precision <= 0){
        stop("precision must be greater than 0")    
    }
    
    if(any(FN_crit > 1 | FN_crit < 0)){
        stop("FN_crit must be a number between 0 and 1")
    } else if(any(FN_crit >= 1 | FN_crit <= 0)){
        warning("FN_crit values equal to 0 or 1 may not provide useful results (over/underconservative)")   
    }
    
    FN_crit_list <- list()
    lockList <- list()
    lockReasonList <- list()
    densityList <- list()
    
    for(FN in 1:length(FN_crit)){
        FN_crit.i <- FN_crit[FN]
        
        ## parameter selection ---
        if(!paramOverride){
            if(is.null(paramFixed)){
                    dNew <- chemSigSelect(data = data, paramList = paramList, ...)[[1]]
                    paramListNew <- names(dNew[-length(dNew)]) # the 'length(dNew)' position is occupied by "Hit" column
                    paramAlpha <- order(paramListNew) ## alphabetical order of inputs (index is left to right as supplied)
            } else if(is.character(paramFixed) & length(paramFixed) < length(paramList)){
                dNew <- data.frame(data[, paramFixed],  # fixed parameters moved to front of data.frame
                    chemSigSelect(data = data, paramList = paramList[!(paramList %in% paramFixed)], ...)[[1]])

                names(dNew)[1:length(paramFixed)] <- paramFixed # add flag to parameter to show it was fixed
                paramListNew <- names(dNew[-length(dNew)])
                paramAlpha <- order(paramListNew) ## alphabetical order of inputs (index is left to right as supplied)
            } else if(!is.character(paramFixed)){
                stop("paramFixed must be a character vector and subset of paramList")
            } else if(length(paramFixed) == length(paramList) & all(paramFixed %in% paramList)){
                dNew <- data[,c(paramFixed, "Hit")]
                paramListNew <- names(dNew[-length(dNew)])
                paramAlpha <- order(paramListNew) ## alphabetical order of inputs (index is left to right as supplied)a
            } else { 
                stop("paramFixed must be a character vector and subset of paramList") 
            }
        } else { ## option to override parameter selection entirely - ONLY look at fixed paramList
            if(!is.null(paramFixed)){
                warning("paramFixed is ignored when paramOverride is TRUE; using paramList instead")
            }
            dNew <- data.frame(data[, paramList], Hit = data$Hit) # keep all data when paramOverride = TRUE
            paramListNew <- paramList
            paramAlpha <- order(paramListNew)
            }
        
    
        ## dataset prep (split by hit/no-hit class, then exclude class)
        nohit <- as.matrix(dNew[!dNew$Hit, -length(dNew)])
        hit <- as.matrix(dNew[dNew$Hit, -length(dNew)])
        allData <- as.matrix(dNew[, -length(dNew)])
        
        

        ## accounting for shift in matrix format (column to vector) when only 1 significant parameter selected
        if(ncol(nohit) == 1) {
            colnames(nohit) <- paramListNew
            colnames(hit) <- paramListNew
            colnames(allData) <- paramListNew
        } else{
            colnames(nohit) <- paramListNew
            colnames(hit) <- paramListNew
            colnames(allData) <- paramListNew
        }
        
        # sample sizes
        nNoHit <- nrow(nohit)
        nHit <- nrow(hit)
        nAllData <- nrow(allData)
        
        ## parameter prep/initialization
        minVal <- matrix(
            do.call(cbind, lapply(as.data.frame(allData), FUN = function(x) {
                min(x, na.rm = TRUE)
                })),
            ncol = length(paramListNew))
        maxVal <- matrix(
            do.call(cbind, lapply(as.data.frame(allData), FUN = function(x) {
                max(x, na.rm = TRUE)
                })),
            ncol = length(paramListNew))
        colnames(minVal) <- paramListNew
        colnames(maxVal) <- paramListNew
        
        # establish "lift" - the size of incremental additions
        lifti <- matrix(
            apply(as.matrix(rbind(minVal, maxVal)), MARGIN = 2, FUN = function(x) {
                   diff(x)/(increment)
               }), 
            byrow = TRUE, nrow = 1) # initial lift to update in loop
        colnames(lifti) <- paramListNew
        
        
        liftf <- lifti ## initializing loop update to lift
        
        ## Limit the number of iterative loops
        estIter <- floor(log10((maxVal - minVal) / (precision * minVal)) / log10(increment)) # taken from user manual
        iter <- estIter # countdown tracker by chemical (lock when iter==0 & FN>FN_crit.i | FP==0 | Cf > maxVal)
        
        if(any(estIter < 1)){
          stop(paste("precision is set too high for one or more analytes:", paste(
                      colnames(estIter)[which(estIter<1)], collapse = ", ")))
        }
        ## setting up data
        
        P <- data.frame(P = (1:100)/100, 
                        do.call(cbind, lapply(as.data.frame(allData), 
                                              quantile, prob = (1:100)/100, type = 7, na.rm = TRUE)))
         
        null = rep(NA, nrow(P))
        E <- data.frame(P = (1:100)/100,
                        TP = null,
                        FN = null,
                        TN = null,
                        FP = null,
                        pFN = null,
                        pFP = null,
                        HR = null,
                        NR = null,
                        PHR = null,
                        PNR = null,
                        FHR = null,
                        FNR = null,
                        OR = null)
        
        ## initial error calculations
        E$TP <- apply(apply(P[2:length(P)], 1, function(x) {
                    apply(hit, 1, function(y) {
                        any(y > x)})}), 2, sum)
        E$FN <- apply(apply(P[2:length(P)], 1, function(x) {
                    apply(hit, 1, function(y) {
                        all(y <= x)})}), 2, sum)
        E$TN <- apply(apply(P[2:length(P)], 1, function(x) {
                    apply(nohit, 1, function(y) {
                        all(y <= x)})}), 2, sum)
        E$FP <- apply(apply(P[2:length(P)], 1, function(x) {
                    apply(nohit, 1, function(y) {
                        any(y > x)})}), 2, sum)
        E$pFN <- with(E, FN / nHit)
        E$pFP <- with(E, FP / nNoHit)
        E$HR <- with(E, 1 - pFN)
        E$NR <- with(E, 1 - pFP)
        E$PHR <- with(E, TP/(TP + FP))
        E$PNR <- with(E, TN/(TN + FN))
        E$FHR <- with(E, 1 - PHR)
        E$FNR <- with(E, 1 - PNR)
        E$OR <- with(E, (TP + TN)/nAllData)
        
        ## Initialize output
        O <- E[1,]
        names(O)[1] <- "FN_crit"
        O[1, 1] <- FN_crit.i
        
        
        ## find initial index for percentile, reduce if necessary to be below FN_crit.i, and stop at index == 1 (minimum)
        index <- which(abs(E$pFN - FN_crit.i) == min(abs(E$pFN - FN_crit.i)))
        
        # take highest percentile among tied values
        if(length(index) > 1) { 
            index <- max(index) 
        } 
        
        repeat{
            if(E$pFN[index] <= FN_crit.i | index == 1){
                break
            } else {
                index <- index - 1
                if(index == 1){
                    warning("percentile index set at minimum; check data manually for FN rates")
                }
            }
        }
        
        ## starting concentrations
        Ci <- matrix(P[index, 2:length(P)], byrow = TRUE, nrow = 1)
        Cf <- matrix(P[index, 2:length(P)], byrow = TRUE, nrow = 1) ## initializing the next step concentrations
        Co <- matrix(P[index, 2:length(P)], byrow = TRUE, nrow = 1) ## initializing output vector
        colnames(Ci) <- paramListNew
        colnames(Cf) <- paramListNew
        colnames(Co) <- paramListNew
        
        
        lock <- matrix(rep(NA, length(paramListNew)), byrow = TRUE, nrow = 1)
        colnames(lock) <- paramListNew

        lockReason <- matrix(rep("", length(paramListNew)), byrow = TRUE, nrow = 1)
        colnames(lockReason) <- paramListNew
        
        ## FP by chemical (initial version - will be recreated each loop)
        FPi <- do.call(cbind, lapply(1:length(Ci), function(col){
                    sum(nohit[,col] > Ci[[col]])
                }))
        colnames(FPi) <- paramListNew
        FPf <- FPi
        
        # START MAIN LOOP ------------------------------------------------------------------------------------------------------------------------
        
        lockN <- 0
        # while(!all(as.logical(lock))){ 
        while(any(is.na(lock))){ 
            # while-loop driven by the lock variable - all chems must be locked to stop the algorithm 
            ## FP by chemical before increment
            
            FPi <- do.call(cbind, lapply(1:length(Ci), function(col){
                        sum(nohit[,col] > Ci[[col]])}))
            
            # Selecting chemical to increase concentration using FP rate
            ## Follows stepwise selection process:
            ## 1. Highest FP (FPi)
            ## 2. Lowest number of iterations conducted (estIter - iter)
            ## 3. Lowest benchmark concentration (Ci)
            chemrank <- matrix(rank(FPi), byrow = TRUE, nrow = 1) # (Step 1); gives non-integer values when ranks tied 
            colnames(chemrank) <- paramListNew
            
            if(max(chemrank) %% 1 != 0){ # if the top rank is tied b/w chems (fraction), then re-evaluate
                maxchemrankindex <- which(chemrank == max(chemrank)) #index for tied chems
                iterRank <- rank(-(estIter[,maxchemrankindex] - iter[,maxchemrankindex])) # negative->least number of iterations goes first
                
                if(max(iterRank) %% 1 != 0){ # if still tied, 
                    maxiterRankindex <- which(iterRank == max(iterRank))
                    set.seed(seed)
                    chemrank[maxchemrankindex][maxiterRankindex] <- length(chemrank) - (rank(unlist(Ci[,maxchemrankindex][maxiterRankindex]), 
                                                                         ties.method = "random") - 1)
                } else {
                    chemrank[maxchemrankindex] <- length(chemrank) # set to max possible rank (Step 2) (maybe arbitrary)
                }
            }
            
            # chem <- paramListNew[!lock][which.max(chemrank[!lock])] # max FP among unlocked chems
            chem <- paramListNew[is.na(lock)][which.max(chemrank[is.na(lock)])] # max FP among unlocked chems

            # increase concentration by lift (updated later)
            Cf[, chem] <- as.numeric(Ci[, chem]) + as.numeric(liftf[, chem])
            # recalculate FP and FN to check if increase was acceptable
            
            # Acceptability is based on three criteria:
            ## 1. pFNf (false negative rate) is less than or equal to the threshold FN_crit.i)
            ## 2. FPf (false positive rate) is greater than 0
            ## 3. New concentration (Cf) is less than or equal to the maximum empirical concentration
            ## ** Note that FPf==0 does not necessarily disqualify a new value; it is accepted if #1,#3 not also triggered.
            
            FPf <- do.call(cbind, lapply(1:length(Cf), function(col){
                            sum(nohit[,col] > Cf[[col]])})) ## recalculate FP
            colnames(FPf) <- paramListNew
            
            TPf <- apply(apply(Cf, 1, function(x) {
                        apply(hit, 1, function(y) {
                            any(y > x)})}), 2, sum)
            FNf <- apply(apply(Cf, 1, function(x) {
                        apply(hit, 1, function(y) {
                            all(y <= x)})}), 2, sum)
            pFNf <- FNf / (TPf + FNf) # recalculate overall pFN
            
            # check acceptability of increase
            if(pFNf > FN_crit.i | 
               FPf[,chem] == 0 |
               (pFNf <= FN_crit.i & FPf[,chem] > 0 & Cf[, chem] > maxVal[, chem])){
                repeat{
                    
                    Cf[, chem] <- as.numeric(Ci[, chem]) + as.numeric(liftf[, chem]) 
                    
                    FPf <- do.call(cbind, lapply(1:length(Cf), function(col){
                                    sum(nohit[, col] > Cf[[col]])}))
                    colnames(FPf) <- colnames(Cf)
                    
                    TPf <- apply(apply(Cf, 1, function(x) {
                                apply(hit, 1, function(y) {
                                    any(y > x)})}), 2, sum)
                    
                    FNf <- apply(apply(Cf, 1, function(x) {
                                apply(hit, 1, function(y) {
                                    all(y <= x)})}), 2, sum)
                    pFNf <- FNf / (TPf + FNf)
                    
                    
                    if((pFNf > FN_crit.i | 
                        FPf[, chem] == 0 |
                        (pFNf <= FN_crit.i & FPf[,chem] > 0 & Cf[, chem] > maxVal[, chem])) & 
                       iter[,chem] == 1){ # any trigger condition met and no more iterations allowed # CD: Changed 0 to 1 to make 4 steps
                        if(FPf[, chem] == 0 & pFNf <= FN_crit.i & Cf[, chem] <= maxVal[,chem]){ # FP as only trigger
                            Ci[, chem] <- Cf[, chem] # Allow final increment when FP is reason for locking
                            lockReason[, chem] <- "FP"
                        } else if (FPf[,chem] > 0 & pFNf > FN_crit.i & Cf[, chem] <= maxVal[, chem]){
                            # no change to Ci for FN lock (unacceptable FN error) 
                            lockReason[, chem] <- "FN"
                            Cf[, chem] <- Ci[, chem]
                        } else if (FPf[,chem] > 0 & pFNf <= FN_crit.i & Cf[, chem] > maxVal[, chem]){
                            # no change to Ci for Max lock (outside empirical data)
                            lockReason[, chem] <- "Max"   
                        } else {
                            lockReason[, chem] <- "Mix" # multiple reasons; no change allowed
                        }
                        lock[, chem] <- lockN + 1
                        lockN <- lockN + 1
                        break # no more iterations allowed; lock chemical and exit loop
                    } else if((pFNf > FN_crit.i | 
                               FPf[, chem] == 0) &
                               iter[, chem] > 1){ # any trigger condition met and iterations still allowed
                        Cf[, chem] <- Ci[, chem] # step back to previous value
                        iter[,chem] <- iter[, chem] - 1 # reduce iterations allowed
                        liftf[, chem] <- liftf[, chem]/increment # shrink lift and repeat loop 
                    } else {
                        Ci[, chem] <- Cf[, chem] # acceptable increase; exit repeat-loop (back to while-loop for next increment)
                        break
                    }
                } ## end repeat-loop for reducing increments
            } else {
                Ci[, chem] <- Cf[, chem]
            }## else condition is when increase was acceptable - no iteration needed (remain in while-loop)
        } # end  while-loop for incrementing up FPM benchmark
        ## --------------------------------------------------------------------------------------------------------------------------
        
        # select empirical concentration closest to calculated Ci
        Co <- as.data.frame(lapply(colnames(Ci), function(col){
                temp <- as.numeric(Ci[, col]) - as.numeric(allData[, col])
                diffPick <- min(temp[temp >= 0])
                return(as.numeric(Ci[, col]) - as.numeric(diffPick))
            }))
        colnames(Co) <- colnames(Ci)
        
        
        # calculate final error calcs and append to final FPM SQBs
        O$TP <- apply(apply(Ci, 1, function(x) {
                    apply(hit, 1, function(y) {
                        any(y > x)})}), 2, sum)
        O$FN <- apply(apply(Ci, 1, function(x) {
                    apply(hit, 1, function(y) {
                        all(y <= x)})}), 2, sum)
        O$TN <- apply(apply(Ci, 1, function(x) {
                    apply(nohit, 1, function(y) {
                        all(y <= x)})}), 2, sum)
        O$FP <- apply(apply(Ci, 1, function(x) {
                    apply(nohit, 1, function(y) {
                        any(y > x)})}), 2, sum)
        O$pFN <- with(O, FN / (TP + FN))
        O$pFP <- with(O, FP / (TN + FP))
        O$HR <- with(O, 1 - pFN)
        O$NR <- with(O, 1 - pFP)
        O$PHR <- with(O, TP/(TP + FP))
        O$PNR <- with(O, TN/(TN + FN))
        O$FHR <- with(O, 1 - PHR)
        O$FNR <- with(O, 1 - PNR)
        O$OR <- with(O, (TP + TN)/nAllData)

        # chemDensity
        dens <- data.frame(
            lapply(paramListNew, function(x) {
                start <- quantile(dNew[, x], p = index/100, na.rm = TRUE)[[1]]
                end <- unlist(Ci[, x])
                return(
                    1 - (end - start)/(max(dNew[, x]) - start)
                )
            }))
        colnames(dens) <- paramListNew
        
        out <- data.frame(Co, signif(O, 3))
        attr(out, "liftf") <- as.data.frame(liftf)
        FN_crit_list[[FN]] <- out
        lockList[[FN]] <- lock
        lockReasonList[[FN]] <- lockReason
        densityList[[FN]] <- dens
    }
    
    # combine results, remove messy rownames, and  return final results
    tmp <- do.call(rbind, FN_crit_list)
    lockTmp <- do.call(rbind, lockList)
    lockReasonTmp <- do.call(rbind, lockReasonList)
    densityTmp <- do.call(rbind, densityList)
    rownames(densityTmp) <- NULL
    rownames(tmp) <- NULL
    rownames(lockTmp) <- NULL
    rownames(lockReasonTmp) <- NULL
    rownames(densityList) <- NULL
    
    if(lockInfo){
        return(list(FPM = tmp, 
                    chemDensity = signif(densityTmp, 3), 
                    lockReason = as.data.frame(lockReasonTmp), 
                    lockOrder = as.data.frame(lockTmp)))
    } else {
        return(list(FPM = tmp, 
                    chemDensity = signif(densityTmp, 3)))
    }
} ## end code

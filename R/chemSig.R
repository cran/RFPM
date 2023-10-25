#' Chemical Variable Selection within the Floating Percentile Model
#' 
#' Determine which chemicals in a dataset have significantly higher concentrations
#' among toxic samples by comparison to non-toxic samples
#' 
#' @param  data data.frame containing, at a minimum, chemical concentrations as columns and a logical \code{Hit} column classifying toxicity
#' @param  paramList character vector of column names of chemical concentration variables in \code{data}
#' @param  testType character string; whether to run parametric or non-parametric tests (default = \code{NULL}). See Details for more information.
#' @param  alpha numeric value between 0 and 1; type-I error rate that, if supplied by the user, will be applied to all tests (default = \code{NULL})
#' @param  alpha.var numeric value between 0 and 1; type-I error rate for testing equal variance assumption (default = \code{0.05})
#' @param  alpha.norm  numeric value between 0 and 1; type-I error rate for testing normality assumption (default = \code{0.05})
#' @param  alpha.test  numeric value between 0 and 1; type-I error rate for testing differences between Hit and No-hit datasets (default = \code{0.05})
#' @param  alternative alternative hypothesis type for equality of central tendency (default = \code{"less"})
#' @param  var.alternative alternative hypothesis type for equal variance test (default = \code{"two.sided"})
#' @param  var.equal logical; whether to assume equal variance (default = \code{NULL})
#' @param  ExcelMode logical; whether to force \code{chemSig} to run like the WA Department of Ecology's Excel-based floating percentile model calculator (default = \code{FALSE})
#' @param  warn logical; whether to generate a warning associated with \code{ExcelMode} (default = \code{TRUE})
#' 
#' @details 
#' \code{chemSig} is called within \code{FPM} via \code{chemSigSelect}, which generates a subset of chemicals to
#' pass into the floating percentile model algorithm. \code{chemSig} only returns a logical vector describing which parameters in \code{paramList} should be selected for
#' benchmark development based on having significantly higher concentrations when \code{Hit == TRUE} than when \code{Hit == FALSE}.
#' The user has the ability to manipulate several of the parameters of the selection algorithm, or they can allow \code{chemSig} to test for
#' assumptions and use appropriate hypothesis tests based on those results. By default, \code{chemSig} will use \code{shapiro.test} to
#' confirm normality, then either \code{var.test} if the data are normal or \code{fligner.test}
#' if the data are non-normal to confirm equal variance. Finally, the function will use \code{t.test} if the data are normal (using the Welch method if unequal variance), 
#' \code{wilcox.test} if non-normal with equal variance, or \code{brunner.munzel.test} if non-normal with unequal variance.
#' 
#' The \code{testType} argument can be one of \code{p}, \code{P}, \code{param}, \code{Param}, \code{parametric}, or \code{Parametric} for parametric
#' test types or \code{non}, \code{Non}, \code{np}, \code{NP}, \code{nonparam}, \code{Nonparam}, \code{non-param}, \code{Non-param}, \code{nonparametric}, \code{Nonparametric}, 
#' \code{non-parametric}, \code{Non-parametric}, or \code{Non-Parametric}.
#' 
#' The user has the option of providing a single \code{alpha} level to apply to all tests (e.g., 0.05) or to specify test-specific alpha levels via \code{alpha.var}, \code{alpha.norm}, and \code{alpha.test}.
#' Note that \code{FPM} by default uses \code{alpha = 0.05} for all tests.
#' 
#' While \code{alternative} and \code{var.alternative} can be adjusted, we strongly recommend that they not be changed from the 
#' default values. For example, changing \code{alternative} from \code{"less"} (default) to \code{"two.sided"} would result in
#' the assumption that chemical concentrations could be significantly higher (as well as lower) when there is no toxicity than when there is, which is inappropriate.
#' The "greater" alternative, which is never appropriate, is not an accepted input for \code{alternative}. Similarly, the
#' assumption of equal variance relates to a \code{"two.sided"} argument, therefore changing the \code{var.alternative} to be
#' \code{"less"} or \code{"greater"} would not be appropriate. 
#' 
#' \code{ExcelMode} assumes \code{testType = "parametric"}, \code{var.equal = TRUE}, \code{alternative = "less"}, and \code{alpha = 0.1}. In
#' actuality, the Excel-based tool uses a one-way ANOVA test to compare two levels of \code{Hit}, which is equivalent to a t-test so long as alpha is adjusted to 0.1.
#' Thus, \code{testType}, \code{alternative}, \code{var.alternative}, and \code{var.equal} are overridden when \code{ExcelMode = TRUE}. 
#' This argument was included for those interested in using 'RFPM' as an alternative to the Excel-based calculator tool to obtain identical benchmark results. Note that 
#' the Excel FPM tool also includes other features which may complicate comparability such as outlier analysis. Outlier analysis is not conducted by RFPM.
#' @return named logical vector
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr all_of
#' @importFrom stats shapiro.test
#' @importFrom stats var.test
#' @importFrom stats fligner.test
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stats ks.test
#' @importFrom lawstat brunner.munzel.test
#' @importFrom dplyr `%>%`
#' @examples
#' paramList = c("Cd", "Cu", "Fe", "Mn", "Ni", "Pb", "Zn")
#' chemSig(h.tristate, paramList, "nonparametric")
#' chemSig(h.tristate, paramList, "parametric")
#' @export
chemSig <- function(data, 
                    paramList, 
                    testType = NULL, 
                    alpha.var = 0.05, 
                    alpha.norm = 0.05, 
                    alpha.test = 0.05, 
                    alpha = NULL, 
                    alternative = "less", 
                    var.alternative = "two.sided", 
                    var.equal = NULL, 
                    warn = TRUE, 
                    ExcelMode = NULL){
    
    if(alternative == "greater"){
        stop("alternative should be 'less' or 'two.sided'; 'greater' should not be used; 'less' is recommended")
    }
    
    if(is.null(ExcelMode)){
        ExcelMode <- FALSE
    } else if(ExcelMode == TRUE){
        testType = "Parametric"
        alpha <- 0.1
        alternative <- "two.sided"
        var.equal <- TRUE

        if(warn){
            message("ExcelMode overrode the following arguments: testType = 'Parametric', alpha = 0.1, alternative = 'two.sided'; var.equal = TRUE")
        }
    }
    
    # alpha argument overrides test-specific alphas
    if(!is.null(alpha)){
       alpha.var <- alpha
       alpha.norm <- alpha
       alpha.test <- alpha
    }
    
    
    if(any((alpha >= 1 | alpha <= 0 | 
       alpha.var >= 1 | alpha.var <= 0 | 
       alpha.norm >= 1 | alpha.norm <= 0 | 
       alpha.test >= 1 | alpha.test <= 0), na.rm = T)) {
        stop("alpha input must be a value >0 and <1")
    }
    
    
    ## message about how melt was previously silenced - supressing any messages still
    data2 <- suppressMessages(tidyr::pivot_longer(data %>% dplyr::select(all_of(c(paramList, "Hit"))), 
                                          cols = all_of(paramList), 
                                          values_drop_na = TRUE, 
                                          values_to = "Result", 
                                          names_to = "Param")
                             )
    
    out <- unlist(
        lapply(split(x = data2, f = data2$Param, drop = TRUE), function(x) {
            if(length(unique(x$Result)) > 1){ 
                if(is.null(testType)){
                        norm <- ifelse(length(unique(x[!x$Hit,]$Result)) == 1 |
                                length(unique(x[x$Hit,]$Result)) == 1, TRUE,
                            ifelse(shapiro.test(x[!x$Hit,]$Result)$p.value < alpha.norm |
                                shapiro.test(x[x$Hit,]$Result)$p.value < alpha.norm, FALSE, TRUE))
                                                                                                                                                                                                           
                        if(norm){
                            if(is.null(var.equal)){
                                var.equal <- var.test(x = x[!x$Hit,]$Result, 
                                            y = x[x$Hit,]$Result, 
                                            alternative = var.alternative)$p.value >= alpha.var 
                            }
                            
                            testSig <- t.test(x = x[!x$Hit,]$Result, 
                                              y = x[x$Hit,]$Result, 
                                              alternative = alternative, paired = FALSE, 
                                              var.equal = var.equal)$p.value < alpha.test
                        } else if(!norm){
                            if(is.null(var.equal)){
                                var.equal <- fligner.test(x = x$Result, g = x$Hit, 
                                    alternative = var.alternative)$p.value >= alpha.var
                            }
                                if(var.equal){
                                    testSig <- wilcox.test(x = x[!x$Hit,]$Result, 
                                        y = x[x$Hit,]$Result, 
                                        alternative = alternative, exact = FALSE)$p.value < alpha.test
                                } else if(!var.equal){
                                    testSig <- suppressWarnings(expr = {
                                        ks.test(x = x[x$Hit,]$Result, 
                                                y = x[!x$Hit,]$Result, 
                                                alternative = alternative, 
                                                exact = FALSE)$p.value < alpha.test
                                        })
                                }
                        }
                        return(testSig)
                        
             } else if(testType %in% c("Parametric", "Param", "P", "Par", 
                                       "parametric", "param", "p", "par")){
                if(is.null(var.equal)){
                    var.equal <- var.test(x = x[!x$Hit,]$Result, 
                        y = x[x$Hit,]$Result, 
                        alternative = var.alternative)$p.value >= alpha.var
                }
                 
                testSig <- t.test(x = x[!x$Hit,]$Result,
                    y = x[x$Hit,]$Result, 
                    alternative = alternative, paired = FALSE, var.equal = var.equal)$p.value < alpha.test
                return(testSig)
                
            } else if(testType %in% c("Non-parametric", "Nonparametric", 
                                      "Non", "Nonpar", "nonpar", 
                                      "NP", "N", "Non-Parametric",
                                      "non-parametric", "nonparametric", "non", 
                                      "nonparam", "np", "n")){
                if(is.null(var.equal)){
                    var.equal <- fligner.test(x = x$Result, g = x$Hit, 
                        alternative = var.alternative)$p.value >= alpha.var
                }
                    if(var.equal){
                        testSig <- wilcox.test(x = x[!x$Hit,]$Result, 
                            y = x[x$Hit,]$Result, 
                            alternative = alternative, exact = FALSE)$p.value < alpha.test
                    } else if(!var.equal){
                        testSig <- lawstat::brunner.munzel.test(x = x[!x$Hit,]$Result, 
                                         y = x[x$Hit,]$Result, 
                                         alternative = alternative)$p.value < alpha.test
                    }
                return(testSig)
                }
            } else{
                return(FALSE) 
            }
        }
        )
    )
    return(out)
}## end code

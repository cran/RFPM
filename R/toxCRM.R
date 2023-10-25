#' Concentration-Response Model, 4-parameter log-logistic curve
#' 
#' Generate dummy toxicity data with strong concentration-response
#' 
#' @param  x numeric vector of chemical concentrations (univariate)
#' @param  max numeric value, asymptotic maximum of model (default = \code{1})
#' @param  min numeric value, asymptotic minimum of model (default = \code{0})
#' @param  steep numeric value, steepness factor for model slope
#' @param  mid numeric value, inflection point concentration for CRM slope
#' @param  eMean numeric value, mean of random-normal error to add to dummy toxicity data (default = \code{0})
#' @param  eSD numeric value, standard deviation of random-normal error to add to dummy toxicity data (default = \code{0})
#' @param  seed numeric value, random seed to set for repeatable random error generating (default = \code{NULL}, i.e., no seed)
#' @details \code{toxCRM} generates dummy toxicity data representing an
#' idealized condition. This approach was used by the authors to simulate data and test \code{FPM} sensitivity
#' and baseline variability. The user must specify all coefficients of the model (though
#' default min/max values are provided) and may add some amount of random error, if desired.
#' Random error adds noise to the dummy toxicity daata, increasing the uncertainty of toxicity predictions using floating percentile model benchmarks. The
#' \code{perfect}, \code{lowNoise}, and \code{highNoise} datasets included in \code{RFPM} were all generated
#' using \code{toxCRM} with increasing levels of \code{eSD}.
#' @return numeric vector
#' @seealso perfect, lowNoise, highNoise
#' @importFrom stats rnorm 
#' @examples
#' concentration = h.northport$Cu
#' toxVals <- toxCRM(concentration, 1, 0, 0.5, 
#'     median(concentration), 0, 0)
#' plot(concentration, toxVals, log="x", main="Perfect estimate")
#' 
#' toxVals_withNoise <- toxCRM(x = concentration,
#'     1, 0, 0.5, median(concentration), 0, 0.1, seed = 1)
#' plot(concentration, toxVals_withNoise, log="x", main="Noisy estimate")
#' @export
toxCRM <- function(x,
                   max = 1, 
                   min = 0, 
                   steep, 
                   mid, 
                   eMean = 0, 
                   eSD = 0, 
                   seed = NULL) {
    if(!is.null(seed)){
        set.seed(seed)
    } else {
        set.seed = NULL
    }
    
    min + 
        ((max - min)/(1 + exp(steep * (log(x) - log(mid))))) + 
        rnorm(length(x), eMean, eSD)
} #end code
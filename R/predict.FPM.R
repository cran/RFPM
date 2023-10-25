#' Predict Toxicity Using the Floating Percentile Model
#' 
#' Use new sediment chemistry data to generate FPM predictions
#' 
#' @param  object FPM class object, created using \code{FPM} (and extracted from the resulting list using \code{.$FPM}, where '.' is the name of the object or FPM call)
#' @param  newdata character vector of column names of chemical concentration variables in \code{data}
#' @param  ... further arguments passed to or from other methods
#' @details There are two things to keep in mind when using predict for 'FPM' objects. Firstly, when \code{FPM} is run, the output object is a list; one of the objects (called "FPM") is the FPM class object that can be used to predict Hits. 
#' 
#' Secondly, unlike other default "predict" methods, \code{predict.FPM} is used strictly to predict toxicity that is not included in the original dataset used to generate \code{fpm}.
#' To predict Hit results for the original data, set \code{hitInfo == TRUE} when running the \code{FPM} function. Note that predicting toxicity for the original dataset may be arbitrary/unnecessary, as the Hit results for the original dataset must be known.
#' 
#' Note that, in order to run \code{predict.FPM}, the \code{newdata} argument must be supplied a \code{data.frame} that includes at least one sample with all of the
#' chemical columns contained in \code{fpm}. Column headers must match exactly.
#' 
#' @return logical
#' @examples
#' # create FPM object with chemical headings that overlap
#' overlap <- intersect(names(h.northport)[1:10], 
#'             names(h.tristate)[7:13])
#' fpm_object = FPM(h.northport, overlap)
#' # run predict on the 'FPM' list item to estimate Hits
#' predict(object = fpm_object$FPM, newdata = h.tristate)
#' @importFrom stats predict
#' @export
predict.FPM <- function(object, newdata, ...){
    tryCatch(expr = {if(!is.data.frame(object)) {object <- object$FPM}},
             error = function(e) stop("'object' must be the list output by FPM or the 'FPM' data.frame extracted from such a list"))
    
    if(!is.data.frame(newdata)){newdata <- as.data.frame(newdata)}
    # Reduce to relevant columns only - matching chemicals
    chems <- names(object)[!names(object) %in% c("FN_crit", "TP", "FN", "TN", "FP", "pFN", "pFP", "sens", "spec", "OR", "FM", "MCC")]
    if(!all(chems %in% names(newdata))) stop("newdata does not include all chemicals needed to predict Hits")
    
    trim <- newdata[chems]
    EF <- data.frame(lapply(chems, function(x) {
        trim[,x]/object[,x]
        }))
    return(apply(EF, 1, function(x) any(x > 1)))
}## end code
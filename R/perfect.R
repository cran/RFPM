#' Simulated Data for Floating Percentile Model, Perfect Toxicity Prediction
#' 
#' Example dataset
#' 
#' @name perfect
#' @docType data
#' @format data.frame containing 147 rows and 11 variables:
#' \describe{
#'   \item{Al}{aluminum concentration, mg/kg dry weight}
#'   \item{As}{arsenic concentration, mg/kg dry weight}
#'   \item{Cu}{copper concentration, mg/kg dry weight}
#'   \item{Cd}{cadmium concentration, mg/kg dry weight}
#'   \item{Cr}{chromium concentration, mg/kg dry weight}
#'   \item{Fe}{iron concentration, mg/kg dry weight}
#'   \item{Pb}{lead concentration, mg/kg dry weight}
#'   \item{Hg}{mercury concentration, mg/kg dry weight}
#'   \item{Ni}{nickel concentration, mg/kg dry weight}
#'   \item{Zn}{zinc concentration, mg/kg dry weight}
#'   \item{Hit}{logical; whether the sample was classified as toxic}
#' }
#' @details \code{perfect} provides a sample of simulated data that were developed to analyze variability and sensitivity of \code{FPM}. 
#' Simulated datasets (n=1000) were originally developed for the analyses; \code{perfect} is a single example realization of the simulation.
#' Simulations were generated using the covariance matrix of \code{h.northport} sediment chemical concentrations and the
#' \code{rmvnorm} function from the \code{splus2R} package (Constantine and Hesterberg 2021). The \code{Hit} values were generated with \code{toxCRM}
#' using simulated concentrations of \code{Cr}, \code{Cu}, \code{Fe}, and \code{Zn}
#' as inputs. The inputs to \code{toxCRM} were based on empirical toxicity ranges (minimum and maximum effect levels) as well as 
#' arbitrarily based on probable effect concentration (PEC) and threshold effect concentration (TEC) values from MacDonald et al (2000). Inflection points were
#' set equal to PECs, and steepnesses were set as the ratio of the PEC to the TEC.
#' Fe inputs were arbitrary chosen to reflect its low toxicity relative to Cr, Cu, and Zn (e.g., high inflection point and relatively low steepness).
#' The median toxicity output (for the 4 chemicals) was compared to a threshold of 75% with \code{Hit == TRUE} assigned to values below that level and
#' \code{Hit == FALSE} assigned to values >=75%.
#' @seealso FPM, toxCRM, h.northport, lowNoise, highNoise
#' @references
#' Constantine W, Hesterberg T. 2021. splus2R: Supplemental S-PLUS Functionality in R. Version 1.3-3 (online). Updated January 30. Available from: https://cran.r-project.org/web/packages/splus2R/index.html.
#' MacDonald DD, Ingersoll CG, Berger TA. 2000. Development and evaluation of consensus-based sediment quality guidelines for freshwater ecosystems. Arch Environ Contam Toxicol 39(5):20-31.
NULL
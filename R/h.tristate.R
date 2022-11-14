#' Tristate Data, \emph{Hyalella azteca} Toxicity
#'
#' Example dataset
#' 
#' @name h.tristate
#' @docType data
#' @format A data frame with 43 rows and 13 variables:
#' \describe{
#'   \item{Sample.ID}{sediment sample ID}
#'   \item{Sample.type}{sediment sample type}
#'   \item{Organism}{organism abbreviation, \code{HA} = \emph{Hyalella azteca}}
#'   \item{Day}{measurement day for toxicity endpoint (post-initiation)}
#'   \item{Endpoint}{toxicity endpoint}
#'   \item{Hit}{logical; whether we classified the sample as toxic (based on effect >25%)}
#'   \item{Cd}{cadmium concentration, mg/kg dry weight}
#'   \item{Cu}{copper concentration, mg/kg dry weight}
#'   \item{Fe}{iron concentration, mg/kg dry weight}
#'   \item{Mn}{manganese concentration, mg/kg dry weight}
#'   \item{Ni}{nickel concentration, mg/kg dry weight}
#'   \item{Pb}{lead concentration, mg/kg dry weight}
#'   \item{Zn}{zinc concentration, mg/kg dry weight}
#' }
#' @details \code{h.trisate} is a field-collected dataset containing sediment chemical concentrations for metals and toxicity classification data for
#' the amphipod \emph{Hyalella azteca}. Data are from freshwater sites in Kansas, Missouri, and Oklahoma (Ingersoll et al 2008). The Hit variable was
#' assigned by the authors of this package.
#' @references
#' Ingersoll CG, MacDonald DD, Besser JM, Brumbaugh WG, Ivey CD, Kemble NE, Kunz JL, May TW, Wang N, Smorong DE. 2008. Sediment chemistry, toxicity, and bioaccumulation data report for the US Environmental Protection Agency - Department of the Interior sampling of metal-contaminated sediment in the Tri-state Mining District in Missouri, Oklahoma, and Kansas. US Geological Survey and MacDonald Environmental Sciences, Ltd., Columbia, MO.
NULL
#' Floating Percentile Model
#' 
#' Floating Percentile Model and supporting functions to inform and improve sediment benchmark development
#' 
#' @details 
#' 'RFPM' is an open-source implementation of the floating percentile model (FPM), which was
#' originally developed by Avocet (2003) using Visual Basic for Applications and distributed to US Pacific Northwest regulatory agencies as a 
#' Microsoft Excel-based tool. 'RFPM' was developed independently by Claire Detering and Brian Church with 
#' support from John Toll and others at Windward Environmental LLC. 
#' 
#' The purpose of the FPM is to generate aquatic toxicity-based sediment quality benchmarks
#' for management of contaminated freshwater sediment sites. These benchmarks are intended to act as classification thresholds,
#' meaning that an exceedance of benchmarks would imply that toxicity (as categorically defined) is likely in the sediment sample.
#' The FPM has been used at sites in the US Pacific Northwest for many years, particularly after being published by the 
#' Washington State Department of Ecology in 2011.
#' 
#' The primary function in 'RFPM' is \code{FPM}, which runs the FPM algorithm on a data.frame object that includes concentrations of
#' chemicals in sediment as well as a logical toxicity classification column called "Hit". Example datasets are provided. The output of \code{FPM}
#' includes a set of sediment quality benchmarks for chemicals with significantly higher concentrations when \code{Hit == TRUE}
#' than when \code{Hit == FALSE}. Plots comparing the \code{Hit == TRUE} and \code{Hit == FALSE} data can also
#' be generated as a diagnostic tool. Supplemental functions (e.g., \code{optimFPM}) can help to optimize \code{FPM} inputs (resulting in more accurate benchmarks) or evaluate
#' the relative importance of each chemical among the FPM benchmarks (i.e., \code{chemVI}).
#' 
#' For 'RFPM', the \code{FPM} algorithm has been changed from the original Avocet (2003) model. Key changes are as follows:
#' 1. A decision tree is implemented to select statistically appropriate hypothesis tests for chemical selection. This can be overridden if the
#' original Excel-based tool method is desired; see \code{?chemSig} for details. The chemical selection and FPM algorithm have been integrated into \code{FPM}, though chemical selection can still be run separately, if desired.
#' 2. The iterative looping of the FPM algorithm over multiple false negative limits was
#' not included in \code{FPM}. We find this functionality of the Excel-based tool to be confusing and unnecessary. 
#' Instead, we believe the results generated for different false negative limits should be independently generated rather than dependent on prior model runs.
#' Thus, \code{FPM} allows for multiple false negative limits to be input in a seqeuential but independent manner, resulting in a data.frame of benchmarks with one 
#' row per false negative limit. In R terminology, \code{FPM} has been vectorized.
#' 3. An optimization function called \code{optimFPM} was developed to optimize the overall reliability of FPM sediment quality benchmarks. Different optimization metrics are provided, and a weight of evidence should be considered when selecting inputs.
#' 4. A function was developed to quickly calculate chemical "variable importance" statistics using \code{chemVI}. These statistics
#' can inform the user about the influence of specific chemicals over the set of FPM benchmarks and, for example,
#' whether certain benchmarks can be ignored without a significant loss of predictive ability.
#' @return bibentry
#' @importFrom utils bibentry
#' @importFrom utils citation
#' @importFrom utils person
#' @references 
#' Avocet. 2003. Development of freshwater sediment quality values for use in Washington State. Phase II report: Development and recommendation of SQVs for freshwater sediments in Washington State. Publication No. 03-09-088. Prepared for Washington Department of Ecology. Avocet Consulting, Kenmore, WA.
#' Ecology. 2011. Development of benthic SQVs for freshwater sediments in Washington, Oregon, and Idaho. Publication no. 11-09-054. Toxics Cleanup Program, Washington State Department of Ecology, Olympia, WA.
#' @export
RFPM <- function(){
    p <- c(
        person("Brian", "Church", role = c("aut", "cre"),
            email = "brianc@windwardenv.com"),
        person("Claire", "Detering", role = "aut", 
            email = "claired@windwardenv.com")
    )
    
    return(bibentry(
        bibtype = "Manual",
        title = "{RFPM}: Floating Percentile Model",
        author = p,
        year = "2022",
        note = "R package version 1.0",
        url = "http://windwardenv.com"))
}

RFPM()
citation("RFPM") ## end code
#' Color Gradient Generator
#' 
#' Generate a gradient of hexadecimal colors from a numeric vector
#' 
#' @param  x numeric vector
#' @param  colors values recognizable by R as colors - text, hexadecimal, numbers, etc. (default = \code{grDevices::heat.colors(10)}).
#' @param  colsteps numeric value, number of unique colors to include in gradient
#' @param  na.rm logical value (default = \code{TRUE}) indicating whether to exclude \code{NA} values in \code{x}
#' @param  ... additional arguments passed to \code{grDevices::colorRamp}
#' @details \code{colorGradient} generates a color gradient based on a numeric input \code{x}.
#' The original version (named 'color.gradient') was written by David Hoop in a 2016 stack.overflow response to question (online). The function was 
#' modified slightly for adaptability of inputs. Note that the full color gradient is used if possible,
#' which can exaggerate small differences in \code{x}. This function is applied by \code{optimFPM} when generating color-based
#' optimization matrix graphics (i.e., when optimizing both \code{alpha} and \code{FN_crit}). 
#' @seealso optimFPM, heat.colors, colorRamp
#' @return character vector
#' @importFrom grDevices heat.colors
#' @importFrom grDevices colorRampPalette 
#' @examples
#' x <- rnorm(n = 100)
#' cols <- colorGradient(x = x, colors = c("red", "white", "blue"))
#' plot(x, col = cols)
#' @export
colorGradient <- function(x, 
                          colors = heat.colors(10), 
                          colsteps = 100,
                          na.rm = TRUE,
                          ...) {
    cols <- colorRampPalette(colors, ...)(colsteps)[findInterval(x, seq(min(x, na.rm), max(x, na.rm), length.out = colsteps))]    
  return(cols)
}## end code
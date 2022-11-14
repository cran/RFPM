## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = F-------------------------------------------------------------
library(RFPM)

## -----------------------------------------------------------------------------
head(h.northport)

## -----------------------------------------------------------------------------
p.northport <- names(h.northport)[1:10] ## all chemical column names
FPM(data = h.northport, paramList = p.northport) ## minimum input - dataset and chemical column names

## -----------------------------------------------------------------------------
plot(chemSigSelect(h.northport[, c("Al", "Cr", "Hit")], paramList = c("Al", "Cr")), type = "boxplot")

## ---- fig.width = 7, fig.height = 6-------------------------------------------
## one-way optimization of the FN Limit - vertical lines show best values based on two metrics
optimFPM(h.northport, p.northport, FN_crit = seq(0.1, 0.9, 0.05), alpha = 0.05)

## two-way optimization of both FN Limit and alpha - black squares show best values based on two metrics
optimFPM(h.northport, p.northport, 
         FN_crit = seq(0.1, 0.9, 0.05), 
         alpha = seq(0.01, 0.2, 0.01))


## ---- fig.width = 7, fig.height = 6-------------------------------------------
cvFPM(h.northport, p.northport, k = 10, FN_crit = seq(0.1, 0.9, 0.05), alpha = 0.05)

## -----------------------------------------------------------------------------

chemVI(h.northport, p.northport)


## -----------------------------------------------------------------------------
FPM(data = h.northport, paramList = p.northport)
FPM(data = h.northport, paramList = c("Cr", "Zn"))

## -----------------------------------------------------------------------------
FPM(h.northport, p.northport, FN_crit = 0.25, alpha = 0.11)$FPM
chemVI(h.northport, p.northport, FN_crit = 0.25, alpha = 0.11)


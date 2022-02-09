################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

## Use command
##  R CMD build --resave-data OptHoldoutSize
## to build package

################################################################################
## Packages and scripts                                                       ##
################################################################################

#' @import matrixStats
#' @import mnormt
#' @import mvtnorm
#' @import mle.tools
#' @import ranger
#' @import stats
#' @import graphics

require("matrixStats")  # Matrix row- and column-wise operations
require("mnormt")       # Multivariate normal distributions
require("mvtnorm")      # Multivariate normal and t distributions
require("ranger")       # Random forests
require("mle.tools")    # Gaussian processes

################################################################################
## Small auxiliary functions                                                  ##
################################################################################

##' Logit
##'
##' @export
##' @name logit
##' @description Logit function: 1/(1+exp(-x))
##' @param x argument
##' @return value of logit(x)
##' @examples
##'
##' # Plot
##' x=seq(-5,5,length=1000)
##' plot(x,logit(x),type="l")
logit=function(x) 1/(1+exp(-x))

##' Logistic
##'
##' @export
##' @name logistic
##' @description Logistic function: -log((1/x)-1)
##' @param x argument
##' @return value of logit(x); na if x is outside (0,1)
##' @examples
##'
##' # Plot
##' x=seq(0,1,length=100)
##' plot(x,logistic(x),type="l")
##'
##' # Logit and logistic are inverses
##' x=seq(-5,5,length=1000)
##' plot(x,logistic(logit(x)),type="l")
logistic=function(x) -log((1/x)-1)





################################################################################
## Data documentation                                                         ##
################################################################################

##' Data for vignette on algorithm comparison
##'
##' @description This object contains data relating to the vignette comparing emulation and parametric algorithms. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"ohs_array"

##' Data for vignette showing general example
##'
##' @description Data for general vignette. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"data_example_simulation"

##' Data for vignette on algorithm comparison
##'
##' @description This object contains data relating to the first plot in the vignette comparing emulation and parametric algorithms. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"ohs_resample"

##' Data for example on asymptotic confidence interval for OHS.
##'
##' @description Data for example for asymptotic confidence interval for OHS. For generation, see example.
##' @docType data
##' @keywords data aspre
"ci_cover_a_yn"

##' Data for example on empirical confidence interval for OHS.
##'
##' @description Data for example for empirical confidence interval for OHS. For generation, see example.
##' @docType data
##' @keywords data aspre
"ci_cover_e_yn"

##' Data for example on asymptotic confidence interval for min cost.
##'
##' @description Data for example for asymptotic confidence interval for min cost. For generation, see example.
##' @docType data
##' @keywords data aspre
"ci_cover_cost_a_yn"

##' Data for example on empirical confidence interval for min cost.
##'
##' @description Data for example for empirical confidence interval for min cost. For generation, see example.
##' @docType data
##' @keywords data aspre
"ci_cover_cost_e_yn"



##' Data for 'next point' demonstration vignette on algorithm comparison using parametric algorithm
##'
##' @description Data containing 'next point selected' information for parametric algorithm in vignette comparing emulation and parametric algorithms. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"data_nextpoint_par"


##' Data for 'next point' demonstration vignette on algorithm comparison using emulation algorithm
##'
##' @description Data containing 'next point selected' information for emulation algorithm in vignette comparing emulation and parametric algorithms. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"data_nextpoint_em"



##' Parametric-based OHS estimation for ASPRE
##'
##' @description This object contains data relating to parametric-based OHS estimation for the ASPRE model. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"aspre_parametric"


##' Emulation-based OHS estimation for ASPRE
##'
##' @description This object contains data relating to emulation-based OHS estimation for the ASPRE model. For generation, see hidden code in vignette, or in pipeline at https://github.com/jamesliley/OptHoldoutSize_pipelines
##' @docType data
##' @keywords data aspre
"aspre_emulation"

##' Parameters of reported ASPRE dataset
##'
##' @description Distribution of covariates for ASPRE dataset; see Rolnik, 2017, NEJM
##' @docType data
##' @keywords data aspre
"params_aspre"



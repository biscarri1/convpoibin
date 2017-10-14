##' Compute the Poisson binomial distribution via convolution
##' 
##' The cdf and pmf for the Poisson binomial distribution computed using
##' convolution methods
##' 
##' \tabular{ll}{ Package: \tab convpoibin\cr Type: \tab Package\cr Version:
##' \tab 1.0\cr Date: \tab 2017-10-03\cr License: \tab GPL-2 \cr LazyLoad: \tab
##' yes\cr }
##' 
##' Two exact and efficient convolution based methods for computing the Poisson
##' binomial cdf and pmf.
##' 
##' @name convpoibin-package
##' @aliases convpoibin-package convpoibin
##' @docType package
##' @author William Biscarri
##' 
##' Maintainer: William Biscarri <wbiscar2@@illinois.edu>
##' @keywords package
##' @useDynLib convpoibin
##' @examples
##' 
##'   p=runif(100)
##'   ppoibin(p,method='direct')
##'   ppoibin(p,method='fft_tree',M=2)
##'   dpoibin(p,method='direct')
##'   dpoibin(p,method='fft_tree',M=2)
##'   
##'   p=runif(1000)
##'   ppoibin(p,method='direct')
##'   ppoibin(p,method='fft_tree',M=4)
##'   dpoibin(p,method='direct')
##'   dpoibin(p,method='fft_tree',M=4)
##' 
NULL
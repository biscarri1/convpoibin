\name{convpoibin-package}
\alias{convpoibin-package}
\alias{convpoibin}
\docType{package}
\title{
Compute the Poisson binomial distribution via convolution

}
\description{
The cdf and pmf for the Poisson binomial distribution computed using convolution methods

}
\details{
\tabular{ll}{
Package: \tab convpoibin\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-9-30\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}

This package implements two exact and efficient convolution based methods for computing the Poisson binomial cdf and pmf.


}
\author{
William Biscarri

Maintainer: William Biscarri <wbiscar2@illinois.edu>

}

\keyword{ package }

\examples{
   p=c(.1,.2,.3,.4,.5)
   pconvpoibin(p)
   pconvpoibin(p,method='fft_tree',M=2)
   dconvpoibin(p)
   dconvpoibin(p,method='fft_tree',M=2)
   p = runif(1000)
   pconvpoibin(p)
   pconvpoibin(p,method='fft_tree',M=8)
   dconvpoibin(p)
   dconvpoibin(p,method='fft_tree',M=8)
}
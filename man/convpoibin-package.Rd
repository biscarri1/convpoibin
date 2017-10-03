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
   kk=0:10
   pp=c(.1,.2,.3,.4,.5)
   ppoibin(kk=kk, pp=pp, method = "DFT-CF",wts=rep(2,5))
   ppoibin(kk=kk, pp=pp, method = "RF",wts=rep(2,5))
   ppoibin(kk=kk, pp=pp, method = "RNA",wts=rep(2,5))
   ppoibin(kk=kk, pp=pp, method = "NA",wts=rep(2,5))
   dpoibin(kk=kk, pp=pp,wts=rep(2,5))
   qpoibin(qq=0:10/10,pp=pp,wts=rep(2,5))
   rpoibin(m=2,pp=pp,wts=rep(2,5))
}
\name{convpoibin}
\alias{ppoibin}
\alias{dpoibin}
\alias{qpoibin}
\alias{rpoibin}
\title{
Compute the Poisson binomial distribution via convolution
}
\description{
The cdf and pmf for the Poisson binomial distribution computed using convolution methods
}
\usage{
pconvpoibin(p,method='direct')
pconvpoibin(p,method='fft_tree',M=2)
dconvpoibin(p,method='direct')
dconvpoibin(p,method='fft_tree',M=2)
}
\arguments{

  \item{p}{ The vector of success probabilities
}
  \item{method}{ "DFT-CF" for the DFT-CF method, "RF" for the recursive formula, "RNA" for the refined normal approximation,
    "NA" for the normal approximation.
}
  \item{qq}{ The values where the quantile function to be evaluated.
}
\details{
See the reference for computational details.
}
}
\value{
      Returns the entire cdf or pmf
}
\references{
Hong, Y. (2012). On computing the distribution function for the Poisson binomial distribution. Computational Statistics & Data Analysis. Available online, http://dx.doi.org/10.1016/j.csda.2012.10.006.
}
\author{

Yili Hong

}


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

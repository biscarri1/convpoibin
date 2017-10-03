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
ppoibin(kk, pp, method = "DFT-CF",wts=NULL)
dpoibin(kk, pp,wts=NULL)
qpoibin(qq, pp,wts=NULL)
rpoibin(m, pp,wts=NULL)
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

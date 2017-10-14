##' Compute the Poisson binomial distribution via convolution
##' 
##' The cdf and pmf for the Poisson binomial distribution computed using
##' convolution methods
##' 
##' See the reference for computational details.
##'
##' @aliases ppoibin dpoibin
##' @param p The vector of success probabilities
##' @param method "direct" for the direct convolution method, "fft_tree" for
##' the binary tree fft convolution method.
##' @param M The number of splits for the tree used in the fft_tree
##' convolution.
##' @return Returns the entire cdf or pmf
##' @export
##' @rdname convpoibin
##' @author
##' 
##' William Biscarri
##' @references test test test
##' @examples
##' 
##'   p=runif(100)
##'   ppoibin(p,method='direct',M=2)
##'   ppoibin(p,method='fft_tree',M=2)
##'   dpoibin(p,method='direct',M=2)
##'   dpoibin(p,method='fft_tree',M=2)
##'   
##'   p=runif(1000)
##'   ppoibin(p,method='direct',M=2)
##'   ppoibin(p,method='fft_tree',M=4)
##'   dpoibin(p,method='direct',M=2)
##'   dpoibin(p,method='fft_tree',M=4)
##'
dpoibin <- function(p,method="direct",M=2){
  
  lp = as.integer(length(p))
  total_size = lp+1
  
  switch(method,
         "direct"={
           density = .C('direct_convolution', as.double(p),lp,as.double(vector("double",total_size)))[[3]]
         },
         "fft_tree" = {
           num_splits = as.integer(M)
           density  = .C('tree_convolution',as.double(p),lp,num_splits,
                         as.double(vector("double",lp+num_splits)))[[4]][1:total_size]
         }
        )
  
  return(density)
}


##' @export
##' @rdname convpoibin
ppoibin <- function(p,method="direct",M=2){
  
  lp = as.integer(length(p))
  total_size = lp+1
  #NOTE: for some reason, saving result and then doing return(result) causes an error
  #however, doing return(cumsum(density)) works, however. Need to understand why.
  switch(method,
         "direct"={
           density = .C('direct_convolution', as.double(p),lp,as.double(vector("double",total_size)))[[3]]
           #result = cumsum(density)
           #density=0
         },
         
         "fft_tree" = {
           num_splits = as.integer(M)
           density  = .C('tree_convolution',as.double(p),lp,num_splits,
                           as.double(vector("double",total_size)))[[4]][1:total_size]
           #result = cumsum(density)
           #density=0
         }
        )
  
  return(cumsum(density))
}

dconvpoibin <- function(p,method="direct",M=2){
  
  lp = as.integer(length(p))
  
  switch(method,
         "direct"={
           
           total_size = lp+1
           
           density = .C('direct_convolution', as.double(p),as.integer(lp),as.double(vector("double",total_size)))[[3]]
           
         },
         
         "fft_tree" = {
                     
           num_splits = as.integer(M)
           
           result = as.double(vector("double",lp+num_splits))
           
           density  = .C('tree_convolution',p,lp,num_splits,result)[[4]][1:(lp+1)]
           
           
           
         }
  )
  
  return(density)
  
}



pconvpoibin <- function(p,method="direct",M=2){
  
  lp = as.integer(length(p))
  
  switch(method,
         "direct"={
           
           total_size = lp+1
           
           density = .C('direct_convolution', as.double(p),as.integer(lp),as.double(vector("double",total_size)))[[3]]
           
           result = cumsum(density)

         },
         
         "fft_tree" = {
           
           num_splits = as.integer(M)
           
           result = as.double(vector("double",lp+num_splits))
           
           density  = .C('tree_convolution',p,lp,num_splits,result)[[4]]
           
           result = cumsum(density[1:(lp+1)])
           
         }
  )
  
  return(result)
  
}
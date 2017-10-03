split_conv = function(p,num_splits){
    
    lp = as.integer(length(p))

    num_splits = as.integer(num_splits)

    result = as.double(vector("double",lp+num_splits))

    density  = .C('splitter',p,lp,num_splits,result)[[4]]
    
    return(density)
    
    }
    
    
    
direct_convolution = function(p){
    
    lp = length(p)
    total_size = lp+1
    density = .C('direct_conv', as.double(p),as.integer(lp),as.double(vector("double",total_size)))[[3]]
    
    return(density)
    
    }
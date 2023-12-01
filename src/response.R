lapply(datasets,function(dataset){
  X <- dataset$X
  y <- dataset$y
  
  W <- matrix(NA,nrow(X),ncol(X))
  
  for (j in 1:ncol(X)){
    x <- X[,j]
    if (length(unique(x)) <= 2 | length(unique(diff(x)))==1){
      next
    }
    
    m <- median(x)
    md <- mad(x,m)
    a <- m-md
    b <- m+md
  
    W[,j] <- x >= a & x <= b
  }
  
  sum(apply(W,1,mean) & y)/sum(y)
  
})


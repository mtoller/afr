abScale <- function(x,a=0,b=1){
  a+(x-min(x))*(b-a)/(max(x)-min(x))
}
benchGauss <- function(x,...){
  z <- c(pcaPP::PCAproj(x,k=1)$scores)
  mu <- mean(z)
  sigma2 <- mean((z-mu)^2)
  
  return(dnorm(mu,mu,sqrt(sigma2)) - dnorm(z,mu,sqrt(sigma2)))
}
benchGaussCopula <- function(X,...){
  n <- nrow(X)
  d <- ncol(X)
  
  skewness <- function(x){
    m <- mean(x)
    (1/(length(x))*sum((x-m)^3)) / (1/length(x)*sum((x-m)^2))^(3/2)
  }
  
  U <- matrix(NA,n,d)
  V <- matrix(NA,n,d)
  W <- matrix(NA,n,d)
  
  for (j in 1:d){
    x <- X[,j]
    
    mu <- mean(x)
    sigma2 <- mean((x-mu)^2)
    
    U[,j] <- pnorm(x,mu,sqrt(sigma2),log.p = TRUE)
    
    V[,j] <- pnorm(-x,-mu,sqrt(sigma2),log.p = TRUE)
    
    if (skewness(x) < 0){
      W[,j] <- U[,j]
    }
    else{
      W[,j] <- V[,j]
    }
  }
  pl <- apply(U,1,function(u)-sum(u))
  pr <- apply(V,1,function(u)-sum(u))
  ps <- apply(W,1,function(u)-sum(u))
  
  return(apply(cbind(pl,pr,ps),1,max))
}
song <- function(x,bw=3,tol=0.5,prop=0.1,...){
  
  m <- median(x)
  
  x <- x - m
  
  #compute loglikelihood
  n <- length(x)
  
  liks <- rep(NaN,n)
  for (t in 1:n){
    liks[t] <- 0
    for (j in 1:n){
      liks[t] <- liks[t] + 1/(n*bw)*dnorm((x[t]-x[j])/bw,0,1)
    }
  }
  
  B_hat <- rep(NA,n)
  sigma2 <- 1
  p <- 0.5
  
  for (t in 1:100){
    #E step
    
    B_hat <- 2*(p)*dnorm(x,0,(sigma2))/((p)*dnorm(x,0,sqrt(sigma2)) + liks)
    B_hat[which(B_hat > 1)] <- 1
    
    
    #M step
    p <- sum(B_hat)/n
    
    sigma2 <- sum(B_hat*x^2)/sum(x)
    
    
  }
  
  if (any(is.na(B_hat))){
    return(rep(FALSE,n))
  }
  #return(1:length(x) %in% rev(order(abScale(B_hat)))[1:prop])
  # return(abScale(B_hat) > tol)
  return(B_hat)
}

uniGaussEM <- function(x,mu_0=mean(x),sigma2_0=var(x),p_0=0.5,factor=1,tol=0.5,prop=0.1,...){
  n <- length(x)
  
  
  mu <- mu_0
  sigma2 <- sigma2_0
  p <- p_0
  
  uni_dens <- factor/(max(x)-min(x))
  B <- rep(NA,n)
  for (t in 1:100){
    B_old <- B
    B <- p*uni_dens/(p*uni_dens + (1-p)*dnorm(x,mu,sqrt(sigma2)))
    
    p_old <- p
    mu_old <- mu
    sigma2_old <- sigma2
    
    
    p <- sum(B)/n
    mu <- sum(B*x)/sum(B)
    sigma2 <- sum(B*x^2)/sum(B)
    
    if (mu_old==mu & sigma2_old == sigma2 & p_old==p){
      break
    }
  }
  #return(1:length(x) %in% rev(order(abScale(B)))[1:prop])
  # return(abScale(B)>tol)
  return(B)
}
uniGaussEmHard <- function(x,mu_0=mean(x),sigma2_0=var(x),p_0=0.5,factor=1,prop=0.1,...){
  n <- length(x)
  
  
  mu <- mu_0
  sigma2 <- sigma2_0
  p <- p_0
  
  uni_dens <- factor/(max(x)-min(x))
  B <- rep(NA,n)
  for (t in 1:100){
    B_old <- B
    B <- dnorm(x,mu,sqrt(sigma2)) <= uni_dens
    
    p_old <- p
    mu_old <- mu
    sigma2_old <- sigma2
    
    
    p <- sum(B)/n
    mu <- mean(x[which(!B)])
    sigma2 <- var(x[which(!B)])
    
    if (mu_old==mu & sigma2_old == sigma2 & p_old==p){
      break
    }
  }
  
  #return(1:length(x) %in% rev(order(abScale(B)))[1:prop])
  return(B)
}

iforest <- function(x,...){
  iforest <- import("pyod.models.iforest")$IForest(n_estimators=as.integer(1000),
                                                   random_state = as.integer(0))
  iforest$fit(x)
  return(iforest$decision_scores_)
}

copod <- function(x,...){
  copod <- import("pyod.models.copod")$COPOD()
  
  copod$fit(x)
  return(copod$decision_scores_)
}

lof <- function(x,...){
  lof <- import("pyod.models.lof")$LOF(n_neighbors = as.integer(5))
  lof$fit(x)
  return(lof$decision_scores_)
}


ocsvm <- function(x,...){
  ocsvm <- import("pyod.models.ocsvm")$OCSVM()
  
  ocsvm$fit(x)
  return(ocsvm$decision_scores_)
}
ocsvmPu <- function(x,a,b,prop,...){
  ocsvm <- import("pyod.models.ocsvm")$OCSVM(contamination = prop/nrow(x))
  known_normal <- x >= a & x <= b
  
  x_train <- as.matrix(x[which(known_normal)])
  x_test <- as.matrix(x[which(!known_normal)])
  
  ocsvm$fit(x_train)
  prediction <- ocsvm$predict(x_test)
  
  B <- rep(NA,nrow(x))
  B[which(known_normal)] <- FALSE
  B[which(!known_normal)] <- as.logical(prediction)
  #return(1:nrow(x) %in% rev(order(ocsvm$decision_scores_))[1:prop])
  return(B)
}

deepSvdd <- function(x,prop,...){
  deep_svdd <- import("pyod.models.deep_svdd")$DeepSVDD(random_state = as.integer(0),verbose = 0)
  deep_svdd$fit(x)
  
  return(deep_svdd$decision_scores_)
}

autoEnc <- function(x,prop,...){
  autoenc <- import("pyod.models.auto_encoder")$AutoEncoder(random_state = as.integer(0),
                                                            hidden_neurons = as.integer(c(1,64,32,32,64,1)),
                                                            verbose = 0)
  autoenc$fit(x)
  
  return(autoenc$decision_scores_)
}

deepSvddPu <- function(x,a,b,prop,...){
  deep_svdd <- import("pyod.models.deep_svdd")$DeepSVDD(contamination = prop/nrow(x),
                                                        random_state = as.integer(0),verbose = 0)
  known_normal <- x >= a & x <= b
  
  x_train <- as.matrix(x[which(known_normal)])
  x_test <- as.matrix(x[which(!known_normal)])
  
  deep_svdd$fit(x_train)
  prediction <- deep_svdd$predict(x_test)
  
  B <- rep(NA,nrow(x))
  B[which(known_normal)] <- FALSE
  B[which(!known_normal)] <- as.logical(prediction)
  #return(1:nrow(x) %in% rev(order(deep_svdd$decision_scores_))[1:prop])
  return(B)
}



prenet <- function(X_train,y_train,X_test,...){
  oldwd <- getwd()
  setwd("semi/prenet/")
  run <- import("prenet_predict")$prenet_predict
  setwd(oldwd)
  
  run(X_train,X_test,y_train,"placeholder",ceiling(sum(y_train)/2),
      batch_size = as.integer(10))[[1]]
}

overlap <- function(X_train,y_train,X_test,...){
  oldwd <- getwd()
  setwd("semi/Overlap-master/")
  run <- import("run")$run_ADSD
  setwd(oldwd)
  
  ratio <- mean(y_train)
  
  run(as.matrix(X_train),as.matrix(y_train),
      X_test,seed = as.integer(0),batch_size = as.integer(10),
      ratio = mean(y_train))
}

# puma <- function(X_train,y_train,X_test,...){
#   oldwd <- getwd()
#   setwd("semi/PU-MIL-AD-main/")
#   run <- import("run_puma")$run
#   setwd(oldwd)
# }
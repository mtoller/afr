abScale <- function(x,a=0,b=1){
  a+(x-min(x))*(b-a)/(max(x)-min(x))
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

iforest <- function(x,prop,...){
  inverted <- FALSE
  if (prop > 0.5*nrow(x)){
    prop <- nrow(x)-prop
    inverted <- TRUE
  }
  iforest <- import("pyod.models.iforest")$IForest(n_estimators=as.integer(1000),
                                                   random_state = as.integer(0),
                                                   contamination = prop/nrow(x))
  iforest$fit(x)
  if (inverted){
    return(!as.logical(iforest$labels_))
  }
  return(as.logical(iforest$labels_))
  #return(1:nrow(x) %in% rev(order(iforest$decision_scores_))[1:prop])
  # return(abScale(iforest$decision_scores_) > 0.5)
}

copod <- function(x,prop,...){
  inverted <- FALSE
  if (prop > 0.5*nrow(x)){
    prop <- nrow(x)-prop
    inverted <- TRUE
  }
  copod <- import("pyod.models.copod")$COPOD(contamination = prop/nrow(x))
  
  copod$fit(x)
  if (inverted){
    return(!as.logical(copod$labels_))
  }
  return(as.logical(copod$labels_))
  #return(1:nrow(x) %in% rev(order(copod$decision_scores_))[1:prop])
  # return(abScale(copod$decision_scores_) > 0.5)
}

lof <- function(x,prop,...){
  inverted <- FALSE
  if (prop > 0.5*nrow(x)){
    prop <- nrow(x)-prop
    inverted <- TRUE
  }
  lof <- import("pyod.models.lof")$LOF(n_neighbors = as.integer(5),contamination = prop/nrow(x))
  
  lof$fit(x)
  if (inverted){
    return(!as.logical(lof$labels_))
  }
  return(as.logical(lof$labels_))
  #return(1:nrow(x) %in% rev(order(lof$decision_scores_))[1:prop])
  # return(abScale(lof$decision_scores_) > 0.5)
}


ocsvm <- function(x,prop,...){
  inverted <- FALSE
  if (prop > 0.5*nrow(x)){
    prop <- nrow(x)-prop
    inverted <- TRUE
  }
  ocsvm <- import("pyod.models.ocsvm")$OCSVM(contamination = prop/nrow(x))
  
  ocsvm$fit(x)
  if (inverted){
    return(!as.logical(ocsvm$labels_))
  }
  return(as.logical(ocsvm$labels_))
  # return(1:nrow(x) %in% rev(order(ocsvm$decision_scores_))[1:prop])
  # return(abScale(ocsvm$decision_scores_) > 0.5)
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
  inverted <- FALSE
  if (prop > 0.5*nrow(x)){
    prop <- nrow(x)-prop
    inverted <- TRUE
  }
  deep_svdd <- import("pyod.models.deep_svdd")$DeepSVDD(contamination = prop/nrow(x),
    random_state = as.integer(0),verbose = 0)
  deep_svdd$fit(x)
  if (inverted){
    return(!as.logical(deep_svdd$labels_))
  }
  return(as.logical(deep_svdd$labels_))
  #return(1:nrow(x) %in% rev(order(deep_svdd$decision_scores_))[1:prop])
  # return(abScale(deep_svdd$decision_scores_) > 0.5)
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
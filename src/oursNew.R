logLikelihood <- function(n,s_x,s_x2,s_B,mu,sigma2,p){
  result <- (n-s_B)*(safeLog(1-p)-safeLog(sqrt(2*pi*sigma2))) + s_B * safeLog(p)
  return(result - (s_x2 - 2*mu*s_x + (n-s_B)*mu^2)/(2*sigma2))
}

checkConstraint <- function(a,b,mu,sigma2,p,P_mean,w,strict=FALSE,eps=1e-8){
  if (strict){
    abs(w^2-(1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2))) -P_mean)^2) <= eps
  }
  else{
    w^2 >= (1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2))) -P_mean)^2
  }
}

optimTarget <- function(a,b,mu,sigma2,p,P_mean,w){
  if (is.nan(p)){
    Inf
  }
  else{
    w^2 - (1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2))) -P_mean)^2
  }
  
}

omega <- function(s_x,s_x2,s_B,n,mu,sigma2s,a,b,info=FALSE){
  if (s_x/(n-s_B) == (a+b)/2){
    stop("There is no constrained solution for s_x/(n-s_B)=(a+b)/2. Try different labels")
  }
  if (info){
    print(paste0("sigma2=",sigma2s))
  }
  sigma2s[which(sigma2s < 0)] <- 0
  sigma2s[which(!is.finite(sigma2s))] <- 0
  
  if (all(sigma2s == 0)){
    stop("You did not pass a positive and finite sigma^2")
    return(list(omega=NaN,sigma2=0))
  }
  
  s <- sqrt(sigma2s)
  
  I <- pnorm(b,mu,s)-pnorm(a,mu,s)
  
  dldmu <- (s_x - (n-s_B)*mu)/sigma2s
  dlds <- (s_x2 - 2 * mu * s_x + (n-s_B)*(mu^2-sigma2s))/s^3
  
  dIdmu <- (exp(-(a-mu)^2/(2*sigma2s))-exp(-(b-mu)^2/(2*sigma2s)))/(sqrt(2*pi)*s)
  dIds <- ((a-mu)*exp(-(a-mu)^2/(2*sigma2s))-(b-mu)*exp(-(b-mu)^2/(2*sigma2s)))/(sqrt(2*pi)*sigma2s)
  if (info){
    
    print(paste0("dldmu=",dldmu))
    print(paste0("dIdmu=",dIdmu))
    print(paste0("dlds=",dlds))
    print(paste0("dIds=",dIds))
    print(paste0("I=",I))
  }

  omega_mu <- dldmu / dIdmu * I
  omega_s <- dlds / dIds * I
  
  if (info){
    print(paste0("omega_mu=",omega_mu))
    print(paste0("omega_s=",omega_s))
  }
  #for the mle sigmas there is exactly one sigma s.t. omega_mu=omega_sigma
  diff <- abs(signif(omega_mu,6) - signif(omega_s,6))
  diff[which(!is.finite(diff))] <- Inf
  
  if (all(diff > 1e-6)){
    return(list(omega=NaN,sigma2=0))
  }
  else if(sum(diff <= 1e-6) > 1){
    warning("More than one omega_mu - omega_s match. Taking smaller one.")
  }
  
  return(list(omega=omega_mu[which.min(diff)],sigma2=sigma2s[which.min(diff)]))
}

sanityCheck <- function(mu,x_bar,x2_bar,a,b){
  alpha <- x2_bar - mu*x_bar + (mu-x_bar)*a
  beta <- x2_bar - mu*x_bar + (mu-x_bar)*b
  m <- 1/2*((a-mu)^2-(b-mu)^2)
  e_power <- exp(-m/beta)
  r <- -alpha/beta * e_power
  inp <- m*(alpha-beta)/beta^2*e_power
  f <- function(x,r)x*exp(x)+r*x
  gamma <- LambertW0(-r*exp(1)) - 1
  
  if (inp == 0 | !is.finite(inp) | !is.finite(r) | is.nan(gamma) | 
      (r < 0 & inp < f(gamma,r))){
    FALSE
  } 
  else{
    TRUE
  } 
}

getParams <- function(mu,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,info=FALSE){
  
  #check for 0 or +- Inf in advance
  if (!sanityCheck(mu,x_bar,x2_bar,a,b)){
    return(list(mu=mu,sigma2=0,p=NaN))
  }
  
  sigma2s <- theLambertSolutionsC(x2_bar,x_bar,mu,a,b)

  if (all(!is.finite(sigma2s))){
    return(list(mu=mu,sigma2=0,p=NaN))
  }
  if (all((sigma2s[which(is.finite(sigma2s))]) <= 0)){
    warning("No positive sigma2")
    return(list(mu=mu,sigma2=0,p=NaN))
  }
  res <- omega(s_x,s_x2,s_B,n,mu,sigma2s,a,b,info=info)
  sigma2 <- res$sigma2
  
  p <- s_B/(n-res$omega)
  if (sigma2 == 0){
    return(list(mu=mu,sigma2=0,p=NaN))
  }
  
  if (!is.finite(p)){
    stop("p is not finite. This should never happen.")
  }
  if (p < 0 || p > 1){
    return(list(mu=mu,sigma2=NaN,p=NaN))
    
  }
  
  return(list(mu=mu,sigma2=sigma2,p=p))
}

getErrorAndParam <- function(mu,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w){
  params <- getParams(mu,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b)
  optimTarget(a,b,params$mu,params$sigma2,params$p,P_mean,w)
}

paramOptimizationBinary <- function(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,epsilon=1e-8){
  asymptote <- (a+b)/2
  sigma2_mle <- 1/(n-s_B)*(s_x2+(n-s_B)*x_bar^2-2*x_bar*s_x)
  if (x_bar == asymptote){
    warning("There is no constrained solution for s_x=(a+b)/2. Returning unconstrained MLE")
    return(list(mu=x_bar,sigma2=sigma2_mle,p=s_B/n))
  }
  
  
  
  if (x_bar < asymptote){
    left <- x_bar - sqrt(sigma2_mle)
    right <- asymptote-epsilon
  }
  else{
    left <- asymptote+epsilon
    right <- x_bar + sqrt(sigma2_mle)
  }
  while(abs(left-right) > epsilon){#binary search to find maximum of constraint function
    mid <- (left+right)/2
    par_minus <- getParams(mid-epsilon,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b)
    par_plus <- getParams(mid+epsilon,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b)
    f_minus <- optimTarget(a,b,par_minus$mu,par_minus$sigma2,par_minus$p,P_mean,w)
    f_plus <- optimTarget(a,b,par_plus$mu,par_plus$sigma2,par_plus$p,P_mean,w)
    
    if (!is.finite(f_minus) | !is.finite(f_plus)){
      warning("Computable range of target function is smaller than epsilon. Returning unconstrained MLE")
      return(list(mu=x_bar,sigma2=sigma2_mle,p=s_B/n))
    }
    
    if (f_plus - f_minus < 0) right <- mid
    else left <- mid
  }
  
  top <- mid
  
  if (x_bar < top){
    left <- x_bar
    right <- top
    increasing <- TRUE
  }
  else{
    right <- x_bar
    left <- top
    increasing <- FALSE
  }
  params <- par_plus
  while(abs(left-right) > epsilon){#binary search to find zero
    mid <- (left+right)/2
    params <- getParams(mid-epsilon,s_x,s_x2,s_B,x_bar,x2_bar,n,a,b)
    f_value <- optimTarget(a,b,params$mu,params$sigma2,params$p,P_mean,w)
    if (increasing){
      if (f_value < 0) left <- mid else right <- mid
    }
    else{
      if (f_value < 0) right <- mid else left <- mid
    }
    
  }
  
  return(params)
}

paramOptimization <- function(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,n_optim_steps=1e+4,...){
  mid <- (a+b)/2
  if (x_bar == mid){
    warning("There is no constrained solution for s_x=(a+b)/2. Returning unconstrained MLE")
  }
  
  if (x_bar < mid){
    lower <- mid - abs(x_bar)*2
    upper <- mid - 1e-7
  }
  else{
    lower <- mid + 1e-7
    upper <- mid + abs(x_bar)*2
  }
  mus <- seq(lower,upper,length=n_optim_steps)
  
  errors <- rep(NaN,n_optim_steps)
  likss <- rep(NaN,n_optim_steps)
  candidate_params <- list()
  for (j in 1:length(mus)){
    print(j)
    params <- getParams(mus[j],s_x,s_x2,s_B,x_bar,x2_bar,n,a,b)
    errors[j] <- optimTarget(a,b,params$mu,params$sigma2,params$p,P_mean,w)
    likss[j] <- logLikelihood(n,s_x,s_x2,s_B,params$mu,params$sigma2,params$p)
    if (j > 1){
      if (sign(errors[j])!=sign(errors[j-1])){
        candidate_params[[length(candidate_params)+1]] <- params
      }
      
    }
  }
  if (all(!is.finite(errors))){
    warning("Optimization failed. Likely too few unique values in data. Returning MLE")
    return(list(mu=x_bar,sigma2=x2_bar - x_bar^2,p=s_B/n))
  }
  if (length(candidate_params)==0){
    warning("No zero crossing in optimization. Returning saddle point instead")
    return(getParams(mus[which.min(abs(errors))],s_x,s_x2,s_B,x_bar,x2_bar,n,a,b))
  }
  liks <- unlist(lapply(candidate_params,function(params)logLikelihood(n,s_x,s_x2,s_B,params$mu,params$sigma2,params$p)))
  plot(mus,errors,type='l',...)
  points(x_bar,optimTarget(a,b,x_bar,1/(n-s_B)*(s_x2+(n-s_B)*x_bar^2-2*x_bar*s_x),s_B/n,P_mean,w),pch=4,col='orange',lwd=3)
  cands <- order(abs(errors))[c(1)]
  points(mus[cands],errors[cands],col="magenta",pch=1,lwd=3)
  lines(rep((a+b)/2,2),c(-1e+6,1e+6),col=rgb(0,0.6,0),lty=3)
  lines(c(-1e+6,1e+6),c(0,0),col='blue',lty=2)
  w_p <- which(errors >= 0)
  polygon(x=c(mus[w_p[1]],mus[w_p[length(w_p)]],rev(mus[w_p])),y=c(errors[w_p[1]],errors[w_p[length(w_p)]],rev(errors[w_p])),
          col = rgb(0.8,0.5,0.6,0.4),border = rgb(0.1,0.1,0.1,0.2))
  
  return(candidate_params[[which.max(liks)]])
}

emAfr <- function(x,a=NULL,b=NULL,alpha=0.05,return_scores=TRUE,ignore_AFR=FALSE,epsilon=1e-8){
  
  if (length(unique(x))<=2){
    warning("anomaly detection on data with less than 3 unique values is meaningless. Returning scaled input")
    return(as.logical(abScale(x)))
  }
  
  if (is.null(a)){
    if (is.null(b)){
      a <- unname(quantile(x,0.25))
      b <- unname(quantile(x,0.75))
    }
    else{
      stop("If 'b' is provided, then 'a' must also be provided")
    }
  }
  if (is.null(b)){
    stop ("If 'a' is provided, then 'b' must also be provided")
  }

  
  #compute quantities that don't change during iterations
  mid <- (a+b)/2
  n <- length(x)
  
  a_cand <- which(x < a | x > b)
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
  
  #initialize B randomly
  B <- rep(FALSE,n)
  # set.seed(0)
  B[a_cand] <- sample(c(FALSE,TRUE),length(a_cand),TRUE,c(0.0,1.0))
  
  Bs <- list()
  
  while (TRUE){
    Bs[[length(Bs)+1]] <- B
    
    if (length(a_cand)==n & sum(B) == n){
      stop("All data points are anomalies and AFR is empty? Bad initialization")
    }
    
    #compute current quantities
    x_bar <- mean(x[which(!B)])
    x2_bar <- mean((x[which(!B)])^2)
    s_x <- sum(x[which(!B)])
    s_x2 <- sum((x[which(!B)])^2)
    s_B <- sum(B)
    
    #compute unconstrained MLE
    mu_mle <- x_bar
    sigma2_mle <- mean((x-mu_mle)^2)
    p_mle <- s_B/n
    
    if (checkConstraint(a,b,mu_mle,sigma2_mle,p_mle,P_mean,w) | ignore_AFR){
      #Case where MLE is good
      mu <- mu_mle
      sigma2 <- sigma2_mle
      p <- p_mle
    }
    else{
      #Case where MLE is not good
      constrained_params <- paramOptimizationBinary(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,epsilon)
      mu <- constrained_params$mu
      sigma2 <- constrained_params$sigma2
      p <- constrained_params$p
      
    }
    
    
    #update B
    B[a_cand] <- FALSE
    if (round(p*n) > 0){
      B[a_cand][order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:round(p*n)]] <- TRUE
    }
    #check convergence
    if (any(unlist(lapply(Bs,function(j)all(j==B))))){
      break
    }
  }
  
  if (return_scores){
    scores <- rep(0,n)
    scores[a_cand] <- dnorm(mu,mu,sqrt(sigma2))-dnorm(x[a_cand],mu,sqrt(sigma2))
    return(scores)
  }
  #else return labels
  return(B)
}

dimensionEmAfr <- function(X,ignore_AFR=FALSE,epsilon=1e-8){
  scores <- matrix(NaN,nrow(X),ncol(X))
  for (j in 1:ncol(X)){
    cat(paste0(j))
    x <- X[,j]
    scores[,j] <- emAfr(x,ignore_AFR=ignore_AFR,epsilon=epsilon)
  }
  cat("\n")
  return(rowSums(scores))
}

randomAfr <- function(x,a=NULL,b=NULL,alpha=0.05,n_Bs=10,ignore_AFR=FALSE,epsilon=1e-8,q1=0.24,q2=0.75){
  if (length(unique(x))<=2){
    warning("anomaly detection on data with less than 3 unique values is meaningless. Returning scaled input")
    return(as.logical(abScale(x)))
  }
  
  if (is.null(a)){
    if (is.null(b)){
      a <- unname(quantile(x,q1))
      b <- unname(quantile(x,q2))
    }
    else{
      stop("If 'b' is provided, then 'a' must also be provided")
    }
  }
  if (is.null(b)){
    stop ("If 'a' is provided, then 'b' must also be provided")
  }
  
  
  #compute quantities that don't change during iterations
  mid <- (a+b)/2
  n <- length(x)
  
  a_cand <- which(x < a | x > b)
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))

  scores <- rep(0,n)
  
  for (j in 1:n_Bs){
    #initialize B randomly
    B <- rep(FALSE,n)
    # set.seed(0)
    B[a_cand] <- sample(c(FALSE,TRUE),length(a_cand),TRUE,c(0.9,0.1))
    while (sum(B) == n-length(a_cand)){ #bad initialization, repeat
      B[a_cand] <- sample(c(FALSE,TRUE),length(a_cand),TRUE,c(0.9,0.1))
    }
    
    #compute current quantities
    x_bar <- mean(x[which(!B)])
    x2_bar <- mean((x[which(!B)])^2)
    s_x <- sum(x[which(!B)])
    s_x2 <- sum((x[which(!B)])^2)
    s_B <- sum(B)
    
    #compute unconstrained MLE
    mu_mle <- x_bar
    sigma2_mle <- mean((x-mu_mle)^2)
    p_mle <- s_B/n
    
    if (checkConstraint(a,b,mu_mle,sigma2_mle,p_mle,P_mean,w) | ignore_AFR){
      #Case where MLE is good
      mu <- mu_mle
      sigma2 <- sigma2_mle
      p <- p_mle
    }
    else{
      #Case where MLE is not good
      constrained_params <- paramOptimizationBinary(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,epsilon)
      mu <- constrained_params$mu
      sigma2 <- constrained_params$sigma2
      p <- constrained_params$p
    }
    
    scores[a_cand] <- scores[a_cand] + dnorm(mu,mu,sqrt(sigma2))-dnorm(x[a_cand],mu,sqrt(sigma2))
  }
  return(scores/n_Bs)
}
dimensionRandomAfr <- function(X,n_estimators=10,ignore_AFR=FALSE,epsilon=1e-8,q1=0.24,q2=0.75){
  scores <- matrix(NaN,nrow(X),ncol(X))
  for (j in 1:ncol(X)){
    cat(paste0(j))
    x <- X[,j]
    scores[,j] <- randomAfr(x,n_Bs = n_estimators,ignore_AFR=ignore_AFR,epsilon=epsilon,q1=q1,q2=q2)
  }
  cat("\n")
  return(rowMeans(scores))
}
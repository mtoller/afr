theOne <- function(x,a,b,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,
                   optim_stop=5e+4,verb=FALSE,...){
  if (verb){
    dev.new()
    plot(x)
  }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  n <- length(x)
  
  a_cand <- which(x < a | x > b)
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
  
  
  if (is.null(B)){
    B <- rep(0,n)
    B[a_cand] <- 1
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  mus <- seq(min(x),max(x),length.out=optim_stop)
  
  
  Bs <- list()
  
  for (i in 1:10){
    
    #Assume constraint not active
    mu <- mean(x[which(!B)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    p <- sum(B)/n
    
    #test if within constraints
    
    a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    if ((a_prob - P_mean)^2 <= w^2){
      if (verb){
        print("Within all constrains, returning standard maximum likelihood estimates")
      }
      
      mu <- mean(x[which(B==0)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      p <- sum(B)/n
    }
    else if ((a_prob_c1 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=0. This should never happen.")
      }
      p <- 0
      B <- rep(FALSE,n)
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      
      
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      
      
      
      p <- P_mean + w
      
      B <- rep(FALSE,n)
      B[a_cand] <- TRUE
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
    }
    else{
      if (verb){
        print("Not within integral constrains")
      }
      
      best_params <- optimizeParams(x,B,a,b,mus,P_mean,w,tol)
      
      mu <- best_params[1]
      sigma2 <- best_params[2]
      p <- best_params[3]
      
    }
    if (verb){
      print(paste0("s_B=",sum(B)))
      print(paste0("x_bar=",mean(x[which(!B)])))
      print(paste0("x2_bar=",mean(x[which(!B)]^2)))
      print(paste0("mu_mle=",mean(x[which(!B)])))
      print(paste0("sigma2_mle=",mean((x[which(!B)]-mean(x[which(!B)]))^2)))
      print(paste0("p_mle=",sum(B)/n))
      print(paste0("mu=",mu))
      print(paste0("sigma2=",sigma2))
      print(paste0("P_mean+w=",P_mean+w))
      print(paste0("p=",p))
      print(paste0("sB/(n-omega)=",sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))))
      print(paste0("omega_mu=",getOmega(x,B,a,b,mu,sigma2,0,info=TRUE)))
      print(paste0("omega_sigma=",getOmega(x,B,a,b,mu,sigma2,1)))
      print(paste0("n-s_B/p=",(n-sum(B)/p)))
      print(paste0("n/(n-omega)=",n/(n- getOmega(x,B,a,b,mu,sigma2))))
      val <- (1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))-P_mean)^2
      print(paste0("Constraint okay? ", val <= w^2, " ",val, "<=",w^2))
      print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
    }
    
    
    
    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
    
    
    B[a_cand] <- flat_prior_const*dnorm(mu,mu,sqrt(sigma2)) >=
      dnorm(x[a_cand],mu,sqrt(sigma2))
    
    B <- as.logical(B)
    
    
    if (all(B==B_old)){
      if (verb){
        print("Converged")
      }
      
      break
    }
    
    if (any(unlist(lapply(Bs,function(j)all(j==B))))){
      if (verb){
        print("circle detected, aborting")
      }
      
      break
    }
    
  }
  if (verb){
    points(which(1:n %in% a_cand & !B),x[which(1:n %in% a_cand & !B)],col='darkblue',pch=3)
    points(which(B),x[which(B)],col='red',pch=4)
    lines(c(1,n),c(b,b),col='cyan',lwd=0.75)
    lines(c(1,n),c(a,a),col='green',lwd=0.75)
    
    legend("topleft",legend = c("Data","Normals outside R","Anomalies","Upper bound","Lower Bound"),
           col=c("black","darkblue","red","cyan","green"),pch=c(1,3,4,NA,NA),lty=c(NA,NA,NA,1,1),lwd=c(1,1,1,0.75,0.75))
    
  }
  
  return(as.logical(B))
  
}
optimizeParams <- function(x,B,a,b,mus,P_mean,w,tol){
  n <- length(x)
  setSigmas <- function(index=1){
    if (is.finite(sigma2s[index])){
      if (sigma2s[index] > tol){
        val <- sum(B)/(n - getOmega(x,B,a,b,mu,sigma2s[index]))
        
        if (!is.finite(val)){
          ps[index] <<- NA
        }
        else{
          if (val > 0 & val < (P_mean+w)){
            ps[index] <<- val
          }
          else{
            ps[index] <<- NA
          }
        }
        
      }
      else{
        ps[index] <<- NA
        sigma2s[index] <<- NA
      }
    }
    else{
      ps[index] <<- NA
      sigma2s[index] <<- NA
    }
  }
  
  best_mu <- NaN
  best_sigma2 <- NaN
  best_e <- Inf
  best_p <- NaN
  
  es <- rep(NA,9)
  ps <- rep(NA,9)
  ps[4:6] <- 0
  ps[7:9] <- P_mean+w
  
  for (mu in mus){
    
    es[1:9] <- NA
    ps[1:3] <- NA
    
    sigma2s <- theLampertSolutions(x,B,mu,a,b)
    setSigmas(1)
    setSigmas(2)
    setSigmas(3)
    sigma2s <- rep(sigma2s,3)
    
    
    for (j in 1:9){
      if (is.finite(ps[j]) & is.finite(sigma2s[j])){
        es[j] <- (1-(1-ps[j]) *
                    (pnorm(b,mu,sqrt(sigma2s[j])) - pnorm(a,mu,sqrt(sigma2s[j]))) - P_mean)^2-w^2
      }
      else{
        es[j] <- Inf
      }
    }
    
    if(all(is.infinite(es))){
      min_e <- Inf
    }
    else{
      min_e <- min(es)
    }
    
    
    if (min_e < best_e){
      best_mu <- mu
      best_sigma2 <- sigma2s[which.min(es)]
      best_e <- min_e
      best_p <- ps[which.min(es)]
      
    }
    
    
    
  }
  return(c(best_mu,best_sigma2,best_p))
}

theLampertSolutions <- function(x,B,mu,a,b,prec=5){
  x_bar <- mean(x[which(!B)])
  x2_bar <- mean(x[which(!B)]^2)
  
  alpha <- x2_bar - mu*x_bar + (mu-x_bar)*a
  beta <- x2_bar - mu*x_bar + (mu-x_bar)*b
  
  t <- 1/beta
  m <- 0.5*((a-mu)^2-(b-mu)^2)
  
  r <- -alpha/beta
  e_power <- exp(-m*t)
  
  inp <- m*(alpha/beta^2-1/beta)*e_power
  
  lamp <- rLambert(inp,r*e_power,prec)
  
  z <- 1/beta + 1/m*lamp
  
  sigma2 <- 1/z
  
  return(sigma2)
  
}

getOmega <- function(x,B,a,b,mu,sigma2,modus=0,info=FALSE){
  n <- length(x)
  
  sum_B <- sum(B)
  sum_x <- sum(x[which(!B)])
  sum_x2 <- sum(x[which(!B)]^2)
  s <- sqrt(sigma2)
  
  I <- pnorm(b,mu,s)-pnorm(a,mu,s)
  
  dldmu <- (sum_x - (n-sum_B)*mu)/sigma2
  dlds <- (sum_x2 - 2 * mu * sum_x + (n-sum_B)*(mu^2-sigma2))/s^3
  
  dIdmu <- (exp(-(a-mu)^2/(2*sigma2))-exp(-(b-mu)^2/(2*sigma2)))/(sqrt(2*pi)*s)
  dIds <- ((a-mu)*exp(-(a-mu)^2/(2*sigma2))-(b-mu)*exp(-(b-mu)^2/(2*sigma2)))/(sqrt(2*pi)*sigma2)
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
  
  if (modus==1){
    return(omega_s)
  }
  
  ##omega_s should always be the same as omega_mu. Switch only used for debugging
  
  return(omega_mu)
  
}
getLogLikelihood <- function(x,B,mu,sigma2,p,a,b,flat_prior_const=1e-8){
  n <- length(x)
  
  s_B <- sum(B)
  s_x <- sum(x[which(!as.logical(B))])
  s_x2 <- sum(x[which(!as.logical(B))]^2)
  
  result <- (n-s_B)*(log(1-p)-log(sqrt(2*pi*sigma2))) + s_B * safeLog(p) + s_B * safeLog(flat_prior_const)
  return(result - (s_x2 - 2*mu*s_x + (n-s_B)*mu^2)/(2*sigma2))
}
safeLog <- function(x){
  if (x == 0){
    0
  }
  else if (x < 0){
    NaN
  }
  else{
    log(x)
  }
  
}

# thePCAOne <- function(X,prop,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,verb=FALSE,...){
thePCAOne <- function(X,prop,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,verb=FALSE,...){
  n <- nrow(X)
  d <- ncol(X)
  
  
  # X <- apply(X,2,scale)
  
  # x <- pcaPP::PCAproj(X)$scores[,1]
  x <- X
  m <- quantile(x,0.5,type=1)
  md <- quantile(abs(x-m),0.5,type=1)
  a <- m-md*0.1
  b <- m+md*0.1
  theOneProp(x,a,b,verb=verb,flat_prior_const = flat_prior_const,optim_stop=1e+3,prop=prop,alpha = alpha)
}

theOneProp <- function(x,a,b,prop,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,
                   optim_stop=1e+3,verb=FALSE,...){
  if (verb){
    dev.new()
    plot(x)
  }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  n <- length(x)
  
  a_cand <- which(x < a | x > b)
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
  
  if (verb){
    print(paste0("P_hat=",P_hat))
    print(paste0("z=",z))
    print(paste0("P_mean=",P_mean))
    print(paste0("w",w))
  }
  
  if (is.null(B)){
    B <- rep(0,n)
    B[a_cand] <- 1
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  mus <- seq(min(x),max(x),length.out=optim_stop)
  
  
  Bs <- list()
  
  for (i in 1:10){
    
    #Assume constraint not active
    mu <- mean(x[which(!B)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    p <- sum(B)/n
    
    #test if within constraints
    
    a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    if ((a_prob - P_mean)^2 <= w^2){
      if (verb){
        print("Within all constrains, returning standard maximum likelihood estimates")
      }
      
      mu <- mean(x[which(B==0)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      p <- sum(B)/n
    }
    else if ((a_prob_c1 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=0. This should never happen.")
      }
      p <- 0
      B <- rep(FALSE,n)
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      
      
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      
      
      
      p <- P_mean + w
      
      B <- rep(FALSE,n)
      B[a_cand] <- TRUE
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
    }
    else{
      if (verb){
        print("Not within integral constrains")
      }
      
      best_params <- optimizeParams(x,B,a,b,mus,P_mean,w,tol)
      
      mu <- best_params[1]
      sigma2 <- best_params[2]
      p <- best_params[3]
      
    }
    if (verb){
      print(paste0("s_B=",sum(B)))
      print(paste0("x_bar=",mean(x[which(!B)])))
      print(paste0("x2_bar=",mean(x[which(!B)]^2)))
      print(paste0("mu_mle=",mean(x[which(!B)])))
      print(paste0("sigma2_mle=",mean((x[which(!B)]-mean(x[which(!B)]))^2)))
      print(paste0("p_mle=",sum(B)/n))
      print(paste0("mu=",mu))
      print(paste0("sigma2=",sigma2))
      print(paste0("P_mean+w=",P_mean+w))
      print(paste0("p=",p))
      print(paste0("sB/(n-omega)=",sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))))
      print(paste0("omega_mu=",getOmega(x,B,a,b,mu,sigma2,0,info=TRUE)))
      print(paste0("omega_sigma=",getOmega(x,B,a,b,mu,sigma2,1)))
      print(paste0("n-s_B/p=",(n-sum(B)/p)))
      print(paste0("n/(n-omega)=",n/(n- getOmega(x,B,a,b,mu,sigma2))))
      val <- (1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))-P_mean)^2
      print(paste0("Constraint okay? ", val <= w^2, " ",val, "<=",w^2))
      print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
    }
    
    
    
    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
    
    
    B[a_cand] <- 1:length(a_cand) %in% order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:prop]
    
    B <- as.logical(B)
    
    
    if (all(B==B_old)){
      if (verb){
        print("Converged")
      }
      
      break
    }
    
    if (any(unlist(lapply(Bs,function(j)all(j==B))))){
      if (verb){
        print("circle detected, aborting")
      }
      
      break
    }
    
  }
  if (verb){
    points(which(1:n %in% a_cand & !B),x[which(1:n %in% a_cand & !B)],col='darkblue',pch=3)
    points(which(B),x[which(B)],col='red',pch=4)
    lines(c(1,n),c(b,b),col='cyan',lwd=0.75)
    lines(c(1,n),c(a,a),col='green',lwd=0.75)
    
    legend("topleft",legend = c("Data","Normals outside R","Anomalies","Upper bound","Lower Bound"),
           col=c("black","darkblue","red","cyan","green"),pch=c(1,3,4,NA,NA),lty=c(NA,NA,NA,1,1),lwd=c(1,1,1,0.75,0.75))
    
  }
  
  return(as.logical(B))
  
}
theExpOne <- function(X,y,prop,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,verb=FALSE,...){
  n <- nrow(X)
  d <- ncol(X)
  
  k <- sum(y)
  # X <- apply(X,2,scale)
  
  # x <- pcaPP::PCAproj(X)$scores[,1]
  x <- X
  
  m <- x[which(!y)[sample(1:k,1)]]
  a <- m-0.0001
  b <- m+0.0001
  theOneProp(x,a,b,verb=verb,flat_prior_const = flat_prior_const,optim_stop=1e+3,prop=prop,alpha = alpha)
}
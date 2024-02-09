theOne <- function(x,a,b,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,
                   optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  if (verb){
    dev.new()
    plot(x,ylim=c(min(c(min(x),a)),max(c(max(x),b))))
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
    #B[a_cand] <- 1
    set.seed(0)
    B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,prob = c(0.05,0.95))
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  mus <- seq(min(x),max(x),length.out=optim_stop)
  
  
  Bs <- list()
  
  for (i in 1:100){
    
    #Assume constraint not active
    mu <- mean(x[which(!B)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    p <- sum(B)/n
    
    
    
    #test if within constraints
    
    a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    
    if (unconstrained){
      a_prob <- P_mean
    }
    
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
      # B <- rep(FALSE,n)
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      
      
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      p <- P_mean + w
      
      # B <- rep(FALSE,n)
      # B[a_cand] <- TRUE
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
    
    
    B[a_cand] <- dnorm(x[a_cand],mu,sqrt(sigma2)) <= flat_prior_const*dnorm(mu,mu,sqrt(sigma2))
    
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
  final_scores <- dnorm(mu,mu,sqrt(sigma2))-dnorm(x,mu,sqrt(sigma2))
  final_scores[which(x >= a & x <= b)] <- 0
  
  return(final_scores)
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
  
  x_bar <- mean(x[which(!B)])
  x2_bar <- mean(x[which(!B)]^2)
  
  for (mu in mus){
    
    es[1:9] <- NA
    ps[1:3] <- NA
    
    
    
    sigma2s <- theLambertSolutionsC(x2_bar,x_bar,mu,a,b)
    setSigmas(1)
    setSigmas(2)
    setSigmas(3)
    # ps[1:9] <- 0
    # setSigmas(4)
    # setSigmas(5)
    # setSigmas(6)
    # setSigmas(7)
    # setSigmas(8)
    # setSigmas(9)
    sigma2s <- rep(sigma2s,3)
    
    
    for (j in 1:9){
      if (is.finite(ps[j]) & is.finite(sigma2s[j])){
        es[j] <- (1-(1-ps[j]) *
                    (pnorm(b,mu,sqrt(sigma2s[j])) - pnorm(a,mu,sqrt(sigma2s[j]))) - P_mean)^2-w^2
        if (is.finite(es[j])){
        }
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
  # print(P_mean)
  # print(w)
  # print(best_e)
  # print(best_mu)
  # print(best_sigma2)
  # print(best_p)
  # stop()
  # if (!(is.finite(best_e))){
  #   return(NaN)
  # }
  # return(best_e)
  return(c(best_mu,best_sigma2,best_p,best_e))
}



theLambertSolutions <- function(x,B,mu,a,b,prec=16){
  x_bar <- mean(x[which(!B)])
  x2_bar <- mean(x[which(!B)]^2)
  
  alpha <- x2_bar - mu*x_bar + (mu-x_bar)*a
  beta <- x2_bar - mu*x_bar + (mu-x_bar)*b
  
  t <- 1/beta
  m <- 0.5*((a-mu)^2-(b-mu)^2)
  
  r <- -alpha/beta
  e_power <- exp(-m*t)
  
  inp <- m*(alpha/beta^2-1/beta)*e_power
  print(paste0("x2_bar=",x2_bar))
  print(paste0("x_bar=",x_bar))
  print(paste0("alpha=",alpha))
  print(paste0("beta=",beta))
  print(paste0("m=",m))
  # print(m*(alpha/beta^2-1/beta))
  # print(-m*t)
  # print(e_power)
  print(inp)
  print(r*e_power)
  
  lamb <- rLambert(inp,r*e_power,prec)
  print(lamb)
  z <- 1/beta + 1/m*lamb
  
  sigma2 <- 1/z
  # sigma2 <- beta*m/(m+beta*lamb)
  print(sigma2)
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
    print(paste0("diff=",abs(omega_mu-omega_s)))
  }
  
  if (modus==1){
    return(omega_s)
  }
  
  ##omega_s should always be the same as omega_mu. Switch only used for debugging
  
  return(omega_mu)
  
}

getOmega2 <- function(s_x,s_x2,s_B,n,a,b,mu,sigma2,modus=0,info=FALSE){
  s <- sqrt(sigma2)
  
  I <- pnorm(b,mu,s)-pnorm(a,mu,s)
  
  dldmu <- (s_x - (n-s_B)*mu)/sigma2
  dlds <- (s_x2 - 2 * mu * s_x + (n-s_B)*(mu^2-sigma2))/s^3
  
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
    print(paste0("diff",abs(omega_mu-omega_s)))
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
thePCAOne <- function(X,B=NULL,alpha=0.05,flat_prior_const=0.05,tol=1e-8,verb=FALSE,...){
  n <- nrow(X)
  d <- ncol(X)
  
  
  X <- apply(X,2,scale)
  
  x <- pcaPP::PCAproj(X)$scores[,1]
  # x <- X
  m <- quantile(x,0.5,type=1)
  md <- quantile(abs(x-m),0.5,type=1)
  a <- m-md*1
  b <- m+md*1
  # theOne(x,a,b,B=B,verb=verb,flat_prior_const = flat_prior_const,optim_stop=1e+3,alpha = alpha,...)
  theOtherOne(x,a,b,B=B,verb=verb,optim_stop=1e+3,alpha = alpha,...)
}

# thePCAOne <- function(X,prop,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,verb=FALSE,...){
theMultiPCAOne <- function(X,B=NULL,alpha=0.05,flat_prior_const=0.05,tol=1e-8,verb=FALSE,...){
  n <- nrow(X)
  d <- ncol(X)
  
  
  X <- apply(X,2,scale)
  
  temp <- pcaPP::PCAproj(X,k=min(c(ncol(X),5)))
  scaled <- temp$scores
  weights <- temp$sdev
  scores <- scaled*0
  for (s in 1:ncol(scaled)){
    x <- scaled[,s]
    m <- quantile(x,0.5,type=1)
    md <- quantile(abs(x-m),0.5,type=1)
    a <- m-md*0.01
    b <- m+md*0.01
    scores[,s] <- theOtherOne(x,a,b,B=B,verb=verb,optim_stop=1e+3,alpha = alpha,...)
    scores[,s] <- scores[,s] * weights[s] / sum(weights)
    #theOtherOne(x,a,b,B=B,verb=verb,optim_stop=1e+3,alpha = alpha,...)
  }
  apply(scores,1,sum)
}

thePCAOneSimple <- function(X,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,verb=FALSE,...){
  n <- nrow(X)
  d <- ncol(X)
  
  
  X <- apply(X,2,scale)
  
  x <- pcaPP::PCAproj(X)$scores[,1]
  # x <- X
  m <- quantile(x,0.5,type=1)
  md <- quantile(abs(x-m),0.5,type=1)
  a <- m-md*0.001
  b <- m+md*0.001
  theSimpleOne(x,a,b,B=B,verb=verb,flat_prior_const = flat_prior_const,optim_stop=1e+3,alpha = alpha,...)
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
theExpOne <- function(x,a,b,B=NULL,alpha=0.05,flat_prior_const=1e-7,tol=1e-8,
                      optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  if (verb){
    dev.new()
    plot(x,ylim=c(min(c(min(x),a)),max(c(max(x),b))))
  }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  x <- (x-min(x))/(max(x)-min(x))
  
  n <- length(x)
  
  a_cand <- which(x < a | x > b)
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
  
  
  if (is.null(B)){
    B <- rep(0,n)
    # B[a_cand] <- 1
    set.seed(0)
    B[a_cand] <- sample(c(0,1),length(a_cand),TRUE)
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  kappas <- seq(0,30,length.out=optim_stop)
  
  
  Bs <- list()
  
  for (i in 1:100){
    
    #Assume constraint not active
    kappa <- 1/mean(x[which(!B)])
    p <- sum(B)/n
    
    
    
    #test if within constraints
    
    a_prob <- 1-(1-p)*(pexp(b,kappa)-pexp(a,kappa))
    a_prob_c1 <- 1-(1-0)*(pexp(b,kappa)-pexp(a,kappa))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pexp(b,kappa)-pnorm(a,kappa))
    
    if (unconstrained){
      a_prob <- P_mean
    }
    
    if ((a_prob - P_mean)^2 <= w^2){
      if (verb){
        print("Within all constrains, returning standard maximum likelihood estimates")
      }
      
      kappa <- 1/mean(x[which(B==0)])
      p <- sum(B)/n
    }
    else if ((a_prob_c1 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=0. This should never happen.")
      }
      p <- 0
      # B <- rep(FALSE,n)
      kappa <- 1/mean(x[which(B==FALSE)])
      
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      p <- P_mean + w
      
      # B <- rep(FALSE,n)
      # B[a_cand] <- TRUE
      kappa <- 1/mean(x[which(B==FALSE)])
    }
    else{
      if (verb){
        print("Not within integral constrains")
      }
      best_params <- optimizeKappa(x,B,a,b,kappas,P_mean,w,tol)
      
      kappa <- best_params[1]
      p <- best_params[2]
      
    }
    if (verb){
      print(paste0("s_B=",sum(B)))
      print(paste0("x_bar=",mean(x[which(!B)])))
      print(paste0("x2_bar=",mean(x[which(!B)]^2)))
      print(paste0("kappa_mle=",1/mean(x[which(!B)])))
      print(paste0("p_mle=",sum(B)/n))
      print(paste0("kappa=",kappa))
      print(paste0("P_mean+w=",P_mean+w))
      print(paste0("p=",p))
      print(paste0("n-s_B/p=",(n-sum(B)/p)))
      val <- (1-(1-p)*(pexp(b,kappa)-pexp(a,kappa))-P_mean)^2
      print(paste0("Constraint okay? ", val <= w^2, " ",val, "<=",w^2))
    }
    
    
    
    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
    
    
    B[a_cand] <- dexp(x[a_cand],kappa) <= flat_prior_const*dexp(0,kappa)
    
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
  
  # return(as.logical(B))
  final_scores <- dexp(0,kappa)-dexp(x,kappa)
  final_scores[which(x >= a & x <= b)] <- 0
  
  return(final_scores)
}

optimizeKappa <- function(x,B,a,b,kappas,P_mean,w,tol){
  n <- length(x)
  
  best_kappa <- NaN
  best_e <- Inf
  best_p <- NaN
  
  es <- rep(NA,3)
  ps <- rep(NA,3)
  ps[2] <- 0
  ps[3] <- P_mean+w
  
  ers <- c()
  for (kappa in kappas){
    
    es[1:3] <- NA
    ps[1] <- NA
    
    if (kappa > tol){
      val <- sum(B)/(n - (sum(B)/kappa - sum(x[which(!B)]))*
                       (exp(-a*kappa)-exp(-b*kappa))/
                       (b*(exp(-b*kappa)-a*exp(-a*kappa))))
      
      if (!is.finite(val)){
        ps[1] <- NA
      }
      else{
        if (val > 0 & val < (P_mean+w)){
          ps[1] <- val
        }
        else{
          ps[1] <- NA
        }
      }
    }
    else{
      ps[1] <- NA
    }
      
    for (j in 1:3){
      if (is.finite(ps[j])){
        es[j] <- (1-(1-ps[j]) * (pexp(b,kappa) - pexp(a,kappa,)) - P_mean)^2-w^2
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
      best_kappa <- kappa
      best_e <- min_e
      best_p <- ps[which.min(es)]
      
    }
    ers <- c(ers,min_e)
  }
  plot.ts(ers)
  return(c(best_kappa,best_p))
}

theSimpleOne <- function(x,a,b,alpha=0.05,tol=1e-8,
                   optim_stop=5e+4,verb=FALSE,...){
  if (verb){
    dev.new()
    plot(x,ylim=c(min(c(min(x),a)),max(c(max(x),b))))
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
  
  
  B <- rep(FALSE,n)

    
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
    # B <- rep(FALSE,n)
    mu <- mean(x[which(B==FALSE)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    
    
  }
  else if ((a_prob_c2 - P_mean)^2 <= w^2){
    if (verb){
      print("Within integral constraints if p=P_mean + w. This should never happen.")
    }
    p <- P_mean + w
    
    # B <- rep(FALSE,n)
    # B[a_cand] <- TRUE
    mu <- mean(x[which(B==FALSE)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
  }
  else{
    if (verb){
      print("Not within integral constrains")
    }
    mus <- seq(min(x),max(x),length.out=optim_stop)
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
  }
    
    
  
  if (verb){
    points(which(1:n %in% a_cand & !B),x[which(1:n %in% a_cand & !B)],col='darkblue',pch=3)
    points(which(B),x[which(B)],col='red',pch=4)
    lines(c(1,n),c(b,b),col='cyan',lwd=0.75)
    lines(c(1,n),c(a,a),col='green',lwd=0.75)
    
    legend("topleft",legend = c("Data","Normals outside R","Anomalies","Upper bound","Lower Bound"),
           col=c("black","darkblue","red","cyan","green"),pch=c(1,3,4,NA,NA),lty=c(NA,NA,NA,1,1),lwd=c(1,1,1,0.75,0.75))
    
  }
  
  final_scores <- dnorm(mu,mu,sqrt(sigma2))-dnorm(x,mu,sqrt(sigma2))
  final_scores[which(x >= a & x <= b)] <- 0
  
  return(final_scores)
}

theOtherOne <- function(x,a,b,B=NULL,alpha=0.05,tol=1e-8,
                   optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  if (verb){
    dev.new()
    plot(x,ylim=c(min(c(min(x),a)),max(c(max(x),b))))
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
    # B[a_cand] <- 1
    set.seed(0)
    B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,c(0.95,0.05))
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  mus <- seq(min(x),max(x),length.out=optim_stop)
  
  
  Bs <- list()
  
  for (i in 1:1000){
    
    #Assume constraint not active
    mu <- mean(x[which(!B)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    p <- sum(B)/n
    
    
    
    #test if within constraints
    
    a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    
    if (unconstrained){
      a_prob <- P_mean
    }
    
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
      # B <- rep(FALSE,n)
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      
      
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      p <- P_mean + w
      
      # B <- rep(FALSE,n)
      # B[a_cand] <- TRUE
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
      # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
    }
    
    
    
    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
    
    
    p_con <- sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))
    
    B[a_cand] <- 0
    
    if (p_con > 1/n){
      B[a_cand][order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:trunc(p_con*n)]] <- 1
    }
    
    
    
    # B[a_cand] <- dnorm(x[a_cand],mu,sqrt(sigma2)) <= flat_prior_const*dnorm(mu,mu,sqrt(sigma2))
    
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
  
  # return(as.logical(B))
  final_scores <- dnorm(mu,mu,sqrt(sigma2))-dnorm(x,mu,sqrt(sigma2))
  final_scores[which(x >= a & x <= b)] <- 0
  
  return(final_scores)
}
theCopulaOne <- function(X,B=NULL,alpha=0.05,tol=1e-8,
                                        optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  # if (verb){
  #   dev.new()
  #   plot(X[,1],ylim=c(min(c(min(X),as[1])),max(c(max(X),bs[1]))))
  # }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  
  n <- nrow(X)
  d <- ncol(X)
  
  U <- matrix(NA,n,d)
  V <- matrix(NA,n,d)
  W <- matrix(NA,n,d)
  
  for (j in 1:d){
    if (d > 50){
      print(paste0(j, " out of ",d))
    }
    
    x <- X[,j]
    
    m <- quantile(x,0.5,type=1)
    md <- quantile(abs(x-m),0.5,type=1)
    if (md == 0){
      md <- 0.01
    }
    a <- m-md*0.01
    b <- m+md*0.01
    
    
    a_cand <- which(x < a | x > b)
    P_hat <- length(a_cand)/n
    
    z <- qnorm(1-alpha/2,0,1)
    
    P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
    w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
    
    
    if (is.null(B)){
      B <- rep(0,n)
      # B[a_cand] <- 1
      set.seed(0)
      B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,c(0.95,0.05))
    }
    B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
    B <- as.logical(B)
    mus <- seq(min(x),max(x),length.out=optim_stop)
    
    
    Bs <- list()
    
    for (i in 1:1000){
      
      #Assume constraint not active
      mu <- mean(x[which(!B)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      p <- sum(B)/n
      
      
      
      #test if within constraints
      
      a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      
      if (unconstrained){
        a_prob <- P_mean
      }
      
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
        # B <- rep(FALSE,n)
        mu <- mean(x[which(B==FALSE)])
        sigma2 <- mean((x[which(!B)]-mu)^2)
        
        
      }
      else if ((a_prob_c2 - P_mean)^2 <= w^2){
        if (verb){
          print("Within integral constraints if p=P_mean + w. This should never happen.")
        }
        p <- P_mean + w
        
        # B <- rep(FALSE,n)
        # B[a_cand] <- TRUE
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
        # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
      }
      
      
      
      B_old <- B
      Bs[[length(Bs)+1]] <- B_old
      
      
      p_con <- sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))
      
      B[a_cand] <- 0
      
      if (p_con > 1/n){
        B[a_cand][order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:trunc(p_con*n)]] <- 1
      }
      
      
      
      # B[a_cand] <- dnorm(x[a_cand],mu,sqrt(sigma2)) <= flat_prior_const*dnorm(mu,mu,sqrt(sigma2))
      
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
    
    # return(as.logical(B))
    U[,j] <- pnorm(x,mu,sqrt(sigma2),log.p = TRUE)
    U[which(x >= a & x <= b),j] <- 0
    
    V[,j] <- pnorm(-x,-mu,sqrt(sigma2),log.p = TRUE)
    V[which(x >= a & x <= b),j] <- 0
    
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

theSimpleCopulaOne <- function(X,B=NULL,alpha=0.05,tol=1e-8,
                         optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  # if (verb){
  #   dev.new()
  #   plot(X[,1],ylim=c(min(c(min(X),as[1])),max(c(max(X),bs[1]))))
  # }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  
  n <- nrow(X)
  d <- ncol(X)
  
  U <- matrix(NA,n,d)
  V <- matrix(NA,n,d)
  W <- matrix(NA,n,d)
  
  for (j in 1:d){
    if (d > 50){
      print(paste0(j, " out of ",d))
    }
    
    x <- X[,j]
    
    m <- quantile(x,0.5,type=1)
    md <- quantile(abs(x-m),0.5,type=1)
    if (md == 0){
      md <- 0.01
    }
    a <- m-md*0.01
    b <- m+md*0.01
    
    
    a_cand <- which(x < a | x > b)
    P_hat <- length(a_cand)/n
    
    z <- qnorm(1-alpha/2,0,1)
    
    P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
    w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
    
    
    B <- rep(FALSE,n)
    
    mus <- seq(min(x),max(x),length.out=optim_stop)
    
    
    #Assume constraint not active
    mu <- mean(x[which(!B)])
    sigma2 <- mean((x[which(!B)]-mu)^2)
    p <- sum(B)/n
      
      
      
    #test if within constraints
      
    a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
    a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      
    if (unconstrained){
      a_prob <- P_mean
    }
      
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
      # B <- rep(FALSE,n)
      mu <- mean(x[which(B==FALSE)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
    }
    else if ((a_prob_c2 - P_mean)^2 <= w^2){
      if (verb){
        print("Within integral constraints if p=P_mean + w. This should never happen.")
      }
      p <- P_mean + w
      
      # B <- rep(FALSE,n)
      # B[a_cand] <- TRUE
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
      # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
    }
    
    U[,j] <- pnorm(x,mu,sqrt(sigma2),log.p = TRUE)
    U[which(x >= a & x <= b),j] <- 0
    
    V[,j] <- pnorm(-x,-mu,sqrt(sigma2),log.p = TRUE)
    V[which(x >= a & x <= b),j] <- 0
    
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
skewness <- function(x){
  m <- mean(x)
  (1/(length(x))*sum((x-m)^3)) / (1/length(x)*sum((x-m)^2))^(3/2)
}
theGeigerCopulaOne <- function(X,B=NULL,alpha=0.05,tol=1e-8,tau=0.5,
                         optim_stop=5e+4,verb=FALSE,unconstrained=FALSE,proportion=NULL,...){
  # if (verb){
  #   dev.new()
  #   plot(X[,1],ylim=c(min(c(min(X),as[1])),max(c(max(X),bs[1]))))
  # }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  
  n <- nrow(X)
  d <- ncol(X)
  
  U <- matrix(NA,n,d)
  V <- matrix(NA,n,d)
  W <- matrix(NA,n,d)
  
  as <- rep(NA,d)
  bs <- rep(NA,d)
  
  a_cand <- matrix(FALSE,n,d)
  
  for (j in 1:d){
    if (d > 50){
      print(paste0(j, " out of ",d))
    }
    
    x <- X[,j]
    
    m <- quantile(x,0.5,type=1)
    md <- quantile(abs(x-m),0.5,type=1)
    if (md == 0){
      md <- 0.01
    }
    as[j] <- m-md*7#*0.01
    bs[j] <- m+md*7#*0.01
    a_cand[,j] <- x < as[j] | x > bs[j]
  }
  
  a_cand <- which(apply(a_cand,1,any))
  P_hat <- length(a_cand)/n
    
  z <- qnorm(1-alpha/2,0,1)
    
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
    
    
  if (is.null(B)){
    B <- rep(0,n)
    # B[a_cand] <- 1
    set.seed(0)
    B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,c(0.95,0.05))
  }
  B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
  B <- as.logical(B)
  
  
    mus <- seq(min(x),max(x),length.out=optim_stop)
    
    
  Bs <- list()
    
  for (i in 1:1000){
    for (j in 1:d){
      a <- as[j]
      b <- bs[j]
      #Assume constraint not active
      mu <- mean(x[which(!B)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      p <- sum(B)/n
      
      
      #test if within constraints
      
      a_prob <- 1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      a_prob_c1 <- 1-(1-0)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      a_prob_c2 <- 1-(1-(P_mean+w))*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))
      
      if (unconstrained){
        a_prob <- P_mean
      }
      
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
        # B <- rep(FALSE,n)
        mu <- mean(x[which(B==FALSE)])
        sigma2 <- mean((x[which(!B)]-mu)^2)
        
        
      }
      else if ((a_prob_c2 - P_mean)^2 <= w^2){
        if (verb){
          print("Within integral constraints if p=P_mean + w. This should never happen.")
        }
        p <- P_mean + w
        
        # B <- rep(FALSE,n)
        # B[a_cand] <- TRUE
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
        # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
      }
      
      U[,j] <- pnorm(x,mu,sqrt(sigma2),log.p = TRUE)
      U[which(x >= a & x <= b),j] <- 0
      
      V[,j] <- pnorm(-x,-mu,sqrt(sigma2),log.p = TRUE)
      V[which(x >= a & x <= b),j] <- 0
      
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
    
    scores <- apply(cbind(pl,pr,ps),1,max)

    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
      
      
    # p_con <- sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))
    # 
    # B[a_cand] <- 0
    #   
    # if (p_con > 1/n){
    #   B[a_cand][order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:trunc(p_con*n)]] <- 1
    # }
      
      
      
    B[a_cand] <- scores[a_cand] >= tau
    
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
    
    # return(as.logical(B))
    
    
  
  
  pl <- apply(U,1,function(u)-sum(u))
  pr <- apply(V,1,function(u)-sum(u))
  ps <- apply(W,1,function(u)-sum(u))
  
  return(apply(cbind(pl,pr,ps),1,max))
  
}
theSimplestOne <- function(X,B=NULL,alpha=0.05,tol=1e-8,
                               optim_stop=5e+4,verb=FALSE,AFR=NULL,...){
  # if (verb){
  #   dev.new()
  #   plot(X[,1],ylim=c(min(c(min(X),as[1])),max(c(max(X),bs[1]))))
  # }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  U <- matrix(NA,n,d)

  as <- rep(NA,d)
  bs <- rep(NA,d)
  
  a_cand <- matrix(FALSE,n,d)
  
  if (is.null(AFR)){
    for (j in 1:d){
      if (d > 50){
        print(paste0(j, " out of ",d))
      }
      
      x <- X[,j]
      
      as[j] <- quantile(x,0.25)
      bs[j] <- quantile(x,0.75)
      a_cand[,j] <- x < as[j] | x > bs[j]
    }
  }
  else{
    
    as <- AFR$as
    bs <- AFR$bs
    
    for (j in 1:d){
      x <- X[,j]
      a_cand[,j] <- x < as[j] | x > bs[j]
    }
  }
  print(as)
  print(bs)
  a_cand <- which(apply(a_cand,1,any))
  if (length(a_cand)==n){
    a_cand <- c()
  }
  ins <- which(!((1:n)%in% a_cand))
  P_hat <- length(a_cand)/n
  
  z <- qnorm(1-alpha/2,0,1)
  
  P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
  w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
  
  if (is.null(B)){
    B <- rep(0,n)
    B[a_cand] <- 1
    # set.seed(0)
    # B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,c(0.5,0.5))
  }
  
  B[ins] <- 0
  B <- as.logical(B)
  
  
  Bs <- list()
  
  for (iter in 1:1000){
    U <- matrix(NA,n,d)
    ps <- rep(NA,d)
    
    for (j in 1:d){
      x <- X[,j]
      mus <- seq(min(x),max(x),length.out=optim_stop)
      
      a <- as[j]
      b <- bs[j]
      #Assume constraint not active
      mu <- mean(x[which(!B)])
      sigma2 <- mean((x[which(!B)]-mu)^2)
      p <- sum(B)/n
      
      if (sigma2 == 0){
        U[,j] <- 0
        next
      }
      
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
        # B <- rep(FALSE,n)
        mu <- mean(x[which(B==FALSE)])
        sigma2 <- mean((x[which(!B)]-mu)^2)
        
        
      }
      else if ((a_prob_c2 - P_mean)^2 <= w^2){
        if (verb){
          print("Within integral constraints if p=P_mean + w. This should never happen.")
        }
        p <- P_mean + w
        
        # B <- rep(FALSE,n)
        # B[a_cand] <- TRUE
        mu <- mean(x[which(B==FALSE)])
        sigma2 <- mean((x[which(!B)]-mu)^2)
      }
      else{
        if (verb){
          print("Not within integral constraints")
        }
        best_params <- optimizeParams(x,B,a,b,mus,P_mean,w,tol)
        mu <- best_params[1]
        sigma2 <- best_params[2]
        p <- best_params[3]
        
      }
      if (verb){
        print(paste0("dimension=",j," out of ",d))
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
        # print(paste0("sB/(n-omega)=",sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))))
        # print(paste0("omega_mu=",getOmega(x,B,a,b,mu,sigma2,0,info=TRUE)))
        # print(paste0("omega_sigma=",getOmega(x,B,a,b,mu,sigma2,1)))
        # print(paste0("n-s_B/p=",(n-sum(B)/p)))
        # print(paste0("n/(n-omega)=",n/(n- getOmega(x,B,a,b,mu,sigma2))))
        # val <- (1-(1-p)*(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))-P_mean)^2
        # print(paste0("Constraint okay? ", val <= w^2, " ",val, "<=",w^2))
        # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
      }
      # if (skewness(x) < 0 ){
      #   U[,j] <- pnorm(x,mu,sqrt(sigma2),log.p = TRUE)
      # }
      # else{
      #   U[,j] <- pnorm(-x,-mu,sqrt(sigma2),log.p = TRUE)
      # }
      
      U[a_cand,j] <- dnorm(x[a_cand],mu,sqrt(sigma2),log = TRUE)
      
      # U[,j] <- pnorm(x,mu,sqrt(sigma2))
      
      # U[which(x >= a & x <= b),j] <- -Inf
      U[ins,j] <- Inf
      ps[j] <- p
    }
    if (verb){
      # print(U)
      print(B)
    }
    
    score <- apply(U,1,sum)
    # score <- rep(NaN,n)
    # for (i in 1:n){
    #   active_row <- U[i,]
    #   score[i] <- mean(apply(U,1,function(u)all(u<=active_row)))
    # }
    # score <- abs(score-0.5)
    # score <- apply(U,1,function(u)-sum(u))
    
    p_all <- mean(ps[which(is.finite(ps))])
    
    if (verb){
      print(paste0("p_all=",p_all))
    }
    
    B_old <- B
    Bs[[length(Bs)+1]] <- B_old
    
    
    # p_con <- sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))
    # 
    print(score)
    B[a_cand] <- 0

    if (p_all > 1/n){
      B[a_cand][order(score[a_cand])[1:trunc(p_all*n)]] <- 1
    }
    
    
    #B[a_cand] <- scores[a_cand] >= tau
    
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
  
  # return(as.logical(B))
  
  if (verb){
    plot(score)
    points(ins,score[ins],col='blue',pch=2)
    print(B)
  }

  return(exp(-score))
  
}

theTrueCopulaOne <- function(X,B=NULL,alpha=0.05,tol=1e-8,
                         optim_stop=5e+4,verb=FALSE,...){
  # if (verb){
  #   dev.new()
  #   plot(X[,1],ylim=c(min(c(min(X),as[1])),max(c(max(X),bs[1]))))
  # }
  
  if (alpha==0 | alpha==1){
    return(rep(FALSE,length(x)))
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  U <- matrix(-Inf,n,d)
  #V <- matrix(NA,n,d)
  #W <- matrix(NA,n,d)
  
  for (j in 1:d){
    if (d > 50){
      print(paste0(j, " out of ",d))
    }
    
    x <- X[,j]
    
    #m <- quantile(x,0.5,type=1)
    #md <- quantile(abs(x-m),0.5,type=1)
    #if (md == 0){
    #  md <- 0.01
    #}
    a <- quantile(x,0.25)#m-md*0.01
    b <- quantile(x,0.75)#m+md*0.01
    
    
    a_cand <- which(x < a | x > b)
    P_hat <- length(a_cand)/n
    
    z <- qnorm(1-alpha/2,0,1)
    
    P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
    w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
    
    if(verb){
      print(paste0("n_a_cand=",length(a_cand)))
      print(paste0("n=",n))
    }
    
    if (is.null(B)){
      B <- rep(0,n)
      # B[a_cand] <- 1
      set.seed(0)
      B[a_cand] <- sample(c(0,1),length(a_cand),TRUE,c(0.05,0.95))
    }
    B[(1:n)[which(!(1:n %in% a_cand))]] <- 0
    B <- as.logical(B)
    mus <- seq(min(x),max(x),length.out=optim_stop)
    
    
    Bs <- list()
    
    for (i in 1:1000){
      
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
        # B <- rep(FALSE,n)
        mu <- mean(x[which(B==FALSE)])
        sigma2 <- mean((x[which(!B)]-mu)^2)
        
        
      }
      else if ((a_prob_c2 - P_mean)^2 <= w^2){
        if (verb){
          print("Within integral constraints if p=P_mean + w. This should never happen.")
        }
        p <- P_mean + w
        
        # B <- rep(FALSE,n)
        # B[a_cand] <- TRUE
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
        # print(paste0("neg. loglik=",-getLogLikelihood(x,B,mu,sigma2,p,a,b,flat_prior_const)))
      }
      
      
      
      B_old <- B
      Bs[[length(Bs)+1]] <- B_old
      
      
      p_con <- sum(B)/(n-getOmega(x,B,a,b,mu,sigma2))
      
      B[a_cand] <- 0
      
      if (p_con > 1/n){
        B[a_cand][order(dnorm(x[a_cand],mu,sqrt(sigma2)))[1:trunc(p_con*n)]] <- 1
      }
      
      
      
      # B[a_cand] <- dnorm(x[a_cand],mu,sqrt(sigma2)) <= flat_prior_const*dnorm(mu,mu,sqrt(sigma2))
      
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
    
    # return(as.logical(B))
    
    if (skewness(x) < 0){
      U[a_cand,j] <- pnorm(x[a_cand],mu,sqrt(sigma2),log.p = TRUE)
    }
    else{
      U[a_cand,j] <- pnorm(-x[a_cand],-mu,sqrt(sigma2),log.p = TRUE)
    }
    
  }
  score <- rep(NaN,n)
  for (i in 1:n){
    active_row <- U[i,]
    score[i] <- mean(apply(U,1,function(u)all(u<=active_row)))
  }
  
  return(score)
  
}


constrainedMLE <- function(x,B,a,b,mu){
  x_bar <- mean(x[which(!B)])
  x2_bar <- mean((x[which(!B)])^2)
  theLambertSolutionsC(x2_bar,x_bar,mu,a,b)
}
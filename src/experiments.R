processArff <- function(file){
  require(farff)
  dataset <- farff::readARFF(file,convert.to.logicals = FALSE)
  return(list(X=as.matrix(dataset[,1:(ncol(dataset)-2)]),y=dataset$outlier=='yes'))
}

evalMetrics <- function(predicted,actual,beta=1,verbose=FALSE){
  
  n <- length(predicted)
  
  tp <- length(which(predicted & actual))
  fp <- length(which(predicted & !actual))
  fn <- length(which(!predicted & actual))
  tn <- length(which(!predicted & !actual))
  
  if (verbose){
    print(c(tp=tp,fp=fp,fn=fn,tn=tn))
  }
  
  accuracy <- (tp+tn)/(tp+fp+tn+fn)
  if (tp == 0){
    precision <- 0
    recall <- 0
    f <- 0
    mcc <- 0
    return(list(tp=tp,fp=fp,fn=fn,tn=tn,
                accuracy=accuracy,precision=precision,recall=recall,f_score=f,mcc=mcc))
  }
  else{
    precision <- tp/(tp+fp)
    recall <- tp/(tp+fn)
    f <- (1+beta^2) * precision * recall / (beta^2 * precision + recall)
  }
  
  if (tn == 0)
  {
    return(list(tp=tp,fp=fp,fn=fn,tn=tn,
                accuracy=accuracy,precision=precision,recall=recall,f_score=f,mcc=0))
  }
  mcc <- (tp*tn - fp*fn)/sqrt(as.numeric(tp+fp)*as.numeric(tp+fn)*as.numeric(tn+fp)*as.numeric(tn+fn))
  
  return(list(tp=tp,fp=fp,fn=fn,tn=tn,accuracy=accuracy,precision=precision,recall=recall,f_score=f,mcc=mcc))
}

macroF1 <- function(predicted,actual){
  (evalMetrics(predicted,actual)$f_score +
     evalMetrics(!predicted,!actual)$f_score)/2
}

phiCoefficient <- function(predicted,actual){
  evalMetrics(predicted,actual)$mcc
}

simulateGaussUniAnomaly <- function(n,p,mu,sigma2,a,b,uni_min,uni_max){
  B <- as.logical(rbinom(n,1,p))
  s_B <- sum(B)
  
  data_n <- rnorm(n-s_B,mu,sqrt(sigma2))
  data_inf <- c()
  
  while (length(data_inf) < s_B){
    anomaly <- runif(1,uni_min,uni_max)
    if (anomaly < a | anomaly > b){
      data_inf <- c(data_inf,anomaly)
    }
  }
  
  
  x <- c(data_inf,data_n)
  
  P <- (1-p)*(1-(pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2)))) + p
  
  
  ordering <- sample(1:n,n)
  B <- ordering %in% (1:length(data_inf))
  
  return(list(x=x[ordering],B=B,P=P))
}

benchmarkExperiment <- function(seed=0){
  
  
  set.seed(seed)
  
  
  detectors <- list(LOF=lof,
                    SVDD=ocsvm,
                    IForest=iforest,
                    DeepSVDD=deepSvdd,
                    COPOD=copod,
                    ours=thePCAOne)
  
  measures <- list(macroF1=macroF1,
                   phiCoefficient=phiCoefficient)
  
  
  
  datasets <- list()
  datasets[[1]] <- processArff("glass.arff")
  datasets[[2]] <- processArff("heart.arff")
  datasets[[3]] <- processArff("hep.arff")
  datasets[[4]] <- processArff("ion.arff")
  datasets[[5]] <- processArff("lymph.arff")
  datasets[[6]] <- processArff("park.arff")
  datasets[[7]] <- processArff("pima.arff")
  datasets[[8]] <- processArff("stamps.arff")
  datasets[[9]] <- processArff("shuttle.arff")
  datasets[[10]] <- processArff("wpbc.arff")
  
  result_table <- matrix(NA,length(datasets),length(detectors)*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures))
  
  
  for (j in 1:length(datasets)){
    print(j)
    
    
    x <- datasets[[j]]$X
    x <- apply(x,2,scale)
    x <- pcaPP::PCAproj(x)$scores[,1]
    
    B_true <- datasets[[j]]$y
    true_a_prop <- sum(B_true)
    
    Bs <- list()
    for (i in 1:length(detectors)){
      print(paste0("Computing ",names(detectors)[i],"..."))
      tryCatch(
        {Bs[[i]] <- detectors[[i]](as.matrix(x),bw=sqrt(sigma2),prop=true_a_prop)},
        error=function(cond){
          print("error")
          Bs[[i]] <<- rep(FALSE,length(B_true))
        }
      )
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)] <- measures[[t]](Bs[[i]],B_true)
      }
      
    }
    print(result_table)
    print(colMeans(result_table,na.rm = TRUE))
    
  }
  
  final_result_table <- matrix(NA,nrow = length(measures),ncol = length(detectors))
  colnames(final_result_table) <- names(detectors)
  
  for (j in 1:length(measures)){
    final_result_table[j,] <- colMeans(result_table[,(1:length(detectors))
                                                    +length(detectors)*(j-1)])
    
  }
  print(colMeans(result_table))
  return(result_table)
}

officeExperiment <- function(seed=0){
  
  set.seed(seed)
  
  detectors <- list(LOF=lof,
                    SVDD=ocsvm,
                    SVDD_PU=ocsvmPu,
                    IForest=iforest,
                    DeepSVDD=deepSvdd,
                    DeepSVDD_Pu=deepSvddPu,
                    COPOD=copod,
                    Song=song,EM=uniGaussEmHard,
                    ours=theOneProp)
  
  measures <- list(macroF1=macroF1,phiCoefficient=phiCoefficient)
  
  result_table <- matrix(NA,1,length(detectors)*2*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures)*2)
  
  for (j in 1:1){
    print(j)
    
    
    
    office <- read.csv("office.csv")
    
    x <- office$time_deviation
    B_true <- office$is_anomaly==1
    a <- -29
    b <- 29
    true_a_prop <- sum(B_true)
    wBninR <- which(x < a | x > b)
    
    Bs <- list()
    for (i in 1:length(detectors)){
      print(paste0("Computing ",names(detectors)[i],"..."))
      Bs[[i]] <- detectors[[i]](as.matrix(x),a=a,b=b,prop=true_a_prop)
      
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)*2] <- measures[[t]](Bs[[i]],B_true)
        result_table[j,i+length(detectors)*(t*2-1)] <- measures[[t]](Bs[[i]][wBninR],B_true[wBninR])
      }
      
      
    }
    
    
    
  }
  
  final_result_table <- matrix(NA,nrow = length(measures)*2,ncol = length(detectors))
  colnames(final_result_table) <- names(detectors)
  
  for (j in 1:length(measures)){
    final_result_table[2*(j-1)+1,] <- result_table[,(1:length(detectors))
                                                   +length(detectors)*(j-1)*2]
    final_result_table[2*j,] <- result_table[,(1:length(detectors))
                                             +length(detectors)*(j*2-1)]
    
  }
  
  row.names(final_result_table) <- unlist(lapply(names(measures),function(nam){
    c(nam,paste0(nam," outside R"))
  }))
  return(final_result_table)
}

simulationExperiment <- function(n_sims=100,seed=0){
  
  set.seed(seed)
  
  detectors <- list(LOF=lof,
                    SVDD=ocsvm,
                    SVDD_PU=ocsvmPu,
                    IForest=iforest,
                    DeepSVDD=deepSvdd,
                    DeepSVDD_Pu=deepSvddPu,
                    COPOD=copod,
                    Song=song,EM=uniGaussEmHard,
                    ours=theOneProp)
  
  measures <- list(macroF1=macroF1,phiCoefficient=phiCoefficient)
  
  result_table <- matrix(NA,n_sims,length(detectors)*2*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures)*2)
  
  for (j in 1:n_sims){
    print(j)
    
    
    
    p <- runif(1,0.005,0.02)
    mu <- runif(1,-5,5)
    sigma2 <- runif(1,0.01,2)
    fact <- runif(1,0.01,3)
    a <- mu-fact*sqrt(sigma2)
    b <- mu+fact*sqrt(sigma2)
    
    sim <- simulateGaussUniAnomaly(1000,p,mu,sigma2,a,b,
                                   uni_min = mu-10*sqrt(sigma2),
                                   uni_max = mu+10*sqrt(sigma2))
    
    x <- sim$x
    B_true <- sim$B
    true_a_prop <- sum(B_true)
    P <- sim$P
    
    wBninR <- which(x < a | x > b)
    
    Bs <- list()
    for (i in 1:length(detectors)){
      print(paste0("Computing ",names(detectors)[i],"..."))
      Bs[[i]] <- detectors[[i]](as.matrix(x),a=a,b=b,bw=sqrt(sigma2),prop=true_a_prop)
      
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)*2] <- measures[[t]](Bs[[i]],B_true)
        result_table[j,i+length(detectors)*(t*2-1)] <- measures[[t]](Bs[[i]][wBninR],B_true[wBninR])
      }
      
      
    }
    
    
    
  }
  
  final_result_table <- matrix(NA,nrow = length(measures)*2,ncol = length(detectors))
  colnames(final_result_table) <- names(detectors)
  
  for (j in 1:length(measures)){
    final_result_table[2*(j-1)+1,] <- colMeans(result_table[,(1:length(detectors))
                                                            +length(detectors)*(j-1)*2])
    final_result_table[2*j,] <- colMeans(result_table[,(1:length(detectors))
                                                      +length(detectors)*(j*2-1)])
    
  }
  
  row.names(final_result_table) <- unlist(lapply(names(measures),function(nam){
    c(nam,paste0(nam," outside R"))
  }))
  print(colMeans(result_table))
  print(apply(result_table,2,sd))
  return(final_result_table)
}

sensitivityExperiment <- function(n_sims=100,n_params=100,seed=0){
  set.seed(seed)
  
  taus <- 10^seq(-1,-10,length.out=n_params)
  
  
  
  measures <- list(macroF1=macroF1)
  
  result_list <- rep(0,length(taus))
  
  for (j in 1:n_sims){
    print(j)
    
    
    
    p <- runif(1,0.005,0.02)
    mu <- runif(1,-5,5)
    sigma2 <- runif(1,0.01,5)
    fact <- runif(1,0.1,3)
    a <- mu-fact*sqrt(sigma2)
    b <- mu+fact*sqrt(sigma2)
    
    sim <- simulateGaussUniAnomaly(1000,p,mu,sigma2,a,b,
                                   uni_min = mu-10*sqrt(sigma2),
                                   uni_max = mu+10*sqrt(sigma2))
    
    x <- sim$x
    B_true <- sim$B
    P <- sim$P
    
    
    for (i in 1:length(taus)){
      tau <- taus[i]
      print(paste0("tau=",tau))
      
      B <- theOne(x,a,b,flat_prior_const=tau,optim_stop = 1e+3)
      result_list[i] <- result_list[i] + macroF1(B,B_true)
    }
    
    
    
  }
  plot(taus,result_list/n_sims,type='l',log='x',ylim=c(0,1))
  
  return(result_list/n_sims)
  
}
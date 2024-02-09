processArff <- function(file){
  require(farff)
  dataset <- farff::readARFF(file,convert.to.logicals = FALSE)
  return(list(X=as.matrix(dataset[,1:(ncol(dataset)-2)]),y=dataset$outlier=='yes'))
}

processNPZ <- function(file){
  if (!exists("np")){
    np <<- import("numpy")
  }
  
  temp <- np$load(file)
  return(list(X=temp$f[["X"]],y=as.logical(temp$f[["y"]])))
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

rocAUC <- function(scores,actual){
  ROCit::rocit(scores,actual)$AUC
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
  
  
  # detectors <- list(LOF=lof,
  #                   SVDD=ocsvm,
  #                   IForest=iforest,
  #                   DeepSVDD=deepSvdd,
  #                   COPOD=copod,
  #                   ours=thePCAOne)
  #detectors <- list(COPOD=copod,IForest=iforest,ours=thePCAOne,ours_copula=theCopulaOne,
                    # ours_simple=thePCAOneSimple,benchGauss=benchGauss,
                    # ours_simple_copula=theSimpleCopulaOne,benchGaussCop=benchGaussCopula)
                    # ours_geiger=theGeigerCopulaOne)
  detectors <- list(LOF=lof,
                    IForest=iforest,
                    COPOD=copod,
                    EM=function(x,...)dimensionEmAfr(x,...,ignore_AFR = TRUE),
                    random=function(X,...)dimensionRandomAfr(X,ignore_AFR=TRUE,...),randomAfr=dimensionRandomAfr,
                    autoEnc=autoEnc,DeepSVDD=deepSvdd)
  
  measures <- list(rocAUC=rocAUC)
  
  
  
  datasets <- list()
  
  #development datasets
  # datasets[[1]] <- processArff("../old_datasets/glass.arff")
  # datasets[[2]] <- processArff("../old_datasets/heart.arff")
  # datasets[[3]] <- processArff("../old_datasets/hep.arff")
  # datasets[[4]] <- processArff("../old_datasets/ion.arff")
  # datasets[[5]] <- processArff("../old_datasets/lymph.arff")
  # datasets[[6]] <- processArff("../old_datasets/park.arff")
  # datasets[[7]] <- processArff("../old_datasets/pima.arff")
  # datasets[[8]] <- processArff("../old_datasets/stamps.arff")
  # datasets[[9]] <- processArff("../old_datasets/shuttle.arff")
  # datasets[[10]] <- processArff("../old_datasets/wpbc.arff")

  
  for (filename in dir("../datasets/")){
    datasets[[length(datasets)+1]] <- processNPZ(paste0("../datasets/",filename))
  }

  result_table <- matrix(NA,length(datasets),length(detectors)*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures))
  
  times_table <- matrix(Inf,length(datasets),length(detectors))
  colnames(times_table) <- names(detectors)
  
  for (j in 1:length(datasets)){
    # print(j)
    # if (!(j %in% c(2:4,10,19,26:29))){
    #   next
    # }
    
    
    x <- datasets[[j]]$X
    #shuffle dataset indices to avoid trivial solutions
    set.seed(seed)
    sam <- sample(1:nrow(x),nrow(x))
    x <- x[sam,]
    
    
    # x <- apply(x,2,scale)
    # x <- pcaPP::PCAproj(x)$scores[,1]
    
    B_true <- datasets[[j]]$y[sam]
    
    Bs <- list()
    for (i in 1:length(detectors)){
      set.seed(seed)
      print(paste0("Computing ",names(detectors)[i],"..."))
      tryCatch(
        {
          t_before <- as.numeric(Sys.time())  
          Bs[[i]] <- detectors[[i]](as.matrix(x))
          t_after <- as.numeric(Sys.time())
        },
        error=function(cond){
          print(paste0("error: ",cond))
          Bs[[i]] <<- rep(0,length(B_true))
        }
      )
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)] <- measures[[t]](Bs[[i]],B_true)
      }
      times_table[j,i] <- t_after-t_before
      
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
  
  print(rowMeans(apply(result_table,1,function(r)rank(-r,TRUE,"last"))))
  print(times_table)
  return(list(results=result_table,times=times_table))
}

semiBenchmarkExperiment <- function(seed=0){
  set.seed(seed)
  
  
  # detectors <- list(LOF=lof,
  #                   SVDD=ocsvm,
  #                   IForest=iforest,
  #                   DeepSVDD=deepSvdd,
  #                   COPOD=copod,
  #                   ours=thePCAOne)
  #detectors <- list(COPOD=copod,IForest=iforest,ours=thePCAOne,ours_copula=theCopulaOne,
  # ours_simple=thePCAOneSimple,benchGauss=benchGauss,
  # ours_simple_copula=theSimpleCopulaOne,benchGaussCop=benchGaussCopula)
  # ours_geiger=theGeigerCopulaOne)
  
  detectors <- list(prenet=prenet)#,overlap=overlap)
  
  measures <- list(rocAUC=rocAUC)
  
  datasets <- list()
  
  
  for (filename in dir("../datasets/")){
    datasets[[length(datasets)+1]] <- processNPZ(paste0("../datasets/",filename))
  }
  
  datasets <- datasets[c(2:4,10,19,26:29)]
  result_table <- matrix(NA,length(datasets),length(detectors)*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures))
  
  times_table <- matrix(Inf,length(datasets),length(detectors))
  colnames(times_table) <- names(detectors)
  
  for (j in 1:length(datasets)){
    
    x <- datasets[[j]]$X
    #shuffle dataset indices to avoid trivial solutions
    set.seed(seed)
    sam <- sample(1:nrow(x),nrow(x))
    x <- x[sam,]
    
    B_true <- datasets[[j]]$y[sam]
    
    #train test split
    splitted <- makeSplit(x,B_true,training_proportion = 0.05)
    X_train <- splitted$X_train
    y_train <- splitted$y_train
    X_test <- splitted$X_test
    y_test <- splitted$y_test
    
    Bs <- list()
    for (i in 1:length(detectors)){
      set.seed(seed)
      print(paste0("Computing ",names(detectors)[i],"..."))
      tryCatch(
        {
          t_before <- as.numeric(Sys.time())  
          Bs[[i]] <- detectors[[i]](X_train,y_train,X_test)
          t_after <- as.numeric(Sys.time())
          },
        error=function(cond){
          print(paste0("error: ",cond))
          Bs[[i]] <<- rep(0,length(y_test))
        }
      )
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)] <- measures[[t]](Bs[[i]],y_test)
      }
      times_table[j,i] <- t_after-t_before
      
    }
    print(result_table)
    print(colMeans(result_table,na.rm = TRUE))
    
  }
  
  # final_result_table <- matrix(NA,nrow = length(measures),ncol = length(detectors))
  # colnames(final_result_table) <- names(detectors)
  # 
  # for (j in 1:length(measures)){
  #   final_result_table[j,] <- colMeans(result_table[,(1:length(detectors))
  #                                                   +length(detectors)*(j-1)])
  #   
  # }
  
  #print(rowMeans(apply(result_table,1,function(r)rank(-r,TRUE,"last"))))
  
  print(times_table)
  return(list(results=result_table,times=times_table))
}

officeExperiment <- function(seed=0){
  
  
  
  detectors <- list(LOF=lof,
                    SVDD=ocsvm,
                    # SVDD_PU=ocsvmPu,
                    IForest=iforest,
                    DeepSVDD=deepSvdd,
                    # DeepSVDD_Pu=deepSvddPu,
                    COPOD=copod,
                    EM=function(x,...)emAfr(x,...,ignore_AFR = TRUE),
                    ours=emAfr,
                    random=function(X,...)randomAfr(X,ignore_AFR=TRUE,...),randomAfr=randomAfr)
  
  measures <- list(rocAUC=rocAUC)
  
  result_table <- matrix(NA,1,length(detectors)*2*length(measures))
  colnames(result_table) <- rep(names(detectors),length(measures)*2)
  
  times_table <- matrix(Inf,1,length(detectors))
  colnames(times_table) <- names(detectors)
  
  for (j in 1:1){
    print(j)
    
    
    
    office <- read.csv("../our_datasets/office.csv")
    
    x <- office$time_deviation
    B_true <- office$is_anomaly==1
    a <- -29
    b <- 29
    true_a_prop <- sum(B_true)
    wBninR <- which(x < a | x > b)
    
    Bs <- list()
    for (i in 1:length(detectors)){
      set.seed(seed)
      print(paste0("Computing ",names(detectors)[i],"..."))
      t_before <- as.numeric(Sys.time())
      Bs[[i]] <- detectors[[i]](as.matrix(x),a=a,b=b)
      t_after <- as.numeric(Sys.time())
      
      for (t in 1:length(measures)){
        result_table[j,i+length(detectors)*(t-1)*2] <- measures[[t]](Bs[[i]],B_true)
        result_table[j,i+length(detectors)*(t*2-1)] <- measures[[t]](Bs[[i]][wBninR],B_true[wBninR])
      }
      times_table[j,i] <- t_after-t_before
      
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
  print(times_table)
  return(final_result_table)
}

officeExperimentNew <- function(){
  set.seed(0)
  office <- read.csv("../our_datasets/office.csv")
  
  x <- office$time_deviation
  X <- matrix(x,ncol=1)
  
  y <- office$is_anomaly==1
  a <- -29
  b <- 29
  
  split10 <- makeSplit(X,y,0.1)
  split20 <- makeSplit(X,y,0.2)
  split50 <- makeSplit(X,y,0.5)
  
  result_table <- matrix(NA,1,7)
  
  print("computing prenet...")
  result_table[1,1] <- rocAUC(prenet(split10$X_train,split10$y_train,split10$X_test),split10$y_test)
  result_table[1,2] <- rocAUC(prenet(split20$X_train,split20$y_train,split20$X_test),split20$y_test)
  result_table[1,3] <- rocAUC(prenet(split50$X_train,split50$y_train,split50$X_test),split50$y_test)
  
  print("computing overlap...")
  result_table[1,4] <- rocAUC(overlap(split10$X_train,split10$y_train,split10$X_test),split10$y_test)
  result_table[1,5] <- rocAUC(overlap(split20$X_train,split20$y_train,split20$X_test),split20$y_test)
  result_table[1,6] <- rocAUC(overlap(split50$X_train,split50$y_train,split50$X_test),split50$y_test)
  
  print("computing CAMLE...")
  result_table[1,7] <- rocAUC(randomAfr(x,a,b),y)
  
  
  return(result_table)
}

simulationExperiment <- function(n_sims=100,seed=0){
  
  set.seed(seed)
  
  detectors <- list(LOF=lof,
                    SVDD=ocsvm,
                    # SVDD_PU=ocsvmPu,
                    IForest=iforest,
                    DeepSVDD=deepSvdd,
                    # DeepSVDD_Pu=deepSvddPu,
                    COPOD=copod,
                    EM=function(x,...)emAfr(x,...,ignore_AFR = TRUE),
                    ours=emAfr,
                    random=function(X,...)randomAfr(X,ignore_AFR=TRUE,...),randomAfr=randomAfr)
  
  measures <- list(rocAUC=rocAUC)
  
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
    P <- sim$P
    
    wBninR <- which(x < a | x > b)
    
    Bs <- list()
    for (i in 1:length(detectors)){
      print(paste0("Computing ",names(detectors)[i],"..."))
      Bs[[i]] <- detectors[[i]](as.matrix(x),a=a,b=b)
      
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

trueSimulationExperiment <- function(n_params=10,n_samples=10,n_data=1000,n_B=10,seed=0){
  set.seed(seed)
  mus <- runif(n_params,-5,5)
  sigmas <- runif(n_params,0.1,2)
  ps <- runif(n_params,0.05,0.95)
  as <- mus-sigmas*0.98
  bs <- mus+sigmas
  
  x <- rep(NaN,n_data)
  
  results <- matrix(NA,n_params*n_samples,15)
  colnames(results) <- c("mu_true","mu_mle","mu_mle_B","mu_c","mu_c_B",
                         "sigma_true","sigma_mle","sigma_mle_B","sigma_c","sigma_c_B",
                         "p_true","p_mle","p_mle_B","p_c","p_c_B")
  for (i in 1:n_params){
    for (j in 1:n_samples){
      
      
      B <- sample(c(FALSE,TRUE),n_data,TRUE,c(1-ps[i],ps[i]))
      s_B <- sum(B)
      while (s_B >= (n_data-1) | s_B <= 1){
        B <- sample(c(FALSE,TRUE),n_data,TRUE,c(1-ps[i],ps[i]))
        s_B <- sum(B)
      }
      
      x[1:(floor(s_B/2))] <- runif(floor(s_B/2),mus[i]-10*sigmas[i],as[i])
      x[(floor(s_B/2)+1):s_B] <- runif(ceiling(s_B/2),bs[i],mus[i]+10*sigmas[i])
      x[(s_B+1):n_data] <- rnorm(n_data-s_B,mus[i],sigmas[i])
      B <- c(rep(TRUE,s_B),rep(FALSE,n_data-s_B))
      a <- as[i]
      b <- bs[i]
      
      n <- length(x)
      
      a_cand <- which(x < a | x > b)
      P_hat <- length(a_cand)/n
      
      alpha <- 0.05
      z <- qnorm(1-alpha/2,0,1)
      P_mean <- 1/(1+z^2/n)*(P_hat+z^2/(2*n))
      w <- z/(1+z^2/n)*sqrt(P_hat*(1-P_hat)/n+z^2/(4*n^2))
      
      x_bar <- mean(x[which(!B)])
      x2_bar <- mean((x[which(!B)])^2)
      s_x <- sum(x[which(!B)])
      s_x2 <- sum((x[which(!B)])^2)
      
      mu_mle <- x_bar
      sigma2_mle <- mean((x[which(!B)]-mu_mle)^2)
      p_mle <- s_B/n
      
      constrained_params <- paramOptimizationBinary(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,1e-8)
      # paramOptimization(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w)
      # invisible(readline(prompt="Press [enter] to continue"))
      mu_c <- constrained_params$mu
      sigma2_c <- constrained_params$sigma2
      p_c <- constrained_params$p
      # 
      # dldmu <- (s_x - (n-s_B)*mu)/sigma2
      # dIdmu <- (exp(-(a-mu)^2/(2*sigma2))-exp(-(b-mu)^2/(2*sigma2)))/(sqrt(2*pi)*sqrt(sigma2))
      # I <- pnorm(b,mu,sqrt(sigma2))-pnorm(a,mu,sqrt(sigma2))
      # 
      # omega <- dldmu/dIdmu*I
      
      mu_mle_B <- 0
      sigma2_mle_B <- 0
      p_mle_B <- 0
      
      mu_c_B <- 0
      sigma2_c_B <- 0
      p_c_B <- 0
      
      for (k in 1:n_B){
        B <- rep(FALSE,n_data)
        B[a_cand] <- sample(c(FALSE,TRUE),length(a_cand),TRUE,c(1-ps[i],ps[i]))
        s_B <- sum(B)
        while (s_B >= length(a_cand) | s_B <= 1){
          B[a_cand] <- sample(c(FALSE,TRUE),length(a_cand),TRUE,c(1-ps[i],ps[i]))
          s_B <- sum(B)
        }
        x_bar <- mean(x[which(!B)])
        x2_bar <- mean((x[which(!B)])^2)
        s_x <- sum(x[which(!B)])
        s_x2 <- sum((x[which(!B)])^2)
        
        mu_mle_B <- mu_mle_B + x_bar
        sigma2_mle_B <-sigma2_mle_B + mean((x[which(!B)]-mu_mle)^2)
        p_mle_B <-p_mle_B + s_B/n
        
        constrained_params <- paramOptimizationBinary(s_x,s_x2,s_B,x_bar,x2_bar,n,a,b,P_mean,w,1e-8)
        
        mu_c_B <- mu_c_B + constrained_params$mu
        sigma2_c_B <- sigma2_c_B + constrained_params$sigma2
        p_c_B <- p_c_B + constrained_params$p
      }
      
      position <- (i-1)*n_samples+j
      results[position,1] <- mus[i]
      results[position,2] <- mu_mle
      results[position,3] <- mu_mle_B/n_B
      results[position,4] <- mu_c
      results[position,5] <- mu_c_B/n_B
      
      results[position,6] <- sigmas[i]
      results[position,7] <- sqrt(sigma2_mle)
      results[position,8] <- sqrt(sigma2_mle_B/n_B)
      results[position,9] <- sqrt(sigma2_c)
      results[position,10] <- sqrt(sigma2_c_B/n_B)
      
      results[position,11] <- ps[i]
      results[position,12] <- p_mle
      results[position,13] <- p_mle_B/n_B
      results[position,14] <- p_c
      results[position,15] <- p_c_B/n_B
      
      # if (position == 45){
      #   plot(x)
      #   points(which(B),x[which(B)],col='red',pch=4)
      #   lines(c(-1000,1000),c(a,a),col='blue')
      #   lines(c(-1000,1000),c(b,b),col='green')
      # }
    }
  }
  results <- as.data.frame(results)
  
  results_final <- rep(NaN,12)
  
  results_final[1] <- mad((results$mu_true - results$mu_mle),na.rm = TRUE)
  results_final[2] <- mad((results$mu_true - results$mu_c),na.rm = TRUE)
  results_final[3] <- mad((results$sigma_true - results$sigma_mle),na.rm = TRUE)
  results_final[4] <- mad((results$sigma_true - results$sigma_c),na.rm = TRUE)
  results_final[5] <- mad((results$p_true - results$p_mle),na.rm = TRUE)
  results_final[6] <- mad((results$p_true - results$p_c),na.rm = TRUE)
  
  results_final[7] <- mad((results$mu_true - results$mu_mle_B),na.rm = TRUE)
  results_final[8] <- mad((results$mu_true - results$mu_c_B),na.rm = TRUE)
  results_final[9] <- mad((results$sigma_true - results$sigma_mle_B),na.rm = TRUE)
  results_final[10] <- mad((results$sigma_true - results$sigma_c_B),na.rm = TRUE)
  results_final[11] <- mad((results$p_true - results$p_mle_B),na.rm = TRUE)
  results_final[12] <- mad((results$p_true - results$p_c_B),na.rm = TRUE)
  
  names(results_final) <- c("mu_mle","mu_c","sigma_mle","sigma_c","p_mle","p_c",
                            "mu_mle_B","mu_c_B","sigma_mle_B","sigma_c_B","p_mle_B","p_c_B")
  
  return(list(table=results,summary=round(results_final,4)))
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
sensitivityStudy <- function(filepath,n_AFR=11){
  dataset <- processNPZ(filepath)
  X <- dataset$X
  y <- dataset$y
  
  offsets <- seq(0.01,0.99,length=n_AFR)
  
  min_a <- 0.0001
  max_b <- 0.999
  
  results <- rep(NaN,n_AFR)
  
  for (i in 1:length(offsets)){
    offset <- offsets[i]
    as <- seq(min_a,max_b-offset,length=n_AFR)
    bs <- seq(min_a+offset,max_b,length=n_AFR)
    
    results[i] <- 0
    for (j in 1:n_AFR){
      set.seed(0)
      results[i] <- results[i] + rocAUC(dimensionRandomAfr(X,q1=as[j],q2=bs[j]),y)
    }
    results[i] <- results[i]/n_AFR
  }
  
  
  val <- rocAUC(dimensionRandomAfr(X),y)
  cairo_pdf(paste0(filepath,"_results.pdf"),4,4)
  plot(offsets,results,type='b',xlab='Offset Δ',ylab='Average AUC-ROC',ylim=c(0,1))
  lines(c(-10,10),c(val,val),col='blue',lty=2)
  dev.off()
  return(list(offsets=offsets,results=results))
}
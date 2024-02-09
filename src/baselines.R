abScale <- function(x,a=0,b=1){
  a+(x-min(x))*(b-a)/(max(x)-min(x))
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
Rcpp::sourceCpp("rLambert_mt_hh.cpp")
library(reticulate)
use_condaenv("anomalyMLE")

source("baselines.R")
# source("ours.R")
source("oursNew.R")
source("experiments.R")
source("utils.R")
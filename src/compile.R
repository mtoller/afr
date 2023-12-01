Rcpp::sourceCpp("rLambert_mt.cpp")
library(reticulate)
use_condaenv("anomalyMLE")

source("baselines.R")
source("ours.R")
source("experiments.R")
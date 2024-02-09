This document explains how to reproduce the results presented in the article.
In the following, it is assumed that Anaconda and R are installed.

Step1: Install python libraries via conda - execute this code in terminal

`conda create --name anomalyMLE --file requirements.txt`

`conda activate anomalyMLE`

`conda install -c conda-forge pyod`
`conda install -c conda-forge tensorflow`


Step2: Install reticulate and link to conda env - execute this code in R launched from this directory.

`install.packages("Rcpp")`
`install.packages("reticulate")`

`source("compile.R")`


Step3: Checkout the external repositories into the folder `./src/semi`.
We include anonymous links to our forks of the original code since we had to modify the interface slightly.

(https://anonymous.4open.science/r/PReNet-8191)

(https://anonymous.4open.science/status/Overlap-D13D)


Step4: Re-run experiments - execute this code in R. Note that it is unavoidable that some results are different from those reported in the paper due to local randomness, e.g., exact package versions, floating point imprecision.

`benchmarkExperiment()`

`officeExperiment()`

`simulationExperiment()`

`sensitivityStudy()`

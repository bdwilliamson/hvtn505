#!/usr/local/bin/Rscript

## MWE of beagle/gizmo problems

library("methods")
library("SuperLearner")
# library("future")
# library("future.apply")
library("e1071")
library("glmnet")
library("xgboost")
library("earth")
library("dplyr")
## only run this if necessary
# devtools::install_github("benkeser/cvma")
# library("cvma")
## only run this if something has changed
# install.packages("HVTN505_2019-4-1.tar.gz", type = "source", repos = NULL)
library("HVTN505")
library("kyotil")
library("argparse")

print("I made it!")
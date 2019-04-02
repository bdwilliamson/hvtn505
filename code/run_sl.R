## run the super learner
## make sure that it is CV.SL, averaged over 10 random starts

## load required libraries and functions
library("SuperLearner")
## only run this if necessary
# devtools::install_github("benkeser/cvma")
library("cvma")
## only run this if something has changed
# install.packages("HVTN505_2019-4-1.tar.gz", type = "source", repos = NULL)
library("HVTN505")

## set up code directory
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  code_dir <- "code/"
} else {
  
}
source(paste0(code_dir, "sl_screens.R")) # set up the screen/algorithm combinations

## ---------------------------------------------------------------------------------
## pre-process the data
## ---------------------------------------------------------------------------------
## read in the full dataset
data("dat.505", package = "HVTN505")
## read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists

library(kyotil)
for (a in var.super$varname) {
  dat.505[[a]]=scale(dat.505[[a]], center=mean(dat.505[[a]][dat.505$trt==1]), scale=sd(dat.505[[a]][dat.505$trt==1]))
  dat.505[[a%.%"_bin"]]=scale(dat.505[[a%.%"_bin"]], center=mean(dat.505[[a%.%"_bin"]][dat.505$trt==1]), scale=sd(dat.505[[a%.%"_bin"]][dat.505$trt==1]))
}
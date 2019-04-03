## run the super learner
## make sure that it is CV.SL, averaged over 10 random starts

## load required libraries and functions
library("SuperLearner")
library("future")
library("future.apply")
library("e1071")
library("glmnet")
library("xgboost")
library("dplyr")
## only run this if necessary
# devtools::install_github("benkeser/cvma")
# library("cvma")
## only run this if something has changed
# install.packages("HVTN505_2019-4-1.tar.gz", type = "source", repos = NULL)
library("HVTN505")
library("kyotil")

## set up code directory
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  code_dir <- "code/"
} else {
  
}
source(paste0(code_dir, "sl_screens.R")) # set up the screen/algorithm combinations
source(paset0(code_dir, "utils.R")) # get CV-AUC for all algs

## ---------------------------------------------------------------------------------
## pre-process the data
## ---------------------------------------------------------------------------------
## read in the full dataset
data("dat.505", package = "HVTN505")
## read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists

## scale vaccine recipients to have mean 0, sd 1 for all vars
for (a in var.super$varname) {
  dat.505[[a]]=scale(dat.505[[a]], center=mean(dat.505[[a]][dat.505$trt==1]), scale=sd(dat.505[[a]][dat.505$trt==1]))
  dat.505[[a%.%"_bin"]]=scale(dat.505[[a%.%"_bin"]], center=mean(dat.505[[a%.%"_bin"]][dat.505$trt==1]), scale=sd(dat.505[[a%.%"_bin"]][dat.505$trt==1]))
}

## set up X, Y for super learning
X_markers <- dat.505 %>% 
  select(var.super$varname, paste0(var.super$varname, "_bin"))
X_exposure <- dat.505 %>% 
  select(age, BMI, bhvrisk)
X <- data.frame(trt = dat.505$trt, X_exposure, X_markers)
weights <- dat.505$wt
Y <- dat.505$case
vaccinees <- cbind.data.frame(Y, weights, X) %>% 
  filter(trt == 1) %>% 
  select(-trt)
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weights
X_vaccine <- vaccinees %>% 
  select(-Y, -weights)

V <- 10
## ---------------------------------------------------------------------------------
## run super learner, with leave-one-out cross-validation and all screens
## do 10 random starts, average over these
## ---------------------------------------------------------------------------------
fits <- vector("list", length = 10)
set.seed(4747)
for (i in 1:10) {
  system.time(fit <- SuperLearner::CV.SuperLearner(Y = Y_vaccine, X = X_vaccine, family = binomial(),
                                                   obsWeights = weights_vaccine,
                                                   SL.library = SL_library, # this comes from sl_screens.R
                                                   method = "method.CC_nloglik",
                                                   cvControl = list(V = V),
                                                   innerCvControl = list(list(V = length(Y_vaccine) - 1)))
              )
  aucs <- get_all_aucs(fit)
  fits[[i]] <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
}

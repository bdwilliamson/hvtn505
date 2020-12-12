#!/usr/local/bin/Rscript

# run the super learner
# make sure that it is CV.SL, averaged over 10 random starts

# load required libraries and functions
library("methods")
library("SuperLearner")
# library("future")
# library("future.apply")
library("e1071")
library("glmnet")
library("xgboost")
library("earth")
library("dplyr")
# only run this if something has changed
# install.packages("HVTN505_2019-4-25.tar.gz", type = "source", repos = NULL)
library("HVTN505")
library("kyotil")
library("argparse")
# only run this if something has changed
# devtools::install_github("bdwilliamson/vimp", upgrade = "never")

# set up code directory
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  code_dir <- "code/"
  data_dir <- "data/"
  # load vimp from user library (to make sure it has correct version)
  library("vimp")
} else {
  code_dir <- ""
  data_dir <- ""
  # load vimp from user library (to make sure it has correct version)
  library("vimp", lib.loc = .libPaths()[2])
}
num_cores <- parallel::detectCores()
print(num_cores)
source(paste0(code_dir, "sl_screens.R")) # set up the screen/algorithm combinations
source(paste0(code_dir, "utils.R")) # get CV-AUC for all algs
source(paste0(code_dir, "measure_auc_ipw.R"))

parser <- ArgumentParser()
parser$add_argument("--aipw", default = "TRUE", 
                    help = "fit the AIPW (TRUE) or IPW (FALSE) estimator")
args <- parser$parse_args()
args$aipw <- as.logical(args$aipw)
if ( !args$aipw ) {
  measure_auc <- measure_auc_ipw
}
# ------------------------------------------------------------------------------
# pre-process the data
# ------------------------------------------------------------------------------
# read in the full dataset
data("dat.505", package = "HVTN505")
# read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists
# note that "var.super" contains individual vars for vaccine-matched antigens,
# and for vaccine-mismatched antigens, has either individual var (if only one)
# or PC1 and/or MDW (only PC1 if cor(PC1, MDW) > 0.9)

# scale vaccine recipients to have mean 0, sd 1 for all vars
for (a in var.super$varname) {
  dat.505[[a]] <- as.vector(
    scale(dat.505[[a]], 
          center = mean(dat.505[[a]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a]][dat.505$trt == 1]))
  )
  dat.505[[a%.%"_bin"]] <- as.vector(
    scale(dat.505[[a%.%"_bin"]], 
          center = mean(dat.505[[a%.%"_bin"]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a%.%"_bin"]][dat.505$trt == 1]))
  )
}
for (a in c("age", "BMI", "bhvrisk")) {
  dat.505[[a]] <- as.vector(
    scale(dat.505[[a]], 
          center = mean(dat.505[[a]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a]][dat.505$trt == 1]))
  )
}

# set up X, Y for super learning
X_markers <- dat.505 %>%
  select(var.super$varname, paste0(var.super$varname, "_bin"))

# only include the following variable sets:
assays <- unique(var.super$assay)
antigens <- unique(var.super$antigen)
# 1. None (baseline variables only)
var_set_none <- rep(FALSE, ncol(X_markers))
this_var_set <- var_set_none
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cat("\n Running noise variables; run ", job_id, "\n")
set.seed(1234)
noise_seeds <- round(runif(50, 1e3, 1e5))
set.seed(noise_seeds[job_id])
X_noise <- as_tibble(replicate(10, rnorm(nrow(X_markers), 0, 1)))

# get the current phase 2 dataset based on the markers of interest
X_markers_varset <- X_noise
X_exposure <- dat.505 %>%
  as_tibble() %>%
  select(age, BMI, bhvrisk)
X <- tibble::tibble(ptid = dat.505$ptid, trt = dat.505$trt,
                    weight = dat.505$wt) %>%
  bind_cols(X_exposure, X_markers_varset)
Y <- tibble(Y = dat.505$case)
vaccinees <- dplyr::bind_cols(Y, X) %>%
  filter(trt == 1) %>%
  select(-trt)
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weight
X_vaccine <- vaccinees %>%
  select(-Y, -weight, -ptid)
# read in the full phase 1 dataset and weights,
# and reorder so that rows match rows of X_vaccine with the remaining rows after
# note that in this case, Z_plus_weights has 2494 rows
# (the number of participants with complete data on age, BMI, bhvrisk)
# and that the number of vaccinees is 1250
Z_plus_weights <- readRDS(file = paste0(data_dir, 
                                        "z_and_weights_for_505_analysis.rds"))
# pull out the participants in the cc cohort who also received the vaccine;
# this matches the rows in vaccinees
all_cc_vaccine <- Z_plus_weights %>%
  filter(ptid %in% vaccinees$ptid, trt == 1)
# pull out the participants who are NOT in the cc cohort and received the vaccine
all_non_cc_vaccine <- Z_plus_weights %>%
  filter(!(ptid %in% vaccinees$ptid), trt == 1)
# put them back together
phase_1_data_vaccine <- dplyr::bind_rows(all_cc_vaccine, all_non_cc_vaccine) %>%
  select(-trt)
Z_vaccine <- phase_1_data_vaccine %>%
  select(-ptid, -weight)
all_ipw_weights_vaccine <- phase_1_data_vaccine %>%
  pull(weight)
C <- (phase_1_data_vaccine$ptid %in% vaccinees$ptid)

V_outer <- 5
V_inner <- length(Y_vaccine) - 1

# get the SL library
# if var_set_none, then don't need screens; otherwise do
sl_lib <- SL_library_noise[!grepl("earth", SL_library_noise)] 
# ------------------------------------------------------------------------------
# run super learner, with leave-one-out cross-validation and all screens
# do 10 random starts, average over these
# ------------------------------------------------------------------------------
# ensure reproducibility
set.seed(4747)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts
fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y_vaccine, 
                           X_mat = X_vaccine, family = "binomial",
                           Z = Z_vaccine, C = C, 
                           z_lib = methods[!grepl("earth", methods)],
                           obsWeights = weights_vaccine,
                           all_weights = all_ipw_weights_vaccine,
                           scale = "logit",
                           sl_lib = sl_lib, # this comes from sl_screens.R
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           vimp = FALSE,
                           mc.cores = num_cores)
saveRDS(fits, paste0("sl_fits_noise_", switch((args$aipw) + 1, "ipw", "aipw"),
                     "_", job_id, ".rds"))
warnings()

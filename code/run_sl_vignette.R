#!/usr/local/bin/Rscript

# run the analysis with a small library for the vignette

# load required libraries and functions
library("methods")
library("SuperLearner")
library("e1071")
library("glmnet")
library("xgboost")
library("earth")
library("dplyr")
# only run this if something has changed
# devtools::install_local("HVTN505_2019-4-25.tar.gz")
library("HVTN505")
# only run this if something has changed
# devtools::install_github("bdwilliamson/vimp", upgrade = "never")
library("vimp")
library("kyotil")
library("argparse")
library("here")

num_cores <- parallel::detectCores()
# print(num_cores)
# read in screens/algos and utility functions (shown above)
source(here("code", "sl_screens.R"))
source(here("code", "utils.R"))

# ---------------------------------------------------------------------------------
# pre-process the data
# ---------------------------------------------------------------------------------
# read in the full dataset
data("dat.505", package = "HVTN505")
# read in the super learner variables
# even if there is a warning message, it still exists
suppressWarnings(data("var.super", package = "HVTN505"))
# note that "var.super" contains individual vars for vaccine-matched antigens,
# and for vaccine-mismatched antigens, has either individual var (if only one)
# or PC1 and/or MDW (only PC1 if cor(PC1, MDW) > 0.9)

# scale vaccine recipients to have mean 0, sd 1 for all vars
# scale vaccine recipients to have mean 0, sd 1 for all vars
for (a in var.super$varname) {
  dat.505[[a]] <- as.vector(
    scale(dat.505[[a]], 
          center = mean(dat.505[[a]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a]][dat.505$trt == 1])))
  dat.505[[a%.%"_bin"]] <- as.vector(
    scale(dat.505[[a%.%"_bin"]], 
          center = mean(dat.505[[a%.%"_bin"]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a%.%"_bin"]][dat.505$trt == 1])))
}
for (a in c("age", "BMI", "bhvrisk")) {
  dat.505[[a]] <- as.vector(
    scale(dat.505[[a]], 
          center = mean(dat.505[[a]][dat.505$trt == 1]), 
          scale = sd(dat.505[[a]][dat.505$trt == 1])))
}

# set up X, Y for super learning
X_markers <- dat.505 %>%
  select(var.super$varname, paste0(var.super$varname, "_bin"))
X_exposure <- dat.505 %>%
  as_tibble() %>% 
  select(age, BMI, bhvrisk)
X <- tibble::tibble(ptid = dat.505$ptid, trt = dat.505$trt,
                    weight = dat.505$wt) %>%
  bind_cols(X_exposure, X_markers)
X_none_all <- tibble::tibble(ptid = dat.505$ptid, trt = dat.505$trt,
                             weight = dat.505$wt) %>%
  bind_cols(X_exposure)
Y <- tibble(Y = dat.505$case)
vaccinees <- dplyr::bind_cols(Y, X) %>%
  filter(trt == 1) %>%
  select(-trt)
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weight
X_vaccine <- vaccinees %>%
  select(-Y, -weight, -ptid)
X_none <- X_none_all %>% 
  filter(trt == 1) %>% 
  select(-trt, -ptid, -weight)
Z_plus_weights <- readRDS(file = here("data", "z_and_weights_for_505_analysis.rds"))
# match the rows in vaccinees to get Z, C
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

# ---------------------------------------------------------------------------------
# run super learner, with leave-one-out cross-validation 
# do 5 random starts, average over these (note that the full analysis uses 10 random starts)
# note that the library is reduced to speed up computation time in this vignette --
# in practice, we recommend using a large library
# ---------------------------------------------------------------------------------
# ensure reproducibility
set.seed(4747)
seeds <- round(runif(5, 1000, 10000)) # average over 5 random starts
fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, 
                           Y = Y_vaccine,
                           X_mat = X_vaccine, 
                           family = "binomial",
                           Z = Z_vaccine, 
                           C = C, z_lib = "SL.glm",
                           obsWeights = weights_vaccine,
                           all_weights = all_ipw_weights_vaccine,
                           scale = "logit",
                           sl_lib = list(SL_library[[1]][c(1, 2)]),
                           method = "method.CC_nloglik",
                           cvControl = 
                             list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = 
                             list(list(V = V_inner)),
                           vimp = FALSE,
                           mc.cores = num_cores
)
saveRDS(fits, here("results", "vignette_sl_fits_varset_11_all.rds"))

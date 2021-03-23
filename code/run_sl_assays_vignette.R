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
Y_vaccine <- tibble(Y = subset(dat.505$case, dat.505$trt == 1))$Y
Z_plus_weights <- readRDS(file = here("data", "z_and_weights_for_505_analysis.rds"))

# only include the following variable sets:
assays <- unique(var.super$assay)
antigens <- unique(var.super$antigen)
# 1. None (baseline variables only)
var_set_none <- rep(FALSE, ncol(X_markers))
# 2. IgG + IgA (all antigens)
var_set_igg_iga <- get_nms_group_all_antigens(X_markers, 
                                              assays = c("IgG", "IgA"),
                                              assays_to_exclude = "IgG3")
# 3. IgG3
var_set_igg3 <- get_nms_group_all_antigens(X_markers, assays = "IgG3")
# 4. T cells (all antigens)
var_set_tcells <- get_nms_group_all_antigens(X_markers, 
                                             assays = c("CD4", "CD8"))
# 5. Fx Ab (all antigens)
var_set_fxab <- get_nms_group_all_antigens(X_markers, 
                                           assays = c("phago", "R2a", "R3a"))
# 6. 1+2+3
var_set_igg_iga_igg3 <- get_nms_group_all_antigens(X_markers, 
                                                   assays = 
                                                     c("IgG", "IgA", "IgG3"))
# 7. 1+2+4
var_set_igg_iga_tcells <- get_nms_group_all_antigens(X_markers,
                                                     assays = 
                                                       c("IgG", "IgA", 
                                                         "CD4", "CD8"),
                                                     assays_to_exclude = 
                                                       "IgG3")
# 8. 1+2+3+4
var_set_igg_iga_igg3_tcells <- get_nms_group_all_antigens(X_markers,
                                                          assays = c("IgG", 
                                                                     "IgA", 
                                                                     "IgG3", 
                                                                     "CD4", 
                                                                     "CD8"))
# 9. 1+2+3+5
var_set_igg_iga_igg3_fxab <- get_nms_group_all_antigens(X_markers,
                                                        assays = c("IgG", 
                                                                   "IgA", 
                                                                   "IgG3", 
                                                                   "phago", 
                                                                   "R2a", 
                                                                   "R3a"))
# 10. 1+4+5
var_set_tcells_fxab <- get_nms_group_all_antigens(X_markers,
                                                  assays = c("CD4", "CD8", 
                                                             "phago", "R2a", 
                                                             "R3a"))
# 11. All
var_set_all <- rep(TRUE, ncol(X_markers))
# 12--14: extra runs to get variable importance
var_set_igg3_fxab <- get_nms_group_all_antigens(X_markers,
                                                assays = c("IgG3", "phago", 
                                                           "R2a", "R3a"))
var_set_igg_iga_tcells_fxab <- get_nms_group_all_antigens(X_markers,
                                                          assays = c("IgG", 
                                                                     "IgA", 
                                                                     "CD4", 
                                                                     "CD8", 
                                                                     "phago", 
                                                                     "R2a", 
                                                                     "R3a"),
                                                          assays_to_exclude = 
                                                            "IgG3")
var_set_igg3_tcells_fxab <- get_nms_group_all_antigens(X_markers,
                                                       assays = c("IgG3", 
                                                                  "CD4", "CD8", 
                                                                  "phago", 
                                                                  "R2a", 
                                                                  "R3a"))

var_set_names <- c("1_baseline_exposure", "2_igg_iga", "3_igg3", 
                   "4_tcells", "5_fxab",
                   "6_igg_iga_igg3", "7_igg_iga_tcells", 
                   "8_igg_iga_igg3_tcells",
                   "9_igg_iga_igg3_fxab", "10_tcells_fxab",
                   "11_all",
                   "12_igg3_fxab", "13_igg_iga_tcells_fxab", 
                   "14_igg3_tcells_fxab")

# set up a matrix of all
var_set_matrix <- rbind(var_set_none, var_set_igg_iga, 
                        var_set_igg3,
                        var_set_tcells, var_set_fxab, 
                        var_set_igg_iga_igg3,
                        var_set_igg_iga_tcells, 
                        var_set_igg_iga_igg3_tcells,
                        var_set_igg_iga_igg3_fxab, 
                        var_set_tcells_fxab,
                        var_set_all, var_set_igg3_fxab,
                        var_set_igg_iga_tcells_fxab, 
                        var_set_igg3_tcells_fxab)
for (i in (1:14)[-11]) {
  job_id <- i
  this_var_set <- var_set_matrix[job_id, ]
  cat("\n Running ", var_set_names[job_id], "\n")
  
  X_markers_varset <- X_markers %>%
    select(names(X_markers)[this_var_set])
  X_exposure <- dat.505 %>%
    as_tibble() %>% 
    select(age, BMI, bhvrisk)
  X <- tibble::tibble(ptid = dat.505$ptid, trt = dat.505$trt, 
                      weight = dat.505$wt) %>% 
    bind_cols(X_exposure, X_markers_varset)
  weights <- dat.505$wt
  Y <- tibble(Y = dat.505$case)
  vaccinees <- dplyr::bind_cols(Y, X) %>%
    filter(trt == 1) %>%
    select(-trt)
  Y_vaccine <- vaccinees$Y
  weights_vaccine <- vaccinees$weight
  X_vaccine <- vaccinees %>%
    select(-Y, -weight, -ptid)
  # match the rows in vaccinees to get Z, C
  all_cc_vaccine <- Z_plus_weights %>%
    filter(ptid %in% vaccinees$ptid, trt == 1)
  # pull out the participants who are NOT in the cc cohort 
  # and received the vaccine
  all_non_cc_vaccine <- Z_plus_weights %>%
    filter(!(ptid %in% vaccinees$ptid), trt == 1)
  # put them back together
  phase_1_data_vaccine <- dplyr::bind_rows(all_cc_vaccine, 
                                           all_non_cc_vaccine) %>%
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
  if (job_id == 1) {
    sl_lib <- methods[1]
  } else {
    sl_lib <- list(SL_library[[1]][c(1, 2)])
  }
  # ---------------------------------------------------------------------------------
  # run super learner, with leave-one-out cross-validation and all screens
  # do 5 random starts, average over these
  # ---------------------------------------------------------------------------------
  # ensure reproducibility
  set.seed(4747)
  seeds <- round(runif(5, 1000, 10000)) # average over five random starts
  fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y_vaccine,
                             X_mat = X_vaccine,
                             family = "binomial",
                             C = C, Z = Z_vaccine, z_lib = "SL.glm",
                             obsWeights = weights_vaccine,
                             all_weights = all_ipw_weights_vaccine,
                             scale = "logit",
                             sl_lib = sl_lib,
                             method = "method.CC_nloglik",
                             cvControl = list(V = V_outer, stratifyCV = TRUE),
                             innerCvControl = list(list(V = V_inner)),
                             vimp = FALSE,
                             mc.cores = num_cores
  )
  # save!
  eval(parse(text = paste0("saveRDS(fits, here('results', 'vignette_sl_fits_varset_", 
                            var_set_names[job_id], ".rds'))")))
}

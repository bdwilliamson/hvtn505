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
measure_auc <- measure_auc_ipw

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
# 2. IgG + IgA (all antigens)
var_set_igg_iga <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG", "IgA"), assays_to_exclude = "IgG3"
)
# 3. IgG3
var_set_igg3 <- get_nms_group_all_antigens(
  X_markers, assays = "IgG3"
)
# 4. T cells (all antigens)
var_set_tcells <- get_nms_group_all_antigens(
  X_markers, assays = c("CD4", "CD8")
)
# 5. Fx Ab (all antigens)
var_set_fxab <- get_nms_group_all_antigens(
  X_markers, assays = c("phago", "R2a", "R3a")
)
# 6. 1+2+3
var_set_igg_iga_igg3 <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG", "IgA", "IgG3")
)
# 7. 1+2+4
var_set_igg_iga_tcells <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG", "IgA", "CD4", "CD8"), assays_to_exclude = "IgG3"
)
# 8. 1+2+3+4
var_set_igg_iga_igg3_tcells <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG", "IgA", "IgG3", "CD4", "CD8")
)
# 9. 1+2+3+5
var_set_igg_iga_igg3_fxab <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG", "IgA", "IgG3", "phago", "R2a", "R3a")
)
# 10. 1+4+5
var_set_tcells_fxab <- get_nms_group_all_antigens(
  X_markers, assays = c("CD4", "CD8", "phago", "R2a", "R3a")
)
# 11. All
var_set_all <- rep(TRUE, ncol(X_markers))
# 12--14: extra runs to get variable importance
var_set_igg3_fxab <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG3", "phago", "R2a", "R3a")
)
var_set_igg_iga_tcells_fxab <- get_nms_group_all_antigens(
  X_markers, 
  assays = c("IgG", "IgA", "CD4", "CD8", "phago", "R2a", "R3a"), 
  assays_to_exclude = "IgG3"
)
var_set_igg3_tcells_fxab <- get_nms_group_all_antigens(
  X_markers, assays = c("IgG3", "CD4", "CD8", "phago", "R2a", "R3a")
)

var_set_names <- c("1_baseline_exposure", "2_igg_iga", 
                   "3_igg3","4_tcells", "5_fxab",
                   "6_igg_iga_igg3", "7_igg_iga_tcells", 
                   "8_igg_iga_igg3_tcells",
                   "9_igg_iga_igg3_fxab", "10_tcells_fxab",
                   "11_all",
                   "12_igg3_fxab", "13_igg_iga_tcells_fxab", 
                   "14_igg3_tcells_fxab")

# set up a matrix of all
var_set_matrix <- rbind(var_set_none, var_set_igg_iga, var_set_igg3, 
                        var_set_tcells, var_set_fxab,
                        var_set_igg_iga_igg3, var_set_igg_iga_tcells, 
                        var_set_igg_iga_igg3_tcells,
                        var_set_igg_iga_igg3_fxab, var_set_tcells_fxab,
                        var_set_all,
                        var_set_igg3_fxab, var_set_igg_iga_tcells_fxab, 
                        var_set_igg3_tcells_fxab)
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_var_set <- var_set_matrix[job_id, ]
cat("\n Running ", var_set_names[job_id], "\n")

# get the current phase 2 dataset based on the markers of interest
X_markers_varset <- X_markers %>%
  select(names(X_markers)[this_var_set])
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
if (job_id == 1) {
  sl_lib <- methods[!grepl("earth", methods)]
} else {
  sl_lib <- SL_library[!grepl("earth", SL_library)] 
}
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
saveRDS(fits, paste0("sl_fits_varset_", var_set_names[job_id], "_ipw.rds"))
warnings()

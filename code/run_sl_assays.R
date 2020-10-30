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
# load vimp from user library (to make sure it has correct version)
library("vimp", lib.loc = .libPaths()[2])

# set up code directory
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  code_dir <- "code/"
} else {
  code_dir <- ""
}
num_cores <- parallel::detectCores()
print(num_cores)
source(paste0(code_dir, "sl_screens.R")) # set up the screen/algorithm combinations
source(paste0(code_dir, "utils.R")) # get CV-AUC for all algs

# ---------------------------------------------------------------------------------
# pre-process the data
# ---------------------------------------------------------------------------------
# read in the full dataset
data("dat.505", package = "HVTN505")
# read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists
# note that "var.super" contains individual vars for vaccine-matched antigens,
# and for vaccine-mismatched antigens, has either individual var (if only one)
# or PC1 and/or MDW (only PC1 if cor(PC1, MDW) > 0.9)

# scale vaccine recipients to have mean 0, sd 1 for all vars
for (a in var.super$varname) {
  dat.505[[a]] <- scale(dat.505[[a]], center = mean(dat.505[[a]][dat.505$trt == 1]), scale = sd(dat.505[[a]][dat.505$trt == 1]))
  dat.505[[a%.%"_bin"]] <- scale(dat.505[[a%.%"_bin"]], center = mean(dat.505[[a%.%"_bin"]][dat.505$trt == 1]), scale = sd(dat.505[[a%.%"_bin"]][dat.505$trt == 1]))
}
for (a in c("age", "BMI", "bhvrisk")) {
  dat.505[[a]] <- scale(dat.505[[a]], center = mean(dat.505[[a]][dat.505$trt == 1]), scale = sd(dat.505[[a]][dat.505$trt == 1]))
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
var_set_igg_iga <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA"), assays_to_exclude = "IgG3")
# 3. IgG3
var_set_igg3 <- get_nms_group_all_antigens(X_markers, assays = "IgG3")
# 4. T cells (all antigens)
var_set_tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
# 5. Fx Ab (all antigens)
var_set_fxab <- get_nms_group_all_antigens(X_markers, assays = c("phago", "R2a", "R3a"))
# 6. 1+2+3
var_set_igg_iga_igg3 <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3"))
# 7. 1+2+4
var_set_igg_iga_tcells <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "CD4", "CD8"), assays_to_exclude = "IgG3")
# 8. 1+2+3+4
var_set_igg_iga_igg3_tcells <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3", "CD4", "CD8"))
# 9. 1+2+3+5
var_set_igg_iga_igg3_fxab <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "IgG3", "phago", "R2a", "R3a"))
# 10. 1+4+5
var_set_tcells_fxab <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8", "phago", "R2a", "R3a"))
# 11. All
var_set_all <- rep(TRUE, ncol(X_markers))
# 12--14: extra runs to get variable importance
var_set_igg3_fxab <- get_nms_group_all_antigens(X_markers, assays = c("IgG3", "phago", "R2a", "R3a"))
var_set_igg_iga_tcells_fxab <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA", "CD4", "CD8", "phago", "R2a", "R3a"), assays_to_exclude = "IgG3")
var_set_igg3_tcells_fxab <- get_nms_group_all_antigens(X_markers, assays = c("IgG3", "CD4", "CD8", "phago", "R2a", "R3a"))

var_set_names <- c("1_baseline_exposure", "2_igg_iga", "3_igg3","4_tcells", "5_fxab",
                   "6_igg_iga_igg3", "7_igg_iga_tcells", "8_igg_iga_igg3_tcells",
                   "9_igg_iga_igg3_fxab", "10_tcells_fxab",
                   "11_all",
                   "12_igg3_fxab", "13_igg_iga_tcells_fxab", "14_igg3_tcells_fxab")

# set up a matrix of all
var_set_matrix <- rbind(var_set_none, var_set_igg_iga, var_set_igg3, var_set_tcells, var_set_fxab,
                        var_set_igg_iga_igg3, var_set_igg_iga_tcells, var_set_igg_iga_igg3_tcells,
                        var_set_igg_iga_igg3_fxab, var_set_tcells_fxab,
                        var_set_all,
                        var_set_igg3_fxab, var_set_igg_iga_tcells_fxab, var_set_igg3_tcells_fxab)
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_var_set <- var_set_matrix[job_id, ]
cat("\n Running ", var_set_names[job_id], "\n")

X_markers_varset <- X_markers %>%
  select(names(X_markers)[this_var_set])

X_exposure <- dat.505 %>%
  select(age, BMI, bhvrisk)
X <- data.frame(trt = dat.505$trt, X_exposure, X_markers_varset)
weights <- dat.505$wt
Y <- dat.505$case
vaccinees <- cbind.data.frame(Y, weights, X) %>%
  filter(trt == 1) %>%
  select(-trt)
Y_vaccine <- vaccinees$Y
weights_vaccine <- vaccinees$weights
X_vaccine <- vaccinees %>%
  select(-Y, -weights)
C <- rep(1, length(Y_vaccine))
Z <- data.frame(Y = Y_vaccine, X_vaccine %>% select(age, BMI, bhvrisk))

V_outer <- 5
V_inner <- length(Y_vaccine) - 1

# get the SL library
# if var_set_none, then don't need screens; otherwise do
if (job_id == 1) {
  sl_lib <- methods[!grepl("earth", methods)]
} else {
  sl_lib <- SL_library[!grepl("earth", SL_library)] # get rid of numerical errors from SL.earth
}
# ---------------------------------------------------------------------------------
# run super learner, with leave-one-out cross-validation and all screens
# do 10 random starts, average over these
# ---------------------------------------------------------------------------------
# ensure reproducibility
set.seed(4747)
seeds <- round(runif(10, 1000, 10000)) # average over 10 random starts
fits <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y_vaccine, X_mat = X_vaccine, family = "binomial",
                           Z = Z, C = C, z_lib = methods[!grepl("earth", methods)],
                           obsWeights = weights_vaccine,
                           scale = "logit",
                           sl_lib = sl_lib, # this comes from sl_screens.R
                           method = "method.CC_nloglik",
                           cvControl = list(V = V_outer, stratifyCV = TRUE),
                           innerCvControl = list(list(V = V_inner)),
                           vimp = FALSE,
                           mc.cores = num_cores)
saveRDS(fits, paste0("sl_fits_varset_", var_set_names[job_id], ".rds"))
warnings()

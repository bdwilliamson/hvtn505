## load in results from individual vimp analysis, make a nice plot

## --------------------------------------------------------------------------------------------------------------------------------------
## set up directories, load required packages
## --------------------------------------------------------------------------------------------------------------------------------------
library("SuperLearner")
library("cvAUC")
library("tidyr")
library("dplyr")
library("cowplot")
# install if anything has changed
#  install.packages("~/Projects/UW/vimp_1.3.0.tar.gz", repos = NULL, type = "source")
library("vimp")
library("HVTN505")
library("kyotil")
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
source("code/auc_plot_assays.R")
source("code/r2_plot_assays.R")
source("code/utils.R")

results_dir <- "results/"
plots_dir <- "plots/"

## --------------------------------------------------------------------------------------------------------------------------------------
## load results objects:
## --------------------------------------------------------------------------------------------------------------------------------------
## baseline vars only
sl_fits_varset_1_baseline_exposure <- readRDS(paste0(results_dir, "sl_fits_varset_1_baseline_exposure.rds"))

## ---------------------------------------------------------------------------------
## pre-process the data
## ---------------------------------------------------------------------------------
## read in the full dataset
data("dat.505", package = "HVTN505")
## read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists
## note that "var.super" contains individual vars for vaccine-matched antigens,
## and for vaccine-mismatched antigens, has either individual var (if only one)
## or PC1 and/or MDW (only PC1 if cor(PC1, MDW) > 0.9)

## scale vaccine recipients to have mean 0, sd 1 for all vars
for (a in var.super$varname) {
  dat.505[[a]] <- scale(dat.505[[a]], center = mean(dat.505[[a]][dat.505$trt == 1]), scale = sd(dat.505[[a]][dat.505$trt == 1]))
  dat.505[[a%.%"_bin"]] <- scale(dat.505[[a%.%"_bin"]], center = mean(dat.505[[a%.%"_bin"]][dat.505$trt == 1]), scale = sd(dat.505[[a%.%"_bin"]][dat.505$trt == 1]))
}
for (a in c("age", "BMI", "bhvrisk")) {
  dat.505[[a]] <- scale(dat.505[[a]], center = mean(dat.505[[a]][dat.505$trt == 1]), scale = sd(dat.505[[a]][dat.505$trt == 1]))
}

## set up X, Y for super learning
X_markers <- dat.505 %>% 
  select(var.super$varname, paste0(var.super$varname, "_bin")) 

## load all individual vars
var_names <- names(X_markers)
for (i in 1:length(var_names)) {
  eval(parse(text = paste0("sl_fit_var_", i, " <- readRDS(paste0(results_dir, 'sl_fits_vimp_', i + 190, '.rds'))")))
}

## --------------------------------------------------------------------------------------------------------------------------------------
## Variable importance analysis
## Presentation: most to least important, within groups defined by full SL importance?
## So Fx ab group on top, IgG + IgA next, IgG3 next, T cells next
## --------------------------------------------------------------------------------------------------------------------------------------
risk_type <- "auc"
for (i in 1:length(var_names)) {
  ## only run if results computation changes
  # eval(parse(text = paste0("vimp_", i, "<- get_cv_vim(full_fit = sl_fit_var_", i, 
  #                          ", reduced_fit = sl_fits_varset_1_baseline_exposure, type = risk_type,
  #                          vimp = TRUE)")))
  # eval(parse(text = paste0("vimp_", i, "_avg <- get_avg_est_ci(vimp_", i, ")")))
  # eval(parse(text = paste0("saveRDS(vimp_", i, "_avg, paste0(results_dir, 'vimp_", i, "_avg.rds'))")))
  eval(parse(text = paste0("vimp_", i, "<- readRDS(paste0(results_dir, 'vimp_", i, "_avg.rds'))")))
}

## make groups of marker variables
fx_ab <- get_nms_group_all_antigens(X_markers, assays =  c("phago", "R2a", "R3a"))
igg_iga <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA"), assays_to_exclude = "IgG3")
igg3 <- get_nms_group_all_antigens(X_markers, assays = "IgG3")
tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))

## make forest plot, with Fx Ab on top, then IgG + IgA, then IgG3, then T cells
vimp_tibble_ind <- tibble(var_name = var_names[1], assay_group = "T cells", 
                          est = vimp_1_avg$est, cil = vimp_1_avg$ci[1],
                          ciu = vimp_1_avg$ci[2])
assay_nms <- c("Fx Ab", "IgG + IgA", "IgG3", "T cells")
for (i in 2:length(var_names)) {
  this_est <- eval(parse(text = paste0("vimp_", i, "_avg")))
  vimp_tibble_ind <- vimp_tibble_ind %>% 
    tibble::add_row(var_name = var_names[i], assay_group = 
                    est = this_est$est, cil = this_est$ci[1],
                    ciu = this_est$ci[2])
}

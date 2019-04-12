## results from full SL analysis
## includes CV-AUC plots

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
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
source("code/auc_plot_assays.R")
source("code/utils.R")

results_dir <- "results/"
plots_dir <- "plots/"

## --------------------------------------------------------------------------------------------------------------------------------------
## load results objects:
## --------------------------------------------------------------------------------------------------------------------------------------
## each is a length 10 list (one for each random start)
var_set_names <- c("1_baseline_exposure", "2_igg_iga", "3_tcells", "4_fxab",
                   "5_igg_iga_tcells", "6_igg_iga_fxab", "7_tcells_fxab",
                   "8_all")
for (i in 1:length(var_set_names)) {
  eval(parse(text = paste0("sl_fits_varset_", var_set_names[i], " <- readRDS(paste0(results_dir, 'sl_fits_varset_', var_set_names[i], '.rds'))")))
}

## check the discrete SL for each fold (baseline exposure only)
lapply(sl_fits_varset_1_baseline_exposure, function(x) x$fit$whichDiscreteSL)

## check the weights
lapply(sl_fits_varset_1_baseline_exposure, function(x) sort(colMeans(x$fit$coef), decreasing=TRUE))

## average the AUCs over the 10 folds, for each
var_set_labels <- c("No markers", "IgG + IgA", "T Cells", "Fx Ab", "IgG + IgA + T Cells",
                    "IgG + IgA + Fx Ab", "T Cells + Fx Ab", "All markers")
for (i in 1:length(var_set_names)) {
  this_name <- paste(unlist(strsplit(var_set_names[i], "_", fixed = TRUE))[-1], collapse = "_")
  eval(parse(text = paste0("all_aucs_i <- as_tibble(do.call(rbind.data.frame, lapply(sl_fits_varset_", var_set_names[i], ", function(x) x$aucs)))")))
  eval(parse(text = paste0("avg_aucs_", var_set_names[i]," <- all_aucs_i %>% 
    group_by(Learner, Screen) %>% 
    summarize(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul)) %>% 
    mutate(assay = this_name, varset_label = var_set_labels[i]) %>% 
    ungroup()")))
}

## combine into a full tibble; add a column to each that is the assay
avg_aucs <- bind_rows(avg_aucs_1_baseline_exposure, avg_aucs_2_igg_iga, avg_aucs_3_tcells,
                     avg_aucs_4_fxab, avg_aucs_5_igg_iga_tcells, avg_aucs_6_igg_iga_fxab,
                     avg_aucs_7_tcells_fxab, avg_aucs_8_all)

## forest plot of AUCs; this one is super nasty, but shows that it works
full_forest_plot_auc <- avg_aucs %>% 
  ggplot(aes(x = AUC, y = factor(paste0(Screen, "_", Learner, "_", assay), levels = paste0(Screen, "_", Learner, "_", assay)[order(AUC)]))) + 
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
  geom_point() +
  xlab("CV-AUC") +
  ylab("Learner")


## --------------------------------------------------------------------------------------------------------------------------------------
## forest plot of CV-AUC for the top learner and SL for each assay combination
## --------------------------------------------------------------------------------------------------------------------------------------
title_font_size <- 18
main_font_size <- 5
fig_width <- fig_height <- 2590
y_title <- 0.96
auc_forest_plot <- auc_plot_assays(avg_aucs, main_font_size)
png(paste0(plots_dir, "cv_auc_forest_plot_sl_plus_top_learner.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
plot_grid(auc_forest_plot$top_learner_nms_plot, auc_forest_plot$top_learner_plot, nrow = 1, align = "h") +
  draw_label("Assay combination", size = title_font_size, x = 0.075, y = y_title) +
  draw_label("Algorithm", size = title_font_size, x = 0.175, y = y_title) +
  draw_label("Screen", size = title_font_size, x = 0.25, y = y_title) +
  draw_label("CV-AUC [95% CI]", size = title_font_size, x = 0.43, y = y_title)
dev.off()

## --------------------------------------------------------------------------------------------------------------------------------------
## Variable importance plot for the different assay combinations:
## All markers yields the full regression fit
## Variable set 1 gives reduced for importance of all markers (group 8)
## Variable set 2 gives reduced for importance of Fx Ab + Tcells (group 7)
## Variable set 3 gives reduced for importance of IgG + IgA + Fx Ab (group 6)
## Variable set 4 gives reduced for importance of IgG + IgA + Tcells (group 5)
## Variable set 5 gives reduced for importance of Fx Ab (group 4)
## Variable set 6 gives reduced for importance of T cells (group 3)
## Variable set 7 gives reduced for importance of IgG + IgA (group 2)
## --------------------------------------------------------------------------------------------------------------------------------------
## load the data
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
  dat.505[[a%.%"_bin"]] <- scale(dat.505[[a%.%"_bin"]], center = mean(dat.505[[a%.%"_bin"]][dat.505$trt == 1]), 
                                 scale = sd(dat.505[[a%.%"_bin"]][dat.505$trt == 1]))
}
for (a in c("age", "BMI", "bhvrisk")) {
  dat.505[[a]] <- scale(dat.505[[a]], center = mean(dat.505[[a]][dat.505$trt == 1]), scale = sd(dat.505[[a]][dat.505$trt == 1]))
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

## do variable importance
## all markers
risk_type <- "r_squared"
# risk_type <- "auc"
vimp_all_markers <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_1_baseline_exposure, type = risk_type)
unlist(lapply(vimp_all_markers, function(x) x$est))
var(unlist(lapply(vimp_all_markers, function(x) x$est)))
mean(unlist(lapply(vimp_all_markers, function(x) x$est)))
median(unlist(lapply(vimp_all_markers, function(x) x$est)))
vimp_all_markers_avg <- get_avg_est_ci(vimp_all_markers)
## T cells + Fx Ab
vimp_tcells_fxab <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_2_igg_iga, type = risk_type)
mean(unlist(lapply(vimp_tcells_fxab, function(x) x$est)))
vimp_tcells_fxab_avg <- get_avg_est_ci(vimp_tcells_fxab)
## IgG + IgA + Fx Ab
vimp_igg_iga_fxab <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_3_tcells, type = risk_type)
mean(unlist(lapply(vimp_igg_iga_fxab, function(x) x$est)))
vimp_igg_iga_fxab_avg <- get_avg_est_ci(vimp_igg_iga_fxab)
## IgG + IgA + T cells
vimp_igg_iga_tcells <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_4_fxab, type = risk_type)
mean(unlist(lapply(vimp_igg_iga_tcells, function(x) x$est)))
vimp_igg_iga_tcells_avg <- get_avg_est_ci(vimp_igg_iga_tcells)
## Fx Ab
vimp_fxab <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_5_igg_iga_tcells, type = risk_type)
mean(unlist(lapply(vimp_fxab, function(x) x$est)))
vimp_fxab_avg <- get_avg_est_ci(vimp_fxab)
## T cells
vimp_tcells <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_6_igg_iga_fxab, type = risk_type)
mean(unlist(lapply(vimp_tcells, function(x) x$est)))
vimp_tcells_avg <- get_avg_est_ci(vimp_tcells)
## IgG + IgA
vimp_igg_iga <- get_cv_vim(full_fit = sl_fits_varset_8_all, reduced_fit = sl_fits_varset_7_tcells_fxab, type = risk_type)
mean(unlist(lapply(vimp_igg_iga, function(x) x$est)))
vimp_igg_iga_avg <- get_avg_est_ci(vimp_igg_iga)

## combine together
vimp_tibble <- tibble(assay_grp = c("All markers", "T Cells + Fx Ab", "IgG + IgA + Fx Ab", 
                                    "IgG + IgA + T Cells", "Fx Ab", "T Cells", "IgG + IgA"),
                      est = c(vimp_all_markers_avg$est, vimp_tcells_fxab_avg$est, vimp_igg_iga_fxab_avg$est,
                              vimp_igg_iga_tcells_avg$est, vimp_fxab_avg$est, vimp_tcells_avg$est,
                              vimp_igg_iga_avg$est),
                      cil = c(vimp_all_markers_avg$ci[1], vimp_tcells_fxab_avg$ci[1], vimp_igg_iga_fxab_avg$ci[1],
                              vimp_igg_iga_tcells_avg$ci[1], vimp_fxab_avg$ci[1], vimp_tcells_avg$ci[1],
                              vimp_igg_iga_avg$ci[1]),
                      ciu = c(vimp_all_markers_avg$ci[2], vimp_tcells_fxab_avg$ci[2], vimp_igg_iga_fxab_avg$ci[2],
                              vimp_igg_iga_tcells_avg$ci[2], vimp_fxab_avg$ci[2], vimp_tcells_avg$ci[2],
                              vimp_igg_iga_avg$ci[2]))

## forest plot of vimp, with labels for the groups
vimp_forest_plot <- vimp_tibble %>% 
  ggplot(aes(x = est, y = factor(assay_grp, levels = assay_grp[order(est)], labels = assay_grp[order(est)]))) +
  geom_point() +
  geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
  ylab("Assay group") +
  xlab(paste0("Variable importance estimate: difference in ", ifelse(risk_type == "r_squared", expression(R^2), "AUC")))
png(paste0(plots_dir, "vimp_forest_plot_", risk_type, ".png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
vimp_forest_plot
dev.off()
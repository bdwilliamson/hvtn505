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
library("HVTN505")
library("kyotil")
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
source("code/plot_assays.R")
source("code/utils.R")
source("code/sl_screens.R")

results_dir <- "results/"
plots_dir <- "plots/"

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## load in the data to get weights
## read in the full dataset
data("dat.505", package = "HVTN505")
## read in the super learner variables
data("var.super", package = "HVTN505") # even if there is a warning message, it still exists

weights_vaccine <- dat.505$wt[dat.505$trt == 1]
## --------------------------------------------------------------------------------------------------------------------------------------
## load results objects:
## --------------------------------------------------------------------------------------------------------------------------------------
## each is a length 10 list (one for each random start)
var_set_names <- c("1_baseline_exposure", "2_igg_iga", "3_igg3","4_tcells", "5_fxab",
                   "6_igg_iga_igg3", "7_igg_iga_tcells", "8_igg_iga_igg3_tcells", 
                   "9_igg_iga_igg3_fxab", "10_tcells_fxab",
                   "11_all")
for (i in 1:length(var_set_names)) {
  eval(parse(text = paste0("sl_fits_varset_", var_set_names[i], " <- readRDS(paste0(results_dir, 'sl_fits_varset_', var_set_names[i], '.rds'))")))
}

## check the discrete SL for each fold (baseline exposure only)
lapply(sl_fits_varset_1_baseline_exposure, function(x) x$fit$whichDiscreteSL)
lapply(sl_fits_varset_11_all, function(x) x$fit$whichDiscreteSL)

## check the weights
lapply(sl_fits_varset_1_baseline_exposure, function(x) sort(colMeans(x$fit$coef), decreasing=TRUE))
lapply(sl_fits_varset_11_all, function(x) sort(colMeans(x$fit$coef), decreasing=TRUE))

## average the AUCs over the 10 folds, for each
var_set_labels <- c("No markers", "IgG + IgA", "IgG3", "T Cells", "Fx Ab", "IgG + IgA + IgG3",
                    "IgG + IgA + T Cells", "IgG + IgA + IgG3 + T Cells",
                    "IgG + IgA + IgG3 + Fx Ab", "T Cells + Fx Ab", "All markers")
for (i in 1:(length(var_set_names))) { 
  this_name <- paste(unlist(strsplit(var_set_names[i], "_", fixed = TRUE))[-1], collapse = "_")
  # eval(parse(text = paste0("all_aucs_i <- as_tibble(do.call(rbind.data.frame, lapply(sl_fits_varset_", var_set_names[i], ", function(x) x$aucs)))")))
  eval(parse(text = paste0("all_aucs_i <- as_tibble(do.call(rbind.data.frame, lapply(sl_fits_varset_", var_set_names[i], ", function(x) get_all_aucs_lst(x, weights = weights_vaccine))))")))
  all_aucs_i <- all_aucs_i %>% 
    filter(!is.na(all_aucs_i$Learner))
  eval(parse(text = paste0("avg_aucs_", var_set_names[i]," <- all_aucs_i %>% 
    group_by(Learner, Screen) %>% 
    summarize(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul)) %>% 
    mutate(assay = this_name, varset_label = var_set_labels[i]) %>% 
    ungroup()")))
}

## combine into a full tibble; add a column to each that is the assay
avg_aucs <- bind_rows(avg_aucs_1_baseline_exposure, avg_aucs_2_igg_iga, avg_aucs_3_igg3,
                      avg_aucs_4_tcells, avg_aucs_5_fxab, avg_aucs_6_igg_iga_igg3,
                      avg_aucs_7_igg_iga_tcells, avg_aucs_8_igg_iga_igg3_tcells,
                      avg_aucs_9_igg_iga_igg3_fxab, avg_aucs_10_tcells_fxab, 
                      avg_aucs_11_all)

## forest plot of AUCs; this one is super nasty, but shows that it works
full_forest_plot_auc <- avg_aucs %>% 
  ggplot(aes(x = AUC, y = factor(paste0(Screen, "_", Learner, "_", assay), levels = paste0(Screen, "_", Learner, "_", assay)[order(AUC)]))) + 
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
  geom_point() +
  xlab("CV-AUC") +
  ylab("Learner")


## --------------------------------------------------------------------------------------------------------------------------------------
## FIGURE 1: forest plot of CV-AUC for the top learner and SL for each assay combination
## --------------------------------------------------------------------------------------------------------------------------------------
title_font_size <- 18
main_font_size <- 5
fig_width <- fig_height <- 2590
y_title <- 0.96
auc_forest_plot <- plot_assays(avg_aucs, type = "auc", main_font_size, main_font_size, sl_only = FALSE, immunoassay = FALSE)
png(paste0(plots_dir, "cv_auc_forest_plot_sl_plus_top_learner.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
plot_grid(auc_forest_plot$top_learner_nms_plot, auc_forest_plot$top_learner_plot, nrow = 1, align = "h") +
  draw_label("Assay combination", size = title_font_size, x = 0.075, y = y_title) +
  draw_label("Algorithm", size = title_font_size, x = 0.175, y = y_title) +
  draw_label("Screen", size = title_font_size, x = 0.25, y = y_title) +
  draw_label("CV-AUC [95% CI]", size = title_font_size, x = 0.43, y = y_title)
dev.off()

## add on immunoassay set
avg_aucs <- avg_aucs %>% 
  mutate(immunoassay_set = get_immunoassay_set(varset_label))

title_font_size <- 26
main_font_size_forest <- 31
main_font_size_lab <- 9.3
fig_width <- fig_height <- 2590
y_title <- 0.945
point_size <- 5
auc_forest_plot <- plot_assays(avg_aucs, type = "auc", main_font_size_forest = main_font_size_forest, 
                                   main_font_size_lab = main_font_size_lab,
                                   sl_only = TRUE, immunoassay = TRUE,
                                   colors = cbbPalette,
                               point_size = point_size)
ggsave(paste0(plots_dir, "cv_auc_forest_plot_sl.png"),
       plot = plot_grid(auc_forest_plot$top_learner_nms_plot, 
                        auc_forest_plot$top_learner_plot, nrow = 1, align = "h") +
         draw_label("Month 7", size = title_font_size, x = 0.14, y = y_title + 0.04, fontface = "bold") +
         draw_label("Marker Set", size = title_font_size, x = 0.1375, y = y_title, fontface = "bold") +
         draw_label("CV-AUC", size = title_font_size, x = 0.4, y = y_title + 0.04, fontface = "bold") +
         draw_label("[95% CI]", size = title_font_size, x = 0.4, y = y_title, fontface = "bold"),
       width = 50, height = 25, units = "cm")

ggsave(paste0(plots_dir, "cv_auc_forest_plot_sl.tiff"),
       plot = plot_grid(auc_forest_plot$top_learner_nms_plot, 
                        auc_forest_plot$top_learner_plot, nrow = 1, align = "h") +
         draw_label("Month 7", size = title_font_size, x = 0.14, y = y_title + 0.04, fontface = "bold") +
         draw_label("Marker Set", size = title_font_size, x = 0.1375, y = y_title, fontface = "bold") +
         draw_label("CV-AUC", size = title_font_size, x = 0.4, y = y_title + 0.04, fontface = "bold") +
         draw_label("[95% CI]", size = title_font_size, x = 0.4, y = y_title, fontface = "bold"),
       width = 50, height = 25, units = "cm")

ggsave(paste0(plots_dir, "cv_auc_forest_plot_sl.pdf"),
       plot = plot_grid(auc_forest_plot$top_learner_nms_plot, 
                        auc_forest_plot$top_learner_plot, nrow = 1, align = "h") +
         draw_label("Month 7", size = title_font_size, x = 0.14, y = y_title + 0.04, fontface = "bold") +
         draw_label("Marker Set", size = title_font_size, x = 0.1375, y = y_title, fontface = "bold") +
         draw_label("CV-AUC", size = title_font_size, x = 0.4, y = y_title + 0.04, fontface = "bold") +
         draw_label("[95% CI]", size = title_font_size, x = 0.4, y = y_title, fontface = "bold"),
       width = 50, height = 25, units = "cm")

## --------------------------------------------------------------------------------------------------------------------------------------
## FIGURE 2: forest plot of CV-R^2 for the top learner and SL for each assay combination
## --------------------------------------------------------------------------------------------------------------------------------------
## get the R-squareds
for (i in 1:length(var_set_names)) {
  this_name <- paste(unlist(strsplit(var_set_names[i], "_", fixed = TRUE))[-1], collapse = "_")
  eval(parse(text = paste0("all_r2s_i <- as_tibble(do.call(rbind.data.frame, lapply(sl_fits_varset_", var_set_names[i], ", function(x) get_all_r2s_lst(x, weights = weights_vaccine))))")))
  all_r2s_i <- all_r2s_i %>% 
    filter(!is.na(all_r2s_i$Learner))
  eval(parse(text = paste0("avg_r2s_", var_set_names[i]," <- all_r2s_i %>% 
    group_by(Learner, Screen) %>% 
    summarize(R2 = mean(R2), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul)) %>% 
    mutate(assay = this_name, varset_label = var_set_labels[i]) %>% 
    ungroup()")))
}
## combine into a full tibble; add a column to each that is the assay
avg_r2s <- bind_rows(avg_r2s_1_baseline_exposure, avg_r2s_2_igg_iga, avg_r2s_3_igg3,
                     avg_r2s_4_tcells, avg_r2s_5_fxab, avg_r2s_6_igg_iga_igg3,
                     avg_r2s_7_igg_iga_tcells, avg_r2s_8_igg_iga_igg3_tcells,
                     avg_r2s_9_igg_iga_igg3_fxab, avg_r2s_10_tcells_fxab,
                     avg_r2s_11_all)

title_font_size <- 18
main_font_size <- 5
fig_width <- fig_height <- 2590
y_title <- 0.96
r2_forest_plot <- plot_assays(avg_r2s, type = "r2", main_font_size, main_font_size,
                              sl_only = FALSE, immunoassay = FALSE)
png(paste0(plots_dir, "cv_r2_forest_plot_sl_plus_top_learner.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
plot_grid(r2_forest_plot$top_learner_nms_plot, r2_forest_plot$top_learner_plot, nrow = 1, align = "h", rel_widths = c(1, 0.55)) +
  draw_label("Assay combination", size = title_font_size, x = 0.075, y = y_title) +
  draw_label("Algorithm", size = title_font_size, x = 0.25, y = y_title) +
  draw_label("Screen", size = title_font_size, x = 0.34, y = y_title) +
  draw_label(expression(paste("CV-", R^2, " [95% CI]", sep = "")), size = title_font_size, x = 0.55, y = y_title)
dev.off()

## add on immunoassay set
avg_r2s <- avg_r2s %>% 
  mutate(immunoassay_set = get_immunoassay_set(varset_label))

title_font_size <- 22
main_font_size_forest <- 15
main_font_size_lab <- 8
fig_width <- fig_height <- 2590
y_title <- 0.96
r2_forest_plot <- plot_assays(avg_r2s, type = "r2", main_font_size_forest = main_font_size_forest, 
                                   main_font_size_lab = main_font_size_lab,
                                   sl_only = TRUE, immunoassay = TRUE,
                                   colors = cbbPalette)
png(paste0(plots_dir, "cv_r2_forest_plot_sl.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
plot_grid(r2_forest_plot$top_learner_nms_plot, 
          r2_forest_plot$top_learner_plot, nrow = 1, align = "h") +
  draw_label("Assay combination", size = title_font_size, x = 0.1, y = y_title) +
  draw_label(expression(paste("CV-", R^2, " [95% CI]", sep = "")), size = title_font_size, x = 0.38, y = y_title)
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

## ----------------------------------------------------------------------------------------------
## FIGURE 3, 4: do variable importance relative to baseline risk vars only
## ----------------------------------------------------------------------------------------------
# risk_type <- "r_squared"
risk_type <- "auc"

## set up the full fits
full_fit_11_all <- sl_fits_varset_11_all
full_fit_10_tcells_fxab <- sl_fits_varset_10_tcells_fxab
full_fit_9_igg_iga_igg3_fxab <- sl_fits_varset_9_igg_iga_igg3_fxab
full_fit_8_igg_iga_igg3_tcells <- sl_fits_varset_8_igg_iga_igg3_tcells
full_fit_7_igg_iga_tcells <- sl_fits_varset_7_igg_iga_tcells
full_fit_6_igg_iga_igg3 <- sl_fits_varset_6_igg_iga_igg3
full_fit_5_fxab <- sl_fits_varset_5_fxab
full_fit_4_tcells <- sl_fits_varset_4_tcells
full_fit_3_igg3 <- sl_fits_varset_3_igg3
full_fit_2_igg_iga <- sl_fits_varset_2_igg_iga

## reduced fit for all
reduced_fit_none <- sl_fits_varset_1_baseline_exposure

## (11) all markers
vimp_all_markers <- get_cv_vim(full_fit = full_fit_11_all,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine,
                               scale = "identity")
vimp_all_markers_avg <- get_avg_est_ci(vimp_all_markers)

## (10) T cells  + Fx Ab
vimp_tcells_fxab <- get_cv_vim(full_fit = full_fit_10_tcells_fxab,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine)
vimp_tcells_fxab_avg <- get_avg_est_ci(vimp_tcells_fxab)

## (9) IgG + IgA + IgG3 + Fx Ab
vimp_igg_iga_igg3_fxab <- get_cv_vim(full_fit = full_fit_9_igg_iga_igg3_fxab,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine)
vimp_igg_iga_igg3_fxab_avg <- get_avg_est_ci(vimp_igg_iga_igg3_fxab)

## (8) IgG + IgA + IgG3 + T cells
vimp_igg_iga_igg3_tcells <- get_cv_vim(full_fit = full_fit_8_igg_iga_igg3_tcells,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine)
vimp_igg_iga_igg3_tcells_avg <- get_avg_est_ci(vimp_igg_iga_igg3_tcells)

## (7) IgG + IgA + T cells
vimp_igg_iga_tcells <- get_cv_vim(full_fit = full_fit_7_igg_iga_tcells,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine)
vimp_igg_iga_tcells_avg <- get_avg_est_ci(vimp_igg_iga_tcells)

## (6) IgG + IgA + IgG3
vimp_igg_iga_igg3 <- get_cv_vim(full_fit = full_fit_6_igg_iga_igg3,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                               weights = weights_vaccine)
vimp_igg_iga_igg3_avg <- get_avg_est_ci(vimp_igg_iga_igg3)

## (5) Fx Ab
vimp_fxab <- get_cv_vim(full_fit = full_fit_5_fxab,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                        weights = weights_vaccine)
vimp_fxab_avg <- get_avg_est_ci(vimp_fxab)

## (4) T cells
vimp_tcells <- get_cv_vim(full_fit = full_fit_4_tcells,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                          weights = weights_vaccine)
vimp_tcells_avg <- get_avg_est_ci(vimp_tcells)

## (3) IgG3
vimp_igg3 <- get_cv_vim(full_fit = full_fit_3_igg3,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                        weights = weights_vaccine)
vimp_igg3_avg <- get_avg_est_ci(vimp_igg3)

## (2) IgG + IgA
vimp_igg_iga <- get_cv_vim(full_fit = full_fit_2_igg_iga,
                               reduced_fit = reduced_fit_none,
                               type = risk_type,
                           weights = weights_vaccine)
vimp_igg_iga_avg <- get_avg_est_ci(vimp_igg_iga)

## combine together
vimp_tibble <- tibble(assay_grp = c("All markers", 
                                    "T Cells + Fx Ab", 
                                    "IgG + IgA + IgG3 + Fx Ab", 
                                    "IgG + IgA + IgG3 + T Cells", 
                                    "IgG + IgA + T Cells",
                                    "IgG + IgA + IgG3",
                                    "Fx Ab", 
                                    "T Cells", 
                                    "IgG3",
                                    "IgG + IgA"),
                      est = c(vimp_all_markers_avg$est, 
                              vimp_tcells_fxab_avg$est, 
                              vimp_igg_iga_igg3_fxab_avg$est,
                              vimp_igg_iga_igg3_tcells_avg$est, 
                              vimp_igg_iga_tcells_avg$est,
                              vimp_igg_iga_igg3_avg$est,
                              vimp_fxab_avg$est, 
                              vimp_tcells_avg$est,
                              vimp_igg3_avg$est,
                              vimp_igg_iga_avg$est),
                      cil = c(vimp_all_markers_avg$ci[1], 
                              vimp_tcells_fxab_avg$ci[1], 
                              vimp_igg_iga_igg3_fxab_avg$ci[1],
                              vimp_igg_iga_igg3_tcells_avg$ci[1], 
                              vimp_igg_iga_tcells_avg$ci[1],
                              vimp_igg_iga_igg3_avg$ci[1],
                              vimp_fxab_avg$ci[1], 
                              vimp_tcells_avg$ci[1],
                              vimp_igg3_avg$ci[1],
                              vimp_igg_iga_avg$ci[1]),
                      ciu = c(vimp_all_markers_avg$ci[2], 
                              vimp_tcells_fxab_avg$ci[2], 
                              vimp_igg_iga_igg3_fxab_avg$ci[2],
                              vimp_igg_iga_igg3_tcells_avg$ci[2], 
                              vimp_igg_iga_tcells_avg$ci[2],
                              vimp_igg_iga_igg3_avg$ci[2],
                              vimp_fxab_avg$ci[2], 
                              vimp_tcells_avg$ci[2],
                              vimp_igg3_avg$ci[2],
                              vimp_igg_iga_avg$ci[2]))
vimp_tibble <- tibble::add_column(vimp_tibble, immunoassay_set = get_immunoassay_set(vimp_tibble$assay_grp))
## save this object for easy loading
saveRDS(vimp_tibble, paste0(results_dir, "vimp_tibble_", risk_type))
vimp_tibble <- readRDS(paste0(results_dir, "vimp_tibble_", risk_type))

title_font_size <- 26
main_font_size_forest <- 31
main_font_size_lab <- 9.3
y_title <- 0.945
point_size <- 5
if (risk_type == "r_squared") {
  lgnd_pos <- c(0, 0.8)
} else {
  lgnd_pos <- c(0.75, 0.8)
}
## forest plot of vimp, with labels for the groups
vimp_forest_plot <- vimp_tibble %>% 
  ggplot(aes(x = est, y = factor(assay_grp, levels = assay_grp[order(est, decreasing = TRUE)], labels = assay_grp[order(est, decreasing = TRUE)]))) +
  geom_errorbarh(aes(xmin = cil, xmax = ciu, color = immunoassay_set), size = point_size/2) +
  geom_point(size = point_size) +
  scale_color_manual(values = cbbPalette[-2]) + # make sure that colors match with other forest plot
  ylab("Month 7 Marker Set") +
  labs(color = "Assay set") +
  xlab(paste0("Variable importance estimate: difference in CV-", ifelse(risk_type == "r_squared", expression(R^2), "AUC"))) +
  theme(legend.position = lgnd_pos, 
        axis.text.y = element_text(size = main_font_size_forest),
        text = element_text(size = main_font_size_forest),
        axis.title = element_text(size = main_font_size_forest), 
        axis.text.x = element_text(size = main_font_size_forest),
        axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0), size = main_font_size_forest),
        plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
ggsave(filename = paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline.png"), 
       plot = vimp_forest_plot,
       width = 50, height = 25, units = "cm")
ggsave(filename = paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline.pdf"), 
       plot = vimp_forest_plot,
       width = 50, height = 25, units = "cm")

## ----------------------------------------------------------------------------------------------
## Look at one high-performing logistic regression model to see ORs, CIs
## first, screen based on a high-performing screen
## these are taken from the full SL, looking over folds; choose
## stepwise with interactions, high-correlation screen
## ----------------------------------------------------------------------------------------------
single_logistic_screen_highcor_x <- screen_highcor_plus_exposure(Y_vaccine, X_vaccine, family = "binomial", obsWeights = weights_vaccine)
X_highcor_single_logistic <- X_vaccine %>% 
  select(names(X_vaccine)[single_logistic_screen_highcor_x])
single_logistic_mod <- SL.step.interaction(Y_vaccine, X_highcor_single_logistic, family = "binomial", obsWeights = weights_vaccine)
single_logistic_mod_summ <- summary(single_logistic_mod$fit$object)
single_logistic_ors <- exp(single_logistic_mod_summ$coefficients[, 1])
single_logistic_ses <- exp(single_logistic_mod_summ$coefficients[, 2])
single_logistic_cis <- exp(cbind(single_logistic_mod_summ$coefficients[, 1] - 1.96*single_logistic_mod_summ$coefficients[, 2],
                                 single_logistic_mod_summ$coefficients[, 1] + 1.96*single_logistic_mod_summ$coefficients[, 2]))
cbind(single_logistic_ors, single_logistic_cis)

# load in results from individual vimp analysis, make a nice plot

# --------------------------------------------------------------------------------------------------------------------------------------
# set up directories, load required packages
# --------------------------------------------------------------------------------------------------------------------------------------
library("SuperLearner")
library("cvAUC")
library("tidyr")
library("dplyr")
library("cowplot")
# install if anything has changed
#  devtools::install_github("bdwilliamson/vimp", upgrade = "never")
library("vimp")
library("HVTN505")
library("kyotil")
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv
source("code/auc_plot_assays.R")
source("code/r2_plot_assays.R")
source("code/utils.R")

results_dir <- "results/"
plots_dir <- "plots/"
# should we compute results, or have they already been computed?
compute_results <- TRUE

# --------------------------------------------------------------------------------------------------------------------------------------
# load results objects:
# --------------------------------------------------------------------------------------------------------------------------------------
# baseline vars only
sl_fits_varset_1_baseline_exposure <- readRDS(paste0(results_dir, "sl_fits_varset_1_baseline_exposure.rds"))

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
Z_plus_weights <- readRDS(file = paste0(data_dir, "z_and_weights_for_505_analysis.rds"))
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
# set up SL library for IPCW
sl_lib_ipcw <- methods[!grepl("SL.glmnet", methods) & !grepl("SL.step.interaction.skinny", methods)]

# load all individual vars
var_names <- names(X_markers)
for (i in 1:length(var_names)) {
  eval(parse(text = paste0("sl_fit_var_", i, " <- readRDS(paste0(results_dir, 'sl_fits_vimp_', i + 190, '.rds'))")))
}

title_font_size <- 22
main_font_size_forest <- 20
main_font_size_lab <- 8
y_title <- 0.96
point_size <- 5
fig_width <- 30
fig_height <- 23
width_mult <- 2
height_mult <- 1.5
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# --------------------------------------------------------------------------------------------------------------------------------------
# Variable importance analysis
# Presentation: most to least important, within groups defined by full SL importance?
# So Fx ab group on top, IgG + IgA next, IgG3 next, T cells next
# --------------------------------------------------------------------------------------------------------------------------------------
risk_type <- "auc"
scale <- "logit"
for (i in 1:length(var_names)) {
  # only run if results computation changes
  if (compute_results) {
    eval(parse(text = paste0("vimp_", i, "<- get_cv_vim(full_fit = sl_fit_var_", i,
                             ", reduced_fit = sl_fits_varset_1_baseline_exposure, type = risk_type,
                             weights = weights_vaccine,
                             scale = scale,
                             C = C, Z = Z,
                             SL.library = sl_lib_ipcw,
                             vimp = TRUE)")))
    eval(parse(text = paste0("vimp_", i, "_avg <- get_avg_est_ci(vimp_", i, ")")))
    eval(parse(text = paste0("saveRDS(vimp_", i, "_avg, paste0(results_dir, 'vimp_", i, "_avg.rds'))")))
  } else {
    eval(parse(text = paste0("vimp_", i, "_avg <- readRDS(paste0(results_dir, 'vimp_", i, "_avg.rds'))")))
  }
}

# make groups of marker variables;
# make outer name of assay (all antigens), inner name of antigen
# outer names
fx_ab <- get_nms_group_all_antigens(X_markers, assays =  c("phago", "R2a", "R3a"))
igg_iga <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA"), assays_to_exclude = "IgG3")
igg3 <- get_nms_group_all_antigens(X_markers, assays = "IgG3")
tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
assay_nms <- rep(NA, length(var_names))
assay_nms[fx_ab] <- "Fx Ab"
assay_nms[igg_iga] <- "IgG + IgA"
assay_nms[igg3] <- "IgG3"
assay_nms[tcells] <- "T cells"
# inner names
all_var_supers <- rbind(var.super, var.super)
any_hiv <- get_nms_group_all_assays(X_markers, "ANYHIV") | all_var_supers$antigen == "ANYHIV"
any_vrc_env <- get_nms_group_all_assays(X_markers, "ANYVRCENV") | all_var_supers$antigen == "ANYVRCENV"
any_vrc_gag <- get_nms_group_all_assays(X_markers, "ANYVRCGAG") | all_var_supers$antigen == "ANYVRCGAG"
any_vrc_nef <- get_nms_group_all_assays(X_markers, "ANYVRCNEF") | all_var_supers$antigen == "ANYVRCNEF"
any_vrc_pol <- get_nms_group_all_assays(X_markers, "ANYVRCPOL") | all_var_supers$antigen == "ANYVRCPOL"
cmv  <- get_nms_group_all_assays(X_markers, "CMV") | all_var_supers$antigen == "CMV"
empty_ad5_vrc  <- get_nms_group_all_assays(X_markers, "EmptyAd5VRC") | all_var_supers$antigen == "EmptyAd5VRC"
gp140  <- get_nms_group_all_assays(X_markers, "gp140") | all_var_supers$antigen == "gp140"
gp41  <- get_nms_group_all_assays(X_markers, "gp41") | all_var_supers$antigen == "gp41"
v3  <- get_nms_group_all_assays(X_markers, "V3") | all_var_supers$antigen == "V3"
gp120  <- get_nms_group_all_assays(X_markers, "gp120") | all_var_supers$antigen == "gp120"
c1  <- get_nms_group_all_assays(X_markers, "C1") | all_var_supers$antigen == "C1"
v1v2  <- get_nms_group_all_assays(X_markers, "V1V2") | all_var_supers$antigen == "V1V2"
c4  <- get_nms_group_all_assays(X_markers, "C4") | all_var_supers$antigen == "C4"
p24  <- get_nms_group_all_assays(X_markers, "p24") | all_var_supers$antigen == "p24"
vrc_nef_b  <- get_nms_group_all_assays(X_markers, "VRCNEFB") | all_var_supers$antigen == "VRCNEFB"
antigen_nms <- rep(NA, length(var_names))
antigen_nms[any_hiv] <- "Any HIV"
antigen_nms[any_vrc_env] <- "Any VRC Env"
antigen_nms[any_vrc_gag] <- "Any VRC Gag"
antigen_nms[any_vrc_nef] <- "Any VRC Nef"
antigen_nms[any_vrc_pol] <- "Any VRC Pol"
antigen_nms[cmv] <- "CMV"
antigen_nms[empty_ad5_vrc] <- "Empty Ad5 VRC"
antigen_nms[gp140] <- "gp140"
antigen_nms[gp41] <- "gp41"
antigen_nms[v3] <- "V3"
antigen_nms[gp120] <- "gp120"
antigen_nms[c1] <- "C1"
antigen_nms[v1v2] <- "V1V2"
antigen_nms[c4] <- "C4"
antigen_nms[p24] <- "p24"
antigen_nms[vrc_nef_b] <- "VRC Nef B"

# make forest plot, with Fx Ab on top, then IgG + IgA, then IgG3, then T cells
# display_var_names removes assay and antigen
vimp_tibble_ind <- tibble(assay_group = assay_nms[1],
                          antigen_group = antigen_nms[1],
                          var_name = var_names[1],
                          est = vimp_1_avg$est, cil = vimp_1_avg$ci[1],
                          ciu = vimp_1_avg$ci[2],
                          p = vimp_1_avg$p,
                          greater_zero = vimp_1_avg$ci[1] > 0,
                          signif_p = vimp_1_avg$p < 0.05)
for (i in 2:length(var_names)) {
  this_est <- eval(parse(text = paste0("vimp_", i, "_avg")))
  vimp_tibble_ind <- vimp_tibble_ind %>%
    tibble::add_row(assay_group = assay_nms[i],
                    antigen_group = antigen_nms[i],
                    var_name = var_names[i],
                    est = this_est$est, cil = this_est$ci[1],
                    ciu = this_est$ci[2],
                    p = this_est$p,
                    greater_zero = this_est$ci[1] > 0,
                    signif_p = this_est$p < 0.05)
}

# --------------------------------------------------------------------------------------------------------------------------------------
# Fx Ab: one plot for each antigen type, make them entirely different plots
# --------------------------------------------------------------------------------------------------------------------------------------
# forest plot of vimp, with labels for the groups
antigen_labs_fxab <- unique((vimp_tibble_ind %>% filter(assay_group == "Fx Ab"))$antigen_group)
vimp_forest_plot_ind_fxab <- assay_antigen_plot_list(vimp_tibble_ind, assay = "Fx Ab",
                                                     antigens = antigen_labs_fxab,
                                                     risk_type = risk_type,
                                                     main_font_size = main_font_size_forest,
                                                     point_size = point_size,
                                                     x_lim = c(-0.2, 0.3),
                                                     cols = cbbPalette[c(1, 3)])
# create the final plots
for (i in 1:length(antigen_labs_fxab)) {
  file_name <- paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_fxab_", antigen_labs_fxab[i], ".png")
  ggsave(file_name, plot = plot_grid(vimp_forest_plot_ind_fxab[[i]], labels = antigen_labs_fxab[i]),
         device = "png",
         height = fig_height, width = fig_width,
         units = "cm", dpi = 300)

}
# --------------------------------------------------------------------------------------------------------------------------------------
# IgG + IgA: one plot for each antigen type
# --------------------------------------------------------------------------------------------------------------------------------------
# forest plot of vimp, with labels for the groups
antigen_labs_igg_iga <- unique((vimp_tibble_ind %>% filter(assay_group == "IgG + IgA"))$antigen_group)
vimp_forest_plot_ind_igg_iga <- assay_antigen_plot_list(vimp_tibble_ind, assay = "IgG + IgA",
                                                     antigens = antigen_labs_igg_iga,
                                                     risk_type = risk_type,
                                                     main_font_size = main_font_size_forest,
                                                     point_size = point_size,
                                                     x_lim = c(-0.2, 0.3),
                                                     cols = cbbPalette[c(1, 3)])
# create the final plots
for (i in 1:length(antigen_labs_igg_iga)) {
  file_name <- paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_igg_iga_", antigen_labs_igg_iga[i], ".png")
  ggsave(file_name, plot = plot_grid(vimp_forest_plot_ind_igg_iga[[i]], labels = antigen_labs_igg_iga[i]),
         device = "png",
         height = fig_height*(height_mult + 0.25), width = fig_width,
         units = "cm", dpi = 300)

}
# --------------------------------------------------------------------------------------------------------------------------------------
# IgG3: one plot for each antigen type
# --------------------------------------------------------------------------------------------------------------------------------------
# forest plot of vimp, with labels for the groups
antigen_labs_igg3 <- unique((vimp_tibble_ind %>% filter(assay_group == "IgG3"))$antigen_group)
vimp_forest_plot_ind_igg3 <- assay_antigen_plot_list(vimp_tibble_ind, assay = "IgG3",
                                                     antigens = antigen_labs_igg3,
                                                     risk_type = risk_type,
                                                     main_font_size = main_font_size_forest,
                                                     point_size = point_size,
                                                     x_lim = c(-0.2, 0.35),
                                                     cols = cbbPalette[c(1, 3)])
# create the final plots
for (i in 1:length(antigen_labs_igg3)) {
  file_name <- paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_igg3_", antigen_labs_igg3[i], ".png")
  ggsave(file_name, plot = plot_grid(vimp_forest_plot_ind_igg3[[i]], labels = antigen_labs_igg3[i]),
         device = "png",
         height = fig_height, width = fig_width,
         units = "cm", dpi = 300)

}
# --------------------------------------------------------------------------------------------------------------------------------------
# T cells: one plot for each antigen type
# --------------------------------------------------------------------------------------------------------------------------------------
# forest plot of vimp, with labels for the groups
antigen_labs_tcells <- unique((vimp_tibble_ind %>% filter(assay_group == "T cells"))$antigen_group)
vimp_forest_plot_ind_tcells <- assay_antigen_plot_list(vimp_tibble_ind, assay = "T cells",
                                                     antigens = antigen_labs_tcells,
                                                     risk_type = risk_type,
                                                     main_font_size = main_font_size_forest,
                                                     point_size = point_size,
                                                     x_lim = c(-0.2, 0.35),
                                                     cols = cbbPalette[c(1, 3)],
                                                     cols2 = cbbPalette[c(1, 4)])
# create the final plots
for (i in 1:length(antigen_labs_tcells)) {
  file_name <- paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_tcells_", antigen_labs_tcells[i], ".png")
  ggsave(file_name, plot = plot_grid(vimp_forest_plot_ind_tcells[[i]], labels = antigen_labs_tcells[i]),
         device = "png",
         height = fig_height, width = fig_width*(width_mult*2/3),
         units = "cm", dpi = 300)

}

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

title_font_size <- 18
main_font_size <- 5
fig_width <- fig_height <- 2590

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

## make groups of marker variables;
## make outer name of assay (all antigens), inner name of antigen
## outer names
fx_ab <- get_nms_group_all_antigens(X_markers, assays =  c("phago", "R2a", "R3a"))
igg_iga <- get_nms_group_all_antigens(X_markers, assays = c("IgG", "IgA"), assays_to_exclude = "IgG3")
igg3 <- get_nms_group_all_antigens(X_markers, assays = "IgG3")
tcells <- get_nms_group_all_antigens(X_markers, assays = c("CD4", "CD8"))
assay_nms <- rep(NA, length(var_names))
assay_nms[fx_ab] <- "Fx Ab"
assay_nms[igg_iga] <- "IgG + IgA"
assay_nms[igg3] <- "IgG3"
assay_nms[tcells] <- "T cells"
## inner names
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

## make forest plot, with Fx Ab on top, then IgG + IgA, then IgG3, then T cells
## display_var_names removes assay and antigen
vimp_tibble_ind <- tibble(assay_group = assay_nms[1], 
                          antigen_group = antigen_nms[1],
                          # var_name = make_nice_variable_name(var_names[1], 
                          #                                    all_var_supers$antigen[1],
                          #                                    all_var_supers$assay[1]),
                          var_name = var_names[1],
                          est = vimp_1_avg$est, cil = vimp_1_avg$ci[1],
                          ciu = vimp_1_avg$ci[2])
for (i in 2:length(var_names)) {
  this_est <- eval(parse(text = paste0("vimp_", i, "_avg")))
  vimp_tibble_ind <- vimp_tibble_ind %>% 
    tibble::add_row(assay_group = assay_nms[i],
                    antigen_group = antigen_nms[i],
                    # var_name = make_nice_variable_name(var_names[i], 
                    #                                    all_var_supers$antigen[i],
                    #                                    all_var_supers$assay[i]), 
                    var_name = var_names[i],
                    est = this_est$est, cil = this_est$ci[1],
                    ciu = this_est$ci[2])
}

## forest plot of vimp, with labels for the groups
vimp_forest_plot_ind_fxab <- vimp_tibble_ind %>% 
  filter(assay_group == "Fx Ab") %>% 
  ggplot(aes(x = est, y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)], labels = var_name[order(est, decreasing = TRUE)]))) +
  geom_point() +
  geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dotted") + 
  ylab("Variable name") +
  xlab(paste0("Variable importance estimate: difference in ", ifelse(risk_type == "r_squared", expression(R^2), "AUC"))) +
  facet_wrap(~antigen_group)
vimp_forest_plot_ind_igg_iga <- vimp_tibble_ind %>% 
  filter(assay_group == "IgG + IgA") %>% 
  ggplot(aes(x = est, y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)], labels = var_name[order(est, decreasing = TRUE)]))) +
  geom_point() +
  geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dotted") + 
  ylab("Variable name") +
  xlab(paste0("Variable importance estimate: difference in ", ifelse(risk_type == "r_squared", expression(R^2), "AUC")))+
  facet_wrap(~antigen_group)
vimp_forest_plot_ind_igg3 <- vimp_tibble_ind %>% 
  filter(assay_group == "IgG3") %>% 
  ggplot(aes(x = est, y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)], labels = var_name[order(est, decreasing = TRUE)]))) +
  geom_point() +
  geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dotted") + 
  ylab("Variable name") +
  xlab(paste0("Variable importance estimate: difference in ", ifelse(risk_type == "r_squared", expression(R^2), "AUC")))+
  facet_wrap(~antigen_group)
vimp_forest_plot_ind_tcells <- vimp_tibble_ind %>% 
  filter(assay_group == "T cells") %>% 
  ggplot(aes(x = est, y = factor(var_name, levels = var_name[order(est, decreasing = TRUE)], labels = var_name[order(est, decreasing = TRUE)]))) +
  geom_point() +
  geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dotted") + 
  ylab("Variable name") +
  xlab(paste0("Variable importance estimate: difference in ", ifelse(risk_type == "r_squared", expression(R^2), "AUC")))+
  facet_wrap(~antigen_group)

# png(paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
# plot_grid(vimp_forest_plot_ind_fxab, vimp_forest_plot_ind_igg_iga,
#           vimp_forest_plot_ind_igg3, vimp_forest_plot_ind_tcells,
#           labels = c("Fx Ab", "IgG + IgA", "IgG3", "T cells"))
# dev.off()
png(paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_fxab.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
vimp_forest_plot_ind_fxab
dev.off()
png(paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_igg_iga.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
vimp_forest_plot_ind_igg_iga
dev.off()
png(paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_igg3.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
vimp_forest_plot_ind_igg3
dev.off()
png(paste0(plots_dir, "vimp_forest_plot_", risk_type, "_rel_to_baseline_ind_tcells.png"), width = 2*fig_width, height = fig_height, units = "px", res = 300)
vimp_forest_plot_ind_tcells
dev.off()

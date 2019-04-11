## results from full SL analysis
## includes CV-AUC plots

## ------------------------------------------------
## set up directories, load required packages
## ------------------------------------------------
library("SuperLearner")
library("cvAUC")
library("tidyr")
library("dplyr")
library("cowplot")
method <- "method.CC_nloglik" # since SuperLearner relies on this to be in GlobalEnv

results_dir <- "results/"

## load results objects:
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

## get SL and the top individual learner/screen combo from each assay
top_learners <- avg_aucs %>% 
  group_by(assay) %>% 
  arrange(desc(AUC), .by_group = TRUE) %>% # arrange in descending order by assay group
  filter((Learner == "SL" & Screen == "All") | row_number() == 1) %>%  # select only the SL row and the top performing learner 
  mutate(learner_nm = make_nice_learner_name(Learner), screen_nm = make_nice_screen_name(Screen)) %>% 
  ungroup() %>% 
  arrange(desc(AUC))

## another forest plot with the groups/learners from the analysis plan
top_learner_plot <- top_learners %>% 
  ggplot(aes(x = AUC, y = factor(paste0(Screen, "_", Learner, "_", assay), levels = paste0(Screen, "_", Learner, "_", assay)[order(AUC)],
                                 labels = paste0(varset_label, " ", learner_nm, " ", screen_nm)))) + 
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
  geom_point() +
  xlab("CV-AUC") +
  ylab("") + 
  xlim(c(0, 1)) +
  theme(axis.text.y = element_blank(),
        plot.margin=unit(c(2,0,0,0),"cm"))
## separate plot with nice names, printed values of the AUCs
top_row <- 0.875
spacing <- 0.06
num_word_space <- 0.1
title_font_size <- 8
main_font_size <- 11
top_learners_labels <- top_learners %>% 
  ungroup() %>% 
  mutate(lab_auc = paste0(round(AUC, 3), " [", round(ci_ll, 3), ", ", round(ci_ul, 3), "]")) %>%
  select(varset_label, learner_nm, screen_nm, lab_auc)
## melt to make a single "value" column
top_learners_labels$var <- 1
top_learners_labs <- melt(top_learners_labels, id.var = "var")
## tack on x, y coordinates
top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1, function(x) which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_auc"))) - 1 + c(0, 0.3, -0.3, 0)[which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_auc")))])
top_learners_labs$y_coord <- rep(as.numeric(rownames(top_learners))[order(top_learners$AUC, decreasing = FALSE)], dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
top_learner_nms_plot <- top_learners_labs %>% 
  ggplot(aes(x = x_coord, y = y_coord, label = value)) +
  geom_text(size = 3.5, hjust = 0, vjust = 0.5) +
  xlim(c(min(top_learners_labs$x_coord) - 1, max(top_learners_labs$x_coord)) + 1) +
  ylim(c(1, max(top_learners_labs$y_coord) + 1)) +
  theme(legend.position="", 
        axis.line=element_blank(),
        axis.text=element_blank(),
        text = element_text(size=3.5),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

plot_grid(top_learner_nms_plot, top_learner_plot, nrow = 1, align = "h") +
  draw_label("Assay combination", size = 1.2*title_font_size, x = 0.05, y = 0.95) +
  draw_label("Algorithm", size = 1.2*title_font_size, x = 0.15, y = 0.95) +
  draw_label("Screen", size = 1.2*title_font_size, x = 0.25, y = 0.95) +
  draw_label("CV-AUC [95% CI]", size = 1.2*title_font_size, x = 0.35, y = 0.95)

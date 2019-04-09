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

## load results object:
## length 10 list (one for each random start)
full_sl_obj <- readRDS(paste0(results_dir, "sl_fits.rds"))

## check the discrete SL for each fold
lapply(full_sl_obj, function(x) x$fit$whichDiscreteSL)

## check the weights
lapply(full_sl_obj, function(x) sort(colMeans(x$fit$coef), decreasing=TRUE))

## average the AUCs over the 10 folds
all_aucs <- as_tibble(do.call(rbind.data.frame, lapply(full_sl_obj, function(x) x$aucs)))
avg_aucs <- all_aucs %>% 
  group_by(Learner, Screen) %>% 
  summarize(AUC = mean(AUC), ci_ll = mean(ci_ll), ci_ul = mean(ci_ul)) %>% 
  ungroup()

## forest plot of AUCs; use better plot names
full_forest_plot_auc <- avg_aucs %>% 
  ggplot(aes(x = AUC, y = factor(paste0(Screen, "_", Learner), levels = paste0(Screen, "_", Learner)[order(AUC)]))) + 
  geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
  geom_point() +
  xlab("CV-AUC") +
  ylab("Learner")

## another forest plot with the groups/learners from the analysis plan

## first, get the SL and top-performing model

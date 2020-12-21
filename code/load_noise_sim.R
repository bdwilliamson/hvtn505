# results from noise simulation

# -------------------------------------------
# load required functions and packages
# -------------------------------------------
library("here")
library("SuperLearner")
library("cvAUC")
library("tidyr")
library("dplyr")
library("ggplot2")
library("cowplot")
theme_set(cowplot::theme_cowplot())
library("data.table")
library("tibble")
# install if anything has changed
#  devtools::install_github("bdwilliamson/vimp", upgrade = "never")
library("vimp")

source(here("code", "noise_screens.R"))

# -------------------------------------------
# load the simulation output
# -------------------------------------------
aipw_output_files <- here("results/noise_sim", 
                      paste0("sl_fits_noise_aipw_", 1:50, ".rds")
                      )
ipw_output_files <- here("results/noise_sim", 
                          paste0("sl_fits_noise_ipw_", 1:50, ".rds")
                         )
read_func <- function(file, method, indx) {
  # read in the file
  lst <- readRDS(file)
  # average AUCs over all 10 
  all_aucs <- as_tibble(rbindlist(lapply(lst, function(x) x$aucs)))
  avg_aucs <- all_aucs %>% 
    group_by(Learner, Screen) %>% 
    summarize(AUC = mean(AUC), se = mean(se), .groups = "drop") %>% 
    mutate(
      ci_ll = plogis(qlogis(AUC) - se * 1 / (AUC + AUC ^ 2) * qnorm(0.975)),
      ci_ul = plogis(qlogis(AUC) + se * 1 / (AUC + AUC ^ 2) * qnorm(0.975)),
      est_type = method, mc_id = indx
      )
  avg_aucs
}
aipw_output_tib <- rbindlist(
  sapply(1:length(aipw_output_files),
         function(i) {
           read_func(aipw_output_files[i], method = "aipw", indx = i)
         },
         simplify = FALSE
         )
  )
ipw_output_tib <- rbindlist(
  sapply(1:length(ipw_output_files),
         function(i) {
           read_func(ipw_output_files[i], method = "ipw", indx = i)
         },
         simplify = FALSE
         )
  )
output_tib <- bind_rows(aipw_output_tib, ipw_output_tib)
# -------------------------------------------
# boxplots with CV-AUC averaged over the 10 
# random starts, for the 50 replicates
# -------------------------------------------
plot_tib <- output_tib %>% 
  mutate(est_fct = factor(est_type, levels = unique(est_type),
                          labels = c("AIPW", "IPW")))

cv_auc_plot_sl <- plot_tib %>% 
  filter(Learner == "SL") %>% 
  ggplot(aes(x = est_fct, y = AUC)) +
  geom_boxplot() +
  geom_point(position = position_dodge2(width = 0.1)) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  ylab("CV-AUC") +
  xlab("Estimator")

ggsave(filename = here("plots", "noise_sim_cv_auc.png"),
       plot = cv_auc_plot_sl,
       width = 10, height = 10, units = "cm")

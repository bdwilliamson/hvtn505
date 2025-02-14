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
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 5,
                    help = "number of replicates for each job")
args <- parser$parse_args()

# -------------------------------------------
# load the simulation output
# -------------------------------------------
num_jobs <- args$nreps_total / args$nreps_per_job * 4
aipw_output_files <- here("results", "noise_sim", 
                      paste0("aipw_", 1:num_jobs, ".rds")
                      )
ipw_output_files <- here("results", "noise_sim", 
                          paste0("ipw_", 1:num_jobs, ".rds")
                         )
read_func <- function(x) {
  tryCatch(readRDS(x), error = function(e) tibble(mc_id = NA, n = NA, est = NA, 
                                                  cil = NA, ciu = NA, 
                                                  test = NA, p_value = NA))
}
aipw_output_tib <- as_tibble(
  rbindlist(
    lapply(as.list(aipw_output_files), read_func)
  )
)
ipw_output_tib <- as_tibble(
  rbindlist(
    lapply(as.list(ipw_output_files), read_func)
  )
)
output_tib <- bind_rows(aipw_output_tib %>% 
                          filter(!is.na(mc_id)) %>% 
                          mutate(est_type = "AIPW"), 
                        ipw_output_tib %>% 
                          filter(!is.na(mc_id)) %>% 
                          mutate(est_type = "IPW"))
# -------------------------------------------
# boxplots with CV-AUC averaged over the 10 
# random starts, for the 50 replicates
# -------------------------------------------
plot_tib <- output_tib %>% 
  mutate(est_fct = factor(est_type, levels = unique(est_type),
                          labels = c("AIPW", "IPW"))) %>% 
  group_by(n, est_fct) %>%
  summarize(mn_est = mean(est), se = sd(est, na.rm = TRUE) / sqrt(1000), 
            .groups = "drop")

cv_auc_plot <- plot_tib %>% 
  ggplot(aes(x = factor(n), y = mn_est, shape = est_fct)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mn_est - 1.96 * se, 
                     ymax = mn_est + 1.96 * se),
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +
  ylab("CV-AUC") +
  xlab("Sample size") +
  labs(shape = "Estimator")

ggsave(filename = here("plots", "noise_sim_cv_auc.png"),
       plot = cv_auc_plot,
       width = 10, height = 10, units = "cm")

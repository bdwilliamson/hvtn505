#!/usr/local/bin/Rscript

# simulation to test AIPW vs IPW estimation of CV-AUC in a two-phase sample

# ------------------------------------
# load required functions and packages
# ------------------------------------
library("SuperLearner")
# only run this if necessary to update package
# devtools::install_github("bdwilliamson/vimp", upgrade = "never")
library("vimp")
library("methods")
library("argparse")
library("xgboost")
library("ranger")
library("gam")
library("dplyr")
library("tibble")
library("data.table")
library("nloptr")
library("glmnet")
library("here")

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  job_id <- 1
  code_prefix <- "code"
  output_prefix <- ""
} else {
  job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  code_prefix <- ""
  output_prefix <- "/fh/fast/huang_y/bwillia2/hvtn505/"
}

source(here(code_prefix, "gen_simdata.R"))
source(here(code_prefix, "run_sim.R"))

# ---------------------------------------------
# pull in command-line arguments,
# set up the simulation
# ---------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "aipw",
                    help = "the name of the simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 5,
                    help = "number of replicates for each job")
args <- parser$parse_args()

print(paste0("Running noise sim using ",
toupper(args$sim_name), "."))

if (!grepl("a", args$sim_name)) {
  source(here(code_prefix, "measure_auc_ipw.R"))
  measure_auc <- measure_auc_ipw
}

# set up static args
V <- 5
SL_V <- 5

# set up dynamic args
ns <- c(200, 500, 1000, 3000)
nreps_per_combo <- args$nreps_total/args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, n = ns)

# get current dynamic args
current_dynamic_args <- param_grid[job_id, ]

# set up SuperLearner library
learner_lib <- c("SL.ranger", "SL.glm", "SL.mean")
#----------------------------------------
# run the simulation nreps_per_job times
#----------------------------------------
current_seed <- current_dynamic_args$mc_id + current_dynamic_args$n + job_id
print(current_seed)
set.seed(current_seed)
system.time(
    sim_output <- sapply(
        1:args$nreps_per_job,
        function(i) {
            run_sim_once(
                iteration = i +
                args$nreps_per_job *
                (current_dynamic_args$mc_id - 1),
                n = current_dynamic_args$n,
                xdim = 50,
                zdim = 2,
                V = V, SL_V = SL_V,
                learner_lib = learner_lib
            )
        }, simplify = FALSE
    )
)
sim_output_tib <- tibble::as_tibble(rbindlist(sim_output))
filename <- paste0(
    output_prefix,
    args$sim_name, "_", job_id, ".rds"
)
if (!dir.exists(output_prefix)) {
    dir.create(output_prefix)
}
saveRDS(
    sim_output_tib,
    file = filename
)

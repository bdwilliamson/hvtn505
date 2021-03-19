# Reproducing the analyses in Niedich et al. (2019)

The code in this folder can be used to reproduce the analyses in [Niedich et al. (2019)](https://doi.org/10.1172/JCI126391). While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments. All analyses for the original paper were done using R version 3.5.3 and `vimp` package version 1.3.0; the current build of the vignette (and the updated code in this repository) uses R version 4.0.2 and `vimp` package version 2.1.8.

The `R` scripts contain the code necessary to run Super Learners and obtain risk predictions based on various sets of baseline and immunoassay variables:

* `gen_analysis_dataset.R` reads in the baseline data from Fong et al. (2018), computed inverse probability of sampling weights (for sampling into the case-cohort sample), and creates a final dataset with the baseline data and weights;
* `sl_screens.R` defines various screening algorithms (each restricts to a maximum of 4 immunoassay variables and forces the baseline variables `age`, `BMI`, and `bhvrisk` into the final model) and learning algorithms for the Super Learner, and finally creates a library consisting of the combinations of these screens and learners;
* `utils.R` contains helper functions for computing estimates of cross-validated AUC, running the cross-validated Super Learner (with nice output formatting), and creating plots;
* `run_sl.R` runs the cross-validated Super Learner based on all variables for each of 10 random seeds;
* `run_sl_assays.R` runs the cross-validated Super Learner for a set of immunoassay variables plus baseline variables and each of 10 random seeds (note that if you are not running this on a Slurm cluster, you will need to edit line 158 to pull in the correct set of variables to run on);
* `run_sl_vimp.R` runs the cross-validated Super Learner for a set of immunoassay variables or an individual immunoassay variable (both with baseline variables) for each of 10 random seeds (note that if you are not running this on a Slurm cluster, you will need to edit line 116 to pull in the correct set of variables to run on);
* `load_full_sl.R` loads the results of the full Super Learner and variable importance analysis using groups;
* `load_individual_vimp_analysis.R` loads the results of the individual assay variable importance analysis;
* `auc_plot_assays.R`, `r2_plot_assays.R`, and `plot_assays.R` provide additional plotting functions;
* other `R` scripts for testing the augmented IPW vs IPW approach (`gen_simdata.R`, `measure_auc_ipw.R`, `noise_screens.R`, `run_sim.R`, `load_noise_sim.R`, `run_sl_assays_sim.R`, `run_sl_assays_test_size.R`)

The bash scripts run the analyses:

* `submit_sl.sh` submits the full Super Learner to the HPC cluster (and calls `run_sl.sh`);
* `submit_sl_assays.sh` submits an array of jobs to the HPC cluster, one job for each immunoassay set (and calls `run_sl_assays.sh`);
* `submit_sl_vimp_assay.sh` submits an array of jobs to the HPC cluster, one job for each immunoassay set necessary for variable importance (that hasn't been run already; this code calls `run_sl_vimp.sh`)
* `submit_sl_groups.sh` submits an array of jobs to the HPC cluster, based on groups of variables for the variable importance analysis (groups that haven't already been run in `submit_sl_assays.sh` or `submit_sl_vimp_assay.sh`; this code calls `run_sl_vimp.sh`);
* `submit_sl_indiv.sh` submits an array of jobs to the HPC cluster, based on invididual variables for the variable importance analysis (this calls `run_sl_vimp.sh`);
* other bash scripts for testing the augmented IPW vs IPW approach (`submit_all_sims.sh`, `submit_sim.sh`, `submit_sl_assays_noise.sh`, `submit_sl-assays_test.sh`).

Again, all analysis code in this repository assumes that you are running on a HPC; in particular, the Fred Hutch `gizmo` cluster. If you are not, you need to modify all of the batch submission scripts `submit_*.sh` to load the appropriate `R` version for your HPC, and/or run the analyses locally using `run_*.sh` (making sure to appropriately edit the `.R` files to take in different command-line arguments for assay or individual variable runs).

## Obtaining an HIV infection risk estimator using all variables

To run the Super Learner on all available variables and obtain estimates of cross-validated AUC, run the following from the command line:
```{bash}
chmod u+x submit_sl.sh
./submit_sl.sh
```  
This submits the full Super Learner analysis to the cluster, requesting 10 cores and allowing the job to run for a maximum of 7 days.

## Obtaining an HIV infection risk estimator using sets of immunoassay variables

To run the Super Learner on each set of immunoassay variables and obtain estimates of cross-validated AUC, run the following from the command line:
```{bash}
chmod u+x submit_sl_assays.sh
./submit_sl_assays.sh
```
This submits a job array with 14 jobs, each requesting 10 cores and running for a maximum of 7 days; these correspond to the following sets of variables (specified in `run_sl_assays.R`):
1. baseline variables only
2. IgG + IgA variables
3. IgG3 variables
4. T cell variables
5. Functional antibody variables
6. Sets 1 + 2 + 3
7. Sets 1 + 2 + 4
8. Sets 1 + 2 + 3 + 4
9. Sets 1 + 2 + 3 + 5
10. Sets 1 + 4 + 5
11. All variables
12. Sets 1 + 3 + 5 (for variable importance)
13. Sets 1 + 2 + 5 (for variable importance)
14. Sets 1 + 3 + 4 + 5 (for variable importance)

## Obtaining an HIV infection risk estimator based on individual immunoassay variables

To run the Super Learner on each individual immunoassay variable and obtain estimates of cross-validated AUC, run the following from the command line:
```{bash}
chmod u+x submit_sl_indiv.sh
./submit_sl_indiv.sh
```
This submits a job array with 419 jobs, each requesting 10 cores and running for a maximum of 7 days. Each job corresponds to studying an individual assay variable (combined with the baseline risk variables).

## Compiling results and creating plots

Once the first two analyses (the full Super Learner and the immunoassay groups) have finished, run `load_full_sl.R` to compile the results and create plots. This can be run either interactively or from the command line using `Rscript load_full_sl.R`

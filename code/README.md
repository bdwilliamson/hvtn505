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

*  

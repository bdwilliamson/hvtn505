## utility functions

## get all AUCs for SL, discrete SL, and individual algorithms
get_all_aucs <- function(sl_fit) {
  # Get the cvAUC of the SuperLearner predictions
  sl_auc <- cvAUC::ci.cvAUC(predictions=sl_fit$SL.predict, labels=sl_fit$Y, folds=sl_fit$folds)
  out <- data.frame(Learner="SL", Screen="All", AUC=sl_auc$cvAUC, ci_ll=sl_auc$ci[1], ci_ul=sl_auc$ci[2])
  
  # Get the cvAUC of the Discrete SuperLearner predictions
  discrete_sl_auc <- cvAUC::ci.cvAUC(predictions=sl_fit$discreteSL.predict, labels=sl_fit$Y, folds=sl_fit$folds)
  out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", AUC=discrete_sl_auc$cvAUC, ci_ll=discrete_sl_auc$ci[1], ci_ul=discrete_sl_auc$ci[2]))
  
  # Get the cvAUC of the individual learners in the library
  get_individual_auc <- function(sl_fit, col) {
    if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
    alg_auc <- cvAUC::ci.cvAUC(predictions = sl_fit$library.predict[, col], labels = sl_fit$Y, folds = sl_fit$folds)
    ## get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
    data.frame(Learner = alg, Screen = screen, AUC = alg_auc$cvAUC, ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2])
  }
  other_aucs <- plyr::ldply(1:ncol(sl_fit$library.predict), function(x) get_individual_auc(sl_fit, x))
  rbind(out, other_aucs)
}

## get the SL and top-performing model for a given Learner + Screen combination


## run CV.SuperLearner for one given random seed
run_cv_sl_once <- function(seed, Y, X_mat, family, obsWeights, sl_lib, method, cvControl, innerCvControl, vimp = FALSE) {
  set.seed(seed)
  fit <- SuperLearner::CV.SuperLearner(Y = Y, X = X_mat, family = family, 
                                                   obsWeights = obsWeights, SL.library = sl_lib, 
                                                   method = method, cvControl = cvControl, 
                                                   innerCvControl = innerCvControl)
  aucs <- get_all_aucs(fit)
  ret_lst <- list(fit = fit, folds = fit$folds, aucs = aucs)
  if (vimp) {
    ret_lst <- list(fit = fit$SL.predict, folds = fit$folds, aucs = aucs)
  }
  return(ret_lst)
}
## run SuperLearner given the set of results from the CV.SL with all markers
## full fit is only the "fit" object from the CV.SL results object
run_reduced_cv_sl_once <- function(seed, full_fit, X_mat, family, obsWeights, sl_lib, method, innerCvControl, vimp = TRUE) {
  ## use the same folds as the CV.SuperLearner
  fold_row_nums <- as.vector(do.call(cbind, full_fit$folds))
  folds_init <- rep(as.numeric(names(full_fit$folds)), each = length(y)/length(full_fit$folds))
  folds_mat <- cbind(fold_row_nums, folds_init)
  folds <- folds_mat[order(folds_mat[, 1]), 2]
  
  ## set the seed, run the SL
  set.seed(seed)
  for (v in 1:)
  
  
}

## get names for multiple assays, all antigens
get_nms_group_all_antigens <- function(X, assays) {
  ## set all vars to be false
  vars <- rep(FALSE, ncol(X))
  ## set vars with assay in name to be true
  ## may be more than one
  for (i in 1:length(assays)) {
    vars[grepl(assays[i], names(X))] <- TRUE
  }
  return(vars)
}

## get the names of a group, for importance
get_nms_group <- function(X, assay, antigen) {
  vars <- rep(TRUE, ncol(X)) ## initially include all vars
  vars[grepl(assay, names(X)) & grepl(antigen, names(X))] <- FALSE # remove the ones in assay x antigen combo, for vimp
  return(vars)
}
## get the names of an individual, for importance
get_nms_ind <- function(X, nm_ind) {
  vars <- rep(TRUE, ncol(X))
  vars[grepl(nm_ind, names(X))] <- FALSE ## remove the one, for vimp
  return(vars)
}

## make nice names for Learner/Screen combos
remove_str <- function(x, str) {
  if (length(x) > 1) {
    return(x[!grepl(str, x)])
  } else {
    return(x)
  }
}
make_nice_learner_name <- function(learners) {
  no_dots <- strsplit(as.character(learners), ".", fixed = TRUE) # split on the dots
  no_sl <- lapply(no_dots, function(x) remove_str(x, "SL")) # remove SL if length is > 1
  no_skinny <- lapply(no_sl, function(x) remove_str(x, "skinny")) # remove "skinny" if it's there
  learner_nms <- unlist(lapply(no_skinny, function(x) paste(x, collapse = "_")))
  return(learner_nms)
}
make_nice_screen_name <- function(screens) {
  no_underscores <- strsplit(as.character(screens), "_", fixed = TRUE)
  no_screen_plus_exposure <- lapply(no_underscores, function(x) x[!grepl("screen", x) & !grepl("plus", x) & !grepl("exposure", x)])
  no_all <- lapply(no_screen_plus_exposure, function(x) remove_str(x, "All"))
  screen_nms <- unlist(lapply(no_all, function(x) paste(x, collapse = "_")))
  return(screen_nms)
}

## get the cv vim for each fold
get_fold_cv_vim <- function(full_fit, reduced_fit, x, type) {
  ## get the outcome, folds
  y <- full_fit[[x]]$fit$Y
  fold_row_nums <- as.vector(do.call(cbind, full_fit[[x]]$fit$folds))
  folds_init <- rep(as.numeric(names(full_fit[[x]]$fit$folds)), each = length(y)/length(full_fit[[x]]$fit$folds))
  folds_mat <- cbind(fold_row_nums, folds_init)
  folds <- folds_mat[order(folds_mat[, 1]), 2]
  
  ## get the full, reduced predictions from the two CV objects
  full_fits <- lapply(as.list(1:length(full_fit[[x]]$fit$folds)), function(i) full_fit[[x]]$fit$SL.predict[folds == i])
  redu_fits <- lapply(as.list(1:length(reduced_fit[[x]]$fit$folds)), function(i) reduced_fit[[x]]$fit$SL.predict[folds == i])
  ys <- lapply(as.list(1:length(full_fit[[x]]$fit$folds)), function(i) full_fit[[x]]$fit$Y[folds == i])
  
  ## variable importance
  vim_est <- cv_vim(Y = y, 
                    f1 = full_fits,
                    f2 = redu_fits,
                    folds = folds,
                    type = type,
                    run_regression = FALSE,
                    alpha = 0.05)
  return(vim_est)
}
## get the CV vim averaged over the 10 folds
get_cv_vim <- function(full_fit, reduced_fit, type) {
  ## get the cv vim for each fold
  all_cv_vims <- lapply(as.list(1:length(full_fit)), get_fold_cv_vim, full_fit = full_fit, reduced_fit = reduced_fit, type = type)
  return(all_cv_vims)
}
## get estimate, CI based on averaging over the 10 random starts
get_avg_est_ci <- function(vimp_lst) {
  ests <- unlist(lapply(vimp_lst, function(x) x$est))
  cis <- do.call(rbind, lapply(vimp_lst, function(x) x$ci))
  est <- mean(ests)
  ci <- colMeans(cis)
  return(list(est = est, ci = ci))
}
## get estimate, CI based on averaging over 10 random starts for individual risk estimators
get_avg_est_ci_individual_risk <- function(vimp_lst) {
  
}

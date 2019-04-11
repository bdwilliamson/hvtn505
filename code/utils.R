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

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
    alg <- alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)]
    screen <- paste0(alg_screen_string[!grepl(".", alg_screen_string, fixed = TRUE)], collapse = "_")
    data.frame(Learner = alg, Screen = screen, AUC = alg_auc$cvAUC, ci_ll = alg_auc$ci[1], ci_ul = alg_auc$ci[2])
  }
  other_aucs <- plyr::ldply(1:ncol(sl_fit$library.predict), function(x) get_individual_auc(sl_fit, x))
  rbind(out, other_aucs)
}

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
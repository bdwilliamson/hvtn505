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
  other_aucs <- plyr::ldply(1:ncol(sl_fit$library.predict), function(col) {
    if(any(is.na(sl_fit$library.predict[,col]))) return(NULL)
    alg_auc <- cvAUC::ci.cvAUC(predictions=sl_fit$library.predict[,col], labels=sl_fit$Y, folds=sl_fit$folds)
    alg <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed=TRUE)[[1]][1]
    screen <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed=TRUE)[[1]][2]
    data.frame(Learner=alg, Screen=screen, AUC=alg_auc$cvAUC, ci_ll=alg_auc$ci[1], ci_ul=alg_auc$ci[2])
  })
  rbind(out, other_aucs)
}

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

## get all R^2s for SL, discrete SL, and individual algorithms
one_r2 <- function(preds, Y) {
  var_y <- mean((Y - mean(Y))^2)
  mse <- mean((Y - preds)^2)
  ic_mse <- (Y - preds)^2 - mse
  ic_var <-  (Y - mean(Y))^2 - var_y
  grad <- matrix(c(1/mse, -1/var_y), nrow = 1)
  ic <- cbind(ic_mse, ic_var)
  se_log_r2 <- sqrt(grad %*% t(ic) %*% ic %*% t(grad))/length(Y)
  est <- 1 - mse/var_y
  ci_low <- 1 - exp(log(mse / var_y) + 1.96 * se_log_r2)
  ci_high <- 1 - exp(log(mse / var_y) - 1.96 * se_log_r2)
  data.frame(r2 = est, cil = ci_low, ciu = ci_high, se = se_log_r2)
}
cv_r2 <- function(preds, Y, folds) {
  V <- length(folds)
  fold_row_nums <- as.vector(do.call(cbind, folds))
  folds_init <- rep(as.numeric(names(folds)), each = length(Y)/length(folds))
  folds_mat <- cbind(fold_row_nums, folds_init)
  folds_numeric <- folds_mat[order(folds_mat[, 1]), 2]
  ests_cis <- do.call(rbind.data.frame, lapply(as.list(1:V), function(v) {
    one_r2(preds = preds[folds_numeric == v], Y[folds_numeric == v])
  }))
  est <- colMeans(ests_cis)[1]
  se <- colMeans(ests_cis)[4]
  ci_low <- 1 - exp(log(1 - est) + 1.96 * se)
  ci_high <- 1 - exp(log(1 - est) - 1.96 * se)
  return(list(r2 = est, ci = c(ci_low, ci_high)))
}
get_all_r2s <- function(sl_fit) {
  # get the CV-R^2 of the SuperLearner predictions
  sl_r2 <- cv_r2(preds = sl_fit$SL.predict, Y = sl_fit$Y, folds = sl_fit$folds)
  out <- data.frame(Learner="SL", Screen="All", R2 = sl_r2$r2, ci_ll = sl_r2$ci[1], ci_ul=sl_r2$ci[2])
  
  # Get the CV-R2 of the Discrete SuperLearner predictions
  discrete_sl_r2 <- cv_r2(preds = sl_fit$discreteSL.predict, Y = sl_fit$Y, folds = sl_fit$folds)
  out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", R2 = discrete_sl_r2$r2, ci_ll = discrete_sl_r2$ci[1], ci_ul = discrete_sl_r2$ci[2]))
  
  # Get the cvr2 of the individual learners in the library
  get_individual_r2 <- function(sl_fit, col) {
    if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
    alg_r2 <- cv_r2(preds = sl_fit$library.predict[, col], Y = sl_fit$Y, folds = sl_fit$folds)
    ## get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
    data.frame(Learner = alg, Screen = screen, R2 = alg_r2$r2, ci_ll = alg_r2$ci[1], ci_ul = alg_r2$ci[2])
  }
  other_r2s <- plyr::ldply(1:ncol(sl_fit$library.predict), function(x) get_individual_r2(sl_fit, x))
  rbind(out, other_r2s)
}

get_all_r2s_lst <- function(sl_fit_lst) {
  # get the CV-R^2 of the SuperLearner predictions
  sl_r2 <- cv_r2(preds = sl_fit_lst$fit$SL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds)
  out <- data.frame(Learner="SL", Screen="All", R2 = sl_r2$r2, ci_ll = sl_r2$ci[1], ci_ul=sl_r2$ci[2])
  
  # Get the CV-R2 of the Discrete SuperLearner predictions
  discrete_sl_r2 <- cv_r2(preds = sl_fit_lst$fit$discreteSL.predict, Y = sl_fit_lst$fit$Y, folds = sl_fit_lst$fit$folds)
  out <- rbind(out, data.frame(Learner="Discrete SL", Screen="All", R2 = discrete_sl_r2$r2, ci_ll = discrete_sl_r2$ci[1], ci_ul = discrete_sl_r2$ci[2]))
  
  # Get the cvr2 of the individual learners in the library
  get_individual_r2 <- function(sl_fit, col) {
    if(any(is.na(sl_fit$library.predict[, col]))) return(NULL)
    alg_r2 <- cv_r2(preds = sl_fit$library.predict[, col], Y = sl_fit$Y, folds = sl_fit$folds)
    ## get the regexp object
    alg_screen_string <- strsplit(colnames(sl_fit$library.predict)[col], "_", fixed = TRUE)[[1]]
    alg <- tail(alg_screen_string[grepl(".", alg_screen_string, fixed = TRUE)], n = 1)
    screen <- paste0(alg_screen_string[!grepl(alg, alg_screen_string, fixed = TRUE)], collapse = "_")
    data.frame(Learner = alg, Screen = screen, R2 = alg_r2$r2, ci_ll = alg_r2$ci[1], ci_ul = alg_r2$ci[2])
  }
  other_r2s <- plyr::ldply(1:ncol(sl_fit_lst$fit$library.predict), function(x) get_individual_r2(sl_fit_lst$fit, x))
  rbind(out, other_r2s)
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
run_reduced_cv_sl_once <- function(seed, Y, X_mat, family, obsWeights, sl_lib, method, innerCvControl, vimp = TRUE) {
  ## pull out the correct set of fitted values
  set.seed(4747)
  seeds <- round(runif(10, 1000, 10000))
  indx <- which(seed == seeds)
  full_fit <- Y[[indx]]$fit
  ## use the same folds as the CV.SuperLearner
  fold_row_nums <- as.vector(do.call(cbind, full_fit$folds))
  folds_init <- rep(as.numeric(names(full_fit$folds)), each = dim(X_mat)[1]/length(full_fit$folds))
  folds_mat <- cbind(fold_row_nums, folds_init)
  folds <- folds_mat[order(folds_mat[, 1]), 2]
  V <- length(unique(folds))
  
  ## set the seed, run the SL
  set.seed(seed)
  preds_lst <- vector("list", length = V)
  for (v in 1:V) {
    ## run an SL of full fit on reduced set of predictors
    inner_sl <- SuperLearner::SuperLearner(Y = full_fit$SL.predict[folds != v], X = X_mat[folds != v, , drop = FALSE],
                                           newX = X_mat[folds == v, , drop = FALSE],
                                           family = family, obsWeights = obsWeights[folds != v], SL.library = sl_lib,
                                           method = method, cvControl = innerCvControl)
    ## get the predicted values on the vth fold
    preds_lst[[v]] <- inner_sl$SL.predict
  } 
  ## make a vector out of the predictions
  preds_mat <- do.call(rbind.data.frame, lapply(preds_lst, function(x) cbind.data.frame(x, row_num = as.numeric(rownames(x)))))
  preds_mat_ordered <- preds_mat[order(preds_mat$row_num), ]
  ## return
  ret_lst <- list(fit = preds_mat_ordered[, 1], folds = full_fit$folds)
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
  get_reduced_fit <- function(i) {
    if (type == "r_squared" & is.null(reduced_fit[[x]]$fit$Y)) {
      tryCatch(reduced_fit[[x]]$fit[folds == i], error = function(e) rep(NA, sum(folds == i)))
    } else {
      reduced_fit[[x]]$fit$SL.predict[folds == i]  
    }
  }
  full_fits <- lapply(as.list(1:length(full_fit[[x]]$fit$folds)), function(i) full_fit[[x]]$fit$SL.predict[folds == i])
  redu_fits <- lapply(as.list(1:length(full_fit[[x]]$fit$folds)), function(i) get_reduced_fit(i))
  ys <- lapply(as.list(1:length(full_fit[[x]]$fit$folds)), function(i) full_fit[[x]]$fit$Y[folds == i])
  
  ## variable importance
  vim_est <- tryCatch(cv_vim(Y = y, 
                    f1 = full_fits,
                    f2 = redu_fits,
                    folds = folds,
                    type = type,
                    run_regression = FALSE,
                    alpha = 0.05), error = function(e) NA)
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
## get the risk estimate, CI based on averaging over the 10 random starts
get_avg_risk_ci <- function(vimp_lst) {
  ests_full <- unlist(lapply(vimp_lst, function(x) x$risk_full))
  ests_redu <- unlist(lapply(vimp_lst, function(x) x$risk_reduced))
  cis_full <- do.call(rbind, lapply(vimp_lst, function(x) x$risk_ci_full))
  cis_redu <- do.call(rbind, lapply(vimp_lst, function(x) x$risk_ci_reduced))
  est_full <- mean(ests_full)
  est_redu <- mean(ests_redu)
  ci_full <- colMeans(cis_full)
  ci_redu <- colMeans(cis_redu)
  return(list(risk_full = est_full, risk_reduced = est_redu,
              risk_ci_full = ci_full, risk_ci_reduced = ci_redu))
}

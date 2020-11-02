# create SL screens, algorithm + screen combinations

# -------------------------------------------------------------------------------------
# SL screens; all models adjust for baseline covariates age, BMI at enrollment, baseline behavioral risk
# -------------------------------------------------------------------------------------
# screen based on logistic regression univariate p-value < level
rank_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
    # logistic regression of outcome on each variable
    listp <- apply(X, 2, function(x, Y, family) {
      summ <- coef(summary(glm(Y ~ x + X$age + X$BMI + X$bhvrisk, family = family, weights = obsWeights)))
        ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
    }, Y = Y, family = family)
    # rank the p-values; give age, BMI, bhvrisk the lowest rank (will always set to TRUE anyways)
    ranked_vars <- rank(listp, ties = "average")
    ranked_vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- 999
    return(ranked_vars)
}
# screen dynamic range: only keep variables with 20th percentile != 80th percentile
screen_dynamic_range_plus_exposure <- function(Y, X, family, obsWeights, id, nVar = 4, ...) {
  # set all vars to false
  vars <- rep(FALSE, ncol(X))

  # keep only those with dynamic range: 20th percentile != 80th percentile
  x_quantiles <- apply(X, 2, function(x) quantile(x, probs = c(0.2, 0.8)))
  vars <- apply(x_quantiles, 2, function(x) round(x[1], 4) != round(x[2], 4))
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>%
    select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE

  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}
# screen dynamic range score: only keep variables with sd(vacinees)/sd(placebo) > 75th percentile
# relies on having var.super loaded in the environment
screen_dynamic_range_score_plus_exposure <- function(Y, X, family, obsWeights, id, var_super = var.super, nVar = 4, ...) {
  # set all to false
  vars <- rep(FALSE, ncol(X))
  # need to apply with the correct label in place of X
  vars_sd_ratio <- ifelse(is.na(var_super$sd.ratio), TRUE, var_super$sd.ratio > quantile(var_super$sd.ratio, probs = c(0.5), na.rm = TRUE))
  vars <- names(X) %in% var_super$varname[vars_sd_ratio] | names(X) %in% paste0(var_super$varname[vars_sd_ratio], "_bin")
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>%
    select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}
# screen based on lasso
screen_glmnet_plus_exposure <- function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, nVar = 4, ...) {
  vars <- screen.glmnet(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...)
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>%
    select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}
# screen based on logistic regression univariate p-value < level
screen_univariate_logistic_pval_plus_exposure <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, nVar = 4, ...) {
    # logistic regression of outcome on each variable
    listp <- apply(X, 2, function(x, Y, family) {
      summ <- coef(summary(glm(Y ~ x + X$age + X$BMI + X$bhvrisk, family = family, weights = obsWeights)))
        ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
    }, Y = Y, family = family)
    vars <- (listp <= minPvalue)
    # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
    vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
    if (sum(vars) < minscreen) {
        warning("number of variables with p value less than minPvalue is less than minscreen")
        vars[rank(listp) <= minscreen] <- TRUE
    }
    # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
    X_initial_screen <- X %>%
      select(names(X)[vars], "age", "BMI", "bhvrisk")
    ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
    vars[vars][ranked_vars > nVar] <- FALSE
    vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
    return(vars)
}
screen_univariate_logistic_pval_plus_exposure_0.01 <- function(Y, X, family, obsWeights, id, minPvalue = 0.01, minscreen = 2, ...) {
  screen_univariate_logistic_pval_plus_exposure(Y, X, family, obsWeights, id, minPvalue, minscreen = 2, ...)
}
screen_univariate_logistic_pval_plus_exposure_0.05 <- function(Y, X, family, obsWeights, id, minPvalue = 0.05, minscreen = 2, ...) {
  screen_univariate_logistic_pval_plus_exposure(Y, X, family, obsWeights, id, minPvalue, minscreen = 2, ...)
}
screen_univariate_logistic_pval_plus_exposure_0.1 <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, ...) {
  screen_univariate_logistic_pval_plus_exposure(Y, X, family, obsWeights, id, minPvalue, minscreen = 2, ...)
}
screen_highcor_plus_exposure <- function(Y, X, family, obsWeights, id, nVar = 4, ...) {
  # set all vars to FALSE
  vars <- rep(FALSE, ncol(X))
  # compute pairwise correlations between all marker vars
  cors <- cor(X, method = "spearman")
  diag(cors) <- NA
  cor_less_0.9 <- (cors <= 0.9)
  # screen out those with r > 0.9
  vars <- apply(cor_less_0.9, 1, function(x) all(x, na.rm = TRUE))
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>%
    select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  # make sure that age, BMI, bhvrisk are true
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# screen to always include the various variable sets
screen_assay_plus_exposure <- function(Y, X, family, obsWeights, id, assays, nVar = 4, ...) {
  # set all vars to be false
  vars <- rep(FALSE, ncol(X))
  # set vars with assay in name to be true
  # may be more than one
  for (i in 1:length(assays)) {
    vars[grepl(assays[i], names(X))] <- TRUE
  }
  # also keep the first three columns of X (correspond to age, BMI, bhvrisk)
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  # keep only a max of nVar immune markers; rank by univariate p-value in a model adjusting for age, BMI, bhvrisk
  X_initial_screen <- X %>%
    select(names(X)[vars], "age", "BMI", "bhvrisk")
  ranked_vars <- rank_univariate_logistic_pval_plus_exposure(Y, X_initial_screen, family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  # set baseline exposure vars to true
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}

# actual screens for antigen combos
# baseline_exposure is the base model
screen_baseline_exposure <- function(Y, X, family, obsWeights, id, ...) {
  # set all vars to false
  vars <- rep(FALSE, ncol(X))
  # set baseline exposure vars to true
  vars[names(X) %in% c("age", "BMI", "bhvrisk")] <- TRUE
  return(vars)
}
screen_igg_iga_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X, family, obsWeights, id, assays = c("IgG", "IgA"))
}
screen_tcells_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X, family, obsWeights, id, assays = c("CD4", "CD8"))
}
screen_fxab_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X, family, obsWeights, id, assays = c("IgG3", "phago", "fcrR2a", "fcrR3a"))
}
screen_igg_iga_tcells_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X,family,  obsWeights, id, assays = c("IgG", "IgA", "CD4", "CD8"))
}
screen_igg_iga_fxab_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X, family, obsWeights, id, assays = c("IgG", "IgA", "IgG3", "phago", "fcrR2a", "fcrR3a"))
}
screen_tcells_fxab_plus_exposure <- function(Y, X, family, obsWeights, id, ...) {
  screen_assay_plus_exposure(Y, X, family, obsWeights, id, assays = c("CD4", "CD8", "IgG3", "phago", "fcrR2a", "fcrR3a"))
}

screens_with_assay_groups <- c("screen_glmnet_plus_exposure", paste0("screen_univariate_logistic_pval_plus_exposure_", c(0.01, 0.05, 0.1)),
             "screen_dynamic_range_plus_exposure", "screen_dynamic_range_score_plus_exposure",
             "screen_baseline_exposure", paste0("screen_", c("igg_iga", "tcells", "fxab", "igg_iga_tcells",
                                                             "igg_iga_fxab", "tcells_fxab"),
                                                "_plus_exposure"),
             "screen_highcor_plus_exposure")
screens <- c("screen_glmnet_plus_exposure", paste0("screen_univariate_logistic_pval_plus_exposure_", c(0.01, 0.05, 0.1)),
                  "screen_dynamic_range_plus_exposure", "screen_dynamic_range_score_plus_exposure",
                  "screen_highcor_plus_exposure")
# -------------------------------------------------------------------------------------
# SL algorithms
# -------------------------------------------------------------------------------------

# --------------------------------------------------------------------------
# define wrappers that are less memory-intensive than the usual SL functions
# --------------------------------------------------------------------------
# skinny glm
SL.glm.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.glm.fit <- SL.glm(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

# skinny glm with interactions
SL.glm.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.glm.fit <- SL.glm.interaction(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, ...)
  SL.glm.fit$fit$object$y <- NULL
  SL.glm.fit$fit$object$model <- NULL
  SL.glm.fit$fit$object$residuals <- NULL
  SL.glm.fit$fit$object$fitted.values <- NULL
  SL.glm.fit$fit$object$effects <- NULL
  SL.glm.fit$fit$object$qr$qr <- NULL
  SL.glm.fit$fit$object$linear.predictors <- NULL
  SL.glm.fit$fit$object$weights <- NULL
  SL.glm.fit$fit$object$prior.weights <- NULL
  SL.glm.fit$fit$object$data <- NULL
  SL.glm.fit$fit$object$family$variance <- NULL
  SL.glm.fit$fit$object$family$dev.resids <- NULL
  SL.glm.fit$fit$object$family$aic <- NULL
  SL.glm.fit$fit$object$family$validmu <- NULL
  SL.glm.fit$fit$object$family$simulate <- NULL
  attr(SL.glm.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.glm.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.glm.fit)
}

# skinny stepwise with interactions
SL.step.interaction.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.step.interaction.fit <- SL.step.interaction(Y = Y, X = X, newX = newX, family = family,
                                                 obsWeights = obsWeights, direction = "forward", ...)
  SL.step.interaction.fit$fit$object$y <- NULL
  SL.step.interaction.fit$fit$object$model <- NULL
  SL.step.interaction.fit$fit$object$residuals <- NULL
  SL.step.interaction.fit$fit$object$fitted.values <- NULL
  SL.step.interaction.fit$fit$object$effects <- NULL
  SL.step.interaction.fit$fit$object$qr$qr <- NULL
  SL.step.interaction.fit$fit$object$linear.predictors <- NULL
  SL.step.interaction.fit$fit$object$weights <- NULL
  SL.step.interaction.fit$fit$object$prior.weights <- NULL
  SL.step.interaction.fit$fit$object$data <- NULL
  SL.step.interaction.fit$fit$object$family$variance <- NULL
  SL.step.interaction.fit$fit$object$family$dev.resids <- NULL
  SL.step.interaction.fit$fit$object$family$aic <- NULL
  SL.step.interaction.fit$fit$object$family$validmu <- NULL
  SL.step.interaction.fit$fit$object$family$simulate <- NULL
  attr(SL.step.interaction.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.interaction.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.interaction.fit)
}

# skinny stepwise (forward)
SL.step.skinny <- function(Y, X, newX, family, obsWeights, ...){
  SL.step.fit <- SL.step(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, direction = "forward", ...)
  SL.step.fit$fit$object$y <- NULL
  SL.step.fit$fit$object$model <- NULL
  SL.step.fit$fit$object$residuals <- NULL
  SL.step.fit$fit$object$fitted.values <- NULL
  SL.step.fit$fit$object$effects <- NULL
  SL.step.fit$fit$object$qr$qr <- NULL
  SL.step.fit$fit$object$linear.predictors <- NULL
  SL.step.fit$fit$object$weights <- NULL
  SL.step.fit$fit$object$prior.weights <- NULL
  SL.step.fit$fit$object$data <- NULL
  SL.step.fit$fit$object$family$variance <- NULL
  SL.step.fit$fit$object$family$dev.resids <- NULL
  SL.step.fit$fit$object$family$aic <- NULL
  SL.step.fit$fit$object$family$validmu <- NULL
  SL.step.fit$fit$object$family$simulate <- NULL
  attr(SL.step.fit$fit$object$terms, ".Environment") <- NULL
  attr(SL.step.fit$fit$object$formula, ".Environment") <- NULL
  return(SL.step.fit)
}

my_SL.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000,
                           max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                           nthread = 1, verbose = 0, save_period = NULL, ...) {
  if (packageVersion("xgboost") < 0.6)
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    if (packageVersion("xgboost") < "1.1.1.1") {
      objective <- "reg:linear"
    }
    else {
      objective <- "reg:squarederror"
    }
    model = xgboost::xgboost(data = xgmat, objective = objective,
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}
# boosted decision stumps
SL.stumpboost <- function(Y, X, newX, family, obsWeights, ...){
  fit <- my_SL.xgboost(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights,
                    max_depth = 1, # so it's only a stump
                    ...)
  return(fit)
}

# naive bayes wrapper
SL.naivebayes <- function(Y, X, newX, family, obsWeights, laplace = 0, ...){
  SuperLearner:::.SL.require("e1071")
  if(family$family == "gaussian"){
    stop("SL.naivebayes only works with binary outcomes")
  }else{
    nb <- naiveBayes(y = Y, x = X, laplace = laplace)
    pred <- predict(nb, newX, type = "raw")[,2]
    out <- list(fit = list(object = nb), pred = pred)
    class(out$fit) <- "SL.naivebayes"
    return(out)
  }
}

# predict method for naive bayes wrapper
predict.SL.naivebayes <- function(object, newdata, ...){
  pred <- predict(object$object, newdata = newdata, type = "raw")[,2]
  return(pred)
}

# methods <- c("SL.glm.skinny", "SL.glm.interaction.skinny", "SL.step.interaction.skinny",
             # "SL.naivebayes", "SL.stumpboost", "SL.glmnet", "SL.earth")
methods <- c("SL.glm.skinny", "SL.glm.interaction.skinny", "SL.step.interaction.skinny",
             "SL.stumpboost", "SL.glmnet", "SL.earth")




# -------------------------------------------------------------------------------------
# Add all alg/screen combinations to global environment, create SL library
# -------------------------------------------------------------------------------------

#' This function takes a super learner method wrapper and a super learner
#' screen wrapper and combines them into a single wrapper and makes that
#' wrapper available in the specified environment. It also makes a predict
#' method available in the specified environment.
#' @param method A super learner method wrapper. See ?SuperLearner::listWrappers(what = "method").
#' @param screen A super learner method wrapper. See ?SuperLearner::listWrappers(what = "screen").
#' @param envir The environment to assign the functions to (default is global environment)
#' @param verbose Print a message with the function names confirming their assignment?
assign_combined_function <- function(method, screen, envir = .GlobalEnv,
                                     verbose = TRUE){
  fn <- eval(parse(text =
          paste0("function(Y, X, newX, obsWeights, family, ...){ \n",
                    "screen_call <- ", screen, "(Y = Y, X = X, newX = newX, obsWeights = obsWeights, family = family, ...) \n",
                    "method_call <- ", method, "(Y = Y, X = X[,screen_call,drop=FALSE], newX = newX[,screen_call,drop = FALSE], obsWeights = obsWeights, family = family, ...) \n",
                    "pred <- method_call$pred \n",
                    "fit <- list(object = method_call$fit$object, which_vars = screen_call) \n",
                    "class(fit) <- paste0('", screen, "', '_', '", method, "') \n",
                    "out <- list(fit = fit, pred = pred) \n",
                    "return(out) \n",
                    "}")))
  fn_name <- paste0(screen,"_",method)
  assign(x = fn_name, value = fn, envir = envir)
  if(verbose){
    message(paste0("Function ", fn_name, " now available in requested environment."))
  }
  if (method == "SL.glmnet") {
      pred_fn <- eval(parse(text =
              paste0("function(object, newdata, ...){ \n",
                        "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                        "pred <- predict(object$object, type = 'response', newx = as.matrix(screen_newdata), s = 'lambda.min', ...) \n",
                        "return(pred) \n",
                     "}")))
    } else if (method == "SL.stumpboost") {
      pred_fn <- eval(parse(text =
              paste0("function(object, newdata, ...){ \n",
                        "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                        "screen_newdata_2 <- matrix(unlist(lapply(screen_newdata, as.numeric)), nrow=nrow(screen_newdata), ncol=ncol(screen_newdata)) \n",
                        "pred <- predict(object$object, newdata = screen_newdata_2, ...) \n",
                        "return(pred) \n",
                     "}")))
    } else if (method == "SL.naivebayes") {
      pred_fn <- eval(parse(text =
              paste0("function(object, newdata, ...){ \n",
                        "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                        'pred <- predict(object$object, newdata = screen_newdata, type = "raw", ...)[,2] \n',
                        "return(pred) \n",
                     "}")))
    } else if (method == "SL.randomForest") {
      pred_fn <- eval(parse(text =
              paste0("function(object, newdata, ...){ \n",
              "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
              "if (object$object$type != 'classification') {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'response')
                                }else {
                                    pred <- predict(object$object, newdata = screen_newdata, type = 'vote')[,
                                        2]
                                }
                                pred",
                     "}")))
    }else {
      pred_fn <- eval(parse(text =
              paste0("function(object, newdata, ...){ \n",
                        "screen_newdata <- newdata[,object$which_vars,drop = FALSE] \n",
                        "pred <- predict(object$object, type = 'response', newdata = screen_newdata, ...) \n",
                        "return(pred) \n",
                     "}")))
    }

  pred_fn_name <- paste0("predict.",screen,"_",method)
  assign(x = pred_fn_name, value = pred_fn, envir = envir)
  if(verbose){
    message(paste0("Function ", pred_fn_name, " now available in requested environment."))
  }
}

# make a data frame of all method/screen combinations needed
screen_method_frame_with_assay_groups <- expand.grid(screen = screens_with_assay_groups, method = methods)
screen_method_frame <- expand.grid(screen = screens, method = methods)

# add to global environment
apply(screen_method_frame_with_assay_groups, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})
apply(screen_method_frame, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})

# create SL library; reference method is glm with baseline exposure vars
SL_library_with_assay_groups <- c(apply(screen_method_frame_with_assay_groups, 1, paste0, collapse = "_"))
SL_library <- c(apply(screen_method_frame, 1, paste0, collapse = "_"))

## create SL screens, algorithm + screen combinations

## -------------------------------------------------------------------------------------
## SL screens; all models adjust for baseline covariates age, BMI at enrollment, baseline behavioral risk
## -------------------------------------------------------------------------------------
## screen dynamic range: only keep variables with 20th percentile != 80th percentile
screen_dynamic_range <- function(Y, X, family, obsWeights, id, ...) {
    # keep only those with dynamic range: 20th percentile != 80th percentile
    x_quantiles <- apply(X, 2, function(x) quantile(x, probs = c(0.2, 0.8)))
    vars <- apply(x_quantiles, 2, function(x) round(x[1], 4) != round(x[2], 4))
    return(vars)
}
## screen dynamic range score: only keep variables with sd(vacinees)/sd(placebo) > 75th percentile
screen_dynamic_range_score <- function(Y, X, family, obsWeights, id, ...) {
  # get dynamic range score; marker is the marker of interest, vaccine is true if vaccinee, false if placebo
  get_dynamic_range_score <- function(marker, vaccine) {
    return(sd(marker[vaccine])/sd(marker[!vaccine]))
  }
  # need to apply with the correct label in place of X
  dynamic_range_scores <- apply(X, 2, function(x) get_dynamic_range_score(x, X))
  vars <- dynamic_range_scores > quantile(dynamic_range_scores, probs = c(0.5))
}
## screen based on lasso 
screen_glmnet <- function(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...) {
    vars <- screen.glmnet(Y, X, family, obsWeights, id, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...)
}
## screen based on logistic regression univariate p-value < level
screen_univariate_logistic_pval <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, ...) {
    ## logistic regression of outcome on each variable
    listp <- apply(X, 2, function(x, Y, family) {
      summ <- coef(summary(glm(Y ~ x, family = family)))
        ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
    }, Y = Y, family = family)
    whichVariable <- (listp <= minPvalue)
    if (sum(whichVariable) < minscreen) {
        warning("number of variables with p value less than minPvalue is less than minscreen")
        whichVariable[rank(listp) <= minscreen] <- TRUE
    }
    return(whichVariable)
}

screens <- c("screen_glmnet", "screen_univariate_logistic_pval",
             "screen_dynamic_range", "screen_dynamic_range_score")

## -------------------------------------------------------------------------------------
## SL algorithms
## -------------------------------------------------------------------------------------

## --------------------------------------------------------------------------
## define wrappers that are less memory-intensive than the usual SL functions
## --------------------------------------------------------------------------
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

## skinny glm with interactions
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

# boosted decision stumps
SL.stumpboost <- function(Y, X, newX, family, obsWeights, ...){
  fit <- SL.xgboost(Y = Y, X = X, newX = newX, family = family, obsWeights = obsWeights, 
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

methods <- c("SL.glm.skinny", "SL.glm.interaction.skinny", "SL.step.interaction.skinny",
             "SL.naivebayes", "SL.stumpboost", "SL.glmnet", "SL.earth")




## -------------------------------------------------------------------------------------
## Add all alg/screen combinations to global environment, create SL library
## -------------------------------------------------------------------------------------

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

## make a data frame of all method/screen combinations needed
screen_method_frame <- expand.grid(screen = screens, method = methods)

## add to global environment
apply(screen_method_frame, 1, function(x) {assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)})

## create SL library
SL_library <- c(apply(screen_method_frame, 1, paste0, collapse = "_"),
"SL.mean")
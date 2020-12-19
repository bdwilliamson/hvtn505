# noise screens -- these don't include baseline exposure variables
noise_rank_univariate_logistic_pval <- function(Y, X, family, 
                                                        obsWeights, id, ...) {
  # logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x, 
                             family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  # rank the p-values
  ranked_vars <- rank(listp, ties = "average")
  return(ranked_vars)
}
# screen dynamic range: only keep variables with 
# 20th percentile != 80th percentile
noise_screen_dynamic_range <- function(Y, X, family, 
                                               obsWeights, id, 
                                               nVar = 4, ...) {
  # set all vars to false
  vars <- rep(FALSE, ncol(X))
  # keep only those with dynamic range: 20th percentile != 80th percentile
  x_quantiles <- apply(X, 2, function(x) quantile(x, probs = c(0.2, 0.8)))
  vars <- apply(x_quantiles, 2, function(x) round(x[1], 4) != round(x[2], 4))
  # keep only a max of nVar immune markers; rank by univariate p-value 
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- noise_rank_univariate_logistic_pval(Y, X_initial_screen, 
                                                     family, obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}
# screen based on lasso
noise_screen_glmnet <- function(Y, X, family, 
                                        obsWeights, id, alpha = 1, 
                                        minscreen = 2, nfolds = 10, 
                                        nlambda = 100, nVar = 4, ...) {
  vars <- screen.glmnet(Y, X, family, obsWeights, id, alpha = 1, 
                        minscreen = 2, nfolds = 10, nlambda = 100, ...)
  # keep only a max of nVar immune markers; rank by univariate p-value 
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- noise_rank_univariate_logistic_pval(Y,  X_initial_screen, 
                                                     family, 
                                                     obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}
# screen based on logistic regression univariate p-value < level
noise_screen_univariate_logistic_pval <- function(Y, X, family, 
                                                          obsWeights, id, 
                                                          minPvalue = 0.1, 
                                                          minscreen = 2, 
                                                          nVar = 4, ...) {
  # logistic regression of outcome on each variable
  listp <- apply(X, 2, function(x, Y, family) {
    summ <- coef(summary(glm(Y ~ x, 
                             family = family, weights = obsWeights)))
    ifelse(dim(summ)[1] > 1, summ[2, 4], 1)
  }, Y = Y, family = family)
  vars <- (listp <= minPvalue)
  if (sum(vars) < minscreen) {
    warning("number of variables with p value less than minPvalue is 
                less than minscreen")
    vars[rank(listp) <= minscreen] <- TRUE
  }
  # keep only a max of nVar immune markers; rank by univariate p-value 
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- noise_rank_univariate_logistic_pval(Y,  X_initial_screen, 
                                                     family, 
                                                     obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}
noise_screen_univariate_logistic_pval_0.01 <- function(Y, X, family, 
                                                               obsWeights, id, 
                                                               minPvalue = 0.01, 
                                                               minscreen = 2, 
                                                               ...) {
  noise_screen_univariate_logistic_pval(Y, X, family, obsWeights, 
                                                id, minPvalue, minscreen = 2, 
                                                ...)
}
noise_screen_univariate_logistic_pval_0.05 <- function(Y, X, family, 
                                                               obsWeights, id, 
                                                               minPvalue = 0.05, 
                                                               minscreen = 2, 
                                                               ...) {
  noise_screen_univariate_logistic_pval(Y, X, family, obsWeights, 
                                                id, minPvalue, 
                                                minscreen = 2, ...)
}
noise_screen_univariate_logistic_pval_0.1 <- function(Y, X, family, 
                                                              obsWeights, id, 
                                                              minPvalue = 0.1, 
                                                              minscreen = 2, 
                                                              ...) {
  noise_screen_univariate_logistic_pval(Y, X, family, obsWeights, 
                                                id, minPvalue, 
                                                minscreen = 2, ...)
}
noise_screen_highcor <- function(Y, X, family, obsWeights, id, 
                                         nVar = 4, ...) {
  # set all vars to FALSE
  vars <- rep(FALSE, ncol(X))
  # compute pairwise correlations between all marker vars
  cors <- cor(X, method = "spearman")
  diag(cors) <- NA
  cor_less_0.9 <- (cors <= 0.9)
  # screen out those with r > 0.9
  vars <- apply(cor_less_0.9, 1, function(x) all(x, na.rm = TRUE))
  # keep only a max of nVar immune markers; rank by univariate p-value 
  X_initial_screen <- X %>%
    select(names(X)[vars])
  ranked_vars <- noise_rank_univariate_logistic_pval(Y,  X_initial_screen, 
                                                     family, 
                                                     obsWeights, id)
  vars[vars][ranked_vars > nVar] <- FALSE
  return(vars)
}

noise_screens <- c(
  "noise_screen_glmnet", 
  paste0("noise_screen_univariate_logistic_pval_", 
         c(0.01, 0.05, 0.1)), 
  "noise_screen_dynamic_range", 
  "noise_screen_highcor"
)
noise_screen_method_frame <- expand.grid(screen = noise_screens, 
                                         method = methods)
apply(noise_screen_method_frame, 1, 
      function(x) {
        assign_combined_function(screen = x[1], method = x[2], verbose = FALSE)
      })
SL_library_noise <- c(apply(noise_screen_method_frame, 1, 
                            paste0, collapse = "_"))
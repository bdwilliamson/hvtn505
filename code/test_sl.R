## test only baseline exposure
system.time(test_sl_baseline_only <- run_cv_sl_once(seed = seeds[1], Y = Y_vaccine, X_mat = X_vaccine, 
                                                    family = "binomial",
                                                    obsWeights = weights_vaccine,
                                                    sl_lib = SL_library[7],
                                                    method = "method.CC_nloglik",
                                                    cvControl = list(V = V_outer),
                                                    innerCvControl = list(list(V = V_inner)),
                                                    vimp = FALSE))

## test all with glm
SL_library_test <- SL_library[grepl("glm.skinny", SL_library)]
system.time(test_sl <- run_cv_sl_once(seed = seeds[1], Y = Y_vaccine, X_mat = X_vaccine, family = "binomial",
                                      obsWeights = weights_vaccine,
                                      sl_lib = SL_library_test, # this comes from sl_screens.R
                                      method = "method.CC_nloglik",
                                      cvControl = list(V = V_outer),
                                      innerCvControl = list(list(V = V_inner)),
                                      vimp = FALSE))

## test only baseline exposure with parallelization
system.time(test_fits_baseline_only <- parallel::mclapply(seeds, FUN = run_cv_sl_once, Y = Y_vaccine, X_mat = X_vaccine, family = "binomial",
                                                          obsWeights = weights_vaccine,
                                                          sl_lib = SL_library[7], 
                                                          method = "method.CC_nloglik",
                                                          cvControl = list(V = V_outer, stratifyCV = TRUE),
                                                          innerCvControl = list(list(V = V_inner)),
                                                          vimp = FALSE,
                                                          mc.cores = num_cores
))

## test with CV on glm
set.seed(seeds[1])
train <- vaccinees %>% 
  group_by(Y) %>% 
  sample_frac(size = 0.8) %>% 
  ungroup()
Y_train <- train$Y
weights_train <- train$weights
X_train <- train %>% 
  select(-Y, -weights)

set.seed(seeds[1])
folds <- caret::createFolds(y = factor(Y_vaccine), k = 5)

## fit glm
glm_coefs_lst <- vector("list", 5)
for (v in 1:5) {
  data_set <- cbind.data.frame(y = Y_vaccine[-unlist(folds[v])], weight = weights_vaccine[-unlist(folds[v])], X_vaccine[-unlist(folds[v]), ])
  glm_coefs_lst[[v]] <- matrix(NA, nrow = 4, ncol = dim(data_set)[1])
  for (i in 1:dim(data_set)[1]) {
    cat("\n v = ", v, "; i = ", i)
    mod <- glm(y ~ age + BMI + bhvrisk, family = "binomial", weights = weight, data = data_set[-i, ])
    glm_coefs_lst[[v]][, i] <- coef(mod)
  }
}
rowMeans(glm_coefs_lst[[1]])
rowMeans(glm_coefs_lst[[2]])
rowMeans(glm_coefs_lst[[3]])
rowMeans(glm_coefs_lst[[4]])
rowMeans(glm_coefs_lst[[5]])

## run it once
system.time(test_small <- run_cv_sl_once(seeds[1], Y = Y_vaccine, X_mat = X_vaccine,
                                         family = "binomial",
                                         obsWeights = weights_vaccine,
                                         sl_lib = SL_library[7], method = "method.CC_nloglik",
                                         cvControl = list(V = V_outer, stratifyCV = TRUE),
                                         innerCvControl = list(list(V = V_inner)), vimp = FALSE))

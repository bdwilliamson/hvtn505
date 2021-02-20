# run the simulation a single time
#   testing the "size" of AIPW vs IPW

run_sim_once <- function(iteration = 1, n = 100, xdim = 50, zdim = 2, V = 5,
SL_V = 5, learner_lib = c("SL.glmnet", "SL.ranger")){
    # generate data
    data_lst <- gen_data(n = n, p = xdim, q = zdim)
    obs_dat <- data_lst$observed_data
    obs_y <- obs_dat$y
    obs_delta <- obs_dat$delta
    obs_x <- obs_dat %>%
        select(-y, -delta)
    names(obs_x) <- paste0("X", 1:ncol(obs_x))
    # note weights are independent of X
    weights <- rep(.25 ^ (-1), length(obs_delta)) 
    # weight_glm <- glm(obs_delta ~ ., family = "binomial", 
    #                   data = obs_x[, 1:(zdim + 1)])
    # weights <- predict(weight_glm, type = "response") ^ (-1)
    # run cross-validated variable importance
    vim_est <- vimp::cv_vim(
        Y = obs_y,
        X = obs_x,
        V = V, stratified = TRUE,
        type = "auc",
        run_regression = TRUE,
        SL.library = learner_lib,
        C = obs_delta,
        Z = c("Y", paste0("X", 1:(zdim + 1))),
        ipc_weights = weights,
        cvControl = list(V = SL_V, stratifyCV = TRUE)
    )
    # return the relevant objects
    output_tib <- tibble::tibble(
        mc_id = iteration, n = n
    ) %>%
        bind_cols(
            tibble::as_tibble(vim_est$mat %>%
            select(-s))
        )
    output_tib
}

# run the simulation a single time
#   testing the "size" of AIPW vs IPW

run_sim_once <- function(iteration = 1, n = 100, xdim = 50, zdim = 2, V = 5, 
                         SL_V = 5, ipc_est_type = "aipw", 
                         learner_lib = c("SL.glmnet", "SL.ranger")){
    # generate data
    data_lst <- gen_data(n = n, p = xdim, q = zdim)
    obs_dat <- data_lst$observed_data
    obs_y <- obs_dat$y
    obs_delta <- obs_dat$delta
    obs_x <- obs_dat %>%
        select(-y, -delta)
    names(obs_x) <- paste0("X", 1:ncol(obs_x))
    # note weights are independent of X
    # weights <- rep(.25 ^ (-1), length(obs_delta)) 
    obs_z <- obs_x[, 1:(zdim + 1)]
    weight_mod <- glm(obs_delta ~ ., family = "binomial", data = obs_z)
    weights <- predict(weight_mod, type = "response") ^ (-1)
    # estimate cross-validated AUC
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
        ipc_est_type = ipc_est_type,
        cvControl = list(V = SL_V, stratifyCV = TRUE)
    )
    pred_se <- vim_est$predictiveness_ci_full[, 2] / qnorm(0.975) - vim_est$est
    p_val <- 1 - pnorm((vim_est$predictiveness_full - 0.5) / sqrt(pred_se))
    pred_tib <- tibble::tibble(
        mc_id = iteration, n = n, est = vim_est$predictiveness_full,
        cil = vim_est$predictiveness_ci_full[, 1], 
        ciu = vim_est$predictiveness_ci_full[, 2],
        test = p_val < 0.05,
        p_value = p_val
    )
    pred_tib
}

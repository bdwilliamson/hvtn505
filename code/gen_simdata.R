# generate data for the "size test" simulation

gen_data <- function(n = 100, p = 50, q = 2) {
    x <- replicate(p, rnorm(n, 0, 1))
    t <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, 0.5)

    # inclusion into subsample
    samp <- rbinom(n, 1, 0.4)

    # ideal dataset
    data_star <- data.frame(y = y, cbind(t, x))
    # observed dataset:
    #   the first two variables are observed on all participants
    x_with_nas <- x
    x_with_nas[samp == 1, (q + 1):p] <- NA
    obs_data <- data.frame(y = y, delta = 1 - samp, cbind(t, x_with_nas))
    list(ideal_data = data_star, observed_data = obs_data)
}

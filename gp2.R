

covmat = function(N, sigma, ell)
{
    # Make a covariance matrix
    C = matrix(nrow=N, ncol=N)
    C = sigma^2*exp(-0.5*(data$t[row(C)] - data$t[col(C)])^2 / ell^2)

    # For numerical stability!
    C = C + 1E-4*diag(N)
    C
}


log_likelihood = function(params, data)
{
    mu = params["mu"]
    C = covmat(length(data$t), params["sigma"], params["ell"])
    y = data$y
    N = length(y)

    # Evaluate log likelihood! Naive calculation
#    logl0 = -0.5*N*log(2*pi) - 0.5*log(det(C))
#                    - 0.5*t(y - mu) %*% solve(C) %*% (y - mu)

    # Better method - specific to this application
    L = t(chol(C))
    logdet = 2*sum(log(diag(L)))
    logl = -0.5*N*log(2*pi) - 0.5*logdet
                - 0.5*t(y - mu) %*% solve(C, y - mu)
#    cat(logl0, logl)
#    cat("\n")

    logl
}

many_log_likelihoods = function(data, num=100)
{
    best = NULL
    best_logl = NULL
    for(i in 1:num)
    {
        params = c(runif(1), runif(1), 1000*runif(1))
        names(params) = c("mu", "sigma", "ell")
        logl = log_likelihood(params, data)
        if(i == 1 || logl > best_logl)
        {
            best = params
            best_logl = logl
        }
    }
    print(best)
    print(best_logl)
    list(best=best, best_logl=best_logl)
}

data = source("data.R")$value


#C = covmat(length(data$t), 0.2, 200)

## Flip the matrix (up <-> down) before making image
#image(apply(C, 2, rev))

time = system.time(many_log_likelihoods(data, 1000))
print(time)


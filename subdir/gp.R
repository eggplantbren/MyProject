

N = 1000

# Make a covariance matrix
C = matrix(nrow=N, ncol=N)

C = exp(-0.5*(row(C) - col(C))^2 / 100^2)

# For numerical stability!
C = C + 1E-4*diag(N)

#pdf("covariance_matrix.pdf")
#image(C)
#dev.off()

U = chol(C)
n = rnorm(1000)
y = t(U) %*% n

plot(y, type="l")


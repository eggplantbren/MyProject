

N = 1000

# Make a covariance matrix
C = matrix(nrow=N, ncol=N)

C = exp(-0.5*(row(C) - col(C))^2 / 100^2)

pdf("covariance_matrix.pdf")
image(C)
dev.off()



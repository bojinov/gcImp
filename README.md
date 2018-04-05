# gImp
An R package for generating multiple imputations using a Gaussian copula
This package generates multiple imputations using Peter Hoff's
sbgcop package. The package can also convert the output to allow for 
seamless integration with the mice package, which can be used to analyze the 
multiple imputed data sets.

# Install using devtools
library(devtools)
devtools::install_github("bojinov/gcImp")
Here is an exmaple: 
# Generate N samples from a multivariate normal
N <- 200
rho <- matrix(0.3, 2, 2)
diag(rho) <- 1
# Compute the Choleski decomposition of rho. 
rho.chol <- chol(rho)
samples.mvn <- matrix(rnorm(N * 2), ncol = 2) %*% rho.chol
# Delete some of the values
p <- rep(0.2, N) # MCAR
# p <- 1/(1 + exp(0.2 - samples.mvn[, 1])) #MAR
R <- sapply(1:N, function(jj) sample(c(NA, 1), size = 1,
                                     prob = c(p[jj], 1 - p[jj])))
# Generate the observed data
samples.mvn[, 2] <- samples.mvn[, 2] * R
out <- gcImp(samples.mvn)
print(out)
# Check the trace plot of the latent correlation 
plot(out$sbgcop.out$C.psamp[1,2,], type = "l")
# Check the trace plot of the mean of the imputed data
plot(colMeans(out$sbgcop.out$Y.imput[,2,]), type = "l")

## Not run:
# Further checks can be performed using the mcmcplots package
library(mcmcplots)
mcmcplots::mcmcplot(colMeans(out$sbgcop.out$Y.imput[,2,]))
mcmcplots::mcmcplot(out$sbgcop.out$C.psamp[1,2,])
## End(Not run)
 
# Convert to a class "mira" and use MICE for analysis
imp <- gc.as.mids(out)
# Stack the imputations and run a linear regression
stacked <- mice::complete(imp, "long")
fit <- lm(V1 ~ V2, data = stacked)
coef(fit)
# Fir separate regressions and combine the output using Rubin's rules
fit <- with(imp, lm(V1 ~ V2))
est <- mice::pool(fit)
summary(est)
# Plot the imputed against the observed
mice::densityplot(imp)

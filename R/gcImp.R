#' Gaussian Copula Imputations
#' 
#' \code{gcImp} generates multiple imputations using Peter Hoff's sbgcop 
#' package. Remember to remove any ID variables from the data matrix input!
#' After generating the imputations, use the \code{gc.as.mire} to convert 
#' the imputations into a "mira" format. We can than use the mice package
#' to carry out all subsequent analysis. 
#' 
#' @param dt A data matrix with missing values. 
#' Remember to remove any ID variables!
#' @param m A positive integer indicating the number of imputations. 
#' Note that m < (nsamp - burn) otherwise there will not be enough samples to 
#' generate the imputations.
#' @param burn A positive integer indicating the burn in. Default value is 300.
#' @param nsamp A positive integer indicating the number of samples from the 
#' MCMC chain. Default value is 1000.
#' @param ... Other arguments to pass to sbgcop.mcmc. 
#' @return An object of class 'gcImp' containing the following components:
#' \item{dt}{The original data set.}
#' \item{m}{The number of imputations.}
#' \item{resp}{The response indicators: a binary matrix of the same dimensions
#' as dt taking a value 1 if the corresponding outcome is observed and 0 if it
#' is missing.}
#' \item{nmis}{The number of missing values per variable.}
#' \item{imp}{A list of length "m" containing the imputed values.}
#' \item{method}{The method used to generated imputations: "Gaussian Copula."}
#' \item{sbgcop.out}{The output from sbgcop.out.}
#' @export
#' @examples
 # Generate N samples from a multivariate normal distribution
#' N <- 200
#' rho <- matrix(0.3, 2, 2)
#' diag(rho) <- 1
#' # Compute the Choleski decomposition of rho. 
#' rho.chol <- chol(rho)
#' # Generate imputations
#' samples.mvn <- matrix(rnorm(N * 2), ncol = 2) %*% rho.chol
#' # Delete some of the values
#' p <- rep(0.2, N)
#' # p <- 1/(1 + exp(0.2 - samples.mvn[, 1]))
#' R <- sapply(1:N, function(jj) sample(c(NA, 1), size = 1,
#'                                      prob = c(p[jj], 1 - p[jj])))
#' # Generate the observed data
#' samples.mvn[, 2] <- samples.mvn[, 2] * R
#' out <- gcImp(samples.mvn)
#' print(out)
#' # Check the trace plot of the latent correlation 
#' plot(out$sbgcop.out$C.psamp[1,2,], type = "l")
#' # Check the trace plot of the mean of the imputed data
#' plot(colMeans(out$sbgcop.out$Y.imput[,2,]), type = "l")
#' \dontrun{# Further checks can be performed using the mcmcplots package
#' library(mcmcplots)
#' mcmcplots::mcmcplot(colMeans(out$sbgcop.out$Y.imput[,2,]))
#' mcmcplots::mcmcplot(out$sbgcop.out$C.psamp[1,2,])
#' } 
#' # Convert to a class "mira" and use MICE for analysis
#' imp <- gc.as.mids(out)
#' # Stack the imputations and run a linear regression
#' stacked <- mice::complete(imp, "long")
#' fit <- lm(V1 ~ V2, data = stacked)
#' coef(fit)
#' # Fir separate regressions and combine the output using Rubin's rules
#' fit <- with(imp, lm(V1 ~ V2))
#' est <- mice::pool(fit)
#' summary(est)
#' # For version of mice > 2.6 we can also plot using the lattice package
#' \dontrun{
#' lattice::densityplot(imp)
#' lattice::bwplot(imp)}
gcImp <- function(dt, m = 20, burn = 300, nsamp = 1000, ...){
  # Sample the posterior
  out <- sbgcop::sbgcop.mcmc(dt, nsamp = nsamp, ...)
  # Extract samples
  imp.data <- out$Y.imp[, , round(seq(burn, nsamp, length = m))]
  # Extract the response indicators
  resp <- 1 - is.na(dt)  
  # Create the output
  gcImp <- list(dt = dt, m = m, resp = resp, nmis = colSums(1 - resp),
                imp = imp.data, method = "Gaussian Copula", sbgcop.out = out,
                iteration = round((nsamp - burn) / m))
  class(gcImp) <- "gcImp"
  return(gcImp)
}

#' @export
print.gcImp <- function(gcImp, ...) {
  cat("Multiply imputed data set")
  cat("\nNumber of multiple imputations: ")
  cat(gcImp$m)
  cat("\nMissing cells per column:")
  print(gcImp$nmis)
  cat("\nImputation Method:")
  print(gcImp$method)
}

#' Convert Gaussian copula to mids
#' 
#' \code{gc.as.mids} converts a gcImp object into a "mira" format. 
#' The mira format can then be used by the mice package to analyze imputed data.
#' 
#' @param out A gcImp calss object, usually the output from \code{gcImp}
#' @return An object of class "mira," see ?mice::mira for help.
#' @export
gc.as.mids <- function(out) {
  imp.data <- 
      lapply(1:ncol(out$dt), 
             function(kk) {
               if (out$nmis[kk] != 0) {
                 mice.style <- out$imp[which(out$resp[, kk] == 0), kk, ]
                 mice.style <- data.frame(mice.style)
                 row.names(mice.style) <- 
                    as.character(which(out$resp[, kk] == 0))
                 colnames(mice.style) <- 
                    as.character(1:out$m)  
               } else {
                mice.style <- list(NULL)
               }
               return(mice.style)
             })
  for (kk in 1:ncol(out$dt)) {
    if (out$nmis[kk] == 0) {
      imp.data[kk] <- list(NULL)
    }
  }
  if (!is.element("data.frame", class(out$dt))) {
    out$data <- data.frame(out$dt)
    colnames(out$data) <- paste("V", as.character(1:ncol(out$data)), sep = "")
  }
  names(imp.data) <- colnames(out$data)
  names(out$nmis) <- colnames(out$data)
  out$imp <- imp.data
  out$call <- "sgbcop(dt)"
  out$method <- rep(out$method, ncol(out$data))
  names(out$method) <- colnames(out$data)
  out$predictorMatrix <- matrix(0, ncol(out$data), ncol(out$data))
  out$predictorMatrix[which(out$nmis != 0), ] <- 1
  diag(out$predictorMatrix) <- 0
  colnames(out$predictorMatrix) <- colnames(out$data)
  rownames(out$predictorMatrix) <- colnames(out$data)
  out$visitSequence <- as.character(1:ncol(out$data))
  out$form <- list(NULL)
  out$post <- list(NULL)
  out$seed <- NA
  out$lastSeedValue <- list(NULL)
  out$chainMean <- list(NULL)
  out$chainVar <- list(NULL)
  out$loggedEvents <- list(NULL)
  out$pad <- list(NULL)
  class(out) <- "mids"
  return(out)
}







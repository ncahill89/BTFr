#-----------------------------------------------------
#' Run MCMC diagnostics - rhat
#'
#' @param mcmc.array an mcmc array from \code{\link{run_core}}
#' @param pars.check the parameter names to be checked
#'
#' @return A message will appear about convergence
#' @export
#' @importFrom coda 'gelman.diag' 'effectiveSize' 'mcmc.list' 'as.mcmc'
gr_diag <- function(mcmc.array,
                    pars.check) {
  R <- rep(NA, length(pars.check))
  p <- 0
  for (parname in pars.check) {
    p <- p + 1 # index
    mcmc.array.temp <- mcmc.array[, , parname]
    mcmc <- coda::mcmc.list()

    for (chain in 1:dim(mcmc.array.temp)[2]) {
      mcmc[[chain]] <- coda::as.mcmc(mcmc.array.temp[, chain])
    }
    r <- coda::gelman.diag(mcmc, autoburnin = FALSE, transform = F)$psrf
    R[p] <- r[, "Point est."]
  }


  names(R) <- pars.check


  if (length(R[R > 1.1]) > 0) {
    cat(paste("Poor/no convergence for:", names(R[R > 1.1]), "(R = ", round(R[R > 1.1], 3), ")", "\n"))
  } else {
    cat(paste0("Rhat looks good, no  convergence issues indicated for checked parameters \n"))
  }

  if (length(R[R > 1.1]) > 0) {
    return(-1)
  } else {
    return(0)
  }
}


#' Run MCMC diagnostics - effective sample size
#'
#' @param mcmc.array an mcmc array from \code{\link{run_core}}
#' @param pars.check the parameter names to be checked
#'
#' @return A message will appear about convergence
#' @export
eff_size <- function(mcmc.array,
                     pars.check) {
  ESS <- rep(NA, length(pars.check))
  p <- 0
  for (parname in pars.check) {
    p <- p + 1 # index
    mcmc.array.temp <- mcmc.array[, , parname]
    mcmc <- mcmc.list()

    for (chain in 1:dim(mcmc.array.temp)[2]) {
      mcmc[[chain]] <- coda::as.mcmc(mcmc.array.temp[, chain])
    }
    es <- coda::effectiveSize(mcmc)
    ESS[p] <- es / (dim(mcmc.array)[1] * dim(mcmc.array)[2])
  }
  names(ESS) <- pars.check
  if (length(ESS[ESS < 0.1]) > 0) {
    cat(paste0("Effective sample size is less than 10% of total iterations for parameter:", " ", names(ESS[ESS < 0.1]), " ", "(", round(ESS[ESS < 0.1], 3) * 100, "%", ")", "\n"))
    cat(paste0("Additional thinning may be required! \n"))
  } else {
    cat(paste0("No apparent autocorrelation issues for checked parameters. \n"))
  }
}

#' Run MCMC diagnostics - monte carlo standard error (mcse)
#'
#' @param mcmc.array an mcmc array from \code{\link{run_core}}
#' @param pars.check the parameter names to be checked
#'
#' @return A message will appear about convergence
#' @export
mcse <- function(mcmc.array,
                 pars.check) {
  MCSE <- rep(NA, length(pars.check))
  p <- 0
  for (parname in pars.check) {
    p <- p + 1 # index
    mcmc.array.temp <- mcmc.array[, , parname]
    mcmc <- mcmc.list()

    for (chain in 1:dim(mcmc.array.temp)[2]) {
      mcmc[[chain]] <- coda::as.mcmc(mcmc.array.temp[, chain])
    }
    es <- coda::effectiveSize(mcmc)
    MCSE[p] <- (stats::sd(mcmc.array.temp) / es) / stats::sd(mcmc.array.temp)
  }
  names(MCSE) <- pars.check
  if (length(MCSE[MCSE > 0.1]) > 0) {
    cat(paste0("The Monte Carlo standard error is greater than 10% of the posterior standard deviation for parameter:", " ", names(MCSE[MCSE > 0.1]), " ", "(", round(MCSE[MCSE > 0.1], 3) * 100, "%", ")", "\n"))
    cat(paste0("Sampling error variation appears too large! \n"))
  } else {
    cat(paste0("The accuracy of the parameter estimation is adequate. \n"))
  }
}

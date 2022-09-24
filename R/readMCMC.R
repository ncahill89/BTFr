#--------------------------------------------------------------------------
#'  Read in JAGS objects
#'
#' @description Read in JAGS objects and constructs \code{mcmc.array}, which is saved to \code{output}.
#' @param ChainIDs Optional: specify which chains to include
#' (to use when you want to exclude a chain that crashed, or which has not finished yet).
#' @param n.samplestot.max Maximum number of posterior samples to save
#' @param use.informative.priors Set to TRUE if using informative priors
#' @param core.run Set to TRUE if running core data
#' @param temp_files temporary JAGS files
#'
#' @return Invisible
#' @export
#'
ConstructMCMCArray <- function(ChainIDs,
                               n.samplestot.max = 15000,
                               core.run = FALSE,
                               use.informative.priors = FALSE,
                               temp_files) {

  # now combine the JAGS files into one mcmc.array
  n.chains <- length(ChainIDs)
  if (n.chains == 1) {
    cat("You need at least two chains!\n")
    return()
  }

  chain <- ifelse(length(ChainIDs) == 1, ChainIDs, ChainIDs[1])
  temp.jags.file <- temp_files[chain]
  load(temp.jags.file)
  n.sim <- dim(mod_upd$BUGSoutput$sims.array)[1]
  n.par <- dim(mod_upd$BUGSoutput$sims.array)[3]
  mcmc.array <- array(NA, c(n.sim, n.chains, n.par))
  dimnames(mcmc.array) <- list(NULL, NULL, names(mod_upd$BUGSoutput$sims.array[
    1,
    1,
  ]))
  for (chain in 1:n.chains) {
    chain_saved <- ifelse(length(ChainIDs) == 1, 1, ChainIDs[chain])
    temp.jags.file <- temp_files[chain_saved]
    load(temp.jags.file)
    mcmc.array[1:n.sim, chain, ] <- mod_upd$BUGSoutput$sims.array[, 1, ]
  }

  if (n.sim > n.samplestot.max) {
    mcmc.array <- mcmc.array[seq(1, n.sample.max, length.out = n.sample.max), , ]
  }

  cat("Run complete!", "\n")

  return(mcmc.array)
}
#----------------------------------------------------------------------
# The End!

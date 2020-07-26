#' Run the core model
#'
#' @param obj An object of class \code{BTF} from \code{\link{run_modern}}
#' @param core_counts.csv A .csv file containing core species counts
#' @param dx Parameter for spacing of spline knots
#' @param ChainNums The number of MCMC chains
#' @param n.iter The number of MCMC iterations
#' @param n.burnin The number of burnin MCMC samples
#' @param n.thin The number of thinning
#' @param run.on.server Set to TRUE for running chains in parallel on a server
#' @param use.informative.priors Set to TRUE if prior information is available
#'
#' @return Output will be saved to \code{output.dir}
#' @export
#' @import R2jags rjags coda
#' @importFrom dplyr "select" "ends_with"
#' @importFrom tidyr "pivot_longer"
#' @examples
#' coremodel<-runcore(core_counts.csv = "core_counts.csv")
run_core<-function( obj,
                    core_counts = NULL,
                    prior_el = NULL,
                    dx=0.2,
                    ChainNums = seq(1,3),
                    n.iter=15000,
                    n.burnin=1000,
                    n.thin=7,
                    validation.run=FALSE,
                    fold=1,
                    use_uniform_prior = FALSE)
{


  dir.create(file.path(getwd(), "temp.JAGSobjects"),showWarnings = FALSE)

  # read in the core data
  if (!is.null(core_counts)) {
    core_dat <- core_counts
  } else core_dat <- BTF::core_counts


  if(!validation.run)
    {
  depth <- core_dat %>% dplyr::pull(Depth)
  core_data_sorted <- sort_core(core_dat = select(core_dat,-"Depth"),
                                species_names = obj$species_names)
  }

  if(validation.run)
  {
    core_data_sorted<-sort_core(core_dat = core_dat,
                                  species_names=obj$species_names)

  }

  #Check if data includes priors
  if(!is.null(prior_el))
  {
  use.informative.priors = TRUE
  prior_emin <- pmax(obj$elevation_min,prior_el$prior_lwr/100)
  prior_emax <- pmin(obj$elevation_max, prior_el$prior_upr/100)
  cat("Running with informative priors")
  }

  #Get other relevant info for the model
  if(validation.run)
  {
    set.seed(3847)
    K <- 10
    folds <- rep(1:K,ceiling(nrow(core_data_sorted$coredata_sorted)/K))
    folds <- folds[sample(1:nrow(core_data_sorted$coredata_sorted))]
    testsamps<-which(folds==fold)
    testsamps

  y0<-core_data_sorted$coredata_sorted[testsamps,]
  }
  if(!validation.run)
  {
    y0<-core_data_sorted$coredata_sorted

  }
  begin0<-obj$data$begin0
  el_mean<-mean(obj$data$x)
  N_count0<-apply(y0,1,sum)
  n0<-nrow(y0)
  m0<-ncol(y0)

  if(is.null(prior_el))
  {
    emin=rep(obj$elevation_min ,n0)
    emax=rep(obj$elevation_max ,n0)
  }

  if(!is.null(prior_el))
  {
    emin=prior_emin
    emax=prior_emax
  }


  #########For Splines
  # This creates the components for the basis function matrix
  xl<-obj$elevation_min
  xr<-obj$elevation_max
  deg=3
  dx<-obj$dx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  n_knots<-length(knots)
  D <- diff(diag(n_knots), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  K <- dim(D)[1]
  Dmat=1
  Delta.hk <- diff(diag(K), diff = Dmat) # difference matrix
  Deltacomb.kh <- t(Delta.hk)%*%solve(Delta.hk%*%t(Delta.hk))

  # Jags model data
  pars=c("p","x0")

  data=list(y=y0,
            n=n0,
            m=m0,
            N_count = N_count0,
            D=D,
            Deltacomb.kh=Deltacomb.kh,
            knots=knots,
            deg=deg,
            n_knots=n_knots,
            begin0=begin0,
            beta0.j= obj$beta0.j,
            beta0_sd = obj$beta0_sd,
            delta0.hj= obj$delta0.hj,
            delta0_sd = obj$delta0_sd,
            alpha0 = obj$alpha0,
            alpha0_sd = obj$alpha0_sd,
            sigma.z = obj$sigma.z,
            emin=emin,
            emax=emax,
            el_mean=el_mean,
            use_uniform_prior = use_uniform_prior)

   for (chainNum in ChainNums){
      cat(paste("Start chain ID ", chainNum), "\n")

      InternalRunCore(chainNum = chainNum,
                      jags_data = data,
                      jags_pars = pars,
                      n.burnin = n.burnin,
                      n.iter = n.iter,
                      n.thin = n.thin)

    }

  data[["depth"]] <- depth
  #Store MCMC output in an array
  get_core_out <- internal_get_core_output(ChainNums = ChainNums, jags_data = data)
  core_out <- list(SWLI = get_core_out$SWLI, mcmc.array=get_core_out$x0.samps)
  #Get parameter diagnostics (for SWLI estimates)
  # get_diagnostics(n=nrow(y0),
  #                 output.dir = output.dir,
  #                 use.informative.priors=use.informative.priors)

  if(validation.run)
  {
    validresults(modern_elevation.csv,fold=fold)
  }

  class(core_out) = 'BTF'

  invisible(core_out)

}

#-----------------------------------------------------
InternalRunCore <- function(#Do MCMC sampling
  ###Do MCMC sampling for one chain
  chainNum, ##<< Chain ID
  jags_data,
  jags_pars,
  n.burnin,
  n.iter,
  n.thin
){
  # set seed before sampling the initial values
  set.seed.chain <- chainNum*209846
  jags.dir <- file.path(getwd(), "temp.JAGSobjects/")
  set.seed(set.seed.chain)
  temp <- rnorm(1)

  #The model for the modern data
cat("
#--------------------------------------------------------------
# Model for BTF
# Niamh Cahill 2020
#--------------------------------------------------------------
    
model{",sep="",append=FALSE, file = file.path("model.txt"), fill = TRUE) 
  
  cat("
  for(i in 1:n)
  {

  # Set up for basis functions where x0 needs to be estimated
  for(k in 1:n_knots)
  {
  J[i,k]<-step(x0[i]-knots[k])
  L[i,k]<-((x0[i]-knots[k])^deg)*J[i,k]
  }

  for(j in begin0:m){
  lambda[i,j] <- 1
  }
  for(j in 1:(begin0-1)){
  spline[i,j]<- alpha0 + beta0.j[j] + inprod(Z0.ih[i,],delta0.hj[,j])

  ##Account for over/under dispersion
  z[i,j] ~ dnorm(spline[i,j],tau_z[j])
  lambda[i,j]<-exp(z[i,j])
  }#End j loop

  y[i,]~dmulti(p[i,],N_count[i])
  lambdaplus[i]<-sum(lambda[i,])

  for(j in 1:m){
  p[i,j]<-lambda[i,j]/lambdaplus[i]
  }#End j loop
  }#End i loop

  #Get basis functions
  B0.ik <- pow(-1,(deg + 1)) * (L %*% t(D))
  Z0.ih <- B0.ik%*%Deltacomb.kh

  #Species specific
  for(j in 1:(begin0-1))
  {
  #Over dispersion parameter
   tau_z[j]<-1/(pow(sigma.z[j],2) +
                pow(beta0_sd[j],2) + 
                pow(alpha0_sd,2) + 
                pow(delta0_sd[j],2))
  
   #tau_z[j]<-1/(pow(sigma.z[j],2))

  }
  ",sep = "",append = TRUE, file = "model.txt",fill = TRUE)
  
  if(jags_data$use_uniform_prior == TRUE)
  {
  cat("
    for(i in 1:n)
    {
    #Prior for x0
    x0[i]~dunif(emin[i],emax[i])
    }
    ", sep = "", append = TRUE, file = "model.txt",fill = TRUE)
  }
  
  if(jags_data$use_uniform_prior == FALSE)
  {
    cat(
    "
    for(i in 1:n)
    {
    #Prior for x0
      x0[i]~dnorm(x0.mean[i],sd.x0[i]^-2)
      x0.mean[i]~dt(el_mean,1,1)T(emin[i],emax[i])
      sd.x0[i] <- 0.1
    }
    ",sep = "", append = TRUE,  file = "model.txt",fill = TRUE)
  }
        
  cat("}"
      ,sep = "", append = TRUE, file = "model.txt",fill = TRUE)
  
  mod <- suppressWarnings(jags(data=jags_data,
            parameters.to.save=jags_pars,
            model.file = "model.txt",
            n.chains=1,
            n.iter=n.iter,
            n.burnin=n.burnin,
            n.thin=n.thin,
            DIC=FALSE,
            jags.seed = set.seed.chain))

  mod_upd <- mod
  save(mod_upd, file=file.path(getwd(), "temp.JAGSobjects", paste0("jags_mod", chainNum, ".Rdata")))

  cat(paste("Hooraah, Chain", chainNum, "has finished!"), "\n")
  return(invisible())
}

internal_get_core_output <- function(ChainNums, jags_data)
{

  mcmc.array <- ConstructMCMCArray(ChainIDs = ChainNums)
  n <- jags_data$n
  
  pars.check<-rep(NA,n)
  for(i in 1:n)
    pars.check[i]<-paste0("x0[",i,"]")
  
  
  # Get gelman diagnostics (Rhat threshold = 1.1)
  # If gelman diagnostic fails then stop!
  gd<-gr_diag(mcmc.array,pars.check = pars.check)
  if(gd==-1)
  {
    cat("WARNING! Convergence issues, check trace plots \n")
  }
  # If gelman diagnostic passes then get other diagnostics  
  if(gd==0)
  {
    eff_size(mcmc.array,pars.check = pars.check)
    mcse(mcmc.array,pars.check = pars.check)
  }
  
  n_samps<-dim(mcmc.array)[1]

  #Get the sorted core data
  Depth<-jags_data$depth
  x0.samps<-array(NA,c(n_samps,n))

  for(i in 1:n)
  {
    parname<-paste0("x0[",i,"]")
    x0.samps[,i]<-mcmc.array[1:n_samps,sample(1,ChainNums),parname]
  }

  SWLI<-apply(x0.samps,2,mean)*100
  SWLI_SD<-apply(x0.samps,2,sd)*100
  lower<- SWLI-2*SWLI_SD
  upper<- SWLI+2*SWLI_SD

  sigma<-(upper-lower)/4
  SWLI_data<-cbind(Depth,SWLI,sigma,lower,upper)

  return(list(SWLI = SWLI_data,x0.samps=x0.samps))
}


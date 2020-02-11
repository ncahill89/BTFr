#-----------------------------------------------------
#' Run the modern (calibration) model
#'
#' @param modern_elevation.csv A .csv file location for modern elevations
#' @param modern_counts.csv A .csv file location for modern counts (to be sorted with \code{\link{sort_modern}})
#' @param dx The elevation interval for spacing the spline knots. Defaults to 0.2.
#' @param ChainNums The number of MCMC chains to run
#' @param n.iter The number of iterations
#' @param n.burnin The number of burnin samples
#' @param n.thin The number of thinning
#' @param run.on.server Set to TRUE if you wish to run the chains in parallel on a server
#'
#' @return Nothing is returned, the relevant data will be saved.
#' @export
#' @import R2jags rjags readr
#'
#' @examples
#' run_modern(modern_elevation.csv = 'elevation.csv',
#' modern_counts.csv = 'counts.csv',
#' n.iter=10,
#' n.thin = 1,
#' n.burnin = 1)

run_modern <- function(modern_elevation.csv = NULL, modern_counts.csv = NULL,
    dx = 0.2, ChainNums = seq(1, 3), n.iter = 40000, n.burnin = 10000,
    n.thin = 15, validation.run = FALSE, fold = 1) {

    dir.create(file.path(getwd(), "temp.JAGSobjects"), showWarnings = FALSE)

    # read in the modern data
    if (!is.null(modern_counts.csv)) {
        modern_dat <- read_csv(modern_counts.csv)
    } else modern_dat <- BTF::modern_counts

    # Get the sorted (by species counts) modern data
    modern_data_sorted <- sort_modern(modern_dat)

    # read in the elevation data
    if (!is.null(modern_elevation.csv)) {
        elevation_dat <- read_csv(modern_elevation.csv)
    } else elevation_dat <- BTF::modern_elevation

    modern_elevation <- elevation_dat$SWLI


    if (validation.run) {
        set.seed(3847)
        K <- 10
        folds <- rep(1:K, ceiling(nrow(modern_data_sorted$moderndat_sorted)/K))
        folds <- folds[sample(1:length(modern_elevation))]
        test_samps <- which(folds == fold)
        test_samps

        y <- modern_data_sorted$moderndat_sorted[-test_samps, ]
        x <- (modern_elevation/100)[-test_samps]

    }

    if (!validation.run) {
        y <- modern_data_sorted$moderndat_sorted
        x <- (modern_elevation/100)
    }

    species_names <- modern_data_sorted$species_names

    # Get min/max elevations (will be used with priors)
    elevation_min <- floor(min(modern_elevation)/10)/10
    elevation_max <- ceiling(max(modern_elevation)/10)/10

    # Get index for the first species (if any) that has all zero counts
    begin0 <- modern_data_sorted$begin0

    # Total species counts
    N_count <- apply(y, 1, sum)

    ###### Regular B Splines Create some basis functions
    res =  bbase(x, xl = elevation_min, xr = elevation_max, dx = dx)  # This creates the basis function matrix
    B.ik <- res$B.ik
    K <- dim(B.ik)[2]

    D = 1
    Delta.hk <- diff(diag(K), diff = D)  # difference matrix
    Deltacomb.kh <- t(Delta.hk) %*% solve(Delta.hk %*% t(Delta.hk))
    Z.ih <- B.ik %*% Deltacomb.kh
    H <- dim(Z.ih)[2]

    # Create for predictions
    grid_size = 50
    SWLI_grid = seq(elevation_min, elevation_max, length = grid_size)
    res_star = bbase(SWLI_grid, xl = elevation_min, xr = elevation_max,
        dx = dx)  # This creates the basis function matrix
    Bstar.ik <- res_star$B
    Zstar.ih <- Bstar.ik %*% Deltacomb.kh

    # Jags model data
    pars = c("p", "beta.j", "sigma.z", "sigma.delta", "delta.hj", "splinestar")

    data = list(y = y, n = nrow(y), m = ncol(y), N_count = N_count, H = H,
        Z.ih = Z.ih, Zstar.ih = Zstar.ih, N_grid = grid_size, begin0 = begin0)


        for (chainNum in ChainNums) {
            cat(paste("Start chain ID ", chainNum), "\n")

            InternalRunOneChain(chainNum = chainNum, jags_data = data,
                jags_pars = pars, n.burnin = n.burnin, n.iter = n.iter,
                n.thin = n.thin)

        }


    # Get model output needed for the core run
    data[["x"]] <- x
    jags_data <- list(data = data, pars = pars, elevation_max = elevation_max,
                      elevation_min = elevation_min, dx = dx, species_names = species_names)

    core_input <- internal_get_core_input(ChainNums = ChainNums, jags_data = jags_data)

    # Update jags_data list
    modern_out <- list(data = data, pars = pars, elevation_max = elevation_max,
        elevation_min = elevation_min, dx = dx, species_names = species_names,
        delta0.hj = core_input$delta0.hj, delta0_sd = core_input$delta0_sd, beta0.j = core_input$beta0.j, beta0_sd = core_input$beta0_sd, sigma.z = core_input$sigma.z, src_dat = core_input$src_dat)

    class(modern_out) = 'BTF'

    invisible(modern_out)

}

#-----------------------------------------------------
InternalRunOneChain <- function(chainNum, jags_data, jags_pars, n.burnin,
    n.iter, n.thin) {
    set.seed.chain <- chainNum * 209846
    jags.dir <- file.path(getwd(), "temp.JAGSobjects/")
    set.seed(set.seed.chain)
    temp <- rnorm(1)

    # The model for the modern data
    modernmodel = "
  model
  {

  for(i in 1:n)
  {
  for(j in begin0:m){
  z[i,j] <- 0
  lambda[i,j]<-exp(z[i,j])
  }
  for(j in 1:(begin0-1)){
  spline[i,j]<-beta.j[j]+inprod(Z.ih[i,],delta.hj[,j])
  z[i,j]~dnorm(spline[i,j],tau.z[j])
  lambda[i,j]<-exp(z[i,j])
  }#End j loop

  y[i,]~dmulti(p[i,],N_count[i])
  lambdaplus[i]<-sum(lambda[i,])
  }#End i loop

  ###Get p's for multinomial
  for(i in 1:n){
  for(j in 1:m){
  p[i,j]<-lambda[i,j]/lambdaplus[i]
  }#End j loop
  }#End i loop


  #####Spline parameters#####
  #Coefficients
  for(j in 1:(begin0-1)){
  for (h in 1:H)
  {
  delta.hj[h,j] ~ dnorm(0, tau.delta)
  }
  }
  #Smoothness
  tau.delta<-pow(sigma.delta,-2)
  sigma.delta~dt(0, 2^-2, 1)T(0,)

  ###Variance parameter###
  for(j in 1:m){
  tau.z[j]<-pow(sigma.z[j],-2)
  sigma.z[j]~dt(0, 2^-2, 1)T(0,)

  ###Intercept (species specific)
  beta.j[j]~dnorm(0,0.01)
  }

  ########Get predictions
  for(i in 1:N_grid){
  for(j in begin0:m){
  zstar[i,j] <- 0
  lambdastar[i,j]<-exp(zstar[i,j])
  }
  for(j in 1:(begin0-1)){
  splinestar[i,j] <- beta.j[j]+inprod(Zstar.ih[i,],delta.hj[,j])
  zstar[i,j] ~ dnorm(splinestar[i,j],tau.z[j])
  lambdastar[i,j] <- exp(zstar[i,j])
  }#End j loop

  }
  }##End model
  "

    mod <- jags(data = jags_data, parameters.to.save = jags_pars, model.file = textConnection(modernmodel),
        n.chains = 1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
        DIC = FALSE, jags.seed = set.seed.chain)

    mod_upd <- mod
    save(mod_upd, file = file.path(getwd(), "temp.JAGSobjects", paste0("jags_mod",
        chainNum, ".Rdata")))

    cat(paste("Hooraah, Chain", chainNum, "has finished!"), "\n")

    return(invisible())
}


internal_get_core_input <- function(ChainNums, jags_data)
{

  mcmc.array <- ConstructMCMCArray(ChainIDs = ChainNums)

  n_samps<-dim(mcmc.array)[1]

  #########For Splines
  # This creates the components for the basis function matrix
  xl<-jags_data$elevation_min
  xr<-jags_data$elevation_max
  begin0 <- jags_data$data$begin0
  deg=3
  dx<-jags_data$dx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  n_knots<-length(knots)
  D <- diff(diag(n_knots), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  K <- dim(D)[1]
  Dmat=1
  Delta.hk <- diff(diag(K), diff = Dmat) # difference matrix
  Deltacomb.kh <- t(Delta.hk)%*%solve(Delta.hk%*%t(Delta.hk))

  ##########Get parameter estimates
  delta.hj_samps<-array(NA,c(n_samps,jags_data$data$H,(begin0-1)))
  beta.j_samps<-sigma.z_samps<-array(NA,c(n_samps,(begin0-1)))

    for(j in 1:(begin0-1))
    {
      for(h in 1:jags_data$data$H)
      {
        parname<-paste0("delta.hj[",h,",",j,"]")
        delta.hj_samps[,h,j]<-mcmc.array[1:n_samps,sample(1,ChainNums),parname]
      }
      parname<-paste0("beta.j[",j,"]")
      beta.j_samps[,j]<-mcmc.array[1:n_samps,sample(1,ChainNums),parname]
    }


  for(j in 1:(begin0-1))
  {
    parname<-paste0("sigma.z[",j,"]")
    sigma.z_samps[,j]<-mcmc.array[1:n_samps,sample(1,ChainNums),parname]
  }


  delta0.hj<-apply(delta.hj_samps,2:3,mean)
  delta0_sd<-apply((apply(delta.hj_samps,2:3,sd)),2,median)

  beta0.j<-apply(beta.j_samps,2,mean)
  beta0_sd<-apply(beta.j_samps,2,sd)

  sigma.z<-apply(sigma.z_samps,2,median)

  # Get model estimates
  # Data
  y <- jags_data$data$y
  n <- nrow(y)
  m <- ncol(y)
  grid_size = 50
  SWLI_grid = seq(jags_data$elevation_min, jags_data$elevation_max, length = grid_size)
  species_names <- jags_data$species_names

  # ---------------------------------------------------- results objects
  p_star <- p_star_all <- spline_star <- spline_star_all <- array(NA,
                                                                  c(n_samps, length(SWLI_grid), m))

  for (i in 1:n_samps) {
    for (j in begin0:m) {
      spline_star_all[i, , j] <- 1
    }
    for (j in 1:(begin0 - 1)) {
      for (k in 1:length(SWLI_grid)) x.index <- seq(1:length(SWLI_grid))
      spline_star_all[i, , j] <- exp(mcmc.array[i, sample(seq(1,
                                                              3), 1), paste0("splinestar[", x.index, ",", j, "]")])
    }
  }

  for (i in 1:n_samps) {
    for (j in 1:m) {
      p_star_all[i, , j] = spline_star_all[i, , j]/apply(spline_star_all[i,
                                                                         , ], 1, sum)

    }
  }

  # Get predicted values
  # ----------------------------------------------------
  pred_pi_median <- apply(p_star_all, 2:3, median)
  pred_pi_high <- apply(p_star_all, 2:3, "quantile", 0.975)
  pred_pi_low <- apply(p_star_all, 2:3, "quantile", 0.025)

  # Plot of predicted output
  # ------------------------------------------------
  df = data.frame(SWLI_grid * 100, pred_pi_median)
  df_low = data.frame(SWLI_grid * 100, pred_pi_low)
  df_high = data.frame(SWLI_grid * 100, pred_pi_high)

  colnames(df) = c("SWLI", species_names)
  colnames(df_low) = c("SWLI", species_names)
  colnames(df_high) = c("SWLI", species_names)


  df_long = df %>% pivot_longer(names_to = "species", values_to = "proportion",
                                -SWLI)
  df_low_long = df_low %>% pivot_longer(names_to = "species", values_to = "proportion_lwr",
                                        -SWLI)
  df_high_long = df_high %>% pivot_longer(names_to = "species", values_to = "proportion_upr",
                                          -SWLI)


  src_dat = df_long %>% dplyr::mutate(proportion_lwr = df_low_long %>%
                                         dplyr::pull(proportion_lwr), proportion_upr = df_high_long %>%
                                         dplyr::pull(proportion_upr))


 return(list(delta0.hj = delta0.hj, delta0_sd = delta0_sd, beta0.j = beta0.j, beta0_sd = beta0_sd, sigma.z = sigma.z, src_dat = src_dat))
}

#' Carry out a validation
#'
#' @param modern_elevation A dataframe of modern elevations
#' @param modern_species A dataframe of modern counts (to be sorted with \code{\link{sort_modern}})
#' @param n_folds Number of folds for CV
#' @param use_uniform_prior change prior on elevation to be uniform
#' @param dx The elevation interval for spacing the spline knots. Defaults to 0.2.
#' @param modern.iter The number of iterations for the modern model
#' @param modern.burnin The number of burnin samples for the modern model
#' @param modern.thin The number of thinning for the modern model
#' @param core.iter The number of iterations for the core model
#' @param core.burnin The number of burnin samples for the core model
#' @param core.thin The number of thinning for the core model
#' @param scale_x Set to TRUE to scale elevation data to have mean 0 and sd 1
#'
#' @return A tibble of validation results
#' @export
#' @importFrom dplyr "full_join"
#' @import tibble
#'
run_valid <- function(modern_elevation = NULL,
                      modern_species = NULL,
                      n_folds = 1,
                      use_uniform_prior = FALSE,
                      dx = 0.2,
                      modern.iter = 40000,
                      modern.burnin = 10000,
                      modern.thin = 15,
                      core.iter = 15000,
                      core.burnin = 1000,
                      core.thin = 7,
                      scale_x = FALSE) {

  # read in the modern data
  if (is.null(modern_species)) {
    modern_species <- BTFr::NJ_modern_species
  }

  # read in the elevation data
  if (is.null(modern_elevation)) {
    elevation_dat <- BTFr::NJ_modern_elevation
  }

  valid_SWLI <- tibble(
    Depth = numeric(),
    SWLI = numeric(),
    sigma = numeric(),
    lower = numeric(),
    True = numeric()
  )
  for (f in 1:n_folds)
  {
    modern_mod <- BTFr::run_modern(modern_elevation,
      modern_species,
      validation.run = TRUE,
      fold = f,
      dx = dx,
      n.iter = modern.iter,
      n.burnin = modern.burnin,
      n.thin = modern.thin,
      scale_x = scale_x
    )


    core_mod <- BTFr::run_core(modern_mod,
      core_species = modern_mod$data$y_test,
      validation.run = TRUE,
      use_uniform_prior = use_uniform_prior,
      n.iter = core.iter,
      n.burnin = core.burnin,
      n.thin = core.thin
    )

    SWLI_res <- as_tibble(core_mod$SWLI)

    if (modern_mod$scale_x) {
      SWLI_res$True <- modern_mod$data$x_test$value * modern_mod$x_scale + modern_mod$x_center
    }
    if (!modern_mod$scale_x) {
      SWLI_res$True <- modern_mod$data$x_test$value * 100
    }
    valid_SWLI <- full_join(valid_SWLI, SWLI_res)
  }

  return(valid_SWLI)
}

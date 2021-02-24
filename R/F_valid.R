#' Carry out a validation
#'
#' @param modern_elevation 
#' @param modern_species 
#' @param n_folds 
#'
#' @return A table of validation results
#' @export
#' @importFrom dplyr "full_join"
#' @import tibble 
#'
#' @examples
#' run_valid()
run_valid <- function(modern_elevation = NULL,
                      modern_species = NULL,
                      n_folds = 1,
                      use_uniform_prior = FALSE,
                      dx = 0.4)
{
  
  # read in the modern data
  if (is.null(modern_species)) {
     modern_species <- BTF::NJ_modern_species
  }
  
  # read in the elevation data
  if (is.null(modern_elevation)) {
    elevation_dat <- BTF::NJ_modern_elevation
  }
  
  valid_SWLI <- tibble(Depth = numeric(),
                       SWLI = numeric(),
                       sigma = numeric(),
                       lower = numeric(),
                       True = numeric())
  for(f in 1:n_folds)
  {
    modern_mod <- run_modern(modern_elevation,
                             modern_species,
                             validation.run = TRUE,
                             fold = f,
                             dx = dx)
    
    
    core_mod <- BTF::run_core(modern_mod,
                              core_species = modern_mod$data$y_test,
                              validation.run = TRUE,
                              use_uniform_prior = use_uniform_prior)
    
  SWLI_res <- as_tibble(core_mod$SWLI)
  SWLI_res$True <- modern_mod$data$x_test$value*modern_mod$x_scale + modern_mod$x_center
    
  valid_SWLI <- full_join(valid_SWLI,SWLI_res)
  
  
  }
  
  return(valid_SWLI)
  
}

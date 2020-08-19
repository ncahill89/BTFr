#' Carry out a 10-fold cross validation
#'
#' @param modern_elevation 
#' @param modern_counts 
#'
#' @return A table of validation results
#' @export
#'
#' @examples
#' run_valid()
run_valid <- function(modern_elevation = NULL,
                      modern_counts = NULL)
{
  
  # read in the modern data
  if (is.null(modern_counts)) {
     modern_counts <- BTF::modern_counts
  }
  
  # read in the elevation data
  if (is.null(modern_elevation)) {
    elevation_dat <- BTF::modern_elevation
  }
  
  valid_SWLI <- tibble(Depth = numeric(),
                       SWLI = numeric(),
                       sigma = numeric(),
                       lower = numeric(),
                       True = numeric())
  for(f in 1:10)
  {
    modern_mod <- run_modern(modern_elevation,
                             modern_counts,
                             validation.run = TRUE,
                             fold = f)
    
    
    core_mod <- BTF::run_core(modern_mod,
                              core_counts = modern_mod$data$y_test,
                              validation.run = TRUE,
                              use_uniform_prior = FALSE)
    
  SWLI_res <- as_tibble(core_mod$SWLI)
  SWLI_res$True <- modern_mod$data$x_test$value*100
    
  valid_SWLI <- full_join(valid_SWLI,SWLI_res)
  }
  
  return(valid_SWLI)
  
}

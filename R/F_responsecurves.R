#' Create Species Response Curves
#'
#' @param modern_mod An object of class \code{BTFr} from \code{\link{run_modern}}
#' @param species_select a vector of species names for which you want to create response curves
#'
#' @return Response curve data files (empirical data and model-based estimates) and species response curve plots
#' @export
#' @import ggplot2 magrittr
#' @importFrom tidyr 'gather'
#' @examples
#' \donttest{
#' test_modern_mod <- run_modern(
#'   modern_elevation = NJ_modern_elevation,
#'   modern_species = NJ_modern_species,
#'   n.iter = 10,
#'   n.burnin = 1,
#'   n.thin = 1
#' )
#' response_curves(test_modern_mod)
#' }
#'
response_curves <- function(modern_mod,
                            species_select = NULL) {

  # Data
  y <- modern_mod$data$y
  n <- nrow(y)
  m <- ncol(y)
  grid_size <- 50
  SWLI_grid <- seq(modern_mod$elevation_min, modern_mod$elevation_max, length = grid_size)
  N_count <- apply(y, 1, sum)
  species_names <- modern_mod$species_names
  begin0 <- modern_mod$data$begin0
  if (is.null(species_select)) species_select <- species_names
  # Get empirical probabilities
  # ----------------------------------------------------

  Pmat <- matrix(NA, n, m)
  for (i in 1:n) {
    Pmat[i, ] <- y[i, ] / N_count[i]
  }

  if (modern_mod$scale_x) {
    empirical_dat <- data.frame((modern_mod$data$x * modern_mod$x_scale) + modern_mod$x_center, Pmat)
  }
  if (!modern_mod$scale_x) {
    empirical_dat <- data.frame(modern_mod$data$x * 100, Pmat)
  }

  colnames(empirical_dat) <- c("SWLI", species_names)

  # Plot of empirical probabilities
  # ------------------------------------------------

  empirical_data_long <- empirical_dat %>%
    tidyr::pivot_longer(
      names_to = "species",
      values_to = "proportion", -SWLI
    ) %>%
    dplyr::filter(species %in% species_select)

  src_dat <- modern_mod$src_dat %>% dplyr::filter(species %in% species_select)

  p <- ggplot(data = empirical_data_long) +
    geom_point(aes(
      x = SWLI, y = proportion,
      colour = "Observed Data"
    ), alpha = 0.2) +
    geom_line(data = src_dat, aes(x = SWLI, y = proportion, colour = "Model Estimates")) +
    geom_ribbon(data = src_dat, aes(x = SWLI, ymin = proportion_lwr, ymax = proportion_upr), alpha = 0.3) +
    scale_colour_manual(name = "", values = c("red", "royalblue2"), guide = guide_legend(override.aes = list(linetype = c(
      "solid",
      "blank"
    ), shape = c(NA, 16)))) +
    theme_minimal() +
    ylim(0,1) +
    ggtitle("Species Response Curves") +
    facet_wrap(~species,
      scales = "free"
    ) +
    theme_classic()


  return(list(
    src_plot = p, src_empirical_dat = empirical_data_long,
    src_model_dat = src_dat
  ))
}

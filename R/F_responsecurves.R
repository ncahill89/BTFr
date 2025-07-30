#' Create Species Response Curves
#'
#' @param modern_mod An object of class \code{BTFr} from \code{\link{run_modern}}
#' @param species_select a vector of species names for which you want to create response curves
#' @param species_order a vector of species names for how you want the species ordered
#'
#' @return Response curve data files (empirical data and model-based estimates) and species response curve plots
#' @export
#' @import ggplot2 magrittr
#' @importFrom tidyr 'gather' 'pivot_longer'
#' @importFrom dplyr 'filter'
#' @examples
#' \donttest{
#' test_modern_mod <- run_modern(
#'   modern_elevation = NJ_modern_elevation,
#'   modern_species = NJ_modern_species,
#'   n.iter = 10,
#'   n.burnin = 1,
#'   n.thin = 1,
#'   parallel = FALSE
#' )
#' response_curves(test_modern_mod)
#' }
#'
response_curves <- function(modern_mod, species_select = NULL, species_order = NULL) {
  # Extract data
  y <- modern_mod$data$y
  n <- nrow(y)
  m <- ncol(y)
  N_count <- rowSums(y)
  species_names <- modern_mod$species_names

  # Define grid
  grid_size <- 50
  SWLI_grid <- seq(modern_mod$elevation_min, modern_mod$elevation_max, length.out = grid_size)

  # If no species selected, use all
  if (is.null(species_select)) {
    species_select <- species_names
  }

  # Normalize counts
  Pmat <- sweep(y, 1, N_count, FUN = "/")

  # Transform x-values
  if (modern_mod$scale_x) {
    SWLI <- (modern_mod$data$x * modern_mod$x_scale) + modern_mod$x_center
  } else {
    SWLI <- modern_mod$data$x * 100
  }

  # Assemble empirical data
  empirical_dat <- data.frame(SWLI, Pmat)
  colnames(empirical_dat) <- c("SWLI", species_names)

  # Pivot to long format
  empirical_data_long <- empirical_dat %>%
    tidyr::pivot_longer(cols = -SWLI, names_to = "species", values_to = "proportion") %>%
    dplyr::filter(species %in% species_select)

  # Filter model estimates
  src_dat <- modern_mod$src_dat %>%
    dplyr::filter(species %in% species_select)

  # Apply custom species order if provided
  if (!is.null(species_order)) {
    empirical_data_long$species <- factor(empirical_data_long$species, levels = species_order)
    src_dat$species <- factor(src_dat$species, levels = species_order)
  }

  # Create the plot
  p <- ggplot(data = empirical_data_long, aes(x = SWLI, y = proportion)) +
    geom_point(aes(color = "Observed"), alpha = 0.3, size = 1.5) +
    geom_line(data = src_dat, aes(x = SWLI, y = proportion, color = "Model"), size = 1) +
    geom_ribbon(
      data = src_dat,
      aes(x = SWLI, ymin = proportion_lwr, ymax = proportion_upr),
      fill = "grey70", alpha = 0.4, inherit.aes = FALSE
    ) +
    facet_wrap(~species, scales = "free_y") +
    scale_color_manual(
      name = "",
      values = c("Observed" = "#0072B2", "Model" = "#D55E00"),
      labels = c("Observed" = "Observed Data", "Model" = "Model Estimates")
    ) +
    labs(
      title = "Species Response Curves",
      subtitle = "Comparing Observed Data and Model Estimates",
      x = "Standardized Water Level Index (SWLI)",
      y = "Proportion"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold")
    )

  return(list(
    src_plot = p,
    src_empirical_dat = empirical_data_long,
    src_model_dat = src_dat
  ))
}

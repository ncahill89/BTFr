#' Sort modern species dataset based on total count size
#'
#' @param modern_dat A modern data object containing species counts
#'
#' @return Returns a modern dataset that has been sorted based on species count size, as well as the sorted species names and the index for the first species that contains all zero counts (if any)
#' @export
#'
sort_modern <- function(modern_dat) {
  # Find and sort the total no. of each species
  species_total <- sort(apply(modern_dat, 2, sum), decreasing = TRUE)

  # Order the species by total count
  species_order <- match(names(species_total), names(modern_dat))

  # Order the modern data based on the species order
  moderndat_sorted <- as.matrix(round(modern_dat[, species_order]))

  # Get the index for the first column of all 0 counts (if there are any)
  if (any(species_total == 0)) {
    begin0 <- min(which(species_total == 0))
  } else {
    begin0 <- length(species_total)
  }

  species_names <- names(species_total)
  return(list(
    moderndat_sorted = moderndat_sorted,
    begin0 = begin0,
    species_names = species_names
  ))
}

#' Sort core dataset to be compatible with the sorted modern data
#'
#' @param core_dat A data object containing core species counts
#' @param species_names A list of species names that have been returned from \code{\link{sort_modern}}
#'
#' @return Returns the sorted core data where the species columns are compatible with the modern data that has been returned from \code{\link{sort_modern}}
#' @export
#'
sort_core <- function(core_dat, species_names) {
  # Match the core species names to the ordered modern species names
  species_match <- match(species_names, names(core_dat))

  # Sort the core data based on the modern species order
  if (any(is.na(species_match))) {
    message("These species appear in the modern but do not appear in the core: \n", paste(species_names[which(is.na(species_match))], collapse = " "), "\nThe model will assume that these species have zero counts.\n", "Check your data and modify if needed, otherwise continue.")
    readline(prompt = "Press [enter] to continue")

    species_add <- species_names[which(is.na(species_match))]
    names_core <- names(core_dat)
    col_add <- matrix(0, nrow = dim(core_dat)[1], ncol = length(species_add))
    core_dat <- cbind(core_dat, col_add)
    names(core_dat) <- c(names_core, species_add)
  }
  species_order <- match(species_names, names(core_dat))
  coredata_sorted <- as.matrix(round(core_dat[, species_order]))

  return(list(coredata_sorted = coredata_sorted))
}

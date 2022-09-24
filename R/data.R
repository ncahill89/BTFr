#' Core species counts from New Jersey, USA.
#'
#' A dataset containing the raw core counts of foraminifera species for field study sites in New Jersey, USA
#'
#' @format A data frame with 69 rows and 24 variables:
#' \describe{
#'   \item{Ab}{count}
#'   \item{Ad}{count}
#'   \item{Ai}{count}
#'   \item{Am}{count}
#'   \item{As}{count}
#'   \item{Calc}{count}
#'   \item{Ea}{count}
#'   \item{Hs}{count}
#'   \item{Jm+Bp}{count}
#'   \item{Mf}{count}
#'   \item{Mp}{count}
#'   \item{PH}{count}
#'   \item{Ph}{count}
#'   \item{Pi}{count}
#'   \item{Pl}{count}
#'   \item{Rs}{count}
#'   \item{TG}{count}
#'   \item{TO}{count}
#'   \item{Tc}{count}
#'   \item{Ti+Sl}{count}
#'   \item{Tq}{count}
#'   \item{Ts}{count}
#'   \item{Tt}{count}
#'   \item{Depth}{core deth, cm}
#' }
#' @source \url{https://cp.copernicus.org/articles/12/525/2016/}
"NJ_core_species"


#' Modern species counts from New Jersey, USA.
#'
#' A dataset containing the raw modern counts of foraminifera species for field study sites in New Jersey, USA
#'
#' @format A data frame with 175 rows and 23 variables:
#' \describe{
#'   \item{Ab}{count}
#'   \item{Ad}{count}
#'   \item{Ai}{count}
#'   \item{Am}{count}
#'   \item{As}{count}
#'   \item{Calc}{count}
#'   \item{Ea}{count}
#'   \item{Hs}{count}
#'   \item{Jm+Bp}{count}
#'   \item{Mf}{count}
#'   \item{Mp}{count}
#'   \item{PH}{count}
#'   \item{Ph}{count}
#'   \item{Pi}{count}
#'   \item{Pl}{count}
#'   \item{Rs}{count}
#'   \item{TG}{count}
#'   \item{TO}{count}
#'   \item{Tc}{count}
#'   \item{Ti+Sl}{count}
#'   \item{Tq}{count}
#'   \item{Ts}{count}
#'   \item{Tt}{count}
#' }
#' @source \url{https://cp.copernicus.org/articles/12/525/2016/}
"NJ_modern_species"

#' Modern elevations from New Jersey, USA.
#'
#' A dataset containing the elevations corresponding to the modern assemblages for field study sites in New Jersey, USA
#'
#' @format A data frame with 175 rows and 3 variables:
#' \describe{
#'   \item{Site}{Name of field site}
#'   \item{SampleID}{ID for the sample}
#'   \item{SWLI}{Standardised water level index (standardised elevation)}
#'
#' }
#' @source \url{https://cp.copernicus.org/articles/12/525/2016/}
"NJ_modern_elevation"

#' Prior elevation ranges for core assemblages
#'
#' A dataset containing upper and lower bounds for core elevations based on a secondary proxy for field study sites in New Jersey, USA
#'
#' @format A data frame with 69 rows and 3 variables:
#' \describe{
#'   \item{Depth}{core depth, cm}
#'   \item{prior_lwr}{elevation lower bound}
#'   \item{prior_upr}{elevarion upper bound}
#'
#' }
#' @source \url{https://cp.copernicus.org/articles/12/525/2016/}
"NJ_priors"



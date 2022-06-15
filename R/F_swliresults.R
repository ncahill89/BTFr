#' Get SWLI estimates and uncertainty
#'
#' @param obj An object of class \code{BTF} from \code{\link{run_core}}
#'
#' @return Table of results and plot
#' @export
#' @importFrom dplyr "select"
#' @examples
#' test_modern_mod <- run_modern(modern_elevation = NJ_modern_elevation,
#'                                 modern_species = NJ_modern_species,
#'                                 n.iter = 10,
#'                                 n.burnin = 1,
#'                                 n.thin = 1)
#' test_core_mod <- run_core(test_modern_mod,
#'                           core_species = NJ_core_species,
#'                           n.iter = 10,
#'                           n.burnin = 1,
#'                           n.thin = 1)
#' swli_results(test_core_mod)
swli_results<-function(obj)
{

 SWLI_dat <- tibble::as_tibble(obj$SWLI)

  p_core<-ggplot(SWLI_dat, aes(y=SWLI,x = Depth))+
    geom_errorbar(aes(ymin = lower,ymax=upper),alpha = 0.5)+
    geom_point(alpha = 0.5)+
    scale_x_reverse() +
    coord_flip()+
    ylab("SWLI")+
    theme(legend.position="none")

  return(list(p_SWLI = p_core, SWLI_dat = SWLI_dat))
}






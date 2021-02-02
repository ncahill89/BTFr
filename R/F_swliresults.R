#' Get SWLI estimates and uncertainty
#'
#' @param obj An object of class \code{BTF} from \code{\link{run_core}}
#'
#' @return Table of results and plot
#' @export
#' @importFrom dplyr "select"
#' @examples
#' swliresults()
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






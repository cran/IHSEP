#' Simulate the child events

#' \code{simchildren} simulates the birth times of all child events
#' spawned from an event relative the birth time of the parent
#' event. This function is to be called by the simulator function for
#' offspring events and is not meant for external use.
#'
#' @param br numerical scalar in [0,1), the branching ratio, or the
#'     expected number of direct children due to an event
#' @param dis, character string, which gives the name of the common
#'     (positive) distribution of the birth times of the child events
#'     relative to the parent event (referred to as the child
#'     birthtime distribution), such as "exp", "gamma", "weibull",
#'     etc.
#' @param par.dis, list, which gives the values of the (named)
#'     parameter(s) of the child birthtime distribution)
#' @return a numeric vector of length giving the birth times of child
#'     events relative to the parent event in ascending order
#' @seealso \code{\link{simoffspring}}
#' @examples
#'   simchildren(br=0.9,dis="exp",par.dis=list(rate=1))

simchildren <- 
  function(br=0.5,dis="exp",par.dis=list(rate=1)){
  N <- rpois(1,br)
  if(N==0)return(numeric(0))
  sort(do.call(paste0("r",dis),c(list(n=N),par.dis)))
}


#' Simulate the offspring events

#' \code{simoffspring} simulates the birth times of offspring events
#' of all generations spawned from an event relative the birth time of
#' the parent event. This function is to be called by the simulator
#' function for Hawkes processes and is not meant for external use.
#'
#' @param br numerical scalar in [0,1), the branching ratio, or the
#'     expected number of direct children due to an event
#' @param dis, character string, which gives the name of the common
#'     (positive) distribution of the birth times of the child events
#'     relative to the parent event (referred to as the child
#'     birthtime distribution), such as "exp", "gamma", "weibull",
#'     etc.
#' @param par.dis, list, which gives the values of the (named)
#'     parameter(s) of the child birthtime distribution)
#' @return a numeric vector of giving the birth times of offspring
#'     events of all generations relative to the parent's birth time
#'     in ascending order
#' @details This function uses recursion, so can break down when the
#'     branching ratio is close to 1, leading to very deep
#'     recursions. In this case, we should use \code{simHawkes1} for
#'     Hawkes process simulation.
#' @seealso \code{\link{simchildren}}
#' 
#' @examples
#'   simoffspring(br=0.9,dis="exp",par.dis=list(rate=1))


simoffspring <- 
    function(br=0.5,dis="exp",par.dis=list(rate=1)){
    tms <- simchildren(br=br,dis=dis,par.dis=par.dis)
    gs <- length(tms)
    if(gs==0)return(numeric(0));
    tmp <- sapply(tms,
                  function(x){
                      c(x,x+simoffspring(br=br,dis=dis,
                                         par.dis=par.dis))
                  })
    sort(unlist(tmp))
}

#' Simulate an (inhomogeneous) Hawkes self-exciting process

#' \code{simHawkes1a} simulates the event times of an inhomogeneous
#' Hawkes process (IHSEP) with background event intensity/rate
#' \eqn{nu(\cdot)\geq 0}, branching ratio \eqn{\eta\in[0,1)}, and
#' offspring birthtime density \eqn{g(\cdot)}, up to a censoring time
#' \eqn{T}.

#' @param nu a function, which gives the background event intensity
#'     function \eqn{\nu(\cdot)}; needs to be a bounded function on
#'     \eqn{[0,T]}.
#' @param cens, a positive scalar, which gives the censoring time.
#' @param nuM, positive scalar, optional argument giving the maximum
#'     of the background intensity function on \eqn{[0,T]}.
#' @param br, scalar in [0,1), giving the branching ratio.
#' @param dis, character string giving the name of the child birthtime
#'     distribution; 'd${dis}' gives the density function \eqn{g(\cdot)}. 
#' @param par.dis, a (named) list giving the values of the parameters of the child birthtime distribution. 
#' @return a vector giving the event times of an inhomogeneous Hawkes process up to the censoring time in ascending order.
simHawkes1a <-
  function(nu=function(x)rep(100,length(x)),cens=1,
           nuM=max(optimize(nu,c(0,cens),maximum=TRUE)$obj,
             nu(0),nu(cens),.Machine$double.eps^0.5)*1.1,
           br=0.5, dis="exp",par.dis=list(rate=1)){
    tms.g0 <- simPois(nu,cens=cens,int.M=nuM)
    tms <- sapply(tms.g0,function(x)c(x,x+simoffspring(br,dis,par.dis)))
    tms <- unlist(tms)
    tms <- tms[tms<=cens]
    sort(tms)
  }

\name{mloglik1e}
\alias{mloglik1e}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minus loglikelihood of an IHSEP model
}
\description{
Calculates the minus loglikelihood of an IHSEP model with given
baseline inensity function \eqn{\nu} and excitation function
\eqn{g(x)=\sum a_i exp(-b_i x)} for event times \code{jtms} on interval
\code{[0,TT]}.   
}
\usage{
mloglik1e(jtms, TT, nuvs, gcoef, InuT)
}
\arguments{
  \item{jtms}{
    A numeric vector, with values sorted in ascending order. Jump times
    to fit the inhomogeneous self-exciting point process model on. 
}
  \item{TT}{
    A scalar. The censoring time, or the terminal time for
    observation. Should be (slightly) greater than the maximum of \code{jtms}.
}
  \item{nuvs}{
    A numeric vector, giving the values of the baseline intensity
    function \eqn{\nu} at the jumptimes \code{jtms}. 
}
  \item{gcoef}{
    A numeric vector (of 2k elements), giving the parameters
    \code{(a_1,...,a_k,b_1,...,b_k)} of the exponential excitation
    function \eqn{g(x)=\sum_{i=1}^k a_i*exp(-b_i*x)}.  
}
  \item{InuT}{
    A numeric value (scalar) giving the integral of \eqn{\nu} on the
    interval \code{[0,TT]}. 
}

}
\details{
  This version of the mloglik function uses external C code to speedup
  the calculations. When given the analytical form of \code{Inu} or a
  quickly calculatable \code{Inu}, it should be (substantially)
  faster than \code{mloglik1a} when calculating the (minus log)
  likelihood when the excitation function is exponential. Otherwise it
  is the same as \code{mloglik0}, \code{mloglik1a}, \code{mloglik1b}. 
}

\value{
  The value of the negative log-liklihood. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
  \code{mloglik0}, \code{mloglik1a} and \code{mloglik1b}
}

\examples{
## simulated data of an IHSEP on [0,1] with baseline intensity function
## nu(t)=200*(2+cos(2*pi*t)) and excitation function
## g(t)=8*exp(-16*t)
data(asep)

## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1e(tms,TT=1,nuvs=200*(2+cos(2*pi*tms)),
          gcoef=8*1:2,
          InuT=integrate(function(x)200*(2+cos(2*pi*x)),0,1)$value)
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
\dontrun{
system.time(est <- optim(c(400,200,2*pi,8,16),
                         function(p){
                           mloglik1e(jtms=tms,TT=1,
                                     nuvs=p[1]+p[2]*cos(p[3]*tms),
                                     gcoef=p[-(1:3)],
                                     InuT=integrate(function(x)p[1]+p[2]*cos(p[3]*x),
                                                    0,1)$value
                                     )
                         },hessian=TRUE,control=list(maxit=5000,trace=TRUE),
                         method="BFGS")
            )
## point estimate by MLE
est$par
## standard error estimates:
diag(solve(est$hessian))^0.5
}
}
\keyword{ point process }
\keyword{ likelihood }

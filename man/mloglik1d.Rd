\name{mloglik1d}
\alias{mloglik1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minus loglikelihood of an IHSEP model
}
\description{
Calculates the minus loglikelihood of an IHSEP model with given
baseline intensity function \eqn{\nu} and excitation function
\eqn{g(x)=\sum a_i exp(-b_i x)} for event times \code{jtms} on interval
\code{[0,TT]}.   
}
\usage{
mloglik1d(jtms, TT, nu, gcoef, Inu)
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
  \item{nu}{
    A (vectorized) function. The baseline intensity function. 
}
  \item{gcoef}{
    A numeric vector (of 2k elements), giving the parameters
    \code{(a_1,...,a_k,b_1,...,b_k)} of the exponential excitation
    function \eqn{g(x)=\sum_{i=1}^k a_i*exp(-b_i*x)}.  
}
  \item{Inu}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the baseline intensity function \eqn{\nu} from 0 to \code{t}. 
}

}
\details{
  This function calculates the minus loglikelihood of the inhomegeneous
  Hawkes model with background intensity function \eqn{\nu(t)} and
  excitation kernel function \eqn{g(t)=\sum_{i=1}^{k} a_i e^{-b_i t}}
  relative to continuous observation of the process from time 0 to time
  \code{TT}. Like \code{mloglik1c}, it takes advantage of the Markovian 
  property of the intensity process and uses external C++ code to speed
  up the computation.  
}

\value{
  The value of the negative log-liklihood. 
}
\author{
  Feng Chen <feng.chen@unsw.edu.au>
}
\seealso{
  \code{mloglik1c}
}

\examples{
## simulated data of an IHSEP on [0,1] with baseline intensity function
## nu(t)=200*(2+cos(2*pi*t)) and excitation function
## g(t)=8*exp(-16*t)
data(asep)

## get the birth times of all generations and sort in ascending order 
tms <- sort(unlist(asep))
## calculate the minus loglikelihood of an SEPP with the true parameters 
mloglik1d(tms,TT=1,nu=function(x)200*(2+cos(2*pi*x)),
          gcoef=8*1:2,
          Inu=function(y)integrate(function(x)200*(2+cos(2*pi*x)),0,y)$value)
## calculate the MLE for the parameter assuming known parametric forms
## of the baseline intensity and excitation functions  
\dontrun{
system.time(est <- optim(c(400,200,2*pi,8,16),
                         function(p){
                           mloglik1d(jtms=tms,TT=1,
                                     nu=function(x)p[1]+p[2]*cos(p[3]*x),
                                     gcoef=p[-(1:3)],
                                     Inu=function(y){
                                      integrate(function(x)p[1]+p[2]*cos(p[3]*x),
                                                0,y)$value
                                     })
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

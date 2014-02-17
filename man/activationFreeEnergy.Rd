\name{activationFreeEnergy}
\alias{activationFreeEnergy}
\title{
Activation Free Energy
}
\description{
Computes the non-adiabatic activation energy ( \eqn{\Delta E^{\ddagger}}{\DeltaEa} ) involved in the calculation of the Landau-Zener and Marcus charge transfer (CT) rates using the following expression:
    \deqn{\Delta E^{\ddagger}=\frac{\left( \lambda + \Delta E^{0} + \Delta E_{Field} \right)^2}{4\lambda}}{%
          \DeltaEa = (\lambda + \DeltaE0 + \DeltaEField)^2/(4\lambda)}

}
\usage{
activationFreeEnergy(lambda, dE0 = 0, dEField = 0)
}
\arguments{
  \item{lambda}{a scalar, a vector, a matrix or an array containing the total reorganization energies (\bold{in eV}) used to calculate the activation energy.}
  \item{dE0}{a scalar, a vector, a matrix or an array containing the site energy differences (\bold{in eV}) used to calculate the activation energy. By default self-exchange reactions are considered (\code{dE0=0}).}
  \item{dEField}{a scalar, a vector, a matrix or an array containing an additional contribution to the site energy differences due to an external electric field (\bold{in eV}). By default no electric field is applied (\code{dEField=0}).}
}
\details{
The arguments of these function can be scalars, vectors, matrices or arrays. Mixing scalar values with vectors, matrices or arrays is allowed but in all other cases the arguments must have the same dimensions and lengths.
}
\value{
Depending on the dimension of the objects passed to the function a scalar, a vector, a matrix or an array containing the activation energies (\bold{in eV}) is returned.
}
\references{
H. Oberhofer and J. Blumberger, \var{Physical Chemistry Chemical Physics}, 14, 13846–13852, \bold{2012}\cr
\cr
B. S. Brunschwig, J. Logan, M. D. Newton and N. Sutin, \var{Journal of the American Chemical Society}, 102, 5798, \bold{1980}\cr
\cr
M. D. Newton and N. Sutin, \var{Annual Review of Physical Chemistry}, 35, 437, \bold{1984}\cr
\cr
M. D. Newton, \var{Chemical Reviews}, 91, 767, \bold{1991}.\cr
\cr
A. Nitzan, Chemical Dynamics in Condensed Phases, Oxford University Press, \bold{2006}
}
\seealso{
\code{\link{Marcus}}, \code{\link{LandauZener}}, \code{\link{dEField}}
}
\examples{
## Produce a map of the non-adiabatic free energy barrier and
## two other maps of the non-adiabatic free energy barrier corrected by
## the deviation between the adiabatic and non-adiabatic potential energy curves at:
##   - the crossing point (dEad = dEna - abs(J))
##   - both the crossing point and the minimas of the curves (dEad = dEna - Delta)
## We assume that the reorganization energies are the following:
lambdaI <- 0.14  # internal reorganization energy in eV
lambdaS <- 36E-3 # external reorganization energy in eV
lambda <- lambdaI + lambdaS

N  <- 301
J  <- seq( 0  , 65,length.out=N)*1E-3 # eV
dE <- seq(-0.5,0.5,length.out=N)      # eV
G  <- expand.grid(J, dE)
J  <- G[,1]
dE <- G[,2]

dEna      <- activationFreeEnergy(lambda = lambda, dE0 = dE, dEField = 0)
dEad      <- dEna - adiabaticCorrection(J = J, lambda = lambda, dE0 = dE, dEField = 0)
dEadCross <- dEna -  J

dEna      <- matrix(dEna     , nrow = N, ncol = N)
dEad      <- matrix(dEad     , nrow = N, ncol = N)
dEadCross <- matrix(dEadCross, nrow = N, ncol = N)

addAxis <- function(bottom = TRUE, left = FALSE, above = FALSE, right = FALSE){
  useless <- lapply(1:4,axis, labels=FALSE)
  if(bottom) axis(1, labels = TRUE)
  if(left  ) axis(2, labels = TRUE)
  if(above ) axis(3, labels = TRUE)
  if(right ) axis(4, labels = TRUE)
  if(bottom) mtext(side=1,line=1.2, text=expression( abs(J)/eV), cex=par("cex"))
  if(left  ) mtext(side=2,line=1.2, text=expression(Delta*E/eV), cex=par("cex"))
  if(right ) mtext(side=4,line=1.2, text=expression(Delta*E/eV), cex=par("cex"))
  box()
}

layout(matrix(1:3, ncol=3))
par(cex=2, lwd=1.5, pty="s", mgp=c(1.1,0.1,0), tck=0.02, mar=rep(0.7,4), oma=rep(2,4))
contour(unique(J), unique(dE), dEna     , levels = seq(-0.1,0.6,0.1),
        xaxt="n", yaxt="n", labcex=2)
addAxis(TRUE, TRUE, FALSE, FALSE)
title(expression(Delta*E[na]), line=1)
contour(unique(J), unique(dE), dEadCross, levels = seq(-0.1,0.6,0.1),
        xaxt="n", yaxt="n", labcex=2)
addAxis(TRUE, FALSE, FALSE, FALSE)
title(expression(Delta*E[na]-J), line=1)
contour(unique(J), unique(dE), dEad     , levels = seq(-0.1,0.6,0.1),
        xaxt="n", yaxt="n", labcex=2)
addAxis(TRUE, FALSE, FALSE, TRUE)
title(expression(Delta*E[na]-Delta), line=1)
}
\keyword{ manip }

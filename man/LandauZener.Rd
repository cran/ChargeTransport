\name{LandauZener}
\alias{LandauZener}
\alias{LandauZenerProbability}
\alias{electronicTransmissionCoefficient}
\alias{adiabaticCorrection}
\title{
Landau-Zener Charge Transfer Rates
}
\description{
Computes charge transfer (CT) rates using the semi-classic Landau-Zener expression obtained by combining transition state theory with the Landau-Zener treatment for non-adiabatic transitions (it includes an adiabatic correction to the Marcus expression):
\deqn{k=\kappa _{el} \nu _{n} \Gamma _{n} exp \left[ - \beta \left( \Delta E^{\ddagger} - \Delta\right) \right]}{%
      k = \kappael * \nun * \Gamman * exp( - \beta * ( \DeltaEa - \Delta ) )}
where \eqn{\kappa _{el}}{\kappael} and \eqn{\Gamma _{n}}{\Gamman} are respectively the thermally averaged electronic transmission coefficient and nuclear tunneling factor. \eqn{\nu _{n}}{\nun} is the frequency of the effective vibrational mode assisting the reaction. \eqn{\Delta E^{\ddagger}}{\DeltaEa} is the non-adiabatic activation energy (see \code{\link{activationFreeEnergy}}) and \eqn{\Delta}{\Delta} is a correction factor relating \eqn{\Delta E^{\ddagger}}{\DeltaEa} to the adiabatic activation energy (\eqn{=\Delta E^{\ddagger} - \Delta}{= \DeltaEa - \Delta}). (\eqn{\beta = 1/k_{B}T}{\beta = 1/kbT} where \eqn{k_{B}}{kb} is the Boltzmann constant and \eqn{T}{T} the temperature).

The electronic transmission coefficient accounting for multiple crossings at the intersection region is given by:
\deqn{\kappa _{el} = \frac{2P_{LZ}}{1+P_{LZ}}}{%
      \kappael = 2PLZ/(1+PLZ)}
\deqn{P_{LZ} = 1 - exp \left( -2 \pi \gamma \right)}{%
      PLZ = 1 - exp(-2\pi\gamma)}
\deqn{2 \pi \gamma = \frac{\pi ^{3/2} \left| J \right| ^2}{h \nu _n \sqrt{ \lambda k_{B}T}}}{%
      2\pi\gamma = ( \pi^{3/2} * |J|^2 ) / ( h \nun * \sqrt( \lambda kb T ) )}
where \eqn{P_{LZ}}{PLZ} is the Landau-Zener transition probability for a single surface crossing event, \eqn{\lambda}{\lambda} is the reorganization energy and \eqn{J}{J} is the electronic coupling. The activation free energy and the adiabatic correction factor are given by:
\deqn{\Delta E^{\ddagger}=\frac{\left( \lambda + \Delta E^{0} + \Delta E_{Field} \right)^2}{4\lambda}}{%
\DeltaEa = (\lambda + \DeltaE0 + \DeltaEField)^2/(4\lambda)}
    \deqn{\Delta=\left|J\right|+\frac{\lambda+\Delta E^{0}+\Delta E_{Field}}{2}-\sqrt{\frac{\left(\lambda+\Delta E^{0}+\Delta E_{Field}\right)^2}{4}+\left|J\right|^2}}{%
\Delta = |J| + (\lambda + \DeltaE0 + \DeltaEField)/2 - \sqrt( ( \lambda+\DeltaE0+\DeltaEField )^2/4 + |J|^2 )}

The parameter \eqn{2 \pi \gamma}{2 \pi \gamma} determines the adiabaticity of the reaction. If \eqn{2 \pi \gamma << 1}{2 \pi \gamma << 1} the reaction is non-adiabatic and the Landau-Zener becomes the Marcus expression. In the opposite limit \eqn{2 \pi \gamma >> 1}{2 \pi \gamma >> 1} the ET is adiabatic and \eqn{P_{LZ}}{PLZ} and \eqn{\kappa _{el}}{\kappael} approach unity. Then, the rate expression is the same as for standard chemical reactions in the classical transition state
approximation:
\deqn{k_{ad} = \nu_{n} exp \left[ - \beta \left( \Delta E^{\ddagger} - \Delta \right) \right]}{%
kad = \nun*exp(-\beta (\DeltaEa - \Delta))}
}
\usage{
LandauZener(J, lambda, nuN, dE0 = 0, dEField = 0, temp = 300, gammaN = 1)

electronicTransmissionCoefficient(J, lambda, nuN, temp = 300)

LandauZenerProbability(J, lambda, nuN, temp = 300)

adiabaticCorrection(J, lambda, dE0 = 0, dEField = 0)
}
\arguments{
  \item{J}{a scalar, a vector, a matrix or an array containing the electronic couplings (\bold{in eV}) used to calculate the CT rates.}
  \item{lambda}{a scalar, a vector, a matrix or an array containing the total reorganization energies (\bold{in eV}) used to calculate the CT rates.}
  \item{nuN}{a scalar giving the frequency of the effective vibrational mode assisting the CT process (\bold{in eV}) used to calculate the CT rates.}
\item{dE0}{a scalar, a vector, a matrix or an array containing the site energy differences (\bold{in eV}) used to calculate the CT rates. By default self-exchange reactions are considered (\code{dE0=0}).}
  \item{dEField}{a scalar, a vector, a matrix or an array containing an additional contribution to the site energy differences due to an external electric field (\bold{in eV}). By default no electric field is applied (\code{dEField=0}).}
  \item{temp}{a scalar giving the temperature (\bold{in Kelvin}) at which to evalute the CT rates. By default CT rates are evaluated at room temperature (temp=300).}
  \item{gammaN}{a scalar giving the nuclear tunneling factor. By default the transition is assumed to occur without any tunneling effect (\code{gammaN=1}).}
}
\details{
The arguments of these function can be scalars, vectors, matrices or arrays. Mixing scalar values with vectors, matrices or arrays is allowed but in all other cases the arguments must have the same dimensions and lengths. Using matrices or arrays is useful to compute simultaneously several charge transfer rates for different pairs of molecules, structures ...
}
\value{
Depending on the dimension of the objects passed to the function a scalar, a vector, a matrix or an array containing the Landau-Zener CT rates (\bold{in s-1}) is returned.
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
\code{\link{energyConversion}}, \code{\link{dEField}}, \code{\link{activationFreeEnergy}}, \code{\link{Marcus}}, \code{\link{MarcusLevichJortner}}, \code{\link{KMC}}
}
\examples{
## Produce a map of the decimal logarithm of the Marcus,
## Marcus-Levich-Jortner and Landau-Zener rate expressions for:
nuN <- 1445 # effective vibrational mode wavenumber in cm-1
lambdaI <- 0.14  # internal reorganization energy in eV
lambdaS <- 36E-3 # external reorganization energy in eV

N  <- 301
J  <- seq( 0  , 65,length.out=N)*1E-3 # eV
dE <- seq(-0.5,0.5,length.out=N)      # eV
G  <- expand.grid(J, dE)
J  <- G[,1]
dE <- G[,2]

kMLJ    <- MarcusLevichJortner(
           J = J, lambdaI = lambdaI, lambdaS = lambdaS,
           hBarW = centimeterMinusOne2electronVolt(nuN), dE0 = dE)
kMarcus <- Marcus(
           J = J, lambda = lambdaI+lambdaS, dE0 = dE)
kLZ     <- LandauZener(
           J = J, lambda = lambdaI+lambdaS,
           nuN = centimeterMinusOne2Hertz(nuN), dE0 = dE)

kMLJ    <- matrix(kMLJ   , nrow = N, ncol = N)
kMarcus <- matrix(kMarcus, nrow = N, ncol = N)
kLZ     <- matrix(kLZ    , nrow = N, ncol = N)

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
contour(unique(J), unique(dE), log10(kMLJ   ),
        zlim = c(1,15), levels = -15:15, xaxt="n", yaxt="n", labcex=3)
addAxis(TRUE, TRUE, FALSE, FALSE)
title("Marcus-Levich-Jortner", line=1)
contour(unique(J), unique(dE), log10(kMarcus),
        zlim = c(1,15), levels = -15:15, xaxt="n", yaxt="n", labcex=3)
addAxis(TRUE, FALSE, FALSE, FALSE)
title("Marcus", line=1)
contour(unique(J), unique(dE), log10(kLZ    ),
        zlim = c(1,15), levels = -15:15, xaxt="n", yaxt="n", labcex=3)
addAxis(TRUE, FALSE, FALSE, TRUE)
title("Landau-Zener", line=1)
}
\keyword{ manip }

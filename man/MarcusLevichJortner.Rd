\name{MarcusLevichJortner}
\alias{MarcusLevichJortner}
\title{
Marcus-Levich-Jortner Charge Transfer Rates
}
\description{
Compute charge transfer (CT) rates using the Marcus-Levich-Jortner expression:
       \deqn{k=\frac{2\pi}{\hbar}J^2\sqrt{\frac{1}{4\pi\lambda_{s} k_{B}T}}\times}{%
             k = 2\pi/h_bar * J^2 * \sqrt(1/4\pi\lambdaskbT)]*}
       \deqn{\sum\limits_{n=0}^{\infty}exp(-S)\frac{S^n}{n!}exp\left(-\frac{\left(\Delta E^{0}+\Delta E_{Field}+\lambda_{s}+n\hbar\omega\right)^{2}}{4\lambda_{s} k_{B}T}\right)}{%
             \sum(n=0,\infty) [ exp(-S) * S^n/n! * exp(-(\DeltaE0+\DeltaEField+\lambdas+nhbar\omega)^2/4\lambdaskbT) ]}
where \eqn{J}{J} and \eqn{\Delta E^{0}=E_{final}-E_{initial}}{\DeltaE0 = Ef - Ei} are respectively the electronic coupling and the site energy difference between the initial and final electronic states involved in the CT reaction, and \eqn{\lambda_{s}}{\lambdas} is the external (intermolecular) reorganization energy. \eqn{S=\lambda_{i}/\hbar\omega}{S = \lambdai / hbar \omega} is the electron-phonon coupling (Huang-Rhys factor), directly related to the internal (intramolecular) reorganization energy ( \eqn{\lambda_{i}}{\lambdai} ) and to the energy of the effective vibrational mode ( \eqn{\hbar\omega}{hbar \omega} ) assisting the CT process. \eqn{\Delta E_{Field}}{dEField} is an additional contribution to the site energy difference due to an external electric field and \eqn{T}{T} is the temperature.}
\usage{
MarcusLevichJortner(J, lambdaI, lambdaS, hBarW,
                    dE0 = 0, dEField = 0, temp = 300, nVib = 50)
}
\arguments{
  \item{J}{a scalar, a vector, a matrix or an array containing the electronic couplings (\bold{in eV}) used to calculate the CT rates.}
  \item{lambdaI}{a scalar, a vector, a matrix or an array containing the internal reorganization energies (\bold{in eV}) used to calculate the CT rates.}
  \item{lambdaS}{a scalar, a vector, a matrix or an array containing the external reorganization energies (\bold{in eV}) used to calculate the CT rates.}
  \item{hBarW}{a scalar giving the energy of the effective vibrational mode assisting the charge transfer process (\bold{in eV}) used to calculate the CT rates.}
  \item{dE0}{a scalar, a vector, a matrix or an array containing the site energy differences (\bold{in eV}) used to calculate the CT rates. By default self-exchange reactions are considered (\code{dE0=0}).}
  \item{dEField}{a scalar, a vector, a matrix or an array containing an additional contribution to the site energy differences due to an external electric field (\bold{in eV}). By default no electric field is applied (\code{dEField=0}).}
  \item{temp}{a scalar giving the temperature (\bold{in Kelvin}) at which to evalute the CT rates. By default CT rates are evaluated at room temperature (temp=300).}
  \item{nVib}{an integer indicating the number of vibrational states used to compute the Franck-Condon factor. By default 50 vibrational states are considered.}
}
\details{
The arguments of these function can be scalars, vectors, matrices or arrays. Mixing scalar values with vectors, matrices or arrays is allowed but in all other cases the arguments must have the same dimensions and lengths. Using matrices or arrays is usefull to compute simultaneously several charge transfer rates for different pairs of molecules, structures ...
}
\value{
Depending on the dimension of the objects passed to the function a scalar, a vector, a matrix or an array containing the Marcus-Levich-Jortner CT rates (\bold{in s-1}) is returned.
}
\references{
R.A. Marcus, \var{Journal of Chemical Physics}, 24:966, \bold{1956}\cr
\cr
R.A. Marcus, \var{Reviews of Modern Physics}, 65:599, \bold{1993}\cr
\cr
S.A. Rice and M.T. Vala, \var{Journal of Chemical Physics}, 42:73, \bold{1965}
}
\seealso{
\code{\link{energyConversion}}, \code{\link{dEField}}, \code{\link{Marcus}}, \code{\link{LandauZener}}, \code{\link{KMC}}
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


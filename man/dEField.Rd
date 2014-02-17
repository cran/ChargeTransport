\name{dEField}
\alias{dEField}
\title{Site Energy Difference due to an Electric Field}
\description{
For a given type of charge carrier, this function computes the site energy differences due to an external electric field ( \eqn{dE_{Field}=-e\overrightarrow{d}.\overrightarrow{F}}{dEField = -ed.F} for a hole and \eqn{dE_{Field}=e\overrightarrow{d}.\overrightarrow{F}}{dEField = ed.F} for an electron )
}
\usage{
dEField(dx, dy, dz, carrier = "e", Fx = 0, Fy = 0, Fz = 1E5)
}
\arguments{
  \item{dx, dy, dz}{scalars, vectors, matrices or arrays containing respectively the x-, y- and z-components of the inter-sites distances (\bold{in Angstrom}).}
  \item{carrier}{a character string indicating the type of charge carrier for which to compute the site energy differences due to the electric field. Can be either "e" or "h", respectively for electron and hole. By default electron transport (\code{carrier="e"}) is considered.}
  \item{Fx, Fy, Fz}{scalar values giving the x-, y- and z-component of the electric field (\bold{in V.cm-1}). By default an electric field of 1E5 V.cm-1 is applied along the z-axis.}
}
\details{
The arguments of these function can be scalars, vectors, matrices or arrays. Mixing scalar values with vectors, matrices or arrays is allowed but in all other cases the arguments must have the same dimensions and lengths.
}
\value{
Depending on the dimension of the objects passed to the function a scalar, a vector, a matrix or an array containing the site energy differences due to the electric field (\bold{in eV}) is returned.
}
\seealso{
\code{\link{energyConversion}}, \code{\link{Marcus}}, \code{\link{MarcusLevichJortner}}, \code{\link{LandauZener}}, \code{\link{KMC}}
}
\examples{
dz <- rnorm(10,3.5,0.3)
Fnorm=1E5
# Electron transport along a 1D-stack of 10 molecules aligned along the z-axis
dEField(dx=0, dy=0, dz=dz, carrier="e", Fx=0, Fy=0, Fz=Fnorm)

# Hole transport along a 1D-stack of 10 molecules aligned along the z-axis
dEField(dx=0, dy=0, dz=dz, carrier="h", Fx=0, Fy=0, Fz=Fnorm)

# Some disorder in the xy-plan
dx <- rnorm(10,0,1)
dy <- rnorm(10,0,1)
dEField(dx=dx, dy=dy, dz=dz, carrier="e", Fx=0, Fy=0, Fz=Fnorm)

# Change the orientation of the electric field
theta = 45
dEField(dx=dx, dy=dy, dz=dz, carrier="h",
        Fx=0, Fy=Fnorm*sin(theta*pi/180), Fz=Fnorm*cos(theta*pi/180))
}
\keyword{ manip }

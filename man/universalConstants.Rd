\name{universalConstants}
\alias{universalConstants}
\docType{data}
\title{universal Constants}
\description{This data set provides various universal constants}
\usage{universalConstants}
\format{
  A data frame containing for each universal constant the following information.
  \describe{
    \item{\code{Quantity}}{a character string describing the constant.}
    \item{\code{Value}}{the value of the constant (numerical).}
    \item{\code{Unit}}{a character string indicating the unit of the constant.}
  }
  Constants such as the speed of light, the Planck or Boltzmann constants are available. 
}
\source{http://www.ebyte.it/library/educards/constants/ConstantsOfPhysicsAndMath.html}

\examples{
# List all the constants
print(universalConstants)

# Data for the speed of light
universalConstants["c",]

# Return the speed of light in m.s-1
universalConstants["c","Value"]

# Return the Planck constant in J.s
universalConstants["h","Value"]

# Return the Planck constant in eV.s
Joule2electronVolt(universalConstants["h","Value"])

}
\keyword{datasets}

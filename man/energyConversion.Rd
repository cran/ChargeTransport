\name{energyConversion}
\alias{energyConversion}
\alias{centimeterMinusOne2electronVolt}
\alias{centimeterMinusOne2Hartree}
\alias{centimeterMinusOne2Hertz}
\alias{centimeterMinusOne2Joule}
\alias{centimeterMinusOne2Kelvin}
\alias{centimeterMinusOne2kiloCaloriePerMole}
\alias{centimeterMinusOne2kiloJoulePerMole}
\alias{electronVolt2centimeterMinusOne}
\alias{electronVolt2Hartree}
\alias{electronVolt2Hertz}        
\alias{electronVolt2Joule}
\alias{electronVolt2Kelvin}        
\alias{electronVolt2kiloCaloriePerMole}
\alias{electronVolt2kiloJoulePerMole}     
\alias{Hartree2centimeterMinusOne}
\alias{Hartree2electronVolt}     
\alias{Hartree2Hertz}
\alias{Hartree2Joule}     
\alias{Hartree2Kelvin}
\alias{Hartree2kiloCaloriePerMole}
\alias{Hartree2kiloJoulePerMole}
\alias{Hertz2centimeterMinusOne}
\alias{Hertz2electronVolt}
\alias{Hertz2Hartree}
\alias{Hertz2Joule}
\alias{Hertz2Kelvin}
\alias{Hertz2kiloCaloriePerMole}
\alias{Hertz2kiloJoulePerMole}
\alias{Joule2centimeterMinusOne}
\alias{Joule2electronVolt}
\alias{Joule2Hartree}
\alias{Joule2Hertz}
\alias{Joule2Kelvin}
\alias{Joule2kiloCaloriePerMole}
\alias{Joule2kiloJoulePerMole}
\alias{Kelvin2centimeterMinusOne}
\alias{Kelvin2electronVolt}
\alias{Kelvin2Hartree}
\alias{Kelvin2Hertz}
\alias{Kelvin2Joule}
\alias{Kelvin2kiloCaloriePerMole}
\alias{Kelvin2kiloJoulePerMole}
\alias{kiloCaloriePerMole2centimeterMinusOne}
\alias{kiloCaloriePerMole2electronVolt}
\alias{kiloCaloriePerMole2Hartree}
\alias{kiloCaloriePerMole2Hertz}
\alias{kiloCaloriePerMole2Joule}
\alias{kiloCaloriePerMole2Kelvin}
\alias{kiloCaloriePerMole2kiloJoulePerMole}
\alias{kiloJoulePerMole2centimeterMinusOne}
\alias{kiloJoulePerMole2electronVolt}
\alias{kiloJoulePerMole2Hartree}
\alias{kiloJoulePerMole2Hertz}
\alias{kiloJoulePerMole2Joule}
\alias{kiloJoulePerMole2Kelvin}
\alias{kiloJoulePerMole2kiloCaloriePerMole}
\title{Energy Conversion}
\description{
Conversion of energy units.
}
\usage{
centimeterMinusOne2electronVolt(x)
centimeterMinusOne2Hartree(x)
centimeterMinusOne2Hertz(x)
centimeterMinusOne2Joule(x)
centimeterMinusOne2Kelvin(x)
centimeterMinusOne2kiloCaloriePerMole(x)
centimeterMinusOne2kiloJoulePerMole(x)

electronVolt2centimeterMinusOne(x)
electronVolt2Hartree(x)
electronVolt2Hertz(x)
electronVolt2Joule(x)
electronVolt2Kelvin(x)
electronVolt2kiloCaloriePerMole(x)
electronVolt2kiloJoulePerMole(x)

Hartree2centimeterMinusOne(x)
Hartree2electronVolt(x)
Hartree2Hertz(x)
Hartree2Joule(x)
Hartree2Kelvin(x)
Hartree2kiloCaloriePerMole(x)
Hartree2kiloJoulePerMole(x)

Hertz2centimeterMinusOne(x)
Hertz2electronVolt(x)
Hertz2Hartree(x)
Hertz2Joule(x)
Hertz2Kelvin(x)
Hertz2kiloCaloriePerMole(x)
Hertz2kiloJoulePerMole(x)

Joule2centimeterMinusOne(x)
Joule2electronVolt(x)
Joule2Hartree(x)
Joule2Hertz(x)
Joule2Kelvin(x)
Joule2kiloCaloriePerMole(x)
Joule2kiloJoulePerMole(x)

Kelvin2centimeterMinusOne(x)
Kelvin2electronVolt(x)
Kelvin2Hartree(x)
Kelvin2Hertz(x)
Kelvin2Joule(x)
Kelvin2kiloCaloriePerMole(x)
Kelvin2kiloJoulePerMole(x)

kiloCaloriePerMole2centimeterMinusOne(x)
kiloCaloriePerMole2electronVolt(x)
kiloCaloriePerMole2Hartree(x)
kiloCaloriePerMole2Hertz(x)
kiloCaloriePerMole2Joule(x)
kiloCaloriePerMole2Kelvin(x)
kiloCaloriePerMole2kiloJoulePerMole(x)

kiloJoulePerMole2centimeterMinusOne(x)
kiloJoulePerMole2electronVolt(x)
kiloJoulePerMole2Hartree(x)
kiloJoulePerMole2Hertz(x)
kiloJoulePerMole2Joule(x)
kiloJoulePerMole2Kelvin(x)
kiloJoulePerMole2kiloCaloriePerMole(x)
}
\arguments{
  \item{x}{an R object containing numerical values to be converted}
}
\details{
These functions use the \code{universalConstants} data set to convert energy values from a given unit to an other. Their are first converted to Joule and then to the desired unit. With care, the functions may be applied also to otherquantities (e.g. temperatures, wavenumbers, etc...)
}
\value{Return \code{x} multiplied by conversion factors.}

\seealso{
\code{\link{universalConstants}}
}
\examples{
# Return the Planck constant in eV.s
Joule2electronVolt(universalConstants["h","Value"])

# Return the Boltzmann constant in eV.K-1
Joule2electronVolt(universalConstants["kb","Value"])

# Convert electron Volt into cm-1 (E = h.c.nu_bar)
electronVolt2centimeterMinusOne(0.2)

# Convert electron Volt into Hz (E = h.nu)
electronVolt2Hertz(0.2)

# Convert Kelvin into electron Volt (E = kb.T)
Kelvin2electronVolt(300)

}
\keyword{ manip }

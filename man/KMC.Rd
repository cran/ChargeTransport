\name{KMC}
\alias{KMC}
\title{Charge Transport Simulation}
\description{Performs a kinetic Monte Carlo (KMC) simulation to propagate a single charge carrier within a given percolation network. Either the Bortz-Kalos-Lebowitz (BKL) or the First Reaction Method (FRM) algorithm can be used to propagate the charge carrier. The function is parallelized over the number of KMC simulations by using the \code{\link{parLapply}} function of the \dQuote{parallel} package.}
\usage{
KMC(cl, con, rates, dx, dy, dz, type = "BKL", nSimu = 10, nHops = 1E7, seed = NULL)
}
\arguments{
  \item{cl}{a cluster object as returned by \code{\link{makeCluster}}.}
  \item{con}{a two-column matrix containing the \dQuote{connectivity} of the percolation network. Each row contains the labels/indexes of a pair of sites within the percolation network.}
  \item{rates}{a matrix containing the charge transfer rates for each pair of sites (columns) and each frame of a molecular dynamics (rows).}
  \item{dx, dy, dz}{matrices containing respectively the x-, y- and z-components of the inter-site distances for each pair of sites (columns) and each frame of a molecular dynamics (rows).}
  \item{type}{a character string specifying the type of KMC simulation to perform. Can be either the \code{FRM} (First Reaction Method) or \code{BKL} (Bortz-Kalos-Lebowitz) algorithm.}
  \item{nSimu}{an integer indicating the number of KMC simulations to be performed (Usually a large number). If possible it has to be a multiple of the number of nodes used for the parallelization.}
  \item{nHops}{an integer indication the number of hopping events to be performed before stopping each KMC simulation.}
  \item{seed}{an integer used to initialize the random number generator. If \code{NULL}, a random number is sampled from \code{\link{.Random.seed}}.}
}
\details{
The labels/indexes of the pair of sites defining the pathways forming the percolation network stored in \code{con} can either be integers or character strings. The dimensions of the different arguments must be consistent: \code{rates}, \code{dx}, \code{dy}, \code{dz}, must have the same dimensions and the number of rows of \code{con} must be equal to the number of columns of \code{rates}, \code{dx}, \code{dy}, \code{dz}. The first dimension (rows) of \code{rates}, \code{dx}, \code{dy}, \code{dz} allow to treat different frames of a molecular dynamics together, assuming that the connectivity of the percolation network is the same for each frame. If the connectivity slightly change during a molecular dynamics, the user can first carefully analysed the trajectory to define a suitable \sQuote{fixed} connectivity (For example, if sites 1 and 2 are neigbours at a point of the dynamics, then, the 1-2 and 2-1 rates have to be computed for all the frames). The \code{\link{dimnames}} of \code{rates} are used to set the \code{\link{dimnames}} of the different matrices/array returned by the function.
}
\value{
A list of class \sQuote{KMC} with the following components:
\item{distx, disty, distz}{matrices containing for each KMC simulation (columns) and each frame (rows) the distance cover by the charge carrier respectively along the x-, y- and z-axis.}
\item{time}{a matrix containing for each KMC simulation (columns) and each frame (rows) the drift time of the charge carrier.}
\item{nhop}{a three-dimensions array containing the number of hopping events that have occured along each pathway of the percolation network (second dimension) for each frame (first dimension) and each KMC simulation (third dimension).}
}
\references{
A. B. Boltz, M. H. Kalos, J. L. Lebowitz, \var{Journal of Computational Physics}, 17:10-18, \bold{1975}\cr
D. T. Gillespie, \var{Journal of Computational Physics} 22:403-434, \bold{1976}
}

\seealso{
\code{\link{Marcus}}, \code{\link{MarcusLevichJortner}}, \code{\link{LandauZener}}%, \code{\link{mobility}}, \code{\link{MSD}}
}
\examples{
########################################################################
## Electron transport within a periodic 1D-stack of 4 molecules:
##
##                        ... 4 | 1 2 3 4 | 1 ...
##
########################################################################
################ Preparation of the percolation network ################
nMols <- 4
## The charge can hop toward the right neigbour
con <- matrix(c(1:nMols,2:nMols,1), ncol=2)
## or the left neigbour
con <- rbind(con,con[,2:1])
## The number of pathways forming the percolation network is then:
nPaths <- 2*nMols
path.names <- paste0("P",1:nPaths)
dimnames(con) <- list(path.names, c("mol.i","mol.f"))
print(con)
########################################################################

############# Preparation of some data for the simulation ##############
## Lets consider a slight deformation of the stack along the z-axis over
## 11 frames of a molecular dynamics
nFrames <- 11
frame.names <- paste0("F", 1:nFrames)
dx <- matrix(0, ncol=nPaths, nrow=nFrames,
             dimnames=list(frame.names, path.names))
dy <- dx
# dz when hoping to the right
dz <- matrix(seq(3.0,4.0,length.out=nMols), ncol=nMols, nrow=nFrames)
# and when hoping to the left
dz <- cbind(dz, -dz)
dimnames(dz) <- dimnames(dx)
## The electronic couplings will slightly decrease
J  <- matrix(seq(10, 110, length.out=nFrames),
             ncol=nMols, nrow=nFrames)*1E-3
J  <- cbind( J, J)
dimnames(J) <- dimnames(dx)
## By symmetry all the sites have the same energy
dE <- 0
## We apply an electric field along the stack (z-axis)
Fn <- 1E5
F <- c(0,0,Fn)
## Lets consider the case of electron transport
carr <- "e"
## This introduce an addition term in the site energy differences.
dEFz <- dEField(dx,dy,dz,carr,F[1],F[2],F[3])
## Lets assum that:
lambda <- 0.14 # eV
########################################################################

############### Calculation of the charge transfer rates ###############
## Using the Marcus expression:
k <- Marcus(J,lambda,dE,dEFz)
format(k, digits=2, scientific=TRUE)
########################################################################

################### Execution of the KMC simulation ####################
cl <- makeCluster(2, outfile="") # Slave nodes output on master stdout
simuFz <- KMC(cl=cl, con=con, rates=k, dx=dx, dy=dy, dz=dz,
              type="BKL", nSimu=2, nHops=1E5)
########################################################################

###### More simulations applying the electric field along x or y #######
F <- c(Fn,0,0)
dEFx <- dEField(dx,dy,dz,carr,F[1],F[2],F[3])
k <- Marcus(J,lambda,dE,dEFx)
simuFx <- KMC(cl=cl, con=con, rates=k, dx=dx, dy=dy, dz=dz,
              type="BKL", nSimu=2, nHops=1E5)
F <- c(0,Fn,0)
dEFy <- dEField(dx,dy,dz,carr,F[1],F[2],F[3])
k <- Marcus(J,lambda,dE,dEFy)
simuFy <- KMC(cl=cl, con=con, rates=k, dx=dx, dy=dy, dz=dz,
              type="BKL", nSimu=2, nHops=1E5)              
########################################################################

################# Calculation of the mobility tensor ###################
driftVelocity <- function(simu){
  V <- c(
         x=1E-8*mean(simu$distx/simu$time),
         y=1E-8*mean(simu$disty/simu$time),
         z=1E-8*mean(simu$distz/simu$time))
  return(V) # Return the average drift velocity in cm.s-1
}
mu.Fx <- -driftVelocity(simuFx)/Fn # The minus is for electron transport
mu.Fy <- -driftVelocity(simuFy)/Fn # The minus is for electron transport
mu.Fz <- -driftVelocity(simuFz)/Fn # The minus is for electron transport
mu <- cbind(mu.Fx, mu.Fy, mu.Fz)
eigen(mu)
########################################################################
}
\keyword{ manip }

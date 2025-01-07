\name{simBMphylo}
\alias{simBMphylo}
\alias{plot.simBMphylo}
\title{Creates a graphical illustration of Brownian motion evolution on a phylogeny}
\usage{
simBMphylo(n, t, sig2, plot=TRUE, ...)
\method{plot}{simBMphylo}(x, ...)
}
\arguments{
	\item{n}{number of taxa to simulate in the output tree.}
	\item{t}{total time for the simulation.}
	\item{sig2}{the rate of evolution under Brownian motion, \eqn{\sigma^2}, or a vector of rates. If the latter the length of the vector must exactly match \code{t}, otherwise the first element of \code{sig2} will just be duplicated \code{t} times.}
	\item{plot}{optional logical value indicating whether or not the simulated object should be plotted.}
	\item{...}{optional arguments to be passed to the \code{plot} method.}
	\item{x}{in \code{plot} method, object of class \code{"simBMphylo"}.}
}
\description{
	Simulates a discrete time phylogeny and Brownian motion trait, and generates a plot.
}
\details{
	The function simulates a discrete-time pure-birth phylogeny (for fixed \emph{N} and \emph{t} using rejection sampling) and then discrete-time Brownian motion on that tree. It then proceeds to generating a plot of the results.
}
\value{
	An object of class \code{"simBMphylo"} or a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bmPlot}}, \code{\link{fastBM}}, \code{\link{pbtree}}
}
\examples{
## set seed
set.seed(777)

## create plot
simBMphylo(n=6,t=100,sig2=0.01)

## reset par to defaults
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
}
\keyword{comparative method}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{continuous character}

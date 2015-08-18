\name{sim.ratebystate}
\alias{sim.ratebystate}
\title{Conduct simulation of state dependent rate variation}
\usage{
sim.ratebystate(tree, sig2x=1, sig2y=1, beta=c(0,1), ...)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{sig2x}{variance of the Brownian process of evolution for \code{x}.}
	\item{sig2y}{variance of the Brownian process of evolution for \code{y} when \code{x-min(x)==1} (for \code{logarithm=FALSE}) or \code{x==0} (for \code{logarithm=TRUE}).}
	\item{beta}{intercept and slope of the relationship between the value of \code{x} and the Brownian rate in \code{y}.}	
	\item{...}{optional arguments which include \code{method} (\code{"by.node"} or \code{"by.branch"} indicating whether to assume the rate varies as a function of the node state or the mean branch state); \code{plot}, a logical value indicating whether or not to plot a phenogram with the branches used for simulation of \code{y} after rescaling by the state of \code{x}; and \code{logarithm}, a logical value indicating whether or not simulate changes in the variance of Brownian evolution for \code{y} as an additive \code{logarithm=FALSE} or multiplicative function of \code{x}. The default is \code{logarithm=FALSE}.}
}
\description{
	This function attempts to simulate two characters under a model in which the rate of evolution for the second (\code{y}) depends on the states for the first \code{x}. See \code{\link{ratebystate}} for more details.
}
\value{
	This function returns a matrix.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastBM}}, \code{\link{ratebystate}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}

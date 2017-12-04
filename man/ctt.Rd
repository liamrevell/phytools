\name{ctt}
\alias{ctt}
\alias{sim.ctt}
\alias{sim.multiCtt}
\title{Generates (or simulates) a 'changes through time' plot from a set of stochastic map character histories}
\usage{
ctt(trees, segments=20, ...)
sim.ctt(tree, Q, anc=NULL, nmaps=100, ...)
sim.multiCtt(tree, Q, anc=NULL, nmaps=100, nsim=100, ...)
}
\arguments{
	\item{trees}{an object of class \code{"multiSimmap"}.}
	\item{segments}{number of segments to break up the history of the tree.}
	\item{tree}{for \code{sim.ctt}, an object of class \code{"phylo"}.}
	\item{Q}{for \code{sim.ctt}, a transition matrix to use for simulation.}
	\item{anc}{ancestral state at the root node for simulation.}
	\item{nmaps}{number of stochastic maps per simulation.}
	\item{nsim}{for \code{sim.multiCtt} only, the number of simulations to run.}
	\item{...}{optional arguments.}
}
\description{
	This function generates a 'changes through time' plot in the style of a lineage-through-time (LTT) plot. It shows the mean rate or the mean number of changes per unit time from a set of stochastic character map trees.
}
\value{
	An object of class \code{"ctt"} or \code{"multiCtt"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ltt}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}

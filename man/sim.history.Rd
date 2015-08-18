\name{sim.history}
\alias{sim.history}
\title{Simulate stochastic character history under some model}
\usage{
sim.history(tree, Q, anc=NULL, nsim=1, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{Q}{a matrix containing the instantaneous transition rates between states. Note that normally this is the \emph{transpose} of the matrix produced by \code{fitDiscrete} in the geiger package or \code{\link{make.simmap}} in phytools; that is to say the transition rate from \code{i -> j} should be given by \code{Q[j,i]}. However, if your matrix is properly conformed (i.e., rows \emph{or} columns sum to 0), then \code{sim.history} will attempt to transpose your matrix correctly & will return an informative message (if \code{message=TRUE}, see below).}
	\item{anc}{an optional value for the state at the root node; if \code{NULL}, a random state will be assigned.}
	\item{nsim}{number of simulations.}
	\item{...}{other optional arguments. Currently only \code{message}, a logical value indicating whether or not to turn on informational messages (defaults to \code{message=TRUE}).}
}
\description{
	This function simulates a stochastic character history for a discretely valued character trait on the tree. The resultant tree is stored as a modified \code{"phylo"} object in stochastic character map (e.g., \code{\link{make.simmap}}) format.
}
\value{
	A modified phylogenetic tree of class \code{"phylo"} (or a modified \code{"multiPhylo"} object, for \code{nsim} > 1) with the following additional elements:
	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
	\item{states}{a vector containing the tip states.}
	\item{node.states}{a matrix containing the states at internal & terminal nodes (according to the dimensions of \code{edge}).}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}, \code{\link{sim.rates}}
}
\keyword{phylogenetics}
\keyword{simulation}

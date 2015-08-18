\name{sim.rates}
\alias{sim.rates}
\title{Brownian simulation with multiple evolutionary rates}
\usage{
sim.rates(tree, sig2, anc=0, nsim=1, internal=F, plot=F)
}
\arguments{
	\item{tree}{is a stochastic map format phylogenetic tree in modified \code{"phylo"} format (e.g., see \code{\link{make.simmap}}).}
	\item{sig2}{a named vector containing the rates for each state; names should be states in \code{mtree}.}
	\item{anc}{optional value for the root state.}
	\item{nsim}{number of simulations.}
	\item{internal}{logical value indicating whether to return states at internal nodes.}
	\item{plot}{logical value indicating whether or not to visual the rate heterogeneity (default value is \code{FALSE}.}
}
\description{
	This function conducts BM simulation on a tree with multiple rates.
}
\value{
	A vector (for \code{nsim=1}) or matrix containing the tip states for the \code{n} species in the tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastBM}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{comparative method}

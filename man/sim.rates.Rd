\name{sim.rates}
\alias{sim.rates}
\alias{multiOU}
\title{Brownian or OU simulation with multiple evolutionary regimes}
\usage{
sim.rates(tree, sig2, anc=0, nsim=1, internal=FALSE, plot=FALSE)
multiOU(tree, alpha, sig2, theta=NULL, a0=NULL, nsim=1, internal=FALSE, ...)
}
\arguments{
	\item{tree}{is a stochastic map format phylogenetic tree in modified \code{"phylo"} format (e.g., see \code{\link{make.simmap}}).}
	\item{sig2}{a named vector containing the rates for each state; names should be states in \code{mtree}.}
	\item{anc}{optional value for the root state.}
	\item{nsim}{number of simulations.}
	\item{internal}{logical value indicating whether to return states at internal nodes.}
	\item{plot}{logical value indicating whether or not to visual the rate heterogeneity (default value is \code{FALSE}.}
	\item{alpha}{single value or vector of values of the OU \eqn{\alpha} parameter.}
	\item{theta}{single value or vector of values of the OU \eqn{\theta} parameter.}
	\item{a0}{optional value of the root state. Defaults to zero.}
	\item{...}{optional arguments.}
}
\description{
	Simulates multi-rate or multi-regime continuous trait evolution on a phylogeny.
}
\details{
	The function \code{sim.rates} conducts BM simulation on a tree with multiple rates.
	
	The function \code{multiOU} conducts multi-regime OU simulations on the tree under a range of conditions. \code{multiOU} uses a difference equation approximation of the OU process.
}
\value{
	A vector (for \code{nsim=1}) or matrix containing the tip states for the \code{n} species in the tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastBM}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{comparative method}
\keyword{continuous character}

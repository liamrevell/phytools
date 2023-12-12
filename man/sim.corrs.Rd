\name{sim.corrs}
\alias{sim.corrs}
\title{Multivariate Brownian simulation with multiple correlations and rates}
\usage{
sim.corrs(tree, vcv, anc=NULL, internal=FALSE)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"} or \code{"simmap"}.}
	\item{vcv}{is a square covariance matrix or named list of matrices (one for each mapped state on the tree).}
	\item{anc}{optional vector of values for the root state.}
	\item{internal}{logical value indicating whether to return states at internal nodes.}
 }
\description{
	Simulates multivariate Brownian motion evolution on a tree with multiple evolutionary correlation/covariance matrices.
}
\details{
	This function conducts BM simulation on a tree with multiple rates and/or multiple evolutionary correlations between characters.
	
	If \code{vcv} is a single matrix, instead of a list of matrices, \code{sim.corrs} will simulate multivariate BM with a single rate matrix.
}
\value{
	A matrix containing the multivariate tip states for the \code{n} species in the tree (and nodes if \code{internal=} \code{TRUE}).
}
\references{
	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastBM}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}, \code{\link{sim.rates}}
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{comparative method}
\keyword{continuous character}

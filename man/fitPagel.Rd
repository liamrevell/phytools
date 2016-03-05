\name{fitPagel}
\alias{fitPagel}
\title{Function to test for correlated evolution of binary traits}
\usage{
fitPagel(tree, x, y, method="fitMk", ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of phenotypic values for a binary trait for the species in \code{tree}.}
	\item{y}{a second binary character for the species in \code{tree}.}
	\item{method}{function to use for optimization (defaults to \code{method="fitMk"}). Other options are \code{"ace"} to use the \code{\link{ace}} function in ape for optimization, or to \code{"fitDiscrete"} (if the geiger package is installed) to use geiger's \code{fitDiscrete} for optimization.}
	\item{...}{optional arguments to be passed to \code{\link{fitMk}}, \code{\link{ace}}, or \code{fitDiscrete}.}
}
\description{
	This function fit's Pagel's (1994) model for the correlated evolution of two binary characters.
}
\value{
	An object of class \code{"fitPagel"} which contains the optimized matrices under an independence & a dependence model, log-likelihoods, a likelihood ratio, and a P-value for the independence model based on a chi-squared test.
}
\references{
	Pagel, M. (1994) Detecting correlated evolution on phylogenies: A general method for the comparative analysis fo discrete characters. \emph{Proceedings of the Royal Society B}, \bold{255}, 37-45.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}


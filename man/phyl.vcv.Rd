\name{phyl.vcv}
\alias{phyl.vcv}
\title{Compute evolutionary VCV matrix for a tree & dataset}
\usage{
phyl.vcv(X, C, lambda)
}
\arguments{
	\item{lambda}{value for \eqn{\lambda} transformation.}
	\item{X}{data matrix.}
	\item{C}{matrix containing the height above the root of each pair of species in the tree. Typically this will have been produced by calling \code{\link{vcv.phylo}}.}
}
\description{
	Primarily an internal function for \code{\link{phyl.pca}}; this can be used to compute the phylogenetic trait variance-covariance matrix given a phylogenetic VCV, lambda, and a data matrix. Should not be confused with \code{\link{vcv.phylo}} in the "ape" package (although one of the objects returned is the output of \code{vcv.phylo}).
}
\details{
	Do not use unless you know what you're doing.
}
\value{
	A list containing three elements, as follows: \code{C}, the matrix \code{vcv.phylo} transformed by \code{lambda}; \code{R}, the among trait variance-covariance matrix for the data in \code{X}; and \code{alpha}, a vector of ancestral states at the root node of the tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}

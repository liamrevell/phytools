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
	Internal function for \code{\link{phyl.pca}}.
}
\details{
	Primarily designed as an internal function for \code{\link{phyl.pca}}; \code{phyl.vcv} can be used to compute the phylogenetic trait variance-covariance matrix given a phylogenetic VCV, \eqn{\lambda}, and a data matrix. 
	
	This function should not be confused with \code{\link{vcv.phylo}} in the \pkg{ape} package (although one of the objects returned is the output of \code{vcv.phylo}).
}
\value{
	A list containing three elements, as follows: \code{C}, the matrix \code{vcv.phylo} transformed by \code{lambda}; \code{R}, the among trait variance-covariance matrix for the data in \code{X}; and \code{alpha}, a vector of ancestral states at the root node of the tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}

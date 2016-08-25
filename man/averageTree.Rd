\name{averageTree}
\alias{averageTree}
\alias{ls.consensus}
\alias{minTreeDist}
\title{Compute an average tree from a set of trees and related operations}
\usage{
averageTree(trees, start=NULL, method="quadratic.path.difference",
	tol=1e-12, quiet=FALSE, ...)
ls.consensus(trees, start=NULL, tol=1e-12, quiet=FALSE, ...)
minTreeDist(tree, trees, method="quadratic.path.difference", ...)
}
\arguments{
	\item{trees}{object of class \code{"multiPhylo"}.}
	\item{tree}{object of class \code{"phylo"}. For \code{minTreeDist} the tree on which to find the edge lengths that minimize the distance to the phylogenies in \code{trees}.}
	\item{start}{starting tree for optimization.}
	\item{method}{distance criterion for minimization. Options are \code{"symmetric.difference"}, \code{"branch.score.difference"}, \code{"path.difference"}, and \code{"quadratic.path.difference"}.}
	\item{tol}{tolerance value for optimization.}
	\item{quiet}{logical value indicating whether to run "quietly" or not.}
	\item{...}{other arguments to be passed internally.}
}
\description{
	These functions compute average trees or consensus trees by various criteria.
}
\value{
	An object of class \code{"phylo"} with edge lengths.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}

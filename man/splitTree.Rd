\name{splitTree}
\alias{splitTree}
\title{Split tree at a point}
\usage{
splitTree(tree, split)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{split}{split encoded as a list with two elements: \code{node}: the node number tipward of the split; and \code{bp}: the position along the branch to break the tree, measured from the rootward end of the edge.}
}
\description{
	Internal function for \code{\link{posterior.evolrate}}.
}
\details{
	This function splits the tree at a given point, and returns the two subtrees as an object of class \code{"multiPhylo"}.
}
\details{
	Probably do not use this unless you can figure out what you are doing.
}
\value{
	Two trees in a list.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

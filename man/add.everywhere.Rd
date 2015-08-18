\name{add.everywhere}
\alias{add.everywhere}
\title{Add tip to all edges in a tree}
\usage{
add.everywhere(tree, tip.name)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{tip.name}{a string containing the name of the tip to add.}
}
\description{
	This function adds a tip to all branches to the tree and returns a a list of trees as an object of class \code{"multiPhylo"}.
}
\value{
	A list of trees as an object of class \code{"multiPhylo"}.  Since the tip can be added to any branch, the length of the list is equal to the number of edges in \code{tree}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\seealso{
	\code{\link{allFurcTrees}}, \code{\link{exhaustiveMP}}
}
\keyword{phylogenetics}
\keyword{utilities}

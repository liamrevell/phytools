\name{modified.Grafen}
\alias{modified.Grafen}
\alias{node.paths}
\title{Computes modified Grafen edge lengths}
\usage{
modified.Grafen(tree, power=2)
node.paths(tree, node)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{power}{power to raise the depths of each node (in nodes).}
	\item{node}{node number for \code{node.paths}.}
}
\description{
	Computes modified Grafen (1989) edge lengths.
}
\details{
	This function computes modified Grafen edge lengths in which the length of the edge is determined not by the number of descendant leaves, but instead by the maximum number of node lengths in the path from the node to any leaf.
	
	\code{node.paths} is used internally by \code{modified.Grafen} and computes the set of paths from a node to all tips descended from that node.
}
\value{
	An object of class \code{"phylo"} with edge lengths.
}
\references{	
	Grafen, A. (1989) The phylogenetic regression. \emph{Philosophical Transactions of the Royal Society of London. Series B. Biological Sciences}, \bold{326}, 119-157.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{compute.brlen}}
}
\keyword{phylogenetics}
\keyword{utilities}

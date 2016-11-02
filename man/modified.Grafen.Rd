\name{modified.Grafen}
\alias{modified.Grafen}
\alias{node.paths}
\title{Computes modified Grafen edge lengths}
\usage{
modified.Grafen(tree,power=2)
node.paths(tree,node)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{power}{power to raise the depths of each node (in nodes).}
	\item{node}{node number for \code{node.paths}.}
}
\description{
	This function computes modified Grafen edge lengths in which the length of the edge is determined not by the number of desecendant leaves, but instead by the maximum number of node lengths in the path from the node to any leaf.
	
	\code{node.paths} is used internally by \code{modified.Grafen} and computes the set of paths from a node to all tips descended from that node.
}
\value{
	An object of class \code{"phylo"} with edge lengths.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

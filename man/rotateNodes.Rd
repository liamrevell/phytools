\name{rotateNodes}
\alias{rotateNodes}
\title{Rotates a node or set of nodes in a phylogenetic tree}
\usage{
rotateNodes(tree, nodes, polytom=c(1,2), ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{nodes}{either a single node number to rotate, a vector of node numbers, or the string \code{"all"}.}
	\item{polytom}{a vector of mode numeric and length two specifying the two clades that should be exchanged in a polytomy (see \code{\link{rotate}}).}
	\item{...}{optional arguments.}
}
\description{
	This function is a wrapper for \code{\link{rotate}} which rotates a set of nodes or all nodes.
}
\details{
	Also addresses the problem that the product of multiple rotations from \code{\link{rotate}} can be non-compliant with the implicit \code{"phylo"} standard because the tip numbers in \code{tree$edge} are not in numerical order \code{1:n} for \code{n} tips.
}
\value{
	An object of class \code{"phylo"} (i.e., a phylogenetic tree).
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

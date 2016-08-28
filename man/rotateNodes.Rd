\name{rotateNodes}
\alias{rotateNodes}
\alias{rotate.multi}
\alias{allRotations}
\title{Rotates a node or set of nodes in a phylogenetic tree}
\usage{
rotateNodes(tree, nodes, polytom=c(1,2), ...)
rotate.multi(tree, node)
allRotations(tree)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{nodes}{either a single node number to rotate, a vector of node numbers, or the string \code{"all"}.}
	\item{polytom}{a vector of mode numeric and length two specifying the two clades that should be exchanged in a polytomy (see \code{\link{rotate}}).}
	\item{node}{a single node to rotate (in the case of \code{rotate.multi}).}
	\item{...}{optional arguments.}
}
\description{
	The function \code{rotateNodes} is a simple wrapper for \code{\link{rotate}} which rotates a set of nodes or all nodes.
	
	The function \code{rotate.multi} finds all possible rotations around a multifurcating node, given by \code{node}. This will be an object of class \code{"multiPhylo"}, assuming that the node specified is indeed a multifurcation.
	
	The function \code{allRotations} computes all possible rotated trees for a given input phylogeny. For a binary tree, this is generally two raised to the power of the number of internal nodes (so a very large number, if \emph{N} is even modest in size).
}
\details{
	All three functions also address the problem that the product of multiple rotations from \code{\link{rotate}} can be non-compliant with the implicit \code{"phylo"} standard because the tip numbers in \code{tree$edge} are not in numerical order \code{1:n} for \code{n} tips.
}
\value{
	An object of class \code{"phylo"} (i.e., a phylogenetic tree), in the case of \code{rotateNodes}, or an object of class \code{"multiPhylo"} for \code{rotate.multi} or \code{allRotations}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

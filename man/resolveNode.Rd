\name{resolveNode}
\alias{resolveNode}
\alias{resolveAllNodes}
\title{Compute all possible resolutions of a node or all nodes in a multifurcating tree}
\usage{
resolveNode(tree,node)
resolveAllNodes(tree)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{node}{for \code{resolveNode}, the node with a polytomy to resolve.}
}
\description{
	Resolves a single multifurcation or all multifurcations in all possible ways.
}
\details{
	This functions resolves a single multifurcation or all multifurcations in a tree in all possible ways. If the input tree has edge lengths, then the resolutions will use internal edges of zero length.

	For \code{resolveNode} applied to a multifurcation with \emph{n} descendants, the number of resolved trees will be equal to the number of possible rooted trees of \emph{n} taxa. (For instance, three for a trifurcation, 15 for a quadrifurcation, and so on.)
	
	For \code{resolveAllNodes} the number of fully resolved trees will be equal to the product of numbers for \code{resolveNode} applied to each multifurcation separately. (For instance, 45 for a tree containing one trifurcation and one quadrifurcation.)
}
\value{
	An object of class \code{"multiPhylo"} - or, if the input tree is already fully resolved, an object of class \code{"phylo"} identical to \code{tree}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{multi2di}}
}
\keyword{phylogenetics}
\keyword{utilities}

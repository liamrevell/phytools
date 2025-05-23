\name{nodeHeights}
\alias{nodeHeights}
\alias{nodeheight}
\title{Compute the heights above the root of each node}
\usage{
nodeHeights(tree, ...)
nodeheight(tree, node, ...)
}
\arguments{
	\item{tree}{a phylogeny as an object of class \code{"phylo"}.}
	\item{node}{for \code{nodeheight}, the node for which we want to compute a height above the root (or including the root edge, for \code{root.edge=TRUE}).}
	\item{...}{optional arguments - presently only \code{root.edge}, a logical value indicating whether or not to include the root edge length in the calculation of node heights.}
}
\description{
	\code{nodeHeights} computes the height above the root for all nodes in the tree. \code{nodeheight} computes the height above the root for a single node.
}
\details{
	The function \code{nodeHeights} also gives a handy way to get the total length of the tree from the root to the highest tip which will be given by \code{max(nodeHeights(tree))}.
	
	Generally speaking, \code{nodeHeights} will be faster if the heights of all or a large proportion of nodes is needed, whereas \code{nodeheight} will be faster if the height of one or a small number of nodes are needed.
}
\value{
	Either a matrix of the same dimensions as \code{tree$edge} containing the height above the root of each node in \code{edge} (for \code{nodeHeights}); or a single positive number (for \code{nodeheight}).
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{vcvPhylo}}
}
\examples{
## load tree
data(vertebrate.tree)

## compute height of all nodes
H<-nodeHeights(vertebrate.tree)
print(H)

## compute total tree depth
max(H)
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{comparative method}

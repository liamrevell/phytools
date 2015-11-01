\name{reroot}
\alias{reroot}
\title{Re-root a tree along an edge}
\usage{
reroot(tree, node.number, position, interactive=FALSE, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{node.number}{ number of the node descending from the target branch in \code{tree$edge} - this can also be a tip in which case the node number is the index number of the tip in \code{tree$tip.label}.}
	\item{position}{position along the target edge at which to re-root the tree.}.
	\item{interactive}{logical value indicating whether to use interactive mode (defaults to \code{interactive=FALSE}).}
	\item{...}{arguments to be passed to \code{plotTree} for \code{interactive=TRUE} only.}
}
\description{
	This function re-roots a phylogenetic tree at an arbitrary position along an edge.
}
\details{
	This function had an error for rootings along edges descended from the root node for phytools<=0.2-47. This should be fixed in the present version. Now uses \code{\link{paste.tree}}, \code{\link{root}}, and \code{\link{splitTree}} internally. Earlier versions also had an error related to node labels. This should be fixed in phytools>=0.4-47.
}
\value{
	A phylogenetic tree in \code{"phylo"} format.
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{splitTree}}, \code{\link{paste.tree}}, \code{\link{root}}
}
\keyword{phylogenetics}
\keyword{utilities}

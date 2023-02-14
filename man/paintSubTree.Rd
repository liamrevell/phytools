\name{paintSubTree}
\alias{paintSubTree}
\alias{paintBranches}
\title{Paint sub-trees with a discrete character}
\usage{
paintSubTree(tree, node, state, anc.state="1", stem=FALSE)
paintBranches(tree, edge, state, anc.state="1")
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"} or a modified object with mapped character traits.}
	\item{node}{an integer specifying the node number tipward of which the function should paint the derived state.}
	\item{edge}{an integer or vector of integers specifying the node or tip numbers of the edges that should be painted in \code{paintBranches}.}
	\item{state}{a string (or numeric value) specifying the state to paint on the tree tipward of \code{node}.}
	\item{anc.state}{the ancestral state to use; will only be applied if there are presently no character values mapped on the tree.}
	\item{stem}{logical or numeric value indicating whether to use the derived state on the stem leading to \code{node} (or not, if \code{stem=FALSE}), or, alternatively, what fraction of the stem should be assigned to the derived clade. Note that for tip clades \code{stem=FALSE} is not allowed.}
}
\description{
	Paints regimes on a tree to create an object of class \code{"simmap"} with mapped regimes.
}
\details{
	These functions map or "paint" arbitrary (i.e., user-specified) discrete character histories on the tree. 
	
	\code{paintSubTree} paints the clade downstream of \code{node} with a particular state; whereas \code{paintBranches} paints only a specified branch.
}
\value{
	An object of class \code{"simmap"} that contains the specified paintings as a mapped discrete character.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}

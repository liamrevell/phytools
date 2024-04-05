\name{averageTree}
\alias{averageTree}
\alias{ls.consensus}
\alias{minTreeDist}
\title{Compute an average tree from a set of trees and related operations}
\usage{
averageTree(trees, start=NULL, method="quadratic.path.difference",
   tol=1e-12, quiet=FALSE, ...)
ls.consensus(trees, start=NULL, tol=1e-12, quiet=FALSE, ...)
minTreeDist(tree, trees, method="quadratic.path.difference", ...)
}
\arguments{
	\item{trees}{object of class \code{"multiPhylo"}.}
	\item{tree}{object of class \code{"phylo"}. For \code{minTreeDist} the tree on which to find the edge lengths that minimize the distance to the phylogenies in \code{trees}.}
	\item{start}{starting tree for optimization.}
	\item{method}{distance criterion for minimization. Options are \code{"symmetric.difference"}, \code{"branch.score.difference"}, \code{"path.difference"}, and \code{"quadratic.path.difference"}.}
	\item{tol}{tolerance value for optimization.}
	\item{quiet}{logical value indicating whether to run "quietly" or not.}
	\item{...}{other arguments to be passed internally.}
}
\description{
	Compute average trees or consensus trees by various criteria.
}
\details{
	The function \code{averageTree} tries to find the (hypothetical) tree topology and branch lengths that has a minimum distance to all the trees in an input set, according to some user-specified tree distance measure.
	
	The function \code{ls.consensus} computes the least-squares consensus tree (Lapointe & Cucumel, 1997) from a set of input trees.
	
	Finally, the function \code{minTreeDist} finds the tree in the input set that is a minimum distance to all the other trees in the set. (This contrasts with \code{averageTree} which can return a tree not in the input set.)
}
\value{
	An object of class \code{"phylo"} with edge lengths.
}
\references{
	Lapointe, F.-J., G. Cucumel (1997) The average consensus procedure: Combination of weighted trees containing identical or overlapping sets of taxa. \emph{Systematic Biology}, \bold{46}, 306-312.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{consensus tree}

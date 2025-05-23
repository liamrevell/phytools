\name{map.to.singleton}
\alias{map.to.singleton}
\alias{plotTree.singletons}
\alias{drop.tip.singleton}
\alias{rootedge.to.singleton}
\title{Converts a tree without singletons to a tree with singleton nodes}
\usage{
map.to.singleton(tree)
plotTree.singletons(tree)
\method{drop.tip}{singleton}(phy, tip, ...)
rootedge.to.singleton(tree)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} (for \code{map.to.singleton},  or a tree with one or more singleton nodes (for \code{plotTree.singletons}, \code{drop.tip.singleton}, and \code{rootedge.to.singleton}).}
	\item{phy}{for \code{drop.tip.singleton}, an object of class \code{"singleton"} or \code{"phylo"}.}
	\item{tip}{for \code{drop.tip.singleton}, a tip label or vector of tip labels.}
	\item{...}{optional arguments for \code{drop.tip.singleton}.}
}
\description{
	\code{map.to.singleton} takes an object of class \code{"simmap"} with a mapped discrete character and converts it to a tree with singleton nodes, in which edge has only one state.
}	
\details{	
	The states for each edge are stored in \code{names(tree$edge.length)}. In a sense this is just an alternative way to use the general structure of the \code{"phylo"} object to encode a tree with a mapped character.	
	
	\code{plotTree.singletons} plots a tree with singleton nodes. Note that \code{\link{plotTree}} and \code{\link[ape]{plot.phylo}} now have no trouble graphing trees with singleton nodes - but they do this by just ignoring the singletons. \code{plotTree.singletons} marks the singletons as nodes on the plotted phylogeny.
	
	\code{drop.tip.singleton} drops tips from the tree leaving ancestral nodes for all remaining tips as singletons. 
	
	Finally, \code{rootedge.to.singleton} converts a tree with a root edge to a tree with a singleton node instead.
}
\value{
	An object of class \code{"phylo"} with singleton nodes. \code{plotTree.singletons} graphs a tree in which the singleton nodes are shown. 
	
	If \code{names(tree$edge.length)!=NULL} \code{plotTree.singletons} will use a different color from \code{\link{palette}} for each mapped state.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{collapse.singles}}, \code{\link[ape]{drop.tip}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{plotting}

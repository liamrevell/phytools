\name{collapse.to.star}
\alias{collapse.to.star}
\title{Collapse a subtree to a star phylogeny}
\usage{
collapse.to.star(tree, node)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{node}{node for the clade to be collapsed.}
}
\description{
	Collapses a subtree to a star.
}
\details{
	If the tree has edge lengths, the function will keep the tips at the same height above the root as in the original tree.
	
	If \code{node} is the global root of the tree a star phylogeny will be created.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{splitTree}}, \code{\link{starTree}}
}
\keyword{phylogenetics}
\keyword{utilities}

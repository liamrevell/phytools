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
	This function collapses a subtree to a star. If the tree has edge lengths, the function will keep the tips at the same height above the root as in the original tree.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{splitTree}}, \code{\link{starTree}}
}
\keyword{phylogenetics}
\keyword{utilities}

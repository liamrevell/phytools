\name{rescaleSimmap}
\alias{rescaleSimmap}
\title{Rescale SIMMAP style tree}
\usage{
rescaleSimmap(tree, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in modified \code{"phylo"} format with a discrete character mapping (e.g., see \code{\link{read.simmap}} or \code{\link{make.simmap}}).}
	\item{...}{other arguments, such as \code{depth}.}
}
\description{
	This function scales a tree with a mapped discrete character to an arbitrary total height, preserving the relative time spent in each state along each edge.
}
\details{
	Replaces \code{rescaleTree} (now \code{rescale.phylo}) in the 'geiger' package for SIMMAP style trees.
}
\value{
	A phylogenetic tree in modified \code{"phylo"} format.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}

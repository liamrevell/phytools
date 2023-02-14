\name{rescale.simmap}
\alias{rescale}
\alias{rescale.simmap}
\alias{rescaleSimmap}
\title{Rescale \code{"simmap"} object}
\usage{
rescale(x, ...)
\method{rescale}{simmap}(x, model="depth", ...)
rescaleSimmap(tree, ...)
}
\arguments{
	\item{x}{object to be rescaled (possibly an object of class \code{"simmap"}).}
	\item{model}{model to use to rescale the tree. Currently the only option is \code{"depth"}.}
	\item{...}{parameter of the model to use in rescaling. Currently the only parameter is \code{depth} for \code{model="depth"}.}
	\item{tree}{for \code{rescaleSimmap}, object of class \code{"simmap"} to be rescaled.}
}
\description{
	Scales a tree with a mapped discrete character (\code{"simmap"} object) to an arbitrary total height, preserving the relative time spent in each state along each edge.
}
\details{
	Replaces \code{rescaleTree} (now \code{rescale.phylo}) in the 'geiger' package for the \code{"simmap"} object class. \code{rescaleSimmap} is now a redundant alias for the method \code{rescale.simmap}.
}
\value{
	An object of class \code{"simmap"}.
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

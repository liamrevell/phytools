\name{untangle}
\alias{untangle}
\title{Attempts to untangle crossing branches for plotting}
\usage{
untangle(tree, method=c("reorder","read.tree"))
}
\arguments{
	\item{tree}{tree as an object of class \code{"phylo"} or \code{"simmap"}.}
	\item{method}{method to use to attempt to untangle branches. \code{method="reorder"} uses two calls of \code{\link{reorder.phylo}} or \code{\link{reorderSimmap}}; \code{method="read.tree"} writes the tree to a text string and then reads it back into memory using \code{\link{read.tree}}.}
}
\description{
	This function attempts to untangle the branches of a tree that are tangled in plotting with \code{\link{plot.phylo}}, \code{\link{plotTree}}, or \code{\link{plotSimmap}}.
}
\value{
	An object of class \code{"phylo"} or \code{"simmap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

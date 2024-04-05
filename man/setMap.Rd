\name{setMap}
\alias{setMap}
\alias{setMap.contMap}
\alias{setMap.densityMap}
\alias{setMap.phyloScattergram}
\alias{setMap.multirateBM_plot}
\title{Set color map for various phylogenetic objects of classes}
\usage{
setMap(x, ...)
\method{setMap}{contMap}(x, ...)
\method{setMap}{densityMap}(x, ...)
\method{setMap}{phyloScattergram}(x, ...)
\method{setMap}{multirateBM_plot}(x, ...)
}
\arguments{
	\item{x}{an object of class \code{"contMap"}, \code{"densityMap"}, \code{"phyloScattergram"}, or \code{"multirateBM_plot"}.}
	\item{...}{arguments to be passed to \code{\link{colorRampPalette}}. Also, the argument \code{invert} which (if \code{invert=TRUE}) will just flip the current color ramp.}
}
\description{
	Changes the color map (ramp) in an object of class \code{"contMap"}, \code{"densityMap"}, \code{"phyloScattergram"}, or \code{"multirateBM_plot"}.
}
\value{
	An object of class \code{"contMap"}, \code{"densityMap"}, \code{"phyloScattergram"}, or \code{"multirateBM_plot"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. (2013) Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
	
	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{multirateBM}}
}
\keyword{phylogenetics}
\keyword{plotting}


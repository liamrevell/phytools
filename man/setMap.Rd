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
	Function to change the color map (ramp) in an object of class \code{"contMap"}, \code{"densityMap"}, \code{"phyloScattergram"}, or \code{"multirateBM_plot"}.
}
\value{
	An object of class \code{"contMap"}, \code{"densityMap"}, \code{"phyloScattergram"}, or \code{"multirateBM_plot"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. 2013. Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{multirateBM}}
}
\keyword{phylogenetics}
\keyword{plotting}


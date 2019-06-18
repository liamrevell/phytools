\name{phylo.to.map}
\alias{phylo.to.map}
\alias{plot.phylo.to.map}
\title{Plot tree with tips linked to geographic coordinates}
\usage{
phylo.to.map(tree, coords, rotate=TRUE, ...)
\method{plot}{phylo.to.map}(x, type=c("phylogram","direct"), ...)
}
\arguments{
	\item{tree}{an object of class "phylo".}
	\item{coords}{a matrix containing the latitude (in column 1) and the longitude of all tip species in the tree. The row names should be the same as \code{tree$tip.label}; however, more than one set of coordinates per species can be supplied by duplicating some row names.}
	\item{rotate}{a logical value indicating whether or not to rotate nodes of the tree to better match longitudinal positions.}
	\item{x}{for \code{plot.phylo.to.map}, an object of class \code{"phylo.to.map"}.}
	\item{type}{a string indicating whether to map the tips of the tree onto a geographic map from a square phylogram (\code{type="phylogram"}) or to project the tree directly onto the map (\code{type="direct"}).}
	\item{...}{various optional arguments. For the function \code{phylo.to.map}, which first creates an object of the special class \code{"phylo.to.map"} and then (optionally) plots this object, arguments include: \code{database} and \code{regions} (see \code{\link{map}}), as well as any arguments that should be passed to \code{plot.phylo.to.map} internally. For \code{phylo.to.map}, optional arguments \code{xlim} and \code{ylim}, which control the plot area for the map; \code{fsize} for the font size of plot labels and \code{ftype} for the font type (following \code{\link{plotSimmap}}; \code{split} which controls the proportion of vertical (or horizontal) space for the tree (first) and map, in a vector; \code{psize} the size of the plotted points on the map - or \code{cex.points}, a vector contain the size of the tip points and geographic coordinate points, respectively; \code{from.tip} a logical value indicating whether to plot the linking lines from the tips (if \code{TRUE}) or from the end of the tip label, the default; \code{colors}, a single value or a vector of colors for the points and the linking lines; \code{pch} a single value or a vector of point types; \code{lwd} and \code{lty} for the linking lines; and \code{pts} a logical value indicating whether or not to plot points at the tips of the tree. \code{mar} and \code{asp} are as in \code{\link{par}}.}
}
\description{
	\code{phylo.to.map} creates an object of class \code{"phylo.to.map"} and (optionally) plots that object.

	\code{plot.phylo.to.map} plots an object of class \code{"phylo.to.map"} in which the tips of the tree point to coordinates on a geographic map.
}
\value{
	\code{phylo.to.map} creates an object of class \code{"phylo.to.map"} and  (if \code{plot=TRUE}) plots a phylogeny projected onto a geographic map.
	
	\code{plot.phylo.to.map} plots on object of class \code{"phylo.to.map"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}


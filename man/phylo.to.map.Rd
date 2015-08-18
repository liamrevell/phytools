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
	\item{coords}{a matrix containing the latitude (in column 1) and the longitude of all tip species in the tree. The row names should be the same as \code{tree$tip.label}.}
	\item{rotate}{a logical value indicating whether or not to rotate nodes of the tree to better match longitudinal positions.}
	\item{x}{for \code{plot.phylo.to.map}, an object of class \code{"phylo.to.map"}.}
	\item{type}{a string indicating whether to map the tips of the tree onto a geographic map from a square phylogram (\code{type="phylogram"}) or to project the tree directly onto the map (\code{type="direct"}).}
	\item{...}{optional arguments. For \code{phylo.to.map}, which creates an object of class \code{"phylo.to.map"} and (optionally) plots that object, optional arguments include: \code{database} and \code{regions} (see \code{\link{map}}); \code{xlim} and \code{ylim}, which control the plot area for the map; \code{fsize} for the font size of plot labels; \code{split} which controls the proportion of vertical space for the tree (first) and map, in a vector; \code{psize} the size of the plotted points on the map; \code{mar} and \code{asp} (see \code{\link{par}}. For \code{plot.phylo.to.map}, the options are the same - excluding \code{database} and \code{regions}.}
}
\description{
	Function plots a tree and tips pointing to coordinates on a global map.
}
\value{
	Plots a phylogeny and map and returns an object of class \code{"phylo.to.map"} invisibly.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}


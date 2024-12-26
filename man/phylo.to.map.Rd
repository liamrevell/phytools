\name{phylo.to.map}
\alias{phylo.to.map}
\alias{plot.phylo.to.map}
\title{Plot tree with tips linked to geographic coordinates}
\usage{
phylo.to.map(tree, coords, rotate=TRUE, ...)
\method{plot}{phylo.to.map}(x, type=c("phylogram","direct"), ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{coords}{a matrix containing the latitude (in column 1) and the longitude of all tip species in the tree. The row names should be the same as \code{tree$tip.label}; however, more than one set of coordinates per species can be supplied by duplicating some row names.}
	\item{rotate}{a logical value indicating whether or not to rotate nodes of the tree to better match longitudinal positions.}
	\item{x}{for \code{plot.phylo.to.map}, an object of class \code{"phylo.to.map"}.}
	\item{type}{a string indicating whether to map the tips of the tree onto a geographic map from a square phylogram (\code{type="phylogram"}) or to project the tree directly onto the map (\code{type="direct"}).}
	\item{...}{various optional arguments. For the function \code{phylo.to.map}, which first creates an object of the special class \code{"phylo.to.map"} and then (optionally) plots this object, arguments include: \code{database} and \code{regions} (see \code{\link[maps]{map}}), as well as any arguments that should be passed to \code{plot.phylo.to.map} internally. For \code{phylo.to.map}, optional arguments \code{xlim} and \code{ylim}, which control the plot area for the map; \code{fsize} for the font size of plot labels and \code{ftype} for the font type (following \code{\link{plotSimmap}}; \code{split} which controls the proportion of vertical (or horizontal) space for the tree (first) and map, in a vector; \code{psize} the size of the plotted points on the map - or \code{cex.points}, a vector contain the size of the tip points and geographic coordinate points, respectively; \code{from.tip} a logical value indicating whether to plot the linking lines from the tips (if \code{TRUE}) or from the end of the tip label, the default; \code{colors}, a single value or a vector of colors for the points and the linking lines; \code{pch} a single value or a vector of point types; \code{lwd} and \code{lty} for the linking lines; and \code{pts} a logical value indicating whether or not to plot points at the tips of the tree. \code{mar} and \code{asp} are as in \code{\link{par}}.}
}
\description{
	Project a phylogeny on a geographic map.
}
\details{
	\code{phylo.to.map} creates an object of class \code{"phylo.to.map"} and (optionally) plots that object.

	\code{plot.phylo.to.map} plots an object of class \code{"phylo.to.map"} in which the tips of the tree point to coordinates on a geographic map.
}
\value{
	\code{phylo.to.map} creates an object of class \code{"phylo.to.map"} and  (if \code{plot=TRUE}) plots a phylogeny projected onto a geographic map.
	
	\code{plot.phylo.to.map} plots on object of class \code{"phylo.to.map"}.
}
\references{
	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\examples{
## generally recommend using higher resolution map
## e.g., from mapdata package
data(tortoise.tree)
data(tortoise.geog)
tortoise.phymap<-phylo.to.map(tortoise.tree,
  tortoise.geog,plot=FALSE,direction="rightwards",
  regions="Ecuador")
plot(tortoise.phymap,direction="rightwards",pts=FALSE,
  xlim=c(-92.25,-89.25),ylim=c(-1.8,0.75),ftype="i",
  fsize=0.8,lty="dashed",map.bg="lightgreen",
  colors="slategrey")
## reset margins
par(mar=c(5.1,4.1,4.1,2.1))
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{biogeography}

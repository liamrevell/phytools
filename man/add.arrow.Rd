\name{add.arrow}
\alias{add.arrow}
\title{Add an arrow pointing to a tip or node on the tree}
\usage{
add.arrow(tree=NULL, tip, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}, \code{"contMap"}, or \code{"densityMap"}. If not supplied, the function will obtain the last plotted phylogeny from the environmental variable \code{last_plot.phylo}.}
	\item{tip}{label of tip or tip or node number. If \code{tree=NULL} then the tip or node number must be supplied.}
	\item{...}{optional arguments to control the shape and size of the arrow including: its length (\code{arrl}) in the units of the plot; the length of the arrowhead (\code{hedl}); the total angle between the wings in the arrowhead (\code{angle}); the line width for the plotted lines (\code{lwd}); the offset from the tip or end of tip label, in character widths (\code{offset}); and the color (\code{col}).}
}
\description{
	This function adds an arrow to a plotted tree.
}
\details{
	This function presently works for radial (\code{type="fan"}) and right facing square phylograms (\code{type="phylogram"}). Trees can be plotted using phytools function \code{plotTree}, \code{plotSimmap}, \code{contMap}, \code{densityMap}, and ape S3 method \code{plot.phylo}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{nodelabels}}
}
\keyword{phylogenetics}
\keyword{utilities}

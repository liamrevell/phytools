\name{phylo.heatmap}
\alias{phylo.heatmap}
\title{Creates a phylogenetic heat map}
\usage{
phylo.heatmap(tree, X, fsize=1, colors=NULL, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix containing data for multiple continuous characters in which \code{rownames} correspond to the tip labels of the tree.}
	\item{fsize}{an integer or vector of length 3 containing the font size for the tip labels, the trait labels, and the legend text. (If a single integer is supplied, then the value will be recycled.}
	\item{colors}{a vector of colors to be passed to \code{\link{image}}. Can be a function call (e.g., \code{heat.colors(n=200)[200:1]).}
	\item{...}{optional arguments.}
}
\description{
	Functions creates a multivariate phylogenetic \code{\link{heatmap}}.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}

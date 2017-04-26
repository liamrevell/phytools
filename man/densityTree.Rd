\name{densityTree}
\alias{densityTree}
\alias{make.transparent}
\title{Plots a posterior sample of trees}
\usage{
densityTree(trees, colors="blue", alpha=NULL, method="plotTree", 
   fix.depth=FALSE, use.edge.length=TRUE, compute.consensus=TRUE, 
   use.gradient=FALSE, show.axis=TRUE, ...)
make.transparent(color, alpha)
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"} or \code{"multiSimmap"}.}
	\item{colors}{a color or a named vector of colors in which names correspond to mapped states in an object of class \code{"multiSimmap"}.}
	\item{alpha}{transparency level for plotted trees which is passed to internally used function, \code{make.transparent}. (\code{0} is fully transparent, which \code{1} is fully opaque.) By default will be one divided by the number of trees.}
	\item{method}{plotting method to be used internally. Can be \code{"plotTree"} or \code{"plotSimmap"}.}
	\item{fix.depth}{logical value indicating whether or not to plot trees with a fixed depth or to permit plotted trees to have different depths.}
	\item{use.edge.length}{logical value indicating whether to use the edge lengths of the input tree. Defaults to \code{use.edge.length=TRUE} unless any input tree edge lengths are \code{NULL}.}
	\item{compute.consensus}{logical value indicating whether or not to use the tip order from a consensus tree. (Defaults to \code{compute.consensus=TRUE} Defaulted to \code{FALSE} in earlier version of this function.)}
	\item{use.gradient}{logical value indicating whether to plot all trees slightly offset using a rainbow color gradient. (Defaults to \code{use.gradient=FALSE}.)}
	\item{show.axis}{logical value indicating whether or not to include a horizontal axis in the plot.}
	\item{...}{arguments to be passed to \code{plotTree} or \code{plotSimmap}. Some may be ignored if they are incompatible with the method.}
	\item{color}{in \code{make.transparent}, the color to render transparent.}
}
\description{
	Functions plots a posterior sample of trees, including with mapped discrete characters.
	
	\code{make.transparent} is used internally and converts a color to transparent with a certain user-specified \code{alpha} level.
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

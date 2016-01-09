\name{density.tree}
\alias{density.tree}
\title{Plots a posterior sample of trees}
\usage{
density.tree(trees, colors="blue", alpha=NULL, method="plotTree", fix.depth=FALSE,
	use.edge.length=TRUE, ...)
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"} or \code{"multiSimmap"}.}
	\item{colors}{a color or a named vector of colors in which names correspond to mapped states in an object of class \code{"multiSimmap"}.}
	\item{alpha}{transparency level for plotted trees. By default will be one divided by the number of trees.}
	\item{method}{plotting method to be used internally. Can be \code{"plotTree"} or \code{"plotSimmap"}.}
	\item{fix.depth}{logical value indicating whether or not to plot trees with a fixed depth or to permit plotted trees to have different depths.}
	\item{use.edge.length}{logical value indicating whether to use the edge lengths of the input tree. Defaults to \code{use.edge.length=TRUE} unless any input tree edge lengths are \code{NULL}.}
	\item{...}{arguments to be passed to \code{plotTree} or \code{plotSimmap}.}
}
\description{
	Functions plots a posterior sample of trees, including with mapped discrete characters.
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

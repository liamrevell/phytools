\name{plotTree.errorbars}
\alias{plotTree.errorbars}
\title{Plot a tree with error bars around divergence dates}
\usage{
plotTree.errorbars(tree, CI, ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{CI}{confidence intervals around internal nodes of the tree, measured in time since the present.}
	\item{...}{optional arguments to be passed to \code{\link{plotTree}}.}
}
\description{
	Plots a tree with error bars around divergence times (nodes).
}
\details{
	The matrix CI show contain (in rows) the lower & upper confidence bounds in time since the present.
	
	Optional arguments specific to the error bar plot include <code>gridlines</code>, <code>bar.lwd</code>,
	<code>cex</code> (for the points plotted at nodes), and <code>bar.col</code>.
}
\value{
	Plots a tree with error bars around internal nodes..
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ltt}}, \code{\link{plotTree}}
}
\keyword{phylogenetics}
\keyword{plotting}

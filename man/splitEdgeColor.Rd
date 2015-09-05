\name{splitEdgeColor}
\alias{splitEdgeColor}
\title{Split edge colors when descendant edges have different mapped states}
\usage{
splitEdgeColor(tree,colors,lwd=2)
}
\arguments{
	\item{tree}{object of class \code{"simmap"}.}
	\item{colors}{named vector of colors to be used for plotting.}
	\item{lwd}{width of the plotted lines.}
}
\description{
	This function splits the vertical line colors to match the daughter edges when the daughters have different states. Mostly to be used with trees generated using \code{\link{paintBranches}} or \code{\link{paintSubTree}}. Also used internally by \code{\link{plotSimmap}}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}

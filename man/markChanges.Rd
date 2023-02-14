\name{markChanges}
\alias{markChanges}
\title{Add marked changes to a plotted tree with mapped discrete character}
\usage{
markChanges(tree, colors=NULL, cex=1, lwd=2, plot=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}.}
	\item{colors}{a named vector of colors used to plot the stochastically mapped character on the tree.}
	\item{cex}{expansion factor for line height.}
	\item{lwd}{line width.}
	\item{plot}{logical value indicating whether the changes should be plotted or not.}
}
\description{
	Adds the reconstructed changes to a plotted tree with a stochastically mapped discrete character.
}
\value{
	This function returns (invisibly) a matrix containing the x & y coordinates of the marked changes on the plotted tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}
}
\keyword{ancestral states}
\keyword{bayesian}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{discrete character}

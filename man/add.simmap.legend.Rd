\name{add.simmap.legend}
\alias{add.simmap.legend}
\title{Add legend to stochastically mapped tree}
\usage{
add.simmap.legend(leg=NULL, colors, prompt=TRUE, vertical=TRUE, ...)
}
\arguments{
	\item{leg}{states for the discrete character in the order of \code{colors}.}
	\item{colors}{colors for the legend in the order of \code{leg}, or, if \code{leg=NULL}, named vector of colors in which \code{names(colors} are the states of the mapped discrete character.}
	\item{prompt}{logical value indicating whether the location of the legend should be obtained interactively (i.e., by clicking in the plotting area).}
	\item{vertical}{logical value indiciating whether to plot the legend vertically (if \code{TRUE}) or horizontally.}
	\item{...}{optional arguments including: \code{x} x-coordinate of the legend (if \code{prompt=FALSE}); \code{y} y-coordinate of the legend; and \code{shape} which can be \code{shape="square"}, the default, or \code{shape="circle"}.}
}
\description{
	This function adds a legend (by default, interactively) to a plotted stochastic character mapped tree.
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

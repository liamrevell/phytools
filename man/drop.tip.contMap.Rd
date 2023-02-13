\name{drop.tip.contMap}
\alias{drop.tip.contMap}
\alias{drop.tip.densityMap}
\alias{keep.tip.contMap}
\alias{keep.tip.densityMap}
\title{Drop tip or tips from an object of class "contMap" or "densityMap"}
\usage{
\method{drop.tip}{contMap}(phy, tip, ...)
\method{drop.tip}{densityMap}(phy, tip, ...)
\method{keep.tip}{contMap}(phy, tip, ...)
\method{keep.tip}{densityMap}(phy, tip, ...)
}
\arguments{
	\item{phy}{an object of class \code{"contMap"} or \code{"densityMap"}.}
	\item{tip}{name or names of species to be dropped or kept.}
	\item{...}{optional arguments to be passed to \code{\link{drop.tip.simmap}}.}
}
\description{
	Drops one or multiple tips from an object of class \code{"contMap"} or \code{"densityMap"}.
}
\details{
	These functions are equivalent to \code{\link{drop.tip}} and \code{\link{keep.tip}} in the \emph{ape} package, but for objects of class \code{"contMap"} and \code{"densityMap"}.

	For more information about objects of class \code{"contMap"} or \code{"densityMap"}, please refer to the documentation pages for \code{\link{contMap}} or \code{\link{densityMap}}, respectively.
}
\value{
	An object of class \code{"contMap"} or \code{"densityMap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{drop.tip}}, \code{\link{drop.tip.simmap}}, \code{\link{keep.tip}}
}
\keyword{phylogenetics}
\keyword{utilities}

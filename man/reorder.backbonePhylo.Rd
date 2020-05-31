\name{reorder.backbonePhylo}
\alias{reorder.backbonePhylo}
\title{Reorders a backbone phylogeny}
\usage{
\method{reorder}{backbonePhylo}(x, order="cladewise", ...)
}
\arguments{
	\item{x}{an object of class \code{"backbonePhylo"}.}
	\item{order}{order. See \code{\link{reorder.phylo}} for possible orderings.}
	\item{...}{optional arguments.}
}
\description{
	Function reorders an object of class \code{"backbonePhylo"}.
}
\value{
	An object of class \code{"backbonePhylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.toBackbone}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{utilities}

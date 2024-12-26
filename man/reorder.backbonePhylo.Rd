\name{reorder.backbonePhylo}
\alias{reorder.backbonePhylo}
\title{Reorders a backbone phylogeny}
\usage{
\method{reorder}{backbonePhylo}(x, order="cladewise", ...)
}
\arguments{
	\item{x}{an object of class \code{"backbonePhylo"}.}
	\item{order}{order. See \code{\link[ape]{reorder.phylo}} for possible orderings.}
	\item{...}{optional arguments.}
}
\description{
	Function reorders an object of class \code{"backbonePhylo"}.
}
\value{
	An object of class \code{"backbonePhylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.toBackbone}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{utilities}

\name{as.multiPhylo}
\alias{as.multiPhylo}
\alias{as.multiPhylo.multiSimmap}
\alias{as.multiPhylo.phylo}
\title{Conversion to object of class \code{"multiPhylo"}}
\usage{
as.multiPhylo(x, ...)
\method{as.multiPhylo}{multiSimmap}(x, ...)
\method{as.multiPhylo}{phylo}(x, ...)
}
\arguments{
	\item{x}{object to be converted to \code{"multiPhylo"}. Presently an object of class \code{"multiSimmap"}, or an object of class \code{"phylo"}. In the latter case an object of class \code{"multiPhylo"} with length \code{1} is generated.}
	\item{...}{optional arguments.}
}
\description{
	This function converts between object classes.
}
\value{
	An object of class \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

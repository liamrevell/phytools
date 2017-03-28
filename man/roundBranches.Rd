\name{roundBranches}
\alias{roundBranches}
\title{Rounds the branch lengths of a tree}
\usage{
roundBranches(tree, digits)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}, \code{"multiPhylo"}, \code{"simmap"}, or \code{"multiSimmap"}.}
	\item{digits}{number of digits for rounding. Passed to \code{\link{round}}.}
}
\description{
	This function rounds the branch lengths of a tree or trees and reconciles any mappings for objects of class \code{"simmap"} or \code{"multiSimmap"}.
}
\value{
	An object of class \code{"phylo"}, \code{"multiPhylo"}, \code{"simmap"}, or \code{"multiSimmap"}, with rounded edge lengths.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

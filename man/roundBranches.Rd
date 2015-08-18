\name{roundBranches}
\alias{roundBranches}
\title{Rounds the branch lengths of a tree}
\usage{
roundBranches(tree, digits)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{digits}{number of digits for rounding. Passed to \code{\link{round}}.}
}
\description{
	This function rounds the branch lengths of a tree or trees and reconciles any mappings (as in \code{\link{read.simmap}}) with the rounded branch lengths.
}
\value{
	A tree with branch lengths, or modified \code{"phylo"} or \code{"multiPhylo"} object with a mapped discrete character.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

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
	Rounds the branch lengths of a phylogenetic tree.

}
\details{
	This function rounds the branch lengths of a tree or trees to a precision indicated by \code{digits}, and reconciles any mappings for objects of class \code{"simmap"} or \code{"multiSimmap"}.
}
\value{
	An object of class \code{"phylo"}, \code{"multiPhylo"}, \code{"simmap"}, or \code{"multiSimmap"}, with rounded edge lengths.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

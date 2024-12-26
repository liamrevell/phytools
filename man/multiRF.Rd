\name{multiRF}
\alias{multiRF}
\title{Computes Robinson-Foulds distance between a set of trees}
\usage{
multiRF(trees,quiet=FALSE,multi2di=FALSE)
}
\arguments{
	\item{trees}{object of class \code{"multiPhylo"} consisting of two or more fully bifurcating, unrooted trees. If trees are rooted, they will be unrooted.}
	\item{quiet}{logical argument indicating whether or not to run quietly. (Defaults to \code{FALSE}.)}
	\item{multi2di}{logical argument indicating whether or not to resolve multifurcating trees. (Defaults to \code{FALSE}.)}
}
\description{
	Computes the Robinson-Foulds (Robinson & Foulds 1981) distance between a set of trees in an object of class \code{"multiPhylo"}.
}
\details{
	Computes the Robinson-Foulds distance between all phylogenies in an object of class \code{"multiPhylo"}. Uses \code{\link[ape]{prop.part}} internally for most of the heavy lifting.
}
\value{
	A matrix containing distances.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Robinson, D. R., Foulds, L. R. (1981) Comparison of phylogenetic trees. \emph{Mathematical Biosciences}, \bold{53}, 131-147.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{phylogeny inference}

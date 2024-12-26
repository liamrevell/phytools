\name{multiC}
\alias{multiC}
\title{Returns a list with phylogenetic VCV matrix for each mapped state}
\usage{
multiC(tree, internal=FALSE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} consisting of a phylogeny with a mapped discrete character.}
	\item{internal}{logical value indicating whether or not internal nodes should be returned.}
}
\description{
	Computes phylogenetic covariance matrices from a \code{"simmap"} object.
}
\details{
	This function takes a modified \code{"phylo"} object as input and returns a set of so-called phylogenetic covariance matrices (e.g., see \code{\link[ape]{vcv.phylo}}) as a list: one for each mapped state.
	
	Used internally by multiple \pkg{phytools} functions, such as \code{\link{brownie.lite}}.
}
\value{
	A list of matrices.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evolvcv.lite}}, \code{\link{read.simmap}}, \code{\link{vcvPhylo}}, \code{\link[ape]{vcv.phylo}}
}
\keyword{phylogenetics}
\keyword{utilities}

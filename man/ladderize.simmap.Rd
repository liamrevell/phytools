\name{ladderize.simmap}
\alias{ladderize.simmap}
\title{Ladderize a tree with a mapped discrete character}
\usage{
ladderize.simmap(tree, right=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"} with a mapped discrete character.}
	\item{right}{a logical specifying how the tree should be ladderized.}
}
\description{
	This function "ladderizes" an object of class \code{"phylo"} with a mapped discrete character. For more information, please see \code{\link{ladderize}}.
}
\value{
	An object of class \code{"phylo"} with a mapped discrete character.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{ladderize}}
}
\keyword{phylogenetics}
\keyword{utilities}

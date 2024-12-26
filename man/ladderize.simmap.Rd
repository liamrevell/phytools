\name{ladderize.simmap}
\alias{ladderize.simmap}
\title{Ladderize a tree with a mapped discrete character}
\usage{
ladderize.simmap(tree, right=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}.}
	\item{right}{a logical specifying how the tree should be ladderized.}
}
\description{
	Ladderizes an object of class \code{"simmap"}.
}
\details{
	This function 'ladderizes' an object of class \code{"simmap"} with a mapped discrete character.
	
	For more information see \code{\link[ape]{ladderize}}.
}
\value{
	A ladderized object of class \code{"simmap"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link[ape]{ladderize}}
}
\keyword{phylogenetics}
\keyword{utilities}

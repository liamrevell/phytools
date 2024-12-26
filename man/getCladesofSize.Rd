\name{getCladesofSize}
\alias{getCladesofSize}
\title{Get all subtrees larger than or equal to a specified size}
\usage{
getCladesofSize(tree, clade.size=2)
}
\arguments{
	\item{tree}{is an object of class \code{"phylo"}.}
	\item{clade.size}{subtree size.}
}
\description{
	This function gets all subtrees that cannot be further subdivided into two reciprocally monophyletic subtrees of size \code{>= clade.size}.
}
\value{
	An object of class \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{extract.clade}}, \code{\link{getDescendants}}
}
\keyword{phylogenetics}
\keyword{utilities}

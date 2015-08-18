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
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{extract.clade}}, \code{\link{getDescendants}}
}
\keyword{phylogenetics}
\keyword{utilities}

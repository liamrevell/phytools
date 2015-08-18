\name{orderMappedEdge}
\alias{orderMappedEdge}
\title{Order the columns of mapped.edge to match across trees}
\usage{
orderMappedEdge(trees, ordering=NULL)
}
\arguments{
	\item{trees}{object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{ordering}{ordering for the columns of \code{$mapped.edge}. If \code{NULL}, then an alphabetical order is assumed. Options are \code{"alphabetical"}, \code{"numerical"}, or any specific ordering of the mapped traits (e.g., \code{c("A","B","C")}.}
}
\description{
	This function takes a modified object of class \code{"multiPhylo"} with a mapped discrete character (e.g., see \code{\link{read.simmap}} and sorts the columns of each \code{tree$mapped.edge} to have the same state ordering. This is handy if we want to, for instance, run \code{brownie.lite} on a set of mapped trees, and then average the fitted parameter values across trees. The function also works for a single tree.
}
\value{
	A modified object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

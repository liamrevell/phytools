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
	Orders the levels of a mapped character to match across trees in a \code{"multiSimmap"} object.
}
\details{
	This function takes a an object of class \code{"multiSimmap"} with a mapped discrete character (e.g., see \code{\link{make.simmap}} and sorts the columns of each \code{tree$mapped.edge} element to have the same state ordering. 
	
	This is handy if we want to, for instance, run \code{brownie.lite} on a set of mapped trees, and then average the fitted parameter values across trees. The function also works for a single tree.
}
\value{
	An object of class \code{"simmap"} or (normally) \code{"multiSimmap"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

\name{drop.tip.simmap}
\alias{drop.tip.simmap}
\alias{keep.tip.simmap}
\alias{extract.clade.simmap}
\title{Drop tips or extract clade from tree with mapped discrete character}
\usage{
\method{drop.tip}{simmap}(phy, tip, ...)
\method{keep.tip}{simmap}(phy, tip, ...)
extract.clade.simmap(tree, node)
}
\arguments{
	\item{phy}{an object of class \code{"simmap"}.}
	\item{tip}{name or names of species to be dropped (or kept).}
	\item{node}{node number for the root node of the clade to be extracted.}
	\item{...}{optional arguments. Currently the logical argument \code{untangle} which if set to \code{TRUE} will call \code{\link{untangle}} before returning the \code{"simmap"} object to the user.}
	\item{tree}{for \code{extract.clade.simmap}, an object of class \code{"simmap"}.}
}
\description{
	This function drops one or multiple tips from an object of class \code{"simmap"}.
}
\details{
	Equivalent to \code{\link{drop.tip}} and \code{\link{keep.tip}} but for a tree with a mapped discrete character.

	\code{extract.clade.simmap} is functionally equivalent to \code{\link{extract.clade}} but preserves discrete character mappings on the tree.
}
\value{
	A modified object of class \code{"phylo"} containing the elements \code{maps} and \code{$mapped.edge} with the time spent in each state along each edge of the tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{drop.tip}}, \code{\link{extract.clade}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{utilities}

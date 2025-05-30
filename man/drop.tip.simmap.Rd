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
	Equivalent to \code{\link[ape]{drop.tip}} and \code{\link[ape]{keep.tip}} but for a tree with a mapped discrete character.

	\code{extract.clade.simmap} is functionally equivalent to \code{\link[ape]{extract.clade}} but preserves discrete character mappings on the tree.
	
	Following \code{\link[ape]{drop.tip}} in \pkg{ape}, the returned tree is always in \code{"cladewise"} order.
}
\value{
	A modified object of class \code{"phylo"} containing the elements \code{maps} and \code{$mapped.edge} with the time spent in each state along each edge of the tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{drop.tip}}, \code{\link[ape]{extract.clade}}, \code{\link[ape]{keep.tip}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{sim.history}}
}
\keyword{phylogenetics}
\keyword{utilities}

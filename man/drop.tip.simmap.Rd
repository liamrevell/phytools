\name{drop.tip.simmap}
\alias{drop.tip.simmap}
\alias{extract.clade.simmap}
\title{Drop tips or extract clade from tree with mapped discrete character}
\usage{
drop.tip.simmap(tree, tip)
extract.clade.simmap(tree, node)
}
\arguments{
	\item{tree}{a modified object of class \code{"phylo"} (see \code{\link{read.simmap}}).}
	\item{tip}{name or names of species to be dropped.}
	\item{node}{node number for the root node of the clade to be extracted.}
}
\description{
	This function drops one or multiple tips from the modified \code{"phylo"} object with a mapped binary or multistate trait (see \code{\link{read.simmap}}) while maintaining the matrix \code{$mapped.edge} and list of mappings by branch \code{maps}. This function is equivalent to \code{\link{drop.tip}} but for a tree with a mapped discrete character.

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

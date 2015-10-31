\name{bind.tip}
\alias{bind.tip}
\title{Attaches a new tip to a tree}
\usage{
bind.tip(tree, tip.label, edge.length=NULL, where=NULL, position=0,
	interactive=FALSE, ...)
}
\arguments{
	\item{tree}{receptor tree.}
	\item{tip.label}{a string containing the species name for the new tip.}
	\item{edge.length}{edge length for the new tip (a scalar).}
	\item{where}{node number to attach new tip. If \code{position>0} then then tip will be attached \emph{below} the specified node. Node numbers can also be tips, in which case the new tip will be added along the terminal edge. To find out the tip number for given species with name \emph{"species"} type: \code{which(tree$tip.label=="species"}.}
	\item{position}{distance \emph{below} node to add tip.}
	\item{interactive}{logical value indicating whether or not the species should be added interactively. (Defaults to \code{FALSE}.)}
	\item{...}{arguments to be passed to \code{plotTree} (for \code{interactive=TRUE}.)}
}
\description{
	Functions adds a new tip to the tree. If the tree is ultrametric and no branch length is specified, then \code{edge.length} is scaled so that the tree remains ultrametric after the new tip is added.
}
\details{
	Wrapper function for 'ape' \code{\link{bind.tree}}. Note that \code{interactive=TRUE} works only for right-facing phylograms.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

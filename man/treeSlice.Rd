\name{treeSlice}
\alias{treeSlice}
\title{Slices the tree at a particular point and returns all subtrees, or the tree rootward of the point}
\usage{
treeSlice(tree, slice, trivial=FALSE, prompt=FALSE, ...)
}
\arguments{
	\item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
	\item{slice}{a real number indicating the height above the root at which to slice the tree.}
	\item{trivial}{a logical value indicating whether or not to return subtrees with a number of tips less than two (default is \code{FALSE}).}
	\item{prompt}{logical value indicating whether or not the height of the slice should be given interactively.}
	\item{...}{for \code{prompt=TRUE}, other arguments to be passed to \code{\link{plotTree}}. In addition, the argument \code{orientation} can be used to specify whether the \code{"tipwards"} subtrees or the \code{"rootwards"} phylogeny are/is to be returned by the function call (using those two argument values, respectively).}
}
\description{
	Cut (or slice) a phylogenetic tree at a particular time point.
}
\details{
	This function slices a tree at a particular height above the root and returns all subtrees or all non-trivial subtrees (i.e., subtrees with more than 1 taxon). Uses \code{\link{extract.clade}} in the "ape" package.
	
	\code{treeSlice} can also be used to crop the terminal fraction of a tree for \code{orientation="rootwards"}.
}
\value{
	An object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{extract.clade}}
}
\keyword{phylogenetics}
\keyword{utilities}

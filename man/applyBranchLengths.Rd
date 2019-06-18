\name{applyBranchLengths}
\alias{applyBranchLengths}
\title{Applies the branch lengths of a reference tree to a target}
\usage{
applyBranchLengths(tree, edge.length)
}
\arguments{
	\item{tree}{target tree.}
	\item{edge.length}{number of digits for rounding. Passed to \code{\link{round}}.}
}
\description{
	This function applies the set of branch lengths from a reference tree to a target tree while reconciling any mappings (as in \code{\link{read.simmap}}) with the new branch lengths.
}
\value{
	An object of class \code{"phylo"} or \code{"simmap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

\name{expand.clade}
\alias{expand.clade}
\alias{plot.expand.clade}
\title{Expands (or contracts) the tip-spacing of a given clade or clades}
\usage{
expand.clade(tree, node, factor=5)
\method{plot}{expand.clade}(x, ...)
}
\arguments{
	\item{tree}{tree an object of class \code{"phylo"} or \code{"simmap"}.}
	\item{node}{node index or vector of node indices.}
	\item{factor}{expansion factor for the tip-spacing of the taxa descended from node or nodes in \code{node}.}
	\item{x}{for \code{plot} method, an object of class \code{"expand.clade"}.}
	\item{...}{optional arguments to be passed to \code{plotTree} or \code{plotSimmap}, depending on the class of \code{x$tree}.}
}
\description{
	The purpose of this function is to compute a custom tip-spacing for users who want to expand or contract the tip-spacing of the descendant taxa from a given node or nodes.
}
\value{
	The function returns an object of class \code{"expand.clade"} which consists of the (possibly re-ordered) tree and a numerical vector with the calculated tip spacing based on the expansion factor specified by the user. This object can be plotted using the S3 \code{plot} method for the object class; or it can be plotted simplying by calling a standard plotting function on the tree & tip spacings.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}

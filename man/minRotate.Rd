\name{minRotate}
\alias{minRotate}
\alias{tipRotate}
\title{Rotates all nodes of the tree to minimize the difference in order with a vector}
\usage{
minRotate(tree, x, ...)
tipRotate(tree, x, ...)
}
\arguments{
	\item{tree}{tree.}
	\item{x}{numeric vector.}
	\item{...}{optional arguments to be used by \code{tipRotate}. Presently optional arguments can be \code{fn}, function to be used to compute the distance between the order of the tip labels in \code{tree} and the numeric vector \code{x} (presently \code{fn=function(x) x^2} by default); \code{methods}, the method or methods of tree traversal (can be \code{"pre"}, \code{"post"}, or \code{c("pre","post")}, for pre-, post-, or both pre- and post-order tree traversal); \code{rotate.multi}, whether to rotate multifurcations in all possible ways using \code{rotate.multi} (defaults to \code{FALSE}); and \code{print}, a logical argument specifying whether to print the search progress or to behave quietly. Only the option \code{print} is available for \code{minRotate}.}
}
\description{
	Rotates all the nodes of the tree to try and minimize the different between the order of the tips and the rank-order of a numeric vector \code{x} or (in the case of \code{tipRotate}) the actual integer vector, \code{x}.
}
\details{
	Both \code{minRotate} and \code{tipRotate} are designed primarily to be used internally by other \pkg{phytools} functions and particularly by \code{\link{phylo.to.map}} (in the case of \code{minRotate}) and by \code{\link{cophylo}} (in the case of \code{tipRotate}).
}
\value{
	A node-rotated object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}

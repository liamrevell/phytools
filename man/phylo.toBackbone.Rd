\name{phylo.toBackbone}
\alias{phylo.toBackbone}
\alias{backbone.toPhylo}
\title{Converts tree to backbone or vice versa}
\usage{
phylo.toBackbone(x, trans, ...)
backbone.toPhylo(x)
}
\arguments{
	\item{x}{an object of class \code{"phylo"} (for the function \code{phylo.toBackbone}), or an object of the special class \code{"backbonePhylo"} (for \code{backbone.toPhylo}).}
	\item{trans}{data frame containing the attributes necessary to translate a backbone tree to an object of class \code{"backbonePhylo"}. The data frame should contain the following variables: \code{tip.label}: the tip labels in the input tree (not all need be included); \code{clade.label}: labels for the unobserved subtrees; \code{N}: number of species in each subtree; and \code{depth}: desired depth of each subtree. \code{depth} for each terminal taxon in \code{x} cannot be greater than the terminal edge length for that taxon.}
	\item{...}{optional arguments.}
}
\description{
	Converts between \code{"phylo"} and \code{"backbonePhylo"}.
}
\value{
	Either an object of class \code{"phylo"} or an object of class \code{"backbonePhylo"}, depending on the method.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plot.backbonePhylo}}
}
\keyword{phylogenetics}
\keyword{plotting}


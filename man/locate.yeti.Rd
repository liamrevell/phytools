\name{locate.yeti}
\alias{locate.yeti}
\title{Locate a cryptic, recently extinct, or missing taxon on a tree}
\usage{
locate.yeti(tree, X, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix with continuous character data.}
	\item{...}{optional arguments including: \code{method} (\code{"ML"} or \code{"REML"}, defaults to \code{"ML"}); \code{search} (\code{"heuristic"} or \code{"exhaustive"}, defaults to \code{"heuristic"}); \code{constraint}, a vector containing the daughter node numbers from \code{tree$edge} for each edge to try; \code{plot} a logical argument specifying whether or not to plot the likelihood profile on edges (defaults to \code{FALSE}); \code{rotate} a logical indicating whether or not to rotate the data based on the input tree; and \code{quiet}, which is logical and has an obvious interpretation.}
}
\description{
	This function uses ML (or REML) to place a recently extinct, cryptic, or missing taxon on an ultrametric (i.e., time-calibrated) phylogeny.
}
\value{
	Optimized tree as an object of class \code{"phylo"}.
}
\references{
	Felsenstein, J. (1981) Maximum likelihood estimation of evolutionary trees from continuous characters. \emph{American Journal of Human Genetics}, 25, 471-492.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J., D. L. Mahler, R. G. Reynolds, and G. J. Slater. (In press) Placing cryptic, recently extinct, or hypothesized taxa into an ultrametric phylogeny using continuous character data: A case study with the lizard \emph{Anolis roosevelti}. \emph{Evolution}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}

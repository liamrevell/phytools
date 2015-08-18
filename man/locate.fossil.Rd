\name{locate.fossil}
\alias{locate.fossil}
\title{Locate a fossil lineage in a tree using continuous characters}
\usage{
locate.fossil(tree, X, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix with continuous character data.}
	\item{...}{optional arguments including \code{time.constraint} which can be a scalar (positive height above the root of the fossil or negative time before present) or a vector (age range of fossil, either positive or negative); \code{edge.constraint}, which is equivalent to \code{constraint} in \code{\link{locate.yeti}}; \code{plot}, \code{rotate}, and \code{quiet}, which have the same interpretation (and defaults) as the equivalent arguments in \code{\link{locate.yeti}}.}
}
\description{
	This function uses ML to place a fossil lineage into a tree using continuous traits.
}
\value{
	Optimized tree as an object of class \code{"phylo"}.
}
\references{
	Felsenstein, J. (1981) Maximum likelihood estimation of evolutionary trees from continuous characters. \emph{American Journal of Human Genetics}, 25, 471-492.

	Felsenstein, J. (2002) Quantitative characters, phylogenies, and morphometrics. In: MacLeod, N. and P. Forey (Eds.) \emph{Morphology, Shape and Phylogeny} (pp. 27-44). Taylor and Francis, London.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J., D. L. Mahler, R. G. Reynolds, and G. J. Slater. (In press) Placing cryptic, recently extinct, or hypothesized taxa into an ultrametric phylogeny using continuous character data: A case study with the lizard \emph{Anolis roosevelti}. \emph{Evolution}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{inference}
\keyword{maximum likelihood}

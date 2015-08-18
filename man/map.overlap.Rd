\name{map.overlap}
\alias{map.overlap}
\title{Proportional overlap between two mapped character histories on a tree}
\usage{
map.overlap(tree1, tree2, tol=1e-6)
}
\arguments{
	\item{tree1}{a modified \code{"phylo"} object (see \code{\link{read.simmap}}.}
	\item{tree2}{a modified \code{"phylo"} object with the same topology and branch lengths, but different map from \code{tree1}.}
	\item{tol}{an optional tolerance value.}
}
\description{
	This function computes the fraction of a stochastic character mapping that is shared between two differently mapped trees.
}
\value{
	A numerical value on the interval 0-1.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}

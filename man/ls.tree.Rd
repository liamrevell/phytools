\name{ls.tree}
\alias{ls.tree}
\title{Least squares branch lengths for a given tree}
\usage{
ls.tree(tree, D)
}
\arguments{
	\item{tree}{phylogeny.}
	\item{D}{distance matrix.}
}
\description{
	Computes the least squares branch lengths conditioned on a topology and distance matrix. Internal function for \code{\link{optim.phylo.ls}}.
}
\details{
	Do not use unless you know what you're doing.
}
\value{
	A tree with branch lengths.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{inference}
\keyword{least squares}

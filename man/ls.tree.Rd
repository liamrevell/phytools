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
	Computes the least squares branch lengths conditioned on a topology and distance matrix.
}
\details{
	Generally intended as a function to be used internally by \code{\link{optim.phylo.ls}}.
}
\value{
	A tree with branch lengths.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{phylogeny inference}
\keyword{least squares}
\keyword{consensus tree}

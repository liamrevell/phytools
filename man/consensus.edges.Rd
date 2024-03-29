\name{consensus.edges}
\alias{consensus.edges}
\title{Compute consensus edges for a tree under some criterion}
\usage{
consensus.edges(trees, method=c("mean.edge","least.squares"), ...)
}
\arguments{
	\item{trees}{object of class \code{"multiPhylo"} in which the trees must have edge lengths. This could be, for instance, a sample from the posterior distribution of trees in a Bayesian analysis.}
	\item{method}{method for computing the edge lengths. Could be the mean of all trees in which the edge is present, or it could be the least-squares edge lengths computed on the mean patristic distance matrices from the input phylogenies in \code{trees}. Note that in the latter case the \pkg{phangorn} function \code{nnls.tree} is used and the option \code{rooted} will be set to \code{is.rooted(tree)} for the consensus tree.}
	\item{...}{optional arguments, the most popular of which is \code{consensus.tree} - a user supplied consensus tree. Another optional argument for \code{method="mean.edge"} is \code{if.absent} which tells the function how to include absent edges in the computation of average edge lengths. Possible values are \code{"zero"} (the default) or \code{"ignore"}.}
}
\description{
	Computes consensus edge lengths from a set of input trees using multiple methods.
}
\value{
	An object of class \code{"phylo"} with edge lengths.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{consensus tree}
\keyword{phylogeny inference}
\keyword{utilities}

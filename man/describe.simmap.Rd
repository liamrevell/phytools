\name{describe.simmap}
\alias{describe.simmap}
\title{Summarizes a stochastic mapped tree or set of trees}
\usage{
describe.simmap(tree, ...)
}
\arguments{
	\item{tree}{a single tree or a set of trees with a mapped discrete character (e.g, see \code{\link{make.simmap}} or \code{\link{read.simmap}}.}
	\item{...}{optional arguments which include: \code{plot}, a logical value indicating whether or not to plot the posterior probabilities at nodes (default is \code{plot=FALSE}); \code{check.equal}, a logical value indicating whether or not to check if all trees are equal using \code{\link{all.equal.phylo}} (default is \code{check.equal=FALSE}); and \code{message}, a logical indicating whether or not to print an informative message to the screen (default is \code{message=TRUE}).}
}
\description{
	This function summarizes the result of one or more stochastic maps.
}
\value{
	A list with the following elements is returned (invisibly if \code{message=TRUE}):
	\item{count}{a matrix containing the number and types of transitions for each tree, if \code{class(tree)="multiPhylo"}.}
	\item{times}{a matrix containg the times spend in each state on each tree.}
	\item{ace}{the posterior probabilities of each node being in each state, if \code{class(tree)="multiPhylo"}.}
	\item{legend}{a vector containing the plot legend, if \code{plot=TRUE}.}
	if \code{class(tree)="phylo"} then the function simply returns the results of \code{\link{countSimmap}} combined with the states at each node of the tree and a matrix containing the total and relative times spent in each state on the tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{comparative method}
\keyword{phylogenetics}
\keyword{utilities}

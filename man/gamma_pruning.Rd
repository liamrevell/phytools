\name{gamma_pruning}
\alias{gamma_pruning}
\title{Internally used function}
\usage{
gamma_pruning(par, nrates=4, tree, x, model=NULL, median=TRUE,
  expm.method="Higham08.b", ...)
}
\arguments{
	\item{par}{function parameters.}
	\item{nrates}{number of discrete rate categories.}
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{data in the form of a binary matrix.}
	\item{model}{design matrix of fitted model.}
	\item{median}{use the median of each rate category.}
	\item{exp.method}{method used for matrix exponentiation.}
	\item{...}{optional arguments.}
}
\description{
	Internally used function to compute the likelihood under a \eqn{\Gamma} model.
}
\value{
	A value of the log-likelihood.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{utilities}

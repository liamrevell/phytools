\name{threshDIC}
\alias{threshDIC}
\title{Deviance Information Criterion from the threshold model}
\usage{
threshDIC(tree, x, mcmc, burnin=NULL, sequence=NULL, method="pD")
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{a named vector containing discrete character states; or a matrix containing the tip species, in rows, and probabilities of being in each state, in columns.}
	\item{mcmc}{list object returned by \code{\link{ancThresh}}.}
	\item{burnin}{number of generations (not samples) to exclude as burn in; if not supplied then 20 percent of generations are excluded.}
	\item{sequence}{assumed ordering of the discrete character state. If not supplied and \code{x} is a vector then numerical-alphabetical order is assumed; if not supplied and \code{x} is a matrix, then the column order of \code{x} is used.}
	\item{method}{method for computing the effective number of parameters (options are \code{"pD"} and \code{"pV"}).}
}
\description{
	This function computes the Deviance Information Criterion from the MCMC object returned by \code{\link{ancThresh}}.
}
\value{
	A vector containing the mean deviance and deviance for the parameter means, the effective number of parameters, and the DIC.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, bold{68}, 743-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}

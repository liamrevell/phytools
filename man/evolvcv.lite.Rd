\name{evolvcv.lite}
\alias{evolvcv.lite}
\title{Likelihood test for a shift in the evolutionary correlation between traits}
\usage{
evolvcv.lite(tree, X, maxit=2000, tol=1e-10)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}.}
	\item{X}{an \code{n} x \code{m} matrix of tip values for \code{m} continuously valued traits in \code{n} species - row names should be species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization - may need to be increased for large trees.}
	\item{tol}{tolerance value for \code{"L-BFGS-B"} optimization.}
}
\description{
	This function takes an object of class \code{"simmap"} with a mapped binary or multistate trait and data for two and only two continuously valued character. It then fits four different evolutionary models: common rates and correlation; different rates, common correlation; different correlations, common rates; no common structure.
}
\value{
	A list with the results summarized for each model.
}
\references{
	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{evol.vcv}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}

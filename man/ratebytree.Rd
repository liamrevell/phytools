\name{ratebytree}
\alias{ratebytree}
\title{Likelihood test for rate variation among trees}
\usage{
ratebytree(trees, x, ...)
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"}.}
	\item{x}{a list of trait vectors for a continuous trait in which the names of each vectors correspond to the tip labels of \code{trees}.}
	\item{...}{optional arguments. Presently, these include the following: \code{model}, the model of continuous trait evolution (options are \code{"BM"}, the default, \code{"OU"}, and \code{"EB"}); \code{tol}, used as a minimum value for the fitting rates, to prevent problems in optimization; \code{trace}, a logical value indicating whether or not to report progress in the optimization; \code{test}, the method for hypothesis testing (options are \code{"chisq"} and \code{"simulation"}); \code{quiet}, a logical value indicating whether or not to run perfectly quietly; and \code{se}, a list of vectors containing the standard errors for each value of \code{x}.}
}
\description{
	This function takes a object of class \code{"multiPhylo"} containing two or more phylogenies (\code{trees}), and a list of trait vectors (\code{x}), and fits two models: one in which the rate (or regime, for models \code{"OU"} and \code{"EB"}) of trait evolution is equal among all trees; and a second in which the rates or regimes can differ between trees.
	
	This model corresponds to an extension the \emph{censored} approach of O'Meara et al. (2006; \emph{Evolution}) and should also be related to the method of Adams (2012; \emph{Systematic Biology}) for comparing rates among traits. See \code{\link{brownie.lite}} for a different implementation of the \emph{noncensored} approach of O'Meara et al. (2006).
}
\details{
	The function also conducts a likelihood-ratio test to compare the two models.
}
\value{
	An object of class \code{"ratebytree"}.
}
\references{
	Adams, D. C. (2012) Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. \emph{Syst. Biol.}, \bold{62}, 181-192.

	O'Meara, B. C., C. Ane, M. J. Sanderson, and P. C. Wainwright. (2006) Testing for different rates of continuous trait evolution using likelihood. \emph{Evolution}, \bold{60}, 922-933.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}

\name{aic.w}
\alias{aic.w}
\title{Computes AIC weights}
\usage{
aic.w(aic)
}
\arguments{
	\item{aic}{vector of AIC values for different fitted models. If the vector has names, these names will be inherited by the vector returned by the function.}
}
\description{
	This function computes AIC weights for a set of fitted models. It returns an object of class \code{"aic.w"} which is just a vector which allows it to be automatically printed with a numerical precision of 8 significant digits.
}
\value{
	A vector of AIC weights.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{statistics}
\keyword{utilities}
\keyword{math}
\keyword{maximum likelihood}
\keyword{information criteria}

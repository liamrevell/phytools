\name{likMlambda}
\alias{likMlambda}
\title{Likelihood for joint lambda}
\usage{
likMlambda(lambda, X, C)
}
\arguments{
	\item{lambda}{scalar, usually on the interval 0,1.}
	\item{X}{data, in a matrix.}
	\item{C}{matrix probably returned by \code{\link{vcv.phylo}}.}
}
\description{
	Computes the likelihood.
}
\details{
	Do not use unless you know what you're doing.
}
\value{
	A scalar.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{math}
\keyword{comparative method}
\keyword{utilities}

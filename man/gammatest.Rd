\name{gammatest}
\alias{gammatest}
\title{Gamma test of Pybus & Harvey (2000)}
\usage{
gammatest(x)
}
\arguments{
	\item{x}{an object of class \code{"ltt"} resulting from a call of the function \code{\link{ltt}}.}
}
\description{
	Conducts \eqn{\gamma}-test of Pybus & Harvey (2000).
}
\value{
	A an object of class \code{"gammatest"} consisting of a list that contains:
	\item{gamma}{a value for the \eqn{\gamma}-statistic.}
	\item{p}{two-tailed P-value for the \eqn{\gamma}-test.}
}
\references{
	Pybus, O. G., and P. H. Harvey (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proc. R. Soc. Lond. B}, \bold{267}, 2267-2272.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ltt}}, \code{\link{mccr}}
}
\examples{
data(darter.tree)
gammatest(ltt(darter.tree,plot=FALSE))
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{diversification}

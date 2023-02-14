\name{threshState}
\alias{threshState}
\title{Computes value for a threshold character from a liability and thresholds}
\usage{
threshState(x, thresholds)
}
\arguments{
	\item{x}{liability.}
	\item{thresholds}{a named vector containing the thresholds.}
}
\description{
	Internal function for \code{\link{ancThresh}}.
}
\details{
	\code{threshState} can also be used to simulate threshold traits.
}
\value{
	A discrete character value.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}, \code{\link{threshDIC}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{discrete character}
\keyword{utilities}

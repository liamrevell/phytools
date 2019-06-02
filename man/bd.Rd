\name{bd}
\alias{bd}
\title{Convert object of class "birthdeath" to raw birth & death rates}
\usage{
bd(x)
}
\arguments{
	\item{x}{object of class \code{"birthdeath"}.}
}
\description{
	This function converts an object of class \code{"birthdeath"} to a vector with the ML birth & death rates.
	
	Note that this is somewhat unecessary as \emph{phytools} now contains functions to fit birth-death and pure-birth diversification models from trees (\code{\link{fit.bd}} and \code{\link{fit.yule}}) that also take into account incomplete sampling fraction.
}
\value{
	A vector.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{birthdeath}}, \code{\link{fit.bd}}, \code{\link{fit.yule}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}

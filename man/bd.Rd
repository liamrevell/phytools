\name{bd}
\alias{bd}
\title{Convert object of class \code{"birthdeath"} to raw birth & death rates}
\usage{
bd(x)
}
\arguments{
	\item{x}{object of class \code{"birthdeath"}.}
}
\description{
	Converts an object of class \code{"birthdeath"} (from \code{\link{birthdeath}}) to a vector with the ML birth & death rates.
}
\details{
	This is now somewhat obsolete as \pkg{phytools} now contains functions to fit birth-death and pure-birth diversification models from trees (\code{\link{fit.bd}} and \code{\link{fit.yule}}) that also take into account incomplete sampling fraction.
}
\value{
	A vector.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{birthdeath}}, \code{\link{fit.bd}}, \code{\link{fit.yule}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}
\keyword{diversification}
\keyword{maximum likelihood}

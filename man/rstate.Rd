\name{rstate}
\alias{rstate}
\title{Pick a random state according to a vector of probabilities}
\usage{
rstate(y)
}
\arguments{
	\item{y}{vector of probabilities. Must have names & should probably add to \code{1.0}.}
}
\description{
	Internal function for \code{\link{make.simmap}}.
}
\details{
	This function picks a random element in a vector according to the probability assigned that element.  It returns the name. Uses \code{\link{rmultinom}}.
}
\value{
	A character or string.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{math}
\keyword{utilities}
\keyword{statistics}
\keyword{discrete character}
\keyword{simulation}

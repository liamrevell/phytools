\name{anceb}
\alias{anceb}
\title{Compute marginal ancestral states}
\usage{
anceb(object, ...)
}
\arguments{
	\item{object}{fitted M\emph{k} model (presently object of class \code{"fitMk"}).}
	\item{...}{optional arguments.}
}
\description{
	Computes marginal ancestral states.
}
\details{
	Still in progress.
}
\value{
	An object of class \code{"anceb"}, which includes a matrix of marginal (empirical Bayes) probabilities.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}

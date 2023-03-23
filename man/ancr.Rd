\name{ancr}
\alias{ancr}
\title{Compute marginal ancestral states}
\usage{
ancr(object, ...)
}
\arguments{
	\item{object}{fitted M\emph{k} model (presently object of class \code{"fitMk"}).}
	\item{...}{optional arguments.}
}
\description{
	Computes marginal ancestral states.
}
\details{
	If the optional argument `tips=TRUE`, then the matrix returned contains empirical Bayes posterior probabilities (marginal scaled likelihoods) for both tips and internal nodes. Otherwise (the default) only node marginal states are returned.
}
\value{
	An object of class \code{"ancr"}, which includes a matrix of marginal (empirical Bayes) probabilities.
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

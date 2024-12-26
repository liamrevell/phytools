\name{likMlambda}
\alias{likMlambda}
\title{Likelihood for joint \eqn{\lambda}}
\usage{
likMlambda(lambda, X, C)
}
\arguments{
	\item{lambda}{scalar, usually on the interval 0,1 (although not required to be).}
	\item{X}{data for various continuous character, in the form of a matrix.}
	\item{C}{\emph{n} x \emph{n} matrix (for \emph{n} taxa) containing the height above the root for each pair of taxa in the tree (e.g., \code{\link[ape]{vcv.phylo}}).}
}
\description{
	Likelihood function for joint estimation of Pagel's \eqn{\lambda} parameter.
}
\details{
	Generally intended to be used internally by other methods that do joint optimization of \eqn{\lambda} (e.g., \code{\link{phyl.pca}}).
}
\value{
	The log-likelihood.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{math}
\keyword{comparative method}
\keyword{utilities}
\keyword{maximum likelihood}

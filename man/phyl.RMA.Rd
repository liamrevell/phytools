\name{phyl.RMA}
\alias{phyl.RMA}
\title{Phylogenetic reduced major axis (RMA) regression}
\usage{
phyl.RMA(x, y, tree, method="BM", lambda=NULL, fixed=FALSE, h0=1.0)
}
\arguments{
	\item{x}{vector with names.}
	\item{y}{vector with names.}
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{method}{method to obtain the correlation structure: can be \code{"BM"} or \code{"lambda"}.}
	\item{lambda}{value of \code{lambda} for fixed \eqn{\lambda}.}
	\item{fixed}{logical value indicating whether or not \eqn{\lambda} should be optimized using likelihood.}
	\item{h0}{null hypothesis for \code{beta}. Defaults to 1.0. Note that a null hypothesis of 0.0 is not allowed.}
}
\description{
	This function performs phylogenetic RMA regression.
}
\details{
	Optionally jointly estimates lambda if \code{method="lambda"}. Likelihood optimization of lambda is performed using \code{\link{optimize}} on the interval (0,1). Note that some statistician think there is \emph{never} a condition in which a reduced-major-axis regression should be used.
}
\value{
	A list with the following components:
	\item{RMA.beta}{a vector of RMA regression coefficients.}
	\item{V}{a VCV matrix for the traits.}
	\item{lambda}{fitted value of lambda (\code{method="lambda"} only).}
	\item{logL}{log-likelihood (\code{method="lambda"} only).}
	\item{test}{a vector containing results for hypothesis tests on \code{beta}.}
	\item{resid}{a vector of residuals for \code{y} given \code{x}.}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phyl.cca}}, \code{\link{phyl.pca}}, \code{\link{phyl.resid}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{statistics}

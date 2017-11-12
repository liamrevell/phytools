\name{phyl.resid}
\alias{phyl.resid}
\title{Phylogenetic size-correction via GLS regression}
\usage{
phyl.resid(tree, x, Y, method="BM")
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{vector containing the single independent variable (e.g., size), or matrix with multiple independent variables in columns.}
	\item{Y}{vector or matrix with one or multiple dependent variables in columns.}
	\item{method}{method to obtain the correlation structure: can be \code{"BM"} or \code{"lambda"}.}
}
\description{
	This function fits one or multiple phylogenetic regressions (depending on the number of columns in \code{Y}) and computes the residuals.  Designed for phylogenetic size correction using GLS regression (e.g., Revell 2009; \emph{Evolution}).
}
\details{
	Optionally fits \eqn{\lambda} for each regression model.  Likelihood optimization of \eqn{\lambda} is performed for \code{method=} \code{"lambda"} using \code{\link{optimize}} on the interval (0,1). This function is theoretically redundant with \code{residuals} applied to a \code{"gls"} object class in which the correlation structure is based on \code{corBrownian} or \code{corPagel}; however some users may find this method simpler, and it provides a good way to cross-check your results & make sure that you are using \code{gls} correctly.
}
\value{
	A list with the following components:
	\item{beta}{a vector or matrix of regression coefficients.}
	\item{resid}{a vector or matrix of residuals for species.}
	\item{lambda}{a vector of lambda values (\code{method="lambda"} only).}
	\item{logL}{a vector of log-likelihoods (\code{method="lambda"} only).}
}
\references{
	Revell, L. J. (2009) Size-correction and principal components for interspecific comparative studies. \emph{Evolution}, \bold{63}, 3258-3268.

	Revell, L. J. (2010) Phylogenetic signal and linear regression on species data. \emph{Methods in Ecology and Evolution}, \bold{1}, 319-329.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phyl.pca}}, \code{\link{gls}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{maximum likelihood}

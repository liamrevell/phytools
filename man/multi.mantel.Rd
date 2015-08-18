\name{multi.mantel}
\alias{multi.mantel}
\title{Multiple matrix regression (partial Mantel test)}
\usage{
multi.mantel(Y, X, nperm=1000)
}
\arguments{
	\item{Y}{single "dependent" square matrix. Can be either a symmetric matrix of class \code{"matrix"} or a distance matrix of class \code{"dist"}.}
	\item{X}{a single independent matrix or multiple independent matrices in a list. As with \code{Y} can be a object of class \code{"matrix"} or class \code{"dist"}.}
	\item{nperm}{number of Mantel permutations.}
}
\description{
	This function conducting a multiple matrix regression (partial Mantel test) and uses Mantel (1967) permutations to test the significance of the model and individual coefficients. It also returns the residual and predicted matrices.
}
\value{
	A list with the following components:
	\item{r.squared}{multiple R-squared.}
	\item{coefficients}{model coefficients, including intercept.}
	\item{tstatistic}{t-statistics for model coefficients.}
	\item{fstatistic}{F-statistic for the overall model.}
	\item{probt}{vector of probabilities, based on permutations, for \code{tstatistic}.}
	\item{probF}{probability of F, based on Mantel permutations.}
	\item{residuals}{matrix of residuals.}
	\item{predicted}{matrix of predicted values.}
}
\references{
	Mantel, N. (1967) The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \bold{27}, 209--220.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}

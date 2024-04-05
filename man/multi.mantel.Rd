\name{multi.mantel}
\alias{multi.mantel}
\title{Multiple matrix regression (partial Mantel test)}
\usage{
multi.mantel(Y, X, nperm=1000)
}
\arguments{
	\item{Y}{single "dependent" square matrix. Can be either a symmetric matrix of class \code{"matrix"} or a distance matrix of class \code{"dist"}.}
	\item{X}{a single independent matrix or multiple independent matrices in a list. As with \code{Y} can be a object of class \code{"matrix"} or class \code{"dist"}, or a list of such objects.}
	\item{nperm}{number of Mantel permutations to be used to compute a P-value of the test.}
}
\description{
	This function conducting a multiple matrix regression (partial Mantel test) and uses Mantel (1967) permutations to test the significance of the model and individual coefficients. It also returns the residual and predicted matrices.
}
\value{
	An object of class \code{"multi.mantel"} consisting of the following elements:
	\item{r.squared}{multiple R-squared.}
	\item{coefficients}{model coefficients, including intercept.}
	\item{tstatistic}{t-statistics for model coefficients.}
	\item{fstatistic}{F-statistic for the overall model.}
	\item{probt}{vector of probabilities, based on permutations, for \code{tstatistic}.}
	\item{probF}{probability of F, based on Mantel permutations.}
	\item{residuals}{matrix of residuals.}
	\item{predicted}{matrix of predicted values.}
	\item{nperm}{the number of permutations used.}
}
\details{
	Printing the object to screen will result in a summary of the analysis similar to \code{summary.lm}, but with p-values derived from Mantel permutations.
	
	Methods \code{residuals} and \code{fitted} can be used to return residual and fitted matrices, respectively.
}	
\references{
	Mantel, N. (1967) The detection of disease clustering and a generalized regression approach. \emph{Cancer Research}, \bold{27}, 209-220.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{distance matrix}

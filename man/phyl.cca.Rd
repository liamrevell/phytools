\name{phyl.cca}
\alias{phyl.cca}
\title{Phylogenetic canonical correlation analysis}
\usage{
phyl.cca(tree, X, Y, lambda=1.0, fixed=TRUE)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{X}{a data matrix with traits in columns.}
	\item{Y}{data matrix with traits in columns, to be correlated with \code{X}.}
	\item{lambda}{optionally, a (fixed) value for \eqn{\lambda}.}
	\item{fixed}{optionally, a logical value indicating whether or not to estimate \eqn{\lambda} using likelihood.}
}
\description{
	Phylogenetic canonical correlation analysis, following Revell & Harrison (2008).
}
\details{
	(Optional) joint optimization of \eqn{\lambda} is performed using \code{\link{optimize}} on the interval (0,1).
}
\value{
	An object of class \code{"phyl.cca"} containing the following elements:
	\item{cor}{canonical correlations.}
	\item{xcoef}{coefficients for the canonical variables for \code{X}.}
	\item{ycoef}{coefficients for the canonical variables for \code{Y}.}
	\item{xscores}{matrix with the canonical scores for \code{X}.}
	\item{yscores}{matrix with the canonical scores for \code{Y}.}
	\item{chisq}{vector of \eqn{\chi^2} values.}
	\item{p}{P-values for the hypothesis test that the \emph{i}th and all subsequent correlations are zero.}
}
\references{
	Revell, L. J., Harrison, A. S. (2008) PCCA: A program for phylogenetic canonical correlation analysis. \emph{Bioinformatics}, \bold{24}, 1018-1020.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phyl.pca}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}

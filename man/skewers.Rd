\name{skewers}
\alias{skewers}
\title{Matrix comparison using the method of random skewers}
\usage{
skewers(X, Y, nsim=100, method=NULL)
}
\arguments{
	\item{X}{covariance matrix.}
	\item{Y}{covariance matrix.}
	\item{nsim}{number of random vectors.}
	\item{method}{method to generate a null distribution of the random skewers correlation between matrices. If \code{method=NULL} then the correlation will be compared to the correlation between random vectors; however this test has type I error substantially above the nominal level for ostensibly random matrices. Other values of \code{method} will be passed as \code{covMethod} to \code{genPositiveDefMat} for a more robust hypothesis test (see below). Recommended values include \code{"unifcorrmat"}.}
}
\description{
	Performs the random skewers matrix comparison method of Cheverud (1996).
}
\details{
	This function performs the random skewers matrix comparison method of Cheverud (1996; also see Cheverud & Marroig 2007 for more details). In addition, it includes a more robust hypothesis test in which random covariance matrices are simulated under a variety of models, and then the mean correlation between response vectors to random skewers are computed.
}
\value{
	A list with the following components:
	\item{r}{mean random skewers correlation.}
	\item{p}{p-value from simulation.}
}
\references{
	Cheverud, J. M. (1996) Quantitative genetic analysis of cranial morphology in the cotton-top (\emph{Saguinus oedipus}) and saddle-back (\emph{S. fuscicollis}) tamarins. \emph{J. Evol. Biol.}, \bold{9}, 5--42.

	Cheverud, J. M. and G. Marroig (2007) Comparing covariance matrices: Random skewers method compared to the common principal components model. \emph{Genetics & Molecular Biology}, \bold{30}, 461--469.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{comparative method}
\keyword{statistics}

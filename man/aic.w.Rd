\name{aic.w}
\alias{aic.w}
\title{Computes Akaike weights}
\usage{
aic.w(aic)
}
\arguments{
	\item{aic}{vector of Akaike Information Criterion (AIC; Akaike, 1974) values for different fitted models. If the vector has names, these names will be inherited by the vector returned by the function.}
}
\description{
	Computes Akaike weights based on a set of AIC values.
}
\details{
	This function computes Akaike weights from a set of AIC values obtained from fitted models.
	
	It returns an object of class \code{"aic.w"} which is just a vector which allows it to be automatically printed with a numerical precision of 8 significant digits.
}
\value{
	A vector of Akaike weights.
}
\references{
	Akaike, H. (1974) A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control}, \bold{19}, 716-723.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{statistics}
\keyword{utilities}
\keyword{math}
\keyword{maximum likelihood}
\keyword{information criteria}

\name{fitThresh}
\alias{fitThresh}
\title{Fits multi-state threshold model using the discrete approximation}
\usage{
fitThresh(tree, x, sequence=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{a named vector.}
	\item{sequence}{assumed ordering of the discrete character state. If not supplied and then alphanumeric order is assumed.}
	\item{...}{optional arguments.}
}
\description{
	This function uses the discrete approximation of Boucher and Demery (2016) to fit a multi-state discrete character threshold model following Felsenstein (2012; Revell 2014).
}
\value{
	An object of class \code{"fitThresh"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.

	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}, \code{\link{bounded_bm}}, \code{\link{fitmultiBM}}, \code{\link{threshBayes}}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{discrete character}

\name{bounded_bm}
\alias{bounded_bm}
\title{Fits bounded and wrapped Brownian motion models}
\usage{
	bounded_bm(tree, x, lims=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{continuous character data in a vector.}
	\item{lims}{bounds. If \code{NULL} will be set to three times the observed range of the trait. (This is designed to be an effectively "unbounded" null model.)}
	\item{...}{optional arguments.}
}
\description{
	Fits a bounded Brownian motion model using the discrete approximation of Boucher & Demery (2016).
}
\details{
	Optional arguments.
}
\value{
	An object of class \code{"bounded_bm"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.
	
	Juhn, M., K. Martinet, F. Boucher, A. Schultz, L. Harmon, K. Burns, and M. Alfaro (In review) A wrapped Brownian motion model for the evolution of color.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{ancestral states}

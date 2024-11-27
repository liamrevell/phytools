\name{fitmultiBM}
\alias{fitmultiBM}
\title{Fits bounded and wrapped Brownian motion models}
\usage{
	fitmultiBM(tree, x, y=NULL, model="ER", ncat=1, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{continuous character data in a vector.}
	\item{y}{optional discrete character, normally in the form of a factor or character vector.}
	\item{model}{discrete character evolution model. (E.g., \code{"ER"}, \code{"SYM"}, \code{"ARD"}.)}
	\item{ncat}{number of rate categories (if hidden rates).}
	\item{...}{optional arguments.}
}
\description{
	Fits a multi-rate discrete-state-dependent Brownian motion model using the discrete approximation of Boucher & Demery (2016).
}
\details{
	Optional arguments.
}
\value{
	An object of class \code{"fitmultiBM"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.
		
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bounded_bm}}, \code{\link{brownie.lite}}, \code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{ancestral states}

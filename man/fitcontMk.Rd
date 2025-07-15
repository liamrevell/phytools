\name{fitcontMk}
\alias{fitcontMk}
\title{Fits binary discrete character evolution model where the transition rate varies as a function of a continuous trait}
\usage{
fitcontMk(tree, x, y, model="sigmoidal", ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{binary discrete character, normally in the form of a factor or character vector.}
	\item{y}{continuous character data in a vector.}
	\item{model}{functional form describing the relationship between the value of the continuous trait and the binary discrete character transition rate. Presently the only option is \code{model="sigmoidal"}.}
	\item{...}{optional arguments.}
}
\description{
	Fits a model where the rate of evolution of a binary trait depends on the value of a continuous character using the discrete approximation of Boucher & Demery (2016).
}
\details{
	\code{fitmultiBM} uses the discrete approximation of Boucher & Demery (2016) to fit a model where the rate of evolution of a binary trait (in \code{x}) depends on the state of a continuous trait (in \code{y}), where the discrete and continuous character are modeled jointly. In some ways this is philosophically equivalent to the binary trait correlation method of Pagel (1994), but where one of the traits is continuously-valued.
	
	Optional arguments include \code{parallel} (whether to distribute matrix exponentiation across nodes, which defaults to \code{parallel=TRUE}), \code{levs} (the number of levels for the discrete approximation, which defaults to \code{levs=100}), and \code{plot_model} (a logical indicating whether or not to graph the model while optimizing (which defaults to \code{plot_model=FALSE}).
}
\value{
	An object of class \code{"fitcontMk"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.
	
	Pagel, M. (1994) Detecting correlated evolution on phylogenies: A general method for the comparative analysis of discrete characters. \emph{Proceedings of the Royal Society B}, \bold{255}, 37-45.
		
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bounded_bm}}, \code{\link{brownie.lite}}, \code{\link{fitMk}}, \code{\link{fitmultiBM}}, \code{\link{fitPagel}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{ancestral states}

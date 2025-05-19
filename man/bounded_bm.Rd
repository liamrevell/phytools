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
	Fits a bounded Brownian motion model using the discrete approximation of Boucher & Demery (2016), or the "wrapped" or circular model of Martinet et al. (In revision).
}
\details{
	Optional arguments include (but are not limited to): \code{parallel}, defaults to \code{parallel=FALSE}; \code{df}, to adjust the degrees of freedom of the model, depending on whether or not the bounds have been set from the observed data; \code{levs}, number of levels for the discretization (defaults to \code{levs=200}); and \code{wrapped}, logical indicating whether to fit the "circular" or wrapped model, instead of the standard bounded model (defaults to \code{wrapped=FALSE}).
}
\value{
	An object of class \code{"bounded_bm"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.
	
	Martinet, K. M., M. S. Juhn, F. C. Boucher, L. J. Harmon, L. J. Revell, S. I. A. Foerster, A. J. Shultz, K. Burns, M. E. Alfaro. (In revision) A wrapped Brownian motion model for traits on a circular scale. 
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}, \code{\link{fitmultiBM}}, \code{\link{fitThresh}}
}
\examples{
\dontrun{
## load tree and data
data(eel.tree)
data(eel.data)

## extract trait of interval
maxTL<-setNames(eel.data$Max_TL_cm,
    rownames(eel.data))

## fit bounded model
eel_bounded<-bounded_bm(eel.tree,maxTL,levs=200,
    parallel=TRUE,lims=range(maxTL))

## fit unbounded model
eel_unbounded<-bounded_bm(eel.tree,maxTL,levs=200,
    parallel=TRUE)

## compare models using AIC
AIC(eel_bounded,eel_unbounded)}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{ancestral states}

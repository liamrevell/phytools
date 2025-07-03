\name{fitmultiBM}
\alias{fitmultiBM}
\alias{fitmultiTrend}
\title{Fits multi-rate (or multi-trend) discrete character dependent continuous trait evolution}
\usage{
fitmultiBM(tree, x, y=NULL, model="ER", ncat=1, ...)
fitmultiTrend(tree, x, y=NULL, model="ER", ncat=1, ...)
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
	Fits a multi-rate (or multi-trend) discrete-state-dependent Brownian motion model using the discrete approximation of Boucher & Demery (2016).
}
\details{
	\code{fitmultiBM} uses the discrete approximation of Boucher & Demery (2016) to fit a model in which the rate of Brownian evolution of a continuous trait (in \code{x}) depends on the state of a discrete character (in \code{y}), where the discrete and continuous character are modeled jointly. If, in fact, the rates of evolution of \code{x} depend on the state of \code{y}, this will lead to more accurate model parameter estimates than when the discrete and continuous characters are separately modeled.
	
	\code{fitmultiTrend} fits the same model, but in which the mean of the Brownian process, and not only its rate, vary as a function of the character \code{y}.
	
	Optional arguments include \code{parallel} (whether to distribute matrix exponentiation across nodes, which defaults to \code{parallel=FALSE}), \code{levs} (the number of levels for the discrete approximation, which defaults to \code{levs=100}), \code{nullmodel} (a logical indicating whether or not to fit a null model without state-dependent continuous trait evolution, which defaults to \code{nullmodel=FALSE}), and \code{plot_model} (a logical indicating whether or not to graph the discretized model prior to fitting, which defaults to \code{plot_model=FALSE}).
}
\value{
	An object of class \code{"fitmultiBM"} or \code{"fitmultiTrend"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.
		
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bounded_bm}}, \code{\link{brownie.lite}}, \code{\link{fitMk}}
}
\examples{
\dontrun{
## load tree and data
data(eel.tree)
data(eel.data)
feed_mode<-setNames(eel.data$feed_mode,
    rownames(eel.data))
lnTL<-setNames(log(eel.data$Max_TL_cm),
    rownames(eel.data))
	
## fit state-dependent model
eel_sdmodel<-fitmultiBM(eel.tree,lnTL,feed_mode,
    model="ARD",parallel=TRUE,levs=100,
    plot_model=TRUE)

## fit null model
eel_nullmodel<-fitmultiBM(eel.tree,lnTL,feed_mode,
    model="ARD",parallel=TRUE,levs=100,
    null_model=TRUE,plot_model=TRUE)

## fit null model with hidden-rates
eel_nullmodel.hrm<-fitmultiBM(eel.tree,lnTL,feed_mode,
    model="ARD",ncat=2,parallel=TRUE,levs=100,
    null_model=TRUE,plot_model=TRUE)

## fit state-dependent model with hidden-rates
eel_sdmodel.hrm<-fitmultiBM(eel.tree,lnTL,feed_mode,
    model="ARD",ncat=2,parallel=TRUE,levs=100,
    plot_model=TRUE)

## compare models
anova(eel_nullmodel,eel_sdmodel,eel_nullmodel.hrm,
    eel_sdmodel.hrm)}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{ancestral states}

\name{ltt}
\alias{ltt}
\alias{ltt.phylo}
\alias{ltt.multiPhylo}
\alias{ltt.simmap}
\alias{ltt.multiSimmap}
\alias{gtt}
\alias{mccr}
\title{Creates lineage-through-time plot (including extinct lineages)}
\usage{
ltt(tree, ...)
\method{ltt}{phylo}(tree, plot=TRUE, drop.extinct=FALSE, log.lineages=TRUE, gamma=TRUE, ...)
\method{ltt}{multiPhylo}(tree, drop.extinct=FALSE, gamma=TRUE, ...)
\method{ltt}{simmap}(tree, plot=TRUE, log.lineages=FALSE, gamma=TRUE, ...)
\method{ltt}{multiSimmap}(tree, gamma=TRUE, ...)
gtt(tree, n=100, ...)
mccr(obj, rho=1, nsim=100, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}, \code{"multiPhylo"}, \code{"simmap"}, or \code{"multiSimmap"}.}
	\item{plot}{a logical value indicating whether or not to create LTT plot.}
	\item{drop.extinct}{logical value indicating whether or not to drop extinct tips from the tree.}
	\item{log.lineages}{logical value indicating whether LTT plot should be on log-linear (default for \code{"phylo"} and \code{"multiPhylo"} objects) or linear-linear scale.}
	\item{gamma}{logical value indicating whether or not to compute \eqn{\gamma} from Pybus & Harvey (2000; \emph{Proc. Roy. Soc. B}).}
	\item{n}{for \code{gtt} the number of time intervals to use to track \eqn{\gamma} through time.}
	\item{obj}{for \code{mccr} an object of class \code{"ltt"}.}
	\item{rho}{for \code{mccr} sampling fraction.}
	\item{nsim}{for \code{mccr} number of simulations to use for the MCCR test.}
	\item{...}{other arguments to be passed to plotting methods. See \code{\link{plot.default}}.}
}
\description{
	Computes and visualizes a lineage through time (LTT) plot, and related measures.
}
\details{
	The function \code{ltt} computes LTT plot with extant and extinct lineages, and optionally conducts \eqn{\gamma}-test of Pybus & Harvey (2000). The object returned by \code{ltt} can be plotted or re-plotted using \code{\link{plot}}.
	
	For the case in which \code{tree} is an object of class \code{"simmap"} or \code{"multiSimmap"} then the object will contain the number of lineages through time (for each tree, in the case of \code{"multiSimmap"} objects) separated by mapped regimes.
	
	The function \code{gtt} computes the value of Pybus & Harvey's \eqn{\gamma} statistic through time by slicing the tree at various points - by default in even intervals from the time above the root at which \emph{N} = 3 to the present day.
	
	The function \code{mccr} performs the MCCR test of Pybus & Harvey (2000) which takes into account incomplete taxon sampling in computing a P-value of the \eqn{\gamma} statistic.

	Although it is calculated here, it's unclear how to interpret the \eqn{\gamma}-statistic if not all the tips in the tree are contemporaneous.
}
\value{
	\code{ltt} returns an object of class \code{"ltt"} which normally includes the following elements:
	\item{times}{a vector of branching times.}
	\item{ltt}{a vector of lineages, or a matrix of lineages in each state over time for objects of class \code{"simmap"} and \code{"multiSimmap"}.}
	\item{gamma}{optionally, a value for the \eqn{\gamma}-statistic.}
	\item{p}{two-tailed P-value for the \eqn{\gamma}-test.}
	If \code{tree} is an object of class \code{"multiPhylo"}, then an object of class \code{"multiLtt"} is returned consisting of a list of object of class \code{"ltt"}.
	
	\code{gtt} returns an object of class \code{"gtt"}.
	
	\code{mccr} returns of object of class \code{"mccr"}.
}
\references{
	Pybus, O. G., and P. H. Harvey (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proc. R. Soc. Lond. B}, \bold{267}, 2267-2272.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{gammatest}}, \code{\link{ltt95}}
}
\examples{
## LTT plots
set.seed(99)
trees<-pbtree(n=100,scale=100,nsim=10)
obj<-ltt(trees,plot=FALSE)
plot(obj,log="y",log.lineages=FALSE,
    bty="l")
title(main="LTT plots for 10 pure-birth trees",
    font.main=3)
tree<-pbtree(b=1,d=0.25,t=4)
obj<-ltt(tree,gamma=FALSE,show.tree=TRUE,
    bty="l")
title(main="LTT plot with superimposed tree",
    font.main=3)
obj
## GTT plot
data(anoletree)
anole.gtt<-gtt(anoletree,n=40)
plot(anole.gtt)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{diversification}

\name{ltt}
\alias{ltt}
\title{Creates lineage-through-time plot (including extinct lineages)}
\usage{
ltt(tree, plot=TRUE, drop.extinct=FALSE, log.lineages=TRUE, gamma=TRUE, ...)
}
\arguments{
	\item{tree}{is a phylogenetic tree in \code{"phylo"} format, or an object of class \code{"multiPhylo"} containing a list of phylogenetic trees.}
	\item{plot}{a logical value indicating whether or not to create LTT plot.}
	\item{drop.extinct}{logical value indicating whether or not to drop extinct tips from the tree.}
	\item{log.lineages}{logical value indicating whether LTT plot should be on log-linear (default) or linear-linear scale.}
	\item{gamma}{logical value indicating whether or not to compute eqn{gamma} from Pybus & Harvey (2000; \emph{Proc. Roy. Soc. B}).}
	\item{...}{other arguments to be passed to plotting methods. See \code{\link{plot.default}}.}
}
\description{
	This function computes LTT plot with extant and extinct lineages, and optionally conducts \eqn{\gamma}-test of Pybus & Harvey (2000). The object returned by \code{ltt} can be plotted or re-plotted using \code{\link{plot}}.
}
\details{
	Although it is calculated here, it's unclear how to interpret the \eqn{\gamma}-statistic if not all the tips in the tree are contemporaneous.
}
\value{
	An object of class \code{"ltt"} which includes the following components:
	\item{times}{a vector of branching times.}
	\item{ltt}{a vector of linages.}
	\item{gamma}{optionally, a value for the gamma-statistic.}
	\item{p}{two-tailed P-value for the gamma-test.}
	If \code{tree} is an object of class \code{"multiPhylo"}, then an object of class \code{"multiLtt"} is returned consisting of a list of object of class \code{"ltt"}.
}
\references{
	Pybus, O. G., and P. H. Harvey (2000) Testing macro-evolutionary models using incomplete molecular phylogenies. \emph{Proc. R. Soc. Lond. B}, \bold{267}, 2267-2272.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{gammatest}}, \code{\link{ltt95}}
}
\examples{
trees<-pbtree(n=100,scale=100,nsim=10)
obj<-ltt(trees,plot=FALSE)
plot(obj,log="y",log.lineages=FALSE,main="lineage through time plots")
tree<-pbtree(b=1,d=0.25,t=4)
obj<-ltt(tree,gamma=FALSE)
obj
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{diversification analysis}

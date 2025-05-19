\name{bmPlot}
\alias{bmPlot}
\alias{sim.absorbing}
\alias{sim.reflective}
\title{Simulates and visualizes continuous trait evolution on a phylogeny}
\usage{
bmPlot(tree, type="BM", anc=0, sig2=1/1000, ngen=1000, ...)
sim.absorbing(tree, x0=0, sig2=1, bounds=c(-Inf,Inf), nsteps=1000, ...)
sim.reflective(tree, x0=0, sig2=1, bounds=c(-Inf,Inf), nsteps=1000, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{type}{the type of plot to create in \code{bmPlot}. See Description.}
	\item{anc}{the ancestral value for the root node in \code{bmPlot}.}
	\item{sig2}{the BM rate (variance of the stochastic diffusion process), \eqn{\sigma^2}.}
	\item{ngen}{number of generations for the simulation: will rescale the tree to this total length in \code{bmPlot}.}
	\item{x0}{the ancestral value of the root node in \code{sim.absorbing}.}
	\item{bounds}{absorbing or reflecting bounds in \code{sim.absorbing}. If \code{bounds=c(-Inf,Inf)} then the simulated process is unbounded.}
	\item{nsteps}{in \code{sim.absorbing}, the number of steps for the time discretization required for forward-in-time simulation.}
	\item{...}{arguments to be passed to different methods.}
}
\description{
	\code{bmPlot} conducts a discrete-time Brownian motion simulation on an input tree, plots the outcome, and returns the tip and internal node states to the user as a named vector. 
	
	\code{sim.absorbing} uses time discretization as well, but to conduct a simulation of either unbounded Brownian motion (if \code{bounds=c(-Inf,Inf)}), reflective bounded Brownian motion (if \code{reflective=TRUE}), or bounded Brownian motion with absorbing bounds.
	
	\code{sim.reflective} is a synonym for \code{sim.absorbing}, but with the optional argument \code{reflective=TRUE} by default.
}
\details{
	\code{bmPlot} will first rescale and round the branch lengths to integer length, if they are not already in integer values. If integer branch lengths are provided, the user should also set \code{ngen} to be equal to the total height of the tree in generations (and thus avoid rescaling). For \code{type="threshold"} the visualization is of the threshold model (Felsenstein 2012), in which the evolving character is liability and the segments of evolution are colored by their value for the threshold trait. If \code{type="threshold"} is used, the function requires at least one addition input: \code{thresholds}, a vector containing the ordered thresholds between states. The user can also provide the colors for plotting in \code{colors}. Note that one more color than threshold should be provided as one threshold implies two states; two thresholds, three states; etc. If no value for \code{colors} is provided, the function will recycle a set of four colors up to the number of times required by \code{thresholds}. Finally, the optional argument \code{return.tree=TRUE} will tell the function to return a list with the tip and note states and an object of class \code{"phylo"} with (for \code{type="threshold"}), the state for the threshold model through time mapped on the branches of the tree in discrete time.
	
	Optional arguments in \code{sim.absorbing} include \code{reflective}, which defaults to \code{reflective=FALSE}. This function uses discretization to conduct a forward-in-time simulation, but does not rescale the tree. \code{sim.reflective} calls \code{sim.absorbing} internally, but with the optional argument \code{reflective=TRUE} by default.
}
\value{
	\code{bmPlot} conducts and plots discrete time Brownian simulation and returns a vector containing the simulated states at internal nodes and tips of the tree. It also returns, by default (although this can be turned off) a tree with the branch lengths in discrete time and with a mapped discrete character (for \code{type="threshold"}).
	
	\code{sim.absorbing} returns a named vector of numeric values containing tip states.
}
\references{
	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bounded_bm}}, \code{\link{fastBM}}, \code{\link{pbtree}}, \code{\link{phenogram}}, \code{\link{threshBayes}}
}
\examples{
## set seed
set.seed(999)

## plot BM simulation on 12 taxon tree
tree<-pbtree(n=12)
par(mfrow=c(1,2),mar=c(5.1,4.1,4.1,0.1))
x<-bmPlot(tree,bty="l")
plotTree(tree,direction="upwards",
    mar=c(5.1,0.1,4.1,1.1),ftype="off")

## reset par to default values
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

## plot simulation of a threshold character
par(mfrow=c(1,2),mar=c(5.1,4.1,4.1,0.1))
tt<-bmPlot(tree,type="threshold",thresholds=c(0,1,2),
    bty="l")
plot(tt$tree,direction="upwards",
    mar=c(5.1,0.1,4.1,1.1),ftype="off",
    colors=setNames(c("black","red","blue"),
    letters[1:3]),lwd=3)

## reset par to default values
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{simulation}
\keyword{continuous character}

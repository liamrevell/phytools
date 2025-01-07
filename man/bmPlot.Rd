\name{bmPlot}
\alias{bmPlot}
\title{Simulates and visualizes discrete-time Brownian evolution on a phylogeny}
\usage{
bmPlot(tree, type="BM", anc=0, sig2=1/1000, ngen=1000, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{type}{the type of plot to create. See Description.}
	\item{anc}{the ancestral value for the root node.}
	\item{sig2}{the BM rate (variance of the Brownian evolution process), \eqn{\sigma^2}.}
	\item{ngen}{number of generations for the simulation: will rescale the tree to this total length.}
	\item{...}{arguments to be passed to different methods.}
}
\description{
	Conducts a discrete-time Brownian motion simulation on an input tree, plots the outcome, and returns the tip and internal node states to the user as a named vector. 
}
\details{
	The function will first rescale and round the branch lengths to integer length, if they are not already in integer values. If integer branch lengths are provided, the user should also set \code{ngen} to be equal to the total height of the tree in generations (and thus avoid rescaling). 
	
	For \code{type="threshold"} the visualization is of the threshold model (Felsenstein 2012), in which the evolving character is liability and the segments of evolution are colored by their value for the threshold trait. If \code{type="threshold"} is used, the function requires at least one addition input: \code{thresholds}, a vector containing the ordered thresholds between states. The user can also provide the colors for plotting in \code{colors}. Note that one more color than threshold should be provided as one threshold implies two states; two thresholds, three states; etc. If no value for \code{colors} is provided, the function will recycle a set of four colors up to the number of times required by \code{thresholds}. Finally, the optional argument \code{return.tree=TRUE} will tell the function to return a list with the tip and note states and an object of class \code{"phylo"} with (for \code{type="threshold"}), the state for the threshold model through time mapped on the branches of the tree in discrete time.
}
\value{
	This function conducts and plots discrete time Brownian simulation and returns a vector containing the simulated states at internal nodes and tips of the tree. 
	
	It also returns, by default (although this can be turned off) a tree with the branch lengths in discrete time and with a mapped discrete character (for \code{type="threshold"}).
}
\references{
	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastBM}}, \code{\link{pbtree}}, \code{\link{phenogram}}, \code{\link{threshBayes}}
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

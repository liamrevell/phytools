\name{phylo.heatmap}
\alias{phylo.heatmap}
\title{Creates a phylogenetic heat map}
\usage{
phylo.heatmap(tree, X, fsize=1, colors=NULL, standardize=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix containing data for multiple continuous characters in which \code{rownames} correspond to the tip labels of the tree.}
	\item{fsize}{an integer or vector of length 3 containing the font size for the tip labels, the trait labels, and the legend text. (If a single integer is supplied, then the value will be recycled.)}
	\item{colors}{a vector of colors to be passed to \code{\link{image}}. Can be a function call (e.g., \code{heat.colors(n=} \code{200)[200:1]}).}
	\item{standardize}{a logical value indicating whether or not to standardize each column of \code{X} to have the same variance & mean prior to analysis.}
	\item{...}{optional arguments. So far these include: \code{legend}, a logical value indicating whether or not to plot a figure legend (defaults to \code{legend=TRUE}); \code{labels}, a logical value indicating whether or not to plot trait labels (defaults to \code{labels=TRUE}); \code{split}, a numeric vector indicating the fraction of the horizontal dimension to use for the tree & heatmap, respectively (defaults to \code{split=c(0.5,0.5)}); \code{xlim}, \code{ylim}, & \code{mar}, defined as in \code{\link{par}}; and \code{ftype}, \code{lwd}, and \code{pts} as defined in \code{\link{plotSimmap}}.}
}
\description{
	Multivariate phylogenetic \code{\link{heatmap}} plot.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load data
data(anoletree)
data(anole.data)

## compute phylogenetic residuals
anole.data<-as.matrix(anole.data)
anole.resids<-cbind(anole.data[,1],
    phyl.resid(anoletree,anole.data[,1,drop=FALSE],
        anole.data[,2:ncol(anole.data)])$resid)
colnames(anole.resids)[1]<-"SVL"

## plot phylogenetic heatmap
phylo.heatmap(anoletree,anole.resids,
    split=c(0.7,0.3),fsize=c(0.4,0.8,0.8),
    standardize=TRUE,pts=FALSE)

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{continuous character}

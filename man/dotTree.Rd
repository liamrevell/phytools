\name{dotTree}
\alias{dotTree}
\alias{dot.legend}
\title{Creates a phylogenetic dot plot}
\usage{
dotTree(tree, x, legend=TRUE, method="plotTree", standardize=FALSE, ...)
dot.legend(x, y, min, max, Ntip, length=5, prompt=FALSE, method="plotTree", 
   ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{vector of trait values; or a matrix. If \code{x} is a vector it must have names that correspond to the tip labels of \code{tree}. If \code{x} is a matrix (and it probably should be a \emph{matrix}, not a data frame) then the row names of the matrix should correspond to the tip labels of the phylogeny. In the case of \code{dot.legend}, the \emph{x} coordinate of the legend.}
	\item{legend}{logical value indicating whether or not a legend should be plotted.}
	\item{method}{tree plotting method to be used internally. Will switch to \code{method="phylogram"} if the number of traits is greater than one. For \code{dot.legend}, it should be the method that was used for the plot.}
	\item{standardize}{a logical value indicating whether or not to standardize \code{x}, or each column of \code{x}, to have a mean of zero & variance of one prior to analysis.}
	\item{y}{\emph{y} coordinate of the legend.}
	\item{min}{minimum value for \code{dot.legend}.}
	\item{max}{maximum value for \code{dot.legend}.}
	\item{Ntip}{number of tips in the plotted tree for \code{dot.legend}.}
	\item{length}{length of legend.}
	\item{prompt}{logical value indicating whether or not to prompt for legend position.}
	\item{...}{optional arguments. In the case of \code{dotTree}, these will be passed to \code{\link{plotTree}} or a different internally used plotting function for \code{method="phylogram"}. See \code{\link{phylo.heatmap}} for more detail on these arguments. Other option for \code{dotTree} also include \code{data.type} (\code{"continuous"} or \code{"discrete"}), \code{colors}, \code{length}, for data type \code{"continuous"} the length of the legend in terms of plotted circles, \code{x.space}, the spacing of the columns in the plotted data matrix, and \code{leg.space}, the spacing of the legend dots (again, for \code{data.type="continuous"} only).}
}
\description{
	Creates a plot in which different sized dots/circles represent different tip values for a quantitative trait.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## log-transform trait data
log.mammal<-log(mammal.data)
## plot dotTree
dotTree(mammal.tree,log.mammal,fsize=0.7,
    standardize=TRUE,length=10)
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{continuous character}
\keyword{discrete character}

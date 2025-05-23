\name{phylomorphospace3d}
\alias{phylomorphospace3d}
\title{Creates three-dimensional phylomorphospace plot}
\usage{
phylomorphospace3d(tree, X, A=NULL, label=TRUE, control=list(), 
   method=c("dynamic","static"), ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{X}{an \code{n} x 3 matrix of tip values for two characters in \code{n} species.}
	\item{A}{an optional \code{m} x 3 matrix (for \code{m} nodes) of values for two traits at internal nodes in the tree - if not supplied, these values will be estimated using \code{\link{anc.ML}}.}
	\item{label}{logical value indicating whether to print tip labels next to terminal nodes in the plot (presently doesn't do anything, but labels can be dropped using \code{control}).}
	\item{control}{a list containing the following optional control parameters: \code{spin}: a logical value indicating whether to animate the plot when created; \code{axes}: a logical indicating whether to plot the axes; \code{box}: a logical value indicating whether to plot in box; \code{simple.axes}: logical value indicating whether to replace \code{box} and \code{axes} with simpler axes; \code{lwd}: line widths; \code{ftype}: font type ("off" turns off labels altogether); \code{col.edge} a vector of colors of length \code{nrow(tree$edge)}.}
	\item{method}{a string either \code{"dynamic"} for a dynamic (animated) plot created using \pkg{rgl}; or \code{"static"} for a flat 3D plot created using \code{scatterplot3d} and base graphics. The latter has the advantage of being very easy to export in standard format.}
	\item{...}{optional arguments to be passed to \code{scatterplot3d}. Most options not available. \code{angle} is an important option that does work here.}
}
\description{
	Creates a phylomorphospace plot in three dimensions.
}
\details{
	This function creates a phylomorphospace plot for three characters using the 3D visualization package, \pkg{rgl} (if available) or statically, by simulating 3D on a flat graphical device.
}
\value{
	This function creates a three dimensional phylomorphospace plot. The function returns a function from \code{spin3d} (for \code{method="dynamic"}); or a series of functions from \code{\link[scatterplot3d]{scatterplot3d}} (for \code{method="static"}).
}
\references{
	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Sidlauskas, B. (2008) Continuous and arrested morphological diversification in sister clades of characiform fishes: A phylomorphospace approach. \emph{Evolution}, \bold{62}, 3135-3156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load data
data(anoletree)
data(anole.data)

## run phylogenetic PCA
anole.pca<-phyl.pca(anoletree,anole.data)

\dontrun{
## create dynamic plot using rgl
phylomorphospace3d(anoletree,scores(anole.pca)[,1:3],
    control=list(spin=FALSE))}

## create static plot
par(cex=0.5) ## adjust font size
phylomorphospace3d(anoletree,scores(anole.pca)[,1:3],
    method="static",angle=-30)
par(cex=1) ## reset to default
}
\seealso{
	\code{\link{fancyTree}}, \code{\link{phenogram}}, \code{\link{phylomorphospace}}
}
\keyword{ancestral states}
\keyword{animation}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{plotting}
\keyword{continuous character}

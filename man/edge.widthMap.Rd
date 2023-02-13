\name{edge.widthMap}
\alias{edge.widthMap}
\alias{plot.edge.widthMap}
\title{Map continuous trait evolution on the tree}
\usage{
edge.widthMap(tree, x, ...)
\method{plot}{edge.widthMap}(x, max.width=0.9, legend="trait value", ...)
}
\arguments{
 	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{a numerical vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}. Or, for \code{plot.edge.widthMap}, an object of class \code{"edge.widthMap"}.}
	\item{max.width}{maximum edge width in plot units.}
	\item{legend}{label for the plot legend.}
	\item{...}{optional arguments - especially for the \code{plot} method. Perhaps the most important of these is \code{min.width}, which defaults to \code{0} but could probably be increased for many datasets and graphical devices. Other arguments are passed internally to \code{\link{plotTree}}.}
}
\description{
	Maps a discrete character onto the edges of the tree using variable edge widths.
}
\value{
	\code{edge.widthMap} returns an object of class \code{"edge.widthMap"}.
	
	\code{plot.edge.widthMap} can be used to plot this object.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{fastAnc}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))
## create "edge.widthMap" object
mammal.ewMap<-edge.widthMap(mammal.tree,ln.bodyMass,
	min.width=0.05)
## plot it
plot(mammal.ewMap,legend="log(body mass)")
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{continuous character}
\keyword{maximum likelihood}

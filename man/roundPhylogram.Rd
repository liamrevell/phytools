\name{roundPhylogram}
\alias{roundPhylogram}
\title{Plot a round phylogram}
\usage{
roundPhylogram(tree, fsize=1.0, ftype="reg", lwd=2, mar=NULL, offset=NULL, 
   direction="rightwards", type="phylogram", xlim=NULL, ylim=NULL, ...)
}
\arguments{
 	\item{tree}{an object of class "phylo" or "multiPhylo" containing one or multiple phylogenies.}
 	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
	\item{mar}{vector containing the margins for the plot to be passed to \code{\link{par}}. If not specified, the default margins are [0.1,0.1,0.1,0.1].}
	\item{offset}{offset for the tip labels.}
	\item{direction}{plotting direction. Only the option \code{direction="rightwards"} is presently supported.}
	\item{type}{plot type. Can be \code{"phylogram"} or \code{"cladogram"}. If \code{type="cladogram"} then the branch lengths are not necessary (and, indeed, are not used).}
	\item{xlim}{x-limits for the plot.}
	\item{ylim}{y-limits for the plot.}
	\item{...}{optional arguments.}
}
\description{
	Function plots one or multiple round phylograms.
}
\details{
	The underscore character \code{"_"} is automatically swapped for a space in tip labels, as in \code{\link{plotSimmap}}.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}, \code{\link{plotTree}}
}
\examples{
data(mammal.tree)
roundPhylogram(compute.brlen(mammal.tree),
    fsize=0.8,ftype="i")
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margin to default
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}

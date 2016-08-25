\name{phenogram}
\alias{phenogram}
\title{Plot phenogram (traitgram)}
\usage{
phenogram(tree, x, fsize=1.0, ftype="reg", colors=NULL, axes=list(),
	add=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}, with or without a mapped discrete character.}
	\item{x}{a vector containing the states at the tips \emph{or} the states at all the tips and the internal nodes of the tree.}
	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{colors}{colors for plotting the mapped character (if available) in \code{tree}. If no character is mapped on the tree, then a single color for all the branches of the tree can be provided.}
	\item{axes}{list of axis dimensions. Items are \code{time} and \code{trait}.}
	\item{add}{optional logical value indicating whether to add to an open plot. If \code{TRUE}, then new axes will not be plotted.}
	\item{...}{optional arguments including \code{xlim}, \code{ylim}, \code{log}, \code{main}, \code{sub}, \code{xlab}, \code{ylab}, \code{asp}, \code{type}, \code{lty}, \code{lwd}, \code{offset}, and \code{digits} are as in \code{\link{plot.default}} or \code{\link{par}}. Note that \code{axes} overrides \code{xlim} and \code{ylim}. \code{spread.labels} is a logical value indicating whether or not to minimize tip label overlap (default is \code{TRUE}); \code{spread.cost} is a numeric vector indicating the relative penalty to be used for label overlap and deviance, respectively (if \code{spread.labels=TRUE}); \code{spread.range} is the range over which to (potentially) spread the labels - note that if labels do not overlap, not all of that range will be used; finally, \code{link} is a numeric value by which to offset the tip labels, linking them to the tips with a dashed line (default is \code{0}, if \code{spread.labels=FALSE}, or 10-percent of the total tree length otherwise). The optional argument \code{offsetFudge} "fudges" the computation of label offset in scaling \code{xlim}. It is 1.37, which is the correct fudge in the Windows R GUI, but this may need to be changed in other systems. \code{hold} indicates whether (or not) the output to the graphical device should be held using \code{\link{dev.hold}} before plotting (defaults to \code{hold=TRUE}). \code{quiet} suppresses some system messages if set to \code{quiet=TRUE}.}
}
\description{
	Function plots a traitgram (Evans et al. 2009), that is, a projection of the phylogenetic tree in a space defined by phenotype (on the \emph{y} axis) and time (on the \emph{x}). If a discrete character is mapped on the tree this will also be plotted.
}
\details{
	For \code{spread.labels=TRUE} numerical optimization is performed to optimize the distribution of the labels vertically, where the solution depends on the vector \code{spread.cost} containing the cost of overlap (first) and the cost of deviation from the vertical position of the tip. Note that because this is done via numerical optimization, plotting may hang briefly while the best solution is found (especially for large trees).
}
\value{
	Plots a traitgram, optionally with a mapped discrete character, and (invisibly) returns a matrix containing the coordinates of the plotted tip labels.
}
\references{
	Evans, M. E. K., Smith, S. A., Flynn, R. S., Donoghue, M. J. (2009) Climate, niche evolution, and diversification of the "bird-cage" evening primroses (Oenothera, sections Anogra and Kleinia). \emph{American Naturalist}, \bold{173}, 225-240.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
\dontrun{
tree<-pbtree(n=20,scale=2)
x<-fastBM(tree)
phenogram(tree,x)
# or, simulate a discrete character history
tree<-sim.history(tree,Q=matrix(c(-1,1,1,-1),2,2),anc="1")
# simulate in which the rate depends on the state
x<-sim.rates(tree,c(1,10))
phenogram(tree,x)
# now use spread.labels
tree<-pbtree(n=40)
x<-fastBM(tree)
phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0))
}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}

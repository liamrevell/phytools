\name{phenogram}
\alias{phenogram}
\title{Plot traitgram (phenogram)}
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
	\item{...}{optional arguments including \code{xlim}, \code{ylim}, \code{log}, \code{main}, \code{sub}, \code{xlab}, \code{ylab}, \code{asp}, \code{type}, \code{lty}, \code{lwd}, \code{offset}, and \code{digits} are as in \code{\link{plot.default}} or \code{\link{par}}. Note that \code{axes} overrides \code{xlim} and \code{ylim}. \code{spread.labels} is a logical value indicating whether or not to minimize tip label overlap (default is \code{TRUE}); \code{spread.cost} is a numeric vector indicating the relative penalty to be used for label overlap and deviance, respectively (if \code{spread.labels=TRUE}); \code{spread.range} is the range over which to (potentially) spread the labels - note that if labels do not overlap, not all of that range will be used; finally, \code{link} is a numeric value by which to offset the tip labels, linking them to the tips with a dashed line (default is \code{0}, if \code{spread.labels=FALSE}, or 10\% of the total tree length otherwise). The optional argument \code{offsetFudge} "fudges" the computation of label offset in scaling \code{xlim}. It is 1.37, which is the correct fudge in the Windows R GUI, but this may need to be changed in other systems. \code{hold} indicates whether (or not) the output to the graphical device should be held using \code{\link{dev.hold}} before plotting (defaults to \code{hold=TRUE}). \code{quiet} suppresses some system messages if set to \code{quiet=TRUE}.}
}
\description{
	Plots a phylogenetic traitgram (Evans et al., 2009).
}
\details{
	Function plots a traitgram (Evans et al. 2009), that is, a projection of the phylogenetic tree in a space defined by phenotype (on the \emph{y} axis) and time (on the \emph{x}). If a discrete character is mapped on the tree this will also be plotted.

	For \code{spread.labels=TRUE} numerical optimization is performed to optimize the distribution of the labels vertically, where the solution depends on the vector \code{spread.cost} containing the cost of overlap (first) and the cost of deviation from the vertical position of the tip. Note that because this is done via numerical optimization, plotting may hang briefly while the best solution is found (especially for large trees).
}
\value{
	Plots a traitgram, optionally with a mapped discrete character, and (invisibly) returns a matrix containing the coordinates of the plotted tip labels.
}
\references{
	Evans, M. E. K., Smith, S. A., Flynn, R. S., Donoghue, M. J. (2009) Climate, niche evolution, and diversification of the "bird-cage" evening primroses (Oenothera, sections Anogra and Kleinia). \emph{American Naturalist}, \bold{173}, 225-240.
	
	Revell, L. J. (2013) Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.

	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
	
	Revell, L. J., K. Schliep, E. Valderrama, and J. E. Richardson (2018) Graphs in phylogenetic comparative analysis: Anscombe's quartet revisited. \emph{Methods in Ecology and Evolution}, \bold{9}, 2145-2154.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)

## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))

## plot traitgram
phenogram(mammal.tree,ln.bodyMass,ftype="i",
    spread.cost=c(1,0),fsize=0.7,color=palette()[4],
    xlab="time (ma)",ylab="log(body mass)",las=1)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{continuous character}

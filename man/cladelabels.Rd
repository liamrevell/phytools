\name{cladelabels}
\alias{cladelabels}
\alias{arc.cladelabels}
\title{Add labels to subtrees of a plotted phylogeny}
\usage{
cladelabels(tree=NULL, text, node, offset=NULL, wing.length=NULL, cex=1,
   orientation="vertical")
arc.cladelabels(tree=NULL, text, node=NULL, ln.offset=1.02,
   lab.offset=1.06, cex=1, orientation="curved", stretch=1, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}. If not supplied, the function will obtain the last plotted phylogeny from the environmental variable \code{last_plot.phylo}.}
	\item{text}{desired clade label text.}
	\item{node}{node number for the most recent common ancestor of members of the clade. For \code{arc.cladelabels} this defaults to \code{NULL} which means that the node of the clade to be labeled should be specified interactively (that is, by clicking on the graphical device).}
	\item{offset}{offset (as a multiplier of character width) for the label. Defaults to \code{offset=1} if \code{tree} is supplied or \code{offset=8} otherwise.}
	\item{wing.length}{length of the wings to add to the top & bottom of the label bar (in character widths).}
	\item{cex}{character expansion factor.}
	\item{orientation}{orientation of the text. Can be \code{orientation = "vertical"} (the default) or \code{"horizontal"}.}
	\item{ln.offset}{line offset (as a function of total tree height) for \code{arc.cladelabels}.}
	\item{lab.offset}{label offset for \code{arc.cladelabels}.}
	\item{stretch}{argument for \code{arc.cladelabels} to be passed to \code{arctext}.}
	\item{...}{optional arguments for \code{arc.cladelabels}.}
}
\description{
	This function adds clade labels to a plotted tree.
}
\details{
	\code{cladelabels} presently works only for rightward facing plotted phylogenies - but no warning will be returned if your tree does not conform to this requirement!
	
	\code{arc.cladelabels} is designed to do a similar thing to \code{cladelabels}, but for plotted fan trees. This function checks to ensure that the most recently plotted tree was plotted with \code{type="fan"} style.
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{nodelabels}}
}
\keyword{phylogenetics}
\keyword{utilities}

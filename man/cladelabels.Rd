\name{cladelabels}
\alias{cladelabels}
\title{Add labels to subtrees of a plotted phylogeny}
\usage{
cladelabels(tree=NULL, text, node, offset=NULL, wing.length=NULL, cex=1,
	orientation="vertical")
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}. If not supplied, the function will obtain the last plotted phylogeny from the environmental variable \code{last_plot.phylo}.}
	\item{text}{desired clade label text.}
	\item{node}{node number for the most recent common ancestor of members of the clade.}
	\item{offset}{offset (as a multiplier of character width) for the label. Defaults to \code{offset=1} if \code{tree} is supplied or \code{offset=8} otherwise.}
	\item{wing.length}{length of the wings to add to the top & bottom of the label bar (in character widths).}
	\item{cex}{character expansion factor.}
	\item{orientation}{orientation of the text. Can be \code{orientation = "vertical"} (the default) or \code{"horizontal"}.}
}
\description{
	This function adds clade labels to a plotted tree.
}
\details{
	This function presently works only for rightward facing plotted phylogenies - but no warning will be returned if your tree does not conform to this requirement!
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{nodelabels}}
}
\keyword{phylogenetics}
\keyword{utilities}

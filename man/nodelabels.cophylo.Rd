\name{nodelabels.cophylo}
\alias{nodelabels.cophylo}
\alias{tiplabels.cophylo}
\alias{edgelabels.cophylo}
\title{Add labels to a plotted "cophylo" object}
\usage{
nodelabels.cophylo(..., which=c("left","right"))
edgelabels.cophylo(..., which=c("left","right"))
tiplabels.cophylo(..., which=c("left","right"))
}
\arguments{
	\item{...}{arguments to be passed to \code{\link{nodelabels}}, \code{\link{edgelabels}}, or \code{\link{tiplabels}}.}
	\item{which}{argument indicated which of the two plotted trees (the \code{"left"} or \code{"right"} tree) to be used.}
}
\description{
	This function adds node, edge, or tip labels to the plotted trees of a \code{"cophylo"} object.
}
\details{
	Note that the order of tips, edges, and nodes may be different in the object of class \code{"cophylo"} than they are in the original input trees, particularly if \code{cophylo(...,rotate=TRUE)} was used.
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{cophylo}}, \code{\link{edgelabels}}, \code{\link{nodelabels}}, \code{\link{tiplabels}}
}
\keyword{phylogenetics}
\keyword{utilities}

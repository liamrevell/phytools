\name{ltt95}
\alias{ltt95}
\alias{plot.ltt95}
\title{Creates a (1-\eqn{\alpha})\% CI for a set of LTTs}
\usage{
ltt95(trees, alpha=0.05, log=FALSE, method=c("lineages","times"),
   mode=c("median","mean"), ...)
\method{plot}{ltt95}(x, ...)
}
\arguments{
	\item{trees}{is an object of class \code{"multiPhylo"} containing a list of phylogenetic trees.}
	\item{alpha}{confidence level. Defaults to \code{alpha=0.05}. \code{alpha=0} will mean that the interval around \emph{all} trees in the set will be plotted.}
	\item{log}{logical value indicating whether or not to plot on the semi-log scale.}
	\item{method}{plot the CI on the number of lineages given time (\code{"lineages"}); or on times given a number of lineages (\code{"times"}).}
	\item{mode}{plot the median or mean LTT.}
	\item{x}{object of class \code{"ltt95"} for plotting method.}
	\item{...}{optional arguments to be used by \code{ltt95} or the plotting method. So far: \code{res} gives the number of time-steps (defaults to \code{res=100}); \code{xaxis} (\code{"standard"}, \code{"negative"}, or \code{"flipped"}) determines the scale (time from the root, time back from the present, or time from the present) of the x-axis of the plot; \code{lend} determines the line end type (as in \code{\link{par}}); \code{shaded} determines whether to plot the (1-\eqn{\alpha})\% CI using dotted lines (if \code{FALSE}) or shading (if \code{TRUE}); and \code{bg} is the background color for shading if \code{shaded=TRUE}.}
}
\description{
	This function computes LTT plots for a set of trees & plots a (1-\eqn{\alpha})\% CI by various methods. (See \code{\link{ltt}} for more details.)
}
\details{
	This function creates a plot and invisibly returns an object of class \code{"ltt95"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ltt}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{diversification}

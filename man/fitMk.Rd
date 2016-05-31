\name{fitMk}
\alias{fitMk}
\alias{plot.fitMk}
\title{Fits Mk model}
\usage{
fitMk(tree, x, model="SYM", fixedQ=NULL, ...)
\method{plot}{fitMk}(x, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names. In the case of \code{plot.fitMk}, an object of class \code{"fitMk"}.}
	\item{model}{model. See \code{make.simmap} or \code{ace} for details.}
	\item{fixedQ}{fixed value of transition matrix \code{Q}, if one is desired.}
	\item{...}{optional arguments, including \code{pi}, the prior distribution at the root node (defaults to \code{pi="equal"}). For \code{plot} method optional arguments include (but may not be limited to): \code{signif}, the number of digits for the rates to be plotted; \code{main}, a character vector of length two with the headings for each subplot; \code{cex.main}, \code{cex.traits}, and \code{cex.rates}, font sizes for the various text elements of the plot; and \code{show.zeros}, a logical argument specifying whether or not to plot arrows with the ML estimated transition rate is not different from zero (with tolerance specified by the optional argument \code{tol}).}
}
\description{
	This function fits a so-called M\emph{k} model for discrete character evolution. It is primarily designed to be used inside of \code{make.simmap}.
}
\value{
	An object of class \code{"fitMk"}. In the case of \code{plot.fitMk}, a plotted M\emph{k} model.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}

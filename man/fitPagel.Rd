\name{fitPagel}
\alias{fitPagel}
\alias{plot.fitPagel}
\title{Function to test for correlated evolution of binary traits}
\usage{
fitPagel(tree, x, y, method="fitMk", model="ARD", dep.var="xy", ...)
\method{plot}{fitPagel}(x, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of phenotypic values for a binary trait for the species in \code{tree}; or a matrix in which the rows of \code{x} give the probability of being in each column state. (The latter option is only supported for \code{method="fitMk"}.) For S3 \code{plot} method, an object of class \code{"fitPagel"}.}
	\item{y}{a second binary character for the species in \code{tree}; or a matrix in which the rows give the probability of being in each column state.}
	\item{method}{function to use for optimization (defaults to \code{method="fitMk"}). Other options are \code{"ace"} to use the \code{\link[ape]{ace}} function in \pkg{ape} for optimization, or to \code{"fitDiscrete"} (if the \pkg{geiger} package is installed) to use \pkg{geiger}'s \code{fitDiscrete} for optimization.}
	\item{model}{model of evolution for the individual characters. Can be \code{model="ER"}, \code{"SYM"} (equivalent to \code{"ER"} in this case),
	and \code{"ARD"}.}
	\item{dep.var}{dependent variable. If \code{dep.var="xy"} than the rate of substitution in \code{x} depends on \code{y} & vice versa. If \code{dep.var="x"} than the substitution rate in \code{x} depends on \code{y}, but not the converse. Finally, if \code{dep.var="y"} than the rate of substitution in \code{y} depends on \code{x}, but not the converse.}
	\item{...}{optional arguments to be passed to \code{\link{fitMk}}, \code{\link[ape]{ace}}, or \code{fitDiscrete}. For \code{plot} method optional arguments include (but may not be limited to): \code{signif}, the number of digits for the rates to be plotted; \code{main}, a character vector of length two with the headings for each subplot; \code{cex.main}, \code{cex.sub}, \code{cex.traits}, and \code{cex.rates}, font sizes for the various text elements of the plot; and \code{lwd.by.rate}, a logical argument specifying whether or not to scale arrow line widths in proportion to the estimated rates.}
}
\description{
	Fits Pagel's (1994) model for the correlated evolution of two binary characters.
}
\details{
	\code{fitPagel} fits both an independent evolution model, as well as Pagel's (1994) binary dependent model, and compares them with a likelihood-ratio test.
	
	\code{plot.fitPagel} plots the fitted models using arrows.
}
\value{
	An object of class \code{"fitPagel"} which contains the optimized matrices under an independence & a dependence model, log-likelihoods, a likelihood ratio, and a P-value for the independence model based on a chi-squared test.
	
	\code{plot.fitPagel} creates a plot showing the different fitted models with arrows.
}
\references{
	Pagel, M. (1994) Detecting correlated evolution on phylogenies: A general method for the comparative analysis of discrete characters. \emph{Proceedings of the Royal Society B}, \bold{255}, 37-45.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{ace}}, \code{\link{fitMk}}, \code{\link{make.simmap}}
}
\examples{
\dontrun{
## load data
data(bonyfish.tree)
data(bonyfish.data)

## extract discrete characters
spawning_mode<-setNames(bonyfish.data$spawning_mode,
    rownames(bonyfish.data))
paternal_care<-setNames(bonyfish.data$paternal_care,
    rownames(bonyfish.data))

## fit correlational model
bonyfish.pagel<-fitPagel(bonyfish.tree,paternal_care,
    spawning_mode)

## test for a difference between models
anova(bonyfish.pagel)

## plot fitted models
plot(bonyfish.pagel,lwd.by.rate=TRUE)

## reset par
par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{discrete character}
\keyword{maximum likelihood}

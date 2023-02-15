\name{plotThresh}
\alias{plotThresh}
\title{Tree plotting with posterior probabilities of ancestral states from the threshold model}
\usage{
plotThresh(tree, x, mcmc, burnin=NULL, piecol, tipcol="input", legend=TRUE, 
   ...)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{a named vector containing discrete character states; or a matrix containing the tip species, in rows, and probabilities of being in each state, in columns.}
	\item{mcmc}{list object returned by \code{\link{ancThresh}}.}
	\item{burnin}{number of generations (not samples) to exclude as burn in; if \code{NULL} then 20\% of generations are excluded as burn-in.}
	\item{piecol}{a named vector containing the colors for the posterior probabilities plotted as pie charts at internal nodes.}
	\item{tipcol}{a string indicating  whether the tip colors should be based on the input data (\code{"input"}) or sampled tip liabilities (\code{"estimated"}). These will only differ if there is uncertainty in the tip states.}
	\item{legend}{logical value or text to be plotted in the legend.}
	\item{...}{other arguments to be passed to \code{\link{plot.phylo}} - \code{label.offset} should be >0 so that tip labels and species names do not overlap.}
}
\description{
	Plots estimated posterior probabilities at nodes under the threshold model.
}
\details{
	This function uses the object returned by \code{\link{ancThresh}} to plot the posterior probabilities of ancestral states under the threshold model.
	
	It is also called internally by \code{\link{ancThresh}}.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, bold{68}, 743-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}, \code{\link{plot.phylo}}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
\keyword{discrete character}

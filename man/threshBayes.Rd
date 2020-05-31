\name{threshBayes}
\alias{threshBayes}
\title{Threshold model using Bayesian MCMC}
\usage{
threshBayes(tree, X, types=NULL, ngen=10000, control=list(), ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix or data frame containing values for a discrete character and a continuous character; or two discrete characters. The row names of \code{X} should be species names.}
	\item{types}{a vector of length \code{ncol(X)} containing the data types for each column of \code{X}, for instance \code{c("discrete","continuous")}.}
	\item{ngen}{a integer indicating the number of generations for the MCMC.}
	\item{control}{a list of control parameters for the MCMC. Control parameters include: \code{sample}, the sampling interval for the MCMC; \code{propvar}, a vector containing (in this order) proposal variances for the two rates (if the type is \code{"discrete"} this will be ignored), the two ancestral states, and the correlation; \code{propliab}, a single proposal variance for the liabilities; \code{pr.mean}, a vector for the mean of the prior probability distributions for each parameter, in the same order as \code{propvar}; \code{pr.liab}, currently ignored; \code{pr.var}, a vector with variances for the prior densities for each parameter, in the same order as \code{pr.mean} - note that for the rates we use an exponential distribution so the first two means are currently ignored; and \code{pr.vliab} currently ignored.}
	\item{...}{other optional arguments. The argument \code{auto.tune} adjusts the proposal distribution to target an proposal acceptance rate (of 0.234, by default).}
}
\description{
	This function uses Bayesian MCMC to fit the quantitative genetics threshold model (Felsenstein 2012) to data for two discrete characters or one discrete and one continuous character.
}
\value{
	This function returns an object of class \code{"threshBayes"} consisting of a list with at least the following two elements: \code{par} a matrix containing the posterior sample for the model parameters (evolutionary rates, ancestral states, and correlation); \code{liab} a matrix containing the posterior sample of the liabilities. For continuous characters, the liabilities are treated as known and so the posterior samples are just the observed values.
}
\details{
	The \code{plot} method for the object class can be used to generate a three panel plot showing the likelihood profile, the mean acceptance rates (using a sliding window), and a profile plot for the correlation coefficient, \emph{r}.

	The \code{density} method for the object can be used to plot a posterior density of the correlation coefficient, \emph{r}. This posterior density is of class \code{"density.threshBayes"} which can also be plotted using an \code{plot} method.
	
	Discrete characters must be binary, but can be coded as any factor.
}	
\references{
	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.Bayes}}, \code{\link{bmPlot}}, \code{\link{evol.rate.mcmc}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
\keyword{continuous character}
\keyword{discrete character}

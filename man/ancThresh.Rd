\name{ancThresh}
\alias{ancThresh}
\title{Ancestral character estimation under the threshold model using Bayesian MCMC}
\usage{
ancThresh(tree, x, ngen=10000, sequence=NULL, method="mcmc", 
   model=c("BM","OU","lambda"), control=list(), ...)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{a named vector containing discrete character states; or a matrix containing the tip species, in rows, and probabilities of being in each state, in columns.}
	\item{ngen}{number of generations to run the MCMC.}
	\item{sequence}{assumed ordering of the discrete character state. If not supplied and \code{x} is a vector then numerical/alphabetical order is assumed; if not supplied and \code{x} is a matrix, then the column order of \code{x} is used.}
	\item{method}{only method currently available is \code{"mcmc"}.}
	\item{model}{model for the evolution of the liability. Options are \code{"BM"} (Brownian motion, the default), \code{"OU"} (Ornstein-Uhlenbeck), or \code{"lambda"} (the lambda model).}
	\item{control}{list containing the following elements: \code{sample}, the sampling interval; \code{propliab} variance of the proposal distribution for liabilities; \code{propthresh} variance on the proposal distribution for the thresholds; \code{propalpha} variance on the proposal distribution for \code{alpha} (for \code{model="OU"}); \code{pr.anc} prior probability distribution on the ancestral states for each node, in a matrix - not all nodes need to be supplied; \code{pr.th} prior density on the thresholds; \code{burnin} number of generations to exclude for burn-in when plotting posterior probabilities on the tree; \code{plot} logical value indicating whether or not to plot the posterior probabilities; \code{print} logical value indicating whether or not to print the state of the MCMC; \code{piecol} colors for the posterior probabilities plotted as pie charts at internal nodes; and \code{tipcol} which indicates whether the tip colors should be based on the input data (\code{"input"}) or sampled tip liabilities (\code{"estimated"}). These will only differ if there is uncertainty in the tip states.}
	\item{...}{additional arguments to be passed to \code{\link{plotThresh}} (called internally).}
}
\description{
	This function uses Bayesian MCMC to estimate ancestral states and thresholds for a discrete character under the threshold model from quantitative genetics (Felsenstein 2012).
}
\details{
	\code{print} and \code{plot} S3 methods are now available for the object class \code{"ancThresh"}.
}
\value{
	This function returns an object of class \code{"ancThresh"} containing the posterior sample from our analysis, although with other components.
}
\references{
	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, bold{68}, 743-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.Bayes}}, \code{\link{threshBayes}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}

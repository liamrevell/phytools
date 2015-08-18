\name{make.simmap}
\alias{make.simmap}
\title{Simulate stochastic character maps on a phylogenetic tree or trees}
\usage{
make.simmap(tree, x, model="SYM", nsim=1, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}, or a list of trees as an object of class \code{"multiPhylo"}.}
	\item{x}{a vector containing the tip states for a discretely valued character, or a matrix containing the prior probabilities of tip states in rows.}
	\item{model}{a character string containing the model - options as in \code{\link{ace}}.}
	\item{nsim}{number of simulations. If \code{tree} is an object of class \code{"multiPhylo"}, then \code{nsim} simulations will be conducted per tree.}
	\item{...}{optional arguments. So far, \code{pi} gives the prior distribution on the root node of the tree - options are \code{"equal"}, \code{"estimated"}, or a vector with the frequencies. If \code{pi="estimated"} then the stationary distribution is estimated by numerically solving \code{pi*Q=0} for \code{pi}, and this is used as a prior on the root. Defaults to \code{pi="equal"} which results in the root node being sampled from the conditional scaled likelihood distribution at the root. \code{message} tells whether or not to print a message containing the rate matrix, \emph{Q} and state frequencies. Defaults to \code{message=TRUE}. For optional argument \code{Q="mcmc"} the mean value of \code{Q} from the posterior sample is printed. \code{tol} gives the tolerance for zero elements in \code{Q}. (Elements less then \code{tol} will be reset to \code{tol}). \code{Q} can be a string (\code{"empirical"} or \code{"mcmc"}), or a fixed value of the transition matrix, \code{Q}. If \code{"empirical"} than a single value of \code{Q}, the most likely value, is used for all simulations. If \code{"mcmc"}, then \code{nsim} values of \code{Q} are first obtained from the posterior distribution for \code{Q} using Bayesian MCMC, then a simulated stochastic character map is generated for each value of \code{Q}. \code{vQ} a single numeric value or a vector containing the (normal) sampling variances for the MCMC. The order of \code{vQ} is assumed to be in the order of the \code{index.matrix} in \code{\link{ace}} for the chosen model. \code{prior} a list containing \code{alpha} and \code{beta} parameters for the gamma prior distribution on the transition rates in \code{Q}. Note that \code{alpha} and \code{beta} can be single values or vectors, if different priors are desired for each value in \code{Q}. As for \code{vQ}, the order of \code{prior} is assumed to be the order of \code{index.matrix} in \code{\link{ace}}. \code{prior} can also be given the optional logical value \code{use.empirical} which tells the function whether or not to give the prior distribution the empirical mean for \code{Q}. If \code{TRUE} then only \code{prior$beta} is used and \code{prior$alpha} is set equal to \code{prior$beta} times the empirical mean of \code{Q}. \code{burnin} and \code{samplefreq} are burn-in and sample frequency for the MCMC, respectively.}
}
\description{
	This function performs stochastic mapping using several methods.

	For \code{Q="empirical"}, it first fits a continuous-time reversible Markov model for the evolution of \code{x} and then simulates stochastic character histories using that model and the tip states on the tree. This is the same procedure that is described in Bollback (2006), except that simulation is performed using a fixed value of the transition matrix, \code{Q}, instead of by sampling \code{Q} from its posterior distribution.

	For \code{Q="mcmc"}, it first samples \code{Q} \code{nsim} times from the posterior probability distribution of \code{Q} using MCMC, then it simulates \code{nsim} stochastic maps conditioned on each sampled value of \code{Q}.

	For \code{Q} set to a matrix, it samples stochastic mappings conditioned on the fixed input matrix.
}
\details{
	Uses code modified from \code{\link{ace}} (by Paradis et al.) to perform Felsenstein's pruning algorithm & compute the likelihood.

	As of phytools>=0.2-33 \code{x} can be a vector of states or a matrix containing the prior probabilities of tip states in rows. In this case the column names of \code{x} should contain the states, and the row names should contain the tip names.

	Note that there was a small (but potentially significant) bug in how node states were simulated by \code{make.simmap} in versions of phytools<=0.2-26. Between phytools 0.2-26 and 0.2-36 there was also a bug for asymmetric models of character change (e.g., \code{model="ARD"}). Finally, between phytools 0.2-33 and phytools 0.2-47 there was an error in use of the conditional likelihoods for the root node, which caused the root node of the tree to be sampled incorrectly. All of these issues should be fixed in the present version.

	\code{Q="mcmc"} and \code{Q} set to a fixed value were introduced to phytools >= 0.2-53. As of the present version of phytools, this method is still somewhat experimental & should be used with caution.

	If \code{tree} is an object of class \code{"multiPhylo"} then \code{nsim} stochastic maps are generated for each input tree.
}
\value{
 	A modified phylogenetic tree of class \code{"phylo"} (or a modified \code{"multiPhylo"} object, for \code{nsim} > 1) with the following additional elements:
 	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
 	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
	\item{Q}{the assumed or sampled value of \code{Q}.}
	\item{logL}{the log-likelihood of the assumed or sampled \code{Q}.}
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

 	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{brownieREML}}, \code{\link{countSimmap}}, \code{\link{describe.simmap}}, \code{\link{evol.vcv}}, \code{\link{plotSimmap}}, \code{\link{read.simmap}}, \code{\link{write.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}

